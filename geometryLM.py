"""
geometryLM.py a fork of geometry3d.py for list mode code


1. 3-d density profile (Dx * Dy * Dz) and 2-d (Dx * Dy) detector.

2. Hexagonal close-packed lattice of Gaussian basis functions.  See
   http://en.wikipedia.org/wiki/Close-packing.

3. Fit density profile that starts at z=Z0
                            _____________________
                           /                    /|
                          /                    / |
                        0/                   Dx  |                
z=Dz..................../..................../................... 
                        |           :        |   |                 
                        |           :        |   |                  
                        |           Dz       |   /                    
                        |           :        |  /                    
                        |.......Dx..:........| / Dy                
                        |           :        |/
z=0 ...................._____________________/.................  <- Detector

Here is a top view of the first two layers of the hexagonal close
packed pattern:

   (0,2,0)     (1,2,0)     (2,2,0)     (2,2,0)     (3,2,0)     (4,2,0)
   (0,1,1)     (1,1,1)     (2,1,1)     (2,1,1)     (3,1,1)     (4,1,1)
         (0,1,0)     (1,1,0)     (2,1,0)     (2,1,0)     (3,1,0)     (4,1,0)
         (0,0,1)     (1,0,1)     (2,0,1)     (2,0,1)     (3,0,1)     (4,0,1)
   (0,0,0)     (1,0,0)     (2,0,0)     (2,0,0)     (3,0,0)     (4,0,0)

See http://en.wikipedia.org/wiki/Spherical_coordinate_system and
http://docs.scipy.org/doc/scipy/reference/sparse.linalg.html
"""
import math, random, numpy, scipy, scipy.io, os, time, scipy.sparse

Dtheta = math.pi/4

exp_cut = -7.0  # Drop terms that are down by e**-7 = 0.0009
def d2oca(c,    # Center: Array of (x,y,z)
         theta, # Angle of incoming ray from zenith in radians
         phi,   # Azimuthal angle of ray in radians
         p      # Point where ray crosses z=0 plane: Array of (x,y,z)
         ):
    """Calculate squared distance between point c and trajectory of
    incoming ray given by (theta,phi,p).
    """
    dx = c-p
    u = numpy.array([math.sin(theta)*math.cos(phi),
                     math.sin(theta)*math.sin(phi),
                     math.cos(theta)
                     ]) # Unit vector along ray
    r = numpy.dot(dx,u) # Radius at point of closest approach
    d = dx-r*u          # difference between poca and c
    return numpy.dot(d,d)
def map_gen():
    """ Return an array index2angles so that given an index k from one
    of Andrew Green's files,
    
    theta,phi = index2angles[k]

    produces the spherical coordinates of the rays.
    """
    ntheta = 70
    Dtheta = 0.39788705*math.pi # Total range
    dtheta = Dtheta/ntheta
    index2angles = []
    j = 0
    for i in xrange(ntheta):
        theta = dtheta*(i+.5)
        Lphi = math.sin(theta)*math.pi*2
        nphi = int(Lphi/dtheta+1)
        dphi = 2*math.pi/nphi
        if False:
            print '%2d   %6.4f   %6.4f   %3d   %4d'%(i,theta/math.pi,
                                                     dphi/math.pi,nphi,j)
        for k in xrange(nphi):
            phi = dphi*(k+.5)
            index2angles.append([theta,phi])
            j += 1
    return numpy.array(index2angles)
def ReadGreen(name,N=-1):
    """ Read one of Andrew Green's files.  Return an Nx5 numpy array
    with one row per bin with entries [x,y,theta,phi,n].

    Takes 7.6 seconds to count 16,000,000 entries and 235 seconds to
    read them.
    """
    import re
    I2As = map_gen()
    colon = re.compile(':')
    pound = re.compile('#')
    if N < 0: #Count entries in file.
        N=0
        File = open(name,'r')
        for line in File.xreadlines():
            if not pound.search(line):
                N += 1
    RV = numpy.empty((N,5))
    i=0
    File = open(name,'r')
    for line in File.xreadlines():
        if pound.search(line):
            x,y =  [float(t)/2.0-.5 for t in pound.split(line)]
        else:
            AngleI,count = [int(t) for t in colon.split(line)]
            AngleI -= 1
            theta,phi = I2As[AngleI]
            RV[i,:] = [x,y,theta,phi,float(count)]
            i += 1
            if i >= N:
                break
    return RV
    
def Reframe(Gamma,  # Array [Gx0, Gx1, Gx2, Gtheta, Gphi]
            In_vecs # N x 5 array. [x,y,theta,phi,n] each row
            ):
    """ Gamma represents the position and orientation of an instrument
    in a lab reference frame.  This function translates the array
    In_vecs from the instrument frame to the lab frame (in the x_c[2]
    = 0 plane).  Letting "a" denote an input vector, "b" an
    intermediate and "c" the vector in the lab frame, I can write the
    steps as follows:

    x_b = A_G * x_a + B_G
    v_b = A_G * F(theta_a,phi_a)   A unit vector on a sphere
    beta = x_b[2]/v_b[2]
    x_c = x_b - beta * v_b
    (theta_c,phi_c) = F_inv(v_b)
    """
    Gtheta,Gphi = Gamma[3:]
    R_theta = numpy.matrix([[math.cos(Gtheta),0,-math.sin(Gtheta)],[0,1.0,0],
                            [math.sin(Gtheta),0,math.cos(Gtheta)]])
    R_phi = numpy.matrix([[math.cos(Gphi),-math.sin(Gphi),0],
                          [math.sin(Gphi),math.cos(Gphi),0],[0,0,1.0]])
    A_G = R_phi*R_theta
    B_G = numpy.array(Gamma[0:3])
    def F(theta,phi):
        return numpy.array([math.cos(phi)*math.sin(theta),
                            math.sin(phi)*math.sin(theta),
                            math.cos(theta)])
    def F_inv(v):
        r = numpy.dot(v,v)
        assert r < 1.0+1.0e-6,'r=%f'%r
        assert r > 1.0-1.0e-6,'r=%f'%r
        phi = math.atan2(v[1],v[0])
        theta = math.acos(v[2])
        return (theta,phi)
    RV = numpy.empty(In_vecs.shape)
    v = numpy.zeros(3)
    for i in xrange(len(RV)):
        v[0:2] = In_vecs[i,0:2]
        x_b = numpy.dot(A_G.A,v) + B_G
        v_b = numpy.dot(A_G,F(In_vecs[i,2],In_vecs[i,3])).A.reshape(-1)
        RV[i,0:2] = (x_b - v_b*(x_b[2]/v_b[2]))[0:2]
        RV[i,2:4] = F_inv(v_b)
        RV[i,4] = In_vecs[i,4]
    return RV
def rebin(data,  # Scipy array Nx5 [x,y,theta,phi,n]
          bs     # Float bin size
          ):
    if bs < .9:
        return data
    dataD = {}
    for row in data:
        qx = bs*int(row[0]/bs) # Quantized x value
        qy = bs*int(row[1]/bs) # Quantized y value
        key = (qx,qy)
        if dataD.has_key(key):
            dataD[key][4] += row[4]
        else:
            dataD[key] = row
    RV = numpy.empty((len(dataD),5))
    keys = dataD.keys()
    for i in xrange(len(keys)):
        RV[i,:] = dataD[keys[i]]
    return RV
def join_data(data_in,      # List of segments
              shifts,       # Array of (x,y) shifts for each segment
              Ts,           # Array of exposure times
              binspace=0.0  # Resolution in cm for rebinning
              ):
    """ Take N_seg segments of data in data_dict and return two arrays
    called dataT and dataB that will be used as traget and background
    data respectively.  The rows in each have the form
    [x,y,theta,phi,n].  For dataT, the data are copied and shifted in
    the (x,y) plane as specified by shifts.
    """
    len_T = 0
    len_B = 0
    N_seg,two = shifts.shape
    assert (two == 2)
    data = [None]*N_seg
    for i in xrange(N_seg):
        data[i] = rebin(data_in[i],binspace)
    for i in xrange(N_seg):
        len_T += len(data[i])
        for j in xrange(N_seg):
            if j != i:
                len_B += len(data[j])
    dataT = numpy.empty((len_T,5))
    dataB = numpy.empty((len_B,5))
    start_T = 0
    start_B = 0
    for i in xrange(N_seg):
        data_i = data[i]
        stop_T = start_T + len(data_i)
        dataT[start_T:stop_T,:] = data_i    # Copy [x,y,theta,phi,n]
        dataT[start_T:stop_T,0:2] += shifts[i,:]
        start_T = stop_T
        tb = Ts.sum()-Ts[i]
        c = Ts[i]/tb
        for j in xrange(N_seg):
            if j == i:
                continue
            data_j = data[j]
            stop_B = start_B + len(data_j)
            dataB[start_B:stop_B,:] = data_j
            dataB[start_B:stop_B,0:2] += shifts[i,:]
            dataB[start_B:stop_B,4] *= c
            start_B = stop_B
    return(dataT,dataB)

class Physics:
    """ A Physics object keeps track of a true density profile and
    how muons interact with that profile.
    """
    def __init__(self,    # Physics
                 Energy,  # Distribution of incident muons
                 density, # Magnitude of density deviation: scalar
                 radius,  # Radius of sphere: scalar
                 center,  # Center of sphere: array
                 Dx,
                 Dy,
                 Dz
                 ):
        """ Create a new instance
        """
        self.Energy = Energy
        self.density = density
        self.radius = radius
        self.center = center
        self.Dx = Dx           
        self.Dy = Dy
        self.Dz = Dz
        return
    def u(self,theta,phi,p):
        dsq = d2oca(self.center,theta,phi,p)
            # Distance squared of point of closest approach to center
        rsq = self.radius**2
        if dsq < rsq:
            LI = 2.0*(rsq-dsq)**.5         # Length of intersection
        else:
            LI = 0
        return (self.Dz/math.cos(theta) + self.density*LI)
    def gen_data(self,    # Physics
                 N_mu     # Number of incident muons to simulate
                 ):
        """ Simulate N_mu incident muons and save transmitted rays in self.rays.
        """
        self.t = float(N_mu) # Explosure time
        nk = 1.0
        self.rays = rays  = []
        for i in xrange(N_mu):
            theta,phi,p = self.Energy.random_ray(self)
            # theta = angle trajectory makes with vert
            # phi azimuthal angle, ie angle from x axis
            # p point where trajectory hits z=0 plane
            integral = self.u(theta,phi,p)
            energy = self.Energy.draw([theta,phi,p])
            if energy > integral:
                rays.append([nk,theta,phi,p])
        return
    def FF(self,    # Physics
           N_mu     # Number of incident muons to simulate
           ):
        """ Simulate N_mu incident muons for flat field.  Save in self
        and return (rays,energies,integrals)
        """
        self.tb = float(N_mu) # Explosure time
        nk = 1.0
        self.FF_rays = rays  = []
        for i in xrange(N_mu):
            theta,phi,p = self.Energy.random_ray(self)
            integral = self.Dz/math.cos(theta)
            energy = self.Energy.draw([theta,phi,p])
            if energy > integral:
                rays.append([nk,theta,phi,p])
        return
class Experiment:
    """ An Experiment object keeps track of detectors and models for
    estimating density.
    """
    def __init__(self,   # Experiment
                 d_0,    # Distance between centers of basis functions
                 Energy, # Distribution of incident muons
                 Dx,
                 Dy,
                 Dz,
                 Z0
                 ):
        """ Create a new instance and initialize centers of basis functions
        """
        self.rho_0 = 1.0
        self.Energy = Energy
        self.Dx = Dx
        self.Dy = Dy
        self.Dz = Dz
        self.Z0 = Z0
        self.d_x = d_0
        self.d_y = self.d_x*3**.5/2
        self.d_z = self.d_x*((2./3.)**.5)
        self.N_x = int(self.Dx/self.d_x)
        self.N_y = int(self.Dy/self.d_y)
        self.N_z = int( (self.Dz-self.Z0)/self.d_z )
        self.NA =numpy.array([self.N_y*self.N_z,self.N_z,1],numpy.int32) 
        # dot(NA,[i_x,i_y,i_z]) gives an address
        self.max = numpy.array([self.N_x,self.N_y,self.N_z],numpy.int32)
        # Might be better to make vectors: self.min, self.max, self.D,
        # self.d and self.N
        self.volume = self.d_x*self.d_y*self.d_z
        # Create functions to calculate offsets of centers
        excess_x = (self.N_x+0.5)*self.d_x - self.Dx
        excess_y = (self.N_y-0.5)*self.d_y + self.d_x - self.Dy
        excess_z = (self.N_z-1)*self.d_z + self.d_x - (self.Dz-self.Z0)
        def xoff(ix,iy,iz):
            """ Calculate shift of center of sphere
            """
            if (iy+iz)%2 == 0:
                return 0.5*self.d_x - excess_x/2
            else:
                return self.d_x - excess_x/2
        def yoff(ix,iy,iz):
            if iz%2 == 0:
                return 0.5*self.d_x- excess_y/2
            else:
                return self.d_x - excess_y/2
        def zoff(ix,iy,iz):
            return Z0 + 0.5*self.d_x - excess_z/2
        # Calculate the centers
        centers = numpy.zeros((self.N_x,self.N_y,self.N_z,3),numpy.float64)
        self.Zs = numpy.zeros(self.N_z)
        for iz in xrange(self.N_z):
            self.Zs[iz] = iz*self.d_z+zoff(0,0,iz)
            for ix in xrange(self.N_x):
                for iy in xrange(self.N_y):
                    centers[ix,iy,iz,:] = [ix*self.d_x+xoff(ix,iy,iz),
                                           iy*self.d_y+yoff(ix,iy,iz),
                                           self.Zs[iz]]
        self.centers_xyz = centers
        self.centers = centers.reshape((-1,3))
        return
    def gen_test(self, # Experiment
                 N_mu,
                 T_name,
                 B_name,
                 GP = None # If None, use ut=Radon*v+ub, else use GP.u()
                 ):
        """ Replacement for GP.gen_data() that puts the trajectories
        on a regular grid and sets the number of muons on each
        trajectory to be proportional to u (either true integral or
        Radon*v).
        """
        import tempfile, amf # Import c implementations for speed
        var2 = 2*(self.s_0*self.d_x)**2 # Basis function range: 2 \sigma^2
        norm = 1/(math.pi*var2)  # normalization for 2-d Gaussian
        T_File = open(T_name,'w')
        B_File = open(B_name,'w')
        name = tempfile.mktemp('C','','/var/lock') # /var/lock is ram
        N = int(math.sqrt(N_mu)) # Use N angles and N positions
        d = math.sqrt(self.Dx*self.Dy/N)
        Xs = numpy.arange(0,self.Dx,d)
        Ys = numpy.arange(0,self.Dy,d)
        var2 = 2*(self.s_0*self.d_x)**2 # Basis function range: 2 \sigma^2
        dtheta = math.sqrt(2*math.pi*(1-math.cos(Dtheta))/N)
        # Azimuthal angle step size is sqrt(solid angle/N)
        N_miss = 0
        N_hit = 0
        i_u = 1 # Dummy value to make uline happy
        Maxlen = -1
        v = self.V
        for theta in numpy.arange(0,Dtheta,dtheta):
            # Aim for squares on sphere about d_theta on a side
            if math.sin(theta) < dtheta/(2*math.pi):
                dphi = 2*math.pi
            else:
                dphi = dtheta/math.sin(theta)
            for phi in numpy.arange(0,2*math.pi,dphi):
                D_xy_u = self.xy_range(theta,phi,var2)
                for x in Xs:
                    for y in Ys:
                        p_0 = numpy.array([x,y,0],numpy.float64)
                        L = amf.uline(self.Zs, D_xy_u, theta, phi, self.d_x,
                             self.d_y, p_0, self.centers_xyz, norm, var2,
                             self.N_x, self.N_y, self.N_z, i_u,name)
                        if L == 0:
                            N_miss += 1
                            continue
                        N_hit += 1
                        if L > Maxlen:
                            Maxlen = L
                            I = numpy.empty(Maxlen,numpy.int32)
                            R = numpy.empty(Maxlen,numpy.float32)
                        amf.readR(L,I,R,name) # Read what uline() wrote
                        os.remove(name)
                        ub = (self.Dz*self.rho_0)/math.cos(theta)
                        nb = self.Energy.F(ub)
                        amf.writeRadon(L,I,R,theta,nb,B_File)
                        if GP == None:
                            ut = amf.dot(v,L,I,R) + ub # Radon*v + ub
                        else:
                            ut = GP.u(theta,phi,p_0)   # True integral
                        nt = self.Energy.F(ut)
                        amf.writeRadon(L,I,R,theta,nt,T_File)
        File = open(T_name+'head','w')
        print >>File,'Nmax: %d\nN_mu: %d\nt: %d'%(Maxlen,N_hit,N_miss+N_hit)
        File = open(B_name+'head','w')
        print >>File,'Nmax: %d\nN_mu: %d\nt: %d'%(Maxlen,N_hit,N_miss+N_hit)
        self.NT_max = Maxlen
        self.T_muons = N_hit
        self.t = float(N_miss+N_hit)
        self.NT_max = Maxlen
        self.B_muons = N_hit
        self.tb = float(N_miss+N_hit)
        self.NB_max = Maxlen
        return
    def index_2_triple(self,j):
        """ Translate an index of centers to an array of integer
        coordinates.
        """
        j_z = j%self.N_z
        k = j/self.N_z
        j_y = k%self.N_y
        j_x = k/self.N_y
        return numpy.array([j_x,j_y,j_z],numpy.int32)
    def triple_2_index(self,ixyz):
        """ Translate an array of integer coordinates to an index of
        centers.
        """
        for i,N in zip(ixyz,[self.N_x,self.N_y,self.N_z]):
            if i< 0 or i >= N:
                return -1
        return numpy.dot(ixyz,self.NA)
    def gen_V(self,    # Experiment
              density, # Magnitude of density deviation: scalar
              radius,  # Radius of sphere: scalar
              center   # Center of sphere: array
              ):
        """ Calculate coefficients V to represent a specified
        spherical deviation from nominal density.  Only used for debugging
        """
        N,three = self.centers.shape
        assert three == 3
        r2 = radius*radius
        self.V = V = numpy.zeros((N,))
        for i in xrange(N):
            ds = center-self.centers[i,:]
            if numpy.dot(ds,ds) < r2:
                V[i] = density*self.volume
        return V

    def TV_prior(self,   # Experiment
                 s_0,    # Scalar: Size of Gaussians, sigma = d_x*s_0
                 alpha,  # Weight of TV regularization
                 beta,   # Smoothing of absolute value
                 RS_x,   # Scalars: RS_* gives the output dimensions of
                 RS_y,   # the resample_matrix
                 RS_z
                 ):
        """Calculate the D matrix for TV regularization:

           T(v) = alpha sqrt(v^T D v*len(v) + beta)

        where alpha and beta are scalars and v_i is a coefficient of a
        Gaussian function that is a component of the density profile.

        Also calculate the resample matrix RSM that maps from vectors
        of coefficients V to a rectangular grid of density values.
        """
        n_centers,three = self.centers.shape
        assert three == 3
        self.s_0 = s_0
        self.alpha = alpha
        self.beta = beta
        n_rho = RS_x*RS_y*RS_z
        name='tmp/priorTV_d%.4g_s%.4g_RSx%d_RSy%d_RSz%d'%(
            self.d_x,s_0,RS_x,RS_y,RS_z,)
        try:
            d = scipy.io.loadmat(name)
            self.D = d['D']
            self.mu_0 = d['V_mean']
            self.RSM = d['RSM']
            return
        except:
            V_mean = numpy.mat(numpy.zeros((1,n_centers),numpy.float64))
            
            ### Build resample matrix to map from V to density field in space
            Range = self.d_x*s_0  # Range of basis functions
            rr2 = Range*Range*2
            norm = 1.0/((2*math.pi)**3 * Range**6)**.5 # Make total mass=1
            RSM = scipy.sparse.lil_matrix((n_rho, n_centers))
            RS_dx = float(self.Dx)/RS_x
            RS_dy = float(self.Dy)/RS_y
            RS_dz = float(self.Dz)/RS_z
            for i_c in xrange(n_centers):
                c = self.centers[i_c,:]
                cx,cy,cz = c
                x_low  = max(0,int((cx-4*Range)/RS_dx))
                x_high = min(RS_x,int((cx+4*Range)/RS_dx))
                y_low  = max(0,int((cy-4*Range)/RS_dy))
                y_high = min(RS_y,int((cy+4*Range)/RS_dy))
                z_low  = max(0,int((cz-4*Range)/RS_dz))
                z_high = min(RS_z,int((cz+4*Range)/RS_dz))
                for ix in xrange(x_low,x_high):
                    x = ix*RS_dx
                    for iy in xrange(y_low,y_high):
                        y = iy*RS_dy
                        for iz in xrange(z_low,z_high):
                            z = iz*RS_dz
                            i_r = iz + RS_z*(iy+RS_y*ix)
                            d = c - numpy.array([x,y,z])
                            exp = -numpy.dot(d,d)/rr2
                            if exp > exp_cut:
                                RSM[i_r,i_c] = norm*math.exp(exp)
                            # i_r is resample index and i_c is center index
            RSM = RSM.tocsr()
            assert RSM.shape == (n_rho,n_centers)
            # Build D
            D = scipy.sparse.lil_matrix((n_centers,n_centers))
            # par_2_list collects separation vectors for the other
            # centers in tetrahedrons above and below a given center.
            # The offsets depend on the parity of the y and z
            # coordinates of the center considered as follows:
            #
            #   Parity       Offsets for tetrahedron above
            #   Y   Z        X   Y   Z
            #   0   0        0   0   +1
            #                -1  0   +1
            #                0   -1  +1
            #   -----------------------
            #   0   1        0   0   +1
            #                +1  0   +1
            #                0   +1  +1
            #   -----------------------
            #   1   0        0   0   +1
            #                +1  0   +1
            #                0   -1  +1
            #   -----------------------
            #   1   1        0   0   +1
            #                -1  0   +1
            #                0   +1  +1
            #   -----------------------
            # associated vals that are above the threshold specified
            # by exp_cut
            par_2_list ={
                (0,0):numpy.array([[0,0],[-1,0],[0,-1]],numpy.int32),
                (0,1):numpy.array([[0,0],[1,0],[0,1]],numpy.int32),
                (1,0):numpy.array([[0,0],[1,0],[0,-1]],numpy.int32),
                (1,1):numpy.array([[0,0],[-1,0],[0,1]],numpy.int32)}
            # This loop augments D with
            #
            #   3  -1  -1  -1
            #  -1   3  -1  -1
            #  -1  -1   3  -1
            #  -1  -1  -1   3
            #
            # placed at coordinates specified by par_2_list
            delta = numpy.empty(3)
            for i in xrange(n_centers):
                ixyz = self.index_2_triple(i)
                par = (ixyz[1]%2,ixyz[2]%2)
                for dz in (-1,1): # Do tetrahedron above and below i
                    List = [i]
                    for dx_dy in par_2_list[par]:
                        delta[0:2] = dx_dy
                        delta[2] = dz
                        jxyz = ixyz+delta
                        j = self.triple_2_index(jxyz)
                        if j >= 0: # triple_2_index returns -1 if out of bounds
                            List.append(j)
                    #If len(List) < 4, this regularization pulls towards 0
                    for j in List:
                        for k in List:
                            if j == k:
                                D[j,k] += 3
                            else:
                                D[j,k] -= 1
            D = D.tocsr()
            scipy.io.savemat(name,{'D':D,'V_mean':V_mean,'RSM':RSM})
            self.D = D
            self.mu_0 = V_mean
            self.RSM = RSM
            return
    def xy_range(self,  # Experiment
                 theta, # ray angle
                 phi,   # ray angle
                 var2   # 2 times variance of basis functions
                 ):
        """ Calculate and return a list of scalar integer offsets that
        are in the closure of the exp_cut neighborhood of an incident
        ray.
        """
        boundary = -(1+math.sqrt(-exp_cut*var2))**2
        p = numpy.zeros((3,)) # Point ray crosses z=0
        def min_perim(d):
            """ Find the minimum distance squared along a perimeter of
            radius d.
            """
            c = numpy.array([self.d_x*d,self.d_y*d,0],numpy.float64) # Test pt
            smallest = d2oca(c,theta,phi,p)
            cs = [c]
            for di in xrange(-d,d):
                cs.append(numpy.array([self.d_x*di,self.d_y*d,0],numpy.float64))
                cs.append(numpy.array([self.d_x*d,self.d_y*di,0],numpy.float64))
            for c in cs:
                rr = d2oca(c,theta,phi,p)
                if rr < smallest:
                    smallest = rr
            return rr
        for d in xrange(1,100):
            rr = min_perim(d)
            if -rr/var2 < boundary:
                break
        #assert d < 25,'d=%d'%d # FixMe why?
        if d > 25:
            print 'd=%d'%d
        D_xy = []
        dxyf = numpy.array([self.d_x,self.d_y,0],numpy.float64)
        for dxi in xrange(-d,d+1):
            for dyi in xrange(-d,d+1):
                dxyi = numpy.array([dxi,dyi,0],numpy.int32)
                c = dxyi*dxyf
                if -d2oca(c,theta,phi,p)/var2 > boundary:
                    D_xy.append(dxyi)
        return numpy.array(D_xy,dtype=numpy.int32,order='C')
    def fminLM(self,            # Experiment
               Radon_T_name,
               Radon_B_name,
               callback=None    #
             ):
        """ List mode variant of fmin.  Reference:
        http://docs.scipy.org/doc/scipy/reference/generated\
        /scipy.optimize.fmin_ncg.html Version 173 used v rather than
        Sigma_0^{-1} v."""
        import amf
        t = float(self.t)
        tb = float(self.tb)
        N_target = self.T_muons
        N_back = self.B_muons
        Max_N = max(self.NT_max,self.NB_max)
        I = numpy.zeros(Max_N,numpy.int32)
        R = numpy.zeros(Max_N,numpy.float32)
        def reg(v):
            Dv = self.D*v
            vDv = numpy.dot(v,Dv)
            T = math.sqrt(vDv*len(v)+self.beta**2)
            if T < 1.0e50:
                print 'T_norm=%9.3g '%(T/self.beta),
                return T,Dv
            vmax = v.max()
            print '\nv.max()=%f\n'%v.max()
            assert vmax < 1.0e50
            return 1.0e50,Dv
            
        def L(v):
            """ Calculate and return the scalar valued objective
            function."""
            t0 = time.time()
            rmsv = math.sqrt(numpy.dot(v,v)/len(v))
            if rmsv > 1e4:
                return rmsv*1e20 # Kludge to prevent overflow in F(k_k)
            T,Dv = reg(v)
            L_ = -self.alpha*T
            File = open(Radon_T_name)
            for k in xrange(N_target):
                n,theta,N_k = amf.readRadon(I,R,File)
                u_k = amf.dot(v,n,I,R);
                u_k += (self.Dz*self.rho_0)/math.cos(theta)
                L_ += N_k*math.log(self.Energy.F(u_k)) # First term
            if not(L_ < 1.0e20 and L_ > -1.0e20):
                print 'L_=%f, u_k=%f, v.max()=%f, R.max()=%f, n=%d\n'%(
                    L_, u_k, v.max(), R.max(), n)
                raise RuntimeError

            File = open(Radon_B_name)
            for j in xrange(N_back):
                n,theta,N_j = amf.readRadon(I,R,File)
                u_j = amf.dot(v,n,I,R);
                u_0 = (self.Dz*self.rho_0)/math.cos(theta)
                u_j += u_0
                F = self.Energy.F(u_j)
                F_0 = self.Energy.F(u_0)
                change = (t/tb)*(N_j*F/F_0) # Next term of return value
                L_ -= change
                if not(L_ < 1.0e20 and L_ > -1.0e20):
                    print '\n\nL_=%f, change=%f, u_j=%f, u_0=%f, F=%f'%(
                        L_, change, u_j, u_0, F)
                    print 'F_0=%f, v.max()=%f, v.min()=%f, R.max()=%f'%(
                        F_0, v.max(), v.min(), R.max())
                    print 'n=%d, Max_N=%d, t=%f, tb=%f, N_j=%f\n'%(
                        n, Max_N, t, tb, N_j)
                    raise RuntimeError
            print '  L in %6.1f ms value: %7.2f'%((time.time()-t0)*1000,
                -L_)
            return (-L_)

        def dL(v):
            """ Calculate and return the gradient of L.  I believe
            that Sigma_0 is a sparse matrix and that rv is an array
            not a matrix."""
            t0 = time.time()
            T,Dv = reg(v)
            dL_ = -self.alpha*Dv/T
            File = open(Radon_T_name)
            for k in xrange(N_target):
                n,theta,N_k = amf.readRadon(I,R,File)
                u_k = amf.dot(v,n,I,R);
                u_k += (self.Dz*self.rho_0)/math.cos(theta)
                F = self.Energy.F(u_k)
                d_F = self.Energy.d_F(u_k)
                c = N_k*d_F/F
                amf.accumulate(dL_,c,n,I,R)

            File = open(Radon_B_name)
            for j in xrange(N_back):
                n,theta,N_j = amf.readRadon(I,R,File)
                u_j = amf.dot(v,n,I,R);
                u_0 = (self.Dz*self.rho_0)/math.cos(theta)
                u_j += u_0
                F_0 = self.Energy.F(u_0)
                d_F = self.Energy.d_F(u_j)
                c = -(t/tb)*(N_j*d_F/F_0)
                amf.accumulate(dL_,c,n,I,R)
            print ' dL in %6.1f ms'%(1000*(time.time()-t0))
            return (-dL_)

        def ddL(v,x):
            """ Calculate and return the Hessian of L at v applied to
            x.  I believe that the return value is an array not a
            matrix."""
            t0 = time.time()
            T,Dv = reg(v)
            Dx = self.D*x
            ddLx_ = -self.alpha*(T*T*Dx-Dv*numpy.dot(v,Dx))/(T*T*T)
            File = open(Radon_T_name)
            for k in xrange(N_target):
                n,theta,N_k = amf.readRadon(I,R,File)
                Rx = amf.dot(x,n,I,R);
                u_k = amf.dot(v,n,I,R);
                u_k += (self.Dz*self.rho_0)/math.cos(theta)
                F = self.Energy.F(u_k)
                d_F = self.Energy.d_F(u_k)
                dd_F = self.Energy.dd_F(u_k)
                c = N_k*(F*dd_F-d_F*d_F)/(F*F)*Rx
                amf.accumulate(ddLx_,c,n,I,R)

            File = open(Radon_B_name)
            for j in xrange(N_back):
                n,theta,N_j = amf.readRadon(I,R,File)
                Rx = amf.dot(x,n,I,R);
                u_j = amf.dot(v,n,I,R);
                u_0 = (self.Dz*self.rho_0)/math.cos(theta)
                u_j += u_0
                c = -(t/tb)*N_j*self.Energy.dd_F(u_j)/self.Energy.F(u_0)*Rx
                amf.accumulate(ddLx_,c,n,I,R)

            print 'ddL in %6.1f ms'%(1000*(time.time()-t0))
            return (-ddLx_)
        import scipy.optimize as OPT
        N_v,three = self.centers.shape
        assert three == 3
        v0 = numpy.zeros((N_v,),numpy.float64)
        #v0 = self.V # Check to see if estimate moves from true
        t0 = time.time()
        #v = OPT.fmin_cg(L, v0, dL, disp=True, callback=callback)
        v = OPT.fmin_ncg(L, v0, dL, fhess_p=ddL, disp=True, callback=callback,
                         avextol=0.1)
        # other options: avextol=1e-8, maxiter=30,
        # disp=True, callback=callback)
        print 'fmin, %f seconds'%(time.time()-t0)
        return v # End of fminLM()

    def writeRadonT(self,   # Experiment
                    GP,     # Physics
                    name    # Name of file to write
                    ):
        """ Open file name and for each detected muon in GP.rays
        (those that went through the target density) write a sparse
        vector that implements the Radon transform for its trajectory."""
        self.NT_max,self.T_muons = self.writeRadon(GP.rays,GP.t,name)
        self.t = GP.t
        return
    def writeRadonB(self,   # Experiment
                    GP,     # Physics
                    name    # Name of file to write
                    ):
        """ Open file name and for each detected muon in GP.FF_rays
        (those that went through the background density) write a
        sparse vector that implements the Radon transform for its
        trajectory."""
        self.NB_max,self.B_muons = self.writeRadon(GP.FF_rays,GP.tb,name)
        self.tb = GP.tb
        return
    def writeRadon(self,      # Experiment
                   rays,      # rays[i]=[theta,phi,p]. p=numpy.array([x,y,0])
                   t,         # Exposure time
                   F_name,    # Name of file to write
                   data=None  # Kludge for alternate format
                   ):
        import tempfile, amf # Import c implementations for speed
        var2 = 2*(self.s_0*self.d_x)**2 # Basis function range: 2 \sigma^2
        norm = 1/(math.pi*var2)  # normalization for 2-d Gaussian
        File = open(F_name,'w')
        name = tempfile.mktemp('C','','/var/lock') # /var/lock is ram
        N_mu = 0  # Count number of Radon vectors written
        Nmax = -1
        i_u = 1   # Dummy value to make uline happy
        if data != None:
            p_0 = numpy.zeros(3)
            rays = data
        for r in rays:
            if data == None:
                nk,theta,phi,p_0 = r
            else:
                x,y,theta,phi,nk = r
                p_0[0:2] = [x,y]
            D_xy_u = self.xy_range(theta,phi,var2)
            N = amf.uline(self.Zs, D_xy_u, theta, phi, self.d_x, self.d_y, p_0,
                          self.centers_xyz, norm, var2, self.N_x, self.N_y,
                          self.N_z, i_u,name)
            if N == 0:
                continue
            if N > Nmax:
                Nmax = N
                I = numpy.empty(N,numpy.int32)
                R = numpy.empty(N,numpy.float32)
            amf.readR(N,I,R,name) # Read what uline() wrote
            os.remove(name)
            amf.writeRadon(N,I,R,theta,nk,File)
            N_mu += 1
        File = open(F_name+'head','w')
        print >>File,'Nmax: %d\nN_mu: %d\nt: %d'%(Nmax,N_mu,t)
        return (Nmax,N_mu)
class Power_E:
    """ Power law energy distribution.  This distribution is
    independent of angle and position.

    if E<scale, Pdf(E)=0
    else, Pdf(E)=scale*E**-p
    """
    def __init__(self,p,scale,theta_min=0,theta_max=Dtheta):
        self.p = p
        self.scale = scale
        self.gen_a = 1-math.cos(theta_min) # For generating random theta
        self.gen_b = 1-math.cos(theta_max) # For generating random theta
        return
    def random_ray(self,  # Power_E
                   EG     # Experiment geometry
                   ):
        t = random.uniform(self.gen_a,self.gen_b)
        theta = math.acos(1-t)  # CDF(theta) = (1-cos(theta))/2
        phi = 2*math.pi*random.random()
        x = EG.Dx*random.random()
        y = EG.Dy*random.random()
        p = numpy.array([x,y,0],numpy.float64)
        return (theta,phi,p)
    def draw(self,  # Power_E
             ray
             ):
        # Random draw from power law distribution. y and x are dummy variables
        y = numpy.random.rand()
        x = y**(1/(1-self.p))
        return self.scale*x
    def F(self, # Power_E
          u     # Scalar energy
          ):
        val = (u/self.scale)**(self.p-1) 
        return min(val,1.0)
    def d_F(self, # Power_E
            u     # Scalar energy
            ):
        s = self.scale
        p = self.p
        us = u*s
        if us > 1.0:
            return us**(-p)*(1-p)/s
        return 0.0
    def dd_F(self, # Power_E
            u     # Scalar energy
            ):
        s = self.scale
        p = self.p
        us = u*s
        if us > 1.0:
            return us**(-(p-1)*p)*(1-p)*(-p)/s**2
        return 0.0
class Exp_E(Power_E):
    """ Exponential energy distribution that is independent of trajectory
    """
    def __init__(self,
                 Lam,
                 theta_min=0,
                 theta_max=Dtheta):  # FixMe: Dtheta is hard coded here
        self.Lam = Lam
        self.gen_a = 1-math.cos(theta_min) # For generating random theta
        self.gen_b = 1-math.cos(theta_max) # For generating random theta
        return
    def draw(self,  # Exp_E
             ray
             ):
        # Random draw from exponential distribution
        return random.expovariate(self.Lam)
    def F(self, # Exp_E
          u     # Scalar energy
          ):
        lam = self.Lam
        val = math.exp(-lam*u)
        return val
    def d_F(self, # Exp_E
            u     # Scalar energy
            ):
        lam = self.Lam
        val = math.exp(-lam*u)
        return -lam*val
    def dd_F(self, # Exp_E
            u     # Scalar energy
            ):
        lam = self.Lam
        val = math.exp(-lam*u)
        return lam*lam*val

class test_E(Power_E):
    """ Fake energy model that maps u to rate with identity
    """
    def __init__(self):
        return
    def F(self, # test_E
          u     # Scalar energy
          ):
        return u
    def d_F(self, # test_E
            u     # Scalar energy
            ):
        return 1.0
    def dd_F(self, # test_E
            u     # Scalar energy
            ):
        return 0.0

if __name__ == '__main__': # Test code
    d_0 = 15.0             # Spacing of Gaussian basis functions
    s_0 = 0.5              # Scale of Gaussian basis functions
    density = -.7          # Density deviation in sphere
    radius = 15.1          # Radius of sphere
    location = numpy.array([70,20,50],numpy.float64)
    Energy = Exp_E(0.05)   # Energy distribution model
    N_mu = 10000           # Number of muons
    N_ux = 9               # Number of x bins in U space
    N_uy = 3
    N_theta = 6            # Number of theta bins in U space
    dev_0 = 0.3            # Variance of prior
    corr_0 = d_0           # Range of density correlations
    RS_x = 12              # Number of x-points for resampling
    RS_y = 6
    RS_z = 6

    T_start = time.time()
    print 'This test takes about ? seconds'
    GE = Experiment(d_0,Energy)
    GP = Physics(Energy,density,radius,location)
    V = GE.gen_V(density,radius,location)
    GP.gen_data(N_mu)
    GP.FF(N_mu)
    GE.fminLM('tmp/tempT','tmp/tempB')
    print 'Test finished in %f seconds'%(time.time()-T_start)

#Local Variables:
#mode:python
#End:
