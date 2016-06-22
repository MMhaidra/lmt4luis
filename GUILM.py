"""
GUILM.py GUI for list mode

Uses:

'python GUILM'                 Start GUI to operate on simulated data
'python GUILM simulation'      Start GUI to operate on simulated data
'python GUILM data'            Start GUI to operate on experimental data
'python GUILM simulation file' Apply commands in file to simulated data
'python GUILM data file'       Apply commands in file to experimental data

"""
import Tkinter, numpy, time, os, tempfile, sys, math
import geometryLM as geom
print '\n'

if len(sys.argv) < 2:
    data_mode = 'simulation'
else:
    if sys.argv[1] == 'simulation':
        data_mode = 'simulation'
    else:
        assert sys.argv[1] == 'data'
        data_mode = 'data'
if data_mode == 'simulation':
    GVD = {'Dx':145,
           'Dy':46,
           'Dz':198,
           'Z0':60,
           'RS_res':2.0,
           'Basis_d':5.0,
           'Basis_s':0.5,
           'alpha':1.0e-5,
           'beta':1.0e4,
           'Sp_d':1.0,
           'Sp_r':20,
           'Sp_x':50,
           'Sp_y':23.0,
           'Sp_z':120,
           'V_y':0.5,
           'Lim_val':1.0,
           'N_mu':1000,
           'N_FF':1000}
    E_index = 2  # Selects model F(x) =n x
if data_mode == 'data':
    X_offset = 10.0
    Y_offset = 31.75
    GVD = {'Dx':288 + X_offset,
           'Dy':48 + 2*Y_offset,
           'Dz':200,
           'Z0':10,
           'RS_res':3.2,
           'Basis_d':5.0,
           'Basis_s':0.5,
           'alpha':1,
           'beta':1.0,
           'V_y':0.5,
           'Lim_val':0.1,
           'XYbin':1.0}
    shifts = numpy.array([[ 0,    0],
                          [45.72, 0],
                          [91.44, 0],
                          [137.16,0]]) # Data from J. Andrew Green
    shifts += [X_offset,Y_offset]
    times = numpy.array([15012.,
                         11187.,
                         13133.,
                         12453.],numpy.float64)
    E_index = 1  # Selects model geom.Exp_E()
GVD['RS_nx'] = int(GVD['Dx']/GVD['RS_res'])
GVD['RS_ny'] = int(GVD['Dy']/GVD['RS_res'])
GVD['RS_nz'] = int(GVD['Dz']/GVD['RS_res'])

##################Incident muon distribution models###################
E_p=2.7
E_scale=5.0
Lam = 0.15      # characterizes exponential particle energy distribution
Energies = [geom.Power_E(E_p,E_scale,theta_min=0,theta_max=1.0),
            geom.Exp_E(Lam=1.0/GVD['Dz']), geom.test_E()
            ]
Energy = Energies[E_index]

############### Flags for preventing use of obsolete data ################
Prior_OK = False
Muons_OK = False
################ A few globals to reduce complaints ##################
GE = None        # Experiment instance
V_0 = None       # Data for plotting simulated density
V_hat = None     # Data for plotting estimated density
Radon_T_name = None
Radon_B_name = None
message = None

############ definition and some instances of alert/help class #########
#From http://mail.python.org/pipermail/python-list/2002-December/174253.html
class HelpDialog(Tkinter.Toplevel):
    def __init__(self,title,text):
        Tkinter.Toplevel.__init__(self)
        self.title(title)
        self.__result = None
        Tkinter.Button(self, command=self.destroy, text='Close ').pack(
          side=Tkinter.TOP)
        Tkinter.Label(self,text=text, justify=Tkinter.LEFT).pack(
          side=Tkinter.LEFT)
def Not_Implemented():
    HelpDialog('Not implemented',
"""I haven't implemented support for this yet.
""")
def PriorHelp():
    HelpDialog('Prior out of date',
"""Sliders that set values for the prior may have changed.
Press the "New Prior" button after closing this alert.
""")
def MuonHelp():
    HelpDialog('Muons out of date',
"""Sliders that affect muon counts may have changed.  Press either the
"Muons" button or one of the "Test" buttons after closing this alert.
""")
def ButtonHelp():
    HelpDialog('Help for the top row of buttons',
"""Quit          Close the application
New Prior  Calculate new prior distribution.  Required if Basis or TV Reg
                 changed.
Muons       Read or make/simulate Radon vectors and single counts stochastically
Test v        Read or make Radon vectors and counts based on coefficients v
Test rho    Read or make Radon vectors and counts based on the "true" density
Estimate   Estimate v_hat based on Radon vectors and counts

Buttons determine energy model that maps u to R, ie, R = F(u):

        Power Law       F(u) = u^{1-p}

        Exponential     F(u) = e^{-\lambda u}

        F(x) = x           F(u) = u    (Only use this model with Test v or
                                                 Test rho)
        
Note that since information calculated for priors and Radon vectors is
cached in files in the directory ./tmp, if you change the code in a
manner that changes such information, you should purge the
corresponding files.""")
def BasisHelp():
    HelpDialog('Help for basis sliders',
"""These sliders set the Gaussian basis functions used
to fit the density perturbation:
                   
d_0  The minimum distance between centers
s_0  The range of each function as fraction of d_0""")
def TVHelp():
    HelpDialog('Help for TV Reg sliders',
"""These sliders adjust the Total Variation Regularization.  They are
both log based 10.

alpha  Scalar coefficient of regularization term
beta   Rounds corner of absolute value function
""")
def NMUHelp():
    HelpDialog('Help for N muon',
"""These sliders set the number of incident muons to
simulate: Left for target and right for flat field.
""")
def SPHelp():
    HelpDialog('Help for Spherical Perturbation sliders',
"""These sliders set:

D   Density of sphere
r   Radius
x   X location of center
y   Y location
z   Z location 
""")
def RebinHelp():
    HelpDialog('Help for V slice sliders',
"""This slider sets the bin size in the xy plane for experimental data
in centimeters.  Since the data have .5 cm resolution to begin with,
they will only get rebinned if this slider is at least 1.0.

The data gets rebinned before the Radon vectors are written when you
press the "Radon Vectors" button.
""")
def VsliceHelp():
    HelpDialog('Help for V slice sliders',
"""These sliders control the plots of the densities in cross sections.
'y' selects the plane whose densities are displayed in terms of the
fraction of the total allowed Y-range, and 'Log range' sets the zrange
for the plot.
""")

############# Initialize pipes for plotting ##########################
plot_data = {}
# Find size of display and calculate plot size and spacing
for line in os.popen('xwininfo -root','r').readlines():
    field = line.split()
    if len(field) < 2:
        continue
    if field[0] == 'Height:':
        X11_height = int(field[1])
    if field[0] == 'Width:':
        X11_width = int(field[1])
P_h = X11_height/3                   # Placement spacing
P_w = min(X11_width/4,int(P_h*1.25)) # Placement spacing
P_H = P_h - 50   # Plot size
P_W = P_w - 15   # Plot size

splot_preface = """
set hidden
set parametric
set contour base
set style data lines
set terminal x11 %d title "%s"
set xlabel 'x'
set ylabel "%s"
set title "%s"
set grid
set cntrparam levels 10
"""
#set cntrparam levels incremental -2, .2, 2
ylabels={'V_0':'z',
         'V_hat':'z'}
titles={'V_0':'$V_0$ = density deviation via coefficients V',
        'V_hat':'$\\\hat V =$mean a-posteriori from ($\\\Sigma_0$,count)'}
surface_name={'V_0':'$V_0$',
         'V_hat':'$\\\hat V$'}
GPD = {} # dict of gnuplot pipes
if data_mode == 'simulation':
    keylist = ('V_0','V_hat')
else:
    keylist = ('V_hat',)
for key in keylist:
    n = len(GPD)
    x = (n%4)*P_w
    y = (n/4)*P_h
    P = os.popen ('gnuplot -geometry %dx%d+%d+%d'%(P_W,P_H,x,y),'w')
    P.write(splot_preface%(n,key,ylabels[key],titles[key]))
    GPD[key] = P
def splot_v(key,    # The name of the pipe to write to and title of surface
            GE,     # An Experiment instance
            v=None, # A one dimensional vector of coefficients            
            tname=None,
            count=None
          ):
    if v != None:
        plot_data[key] = v
    else:
        if plot_data.has_key(key):
            v = plot_data[key]     # For plotting new view when y slider changes
        else:
            return
    if tname == None:
        tname = tempfile.mktemp('',key+'.')
    iy = min(int(GVD['V_y']*GVD['RS_ny']),GVD['RS_ny']-1)
    rho = numpy.asarray(GE.RSM*v)
    RS_nx,RS_ny,RS_nz = (GVD['RS_nx'],GVD['RS_ny'],GVD['RS_nz'])
    rho = rho.reshape((RS_nx,RS_ny,RS_nz))
    rho = rho[:,iy,:]
    rho = rho.reshape((RS_nx,RS_nz))
    t_file = open(tname,'w')
    d_x = float(GE.Dx)/RS_nx
    d_z = float(GE.Dz)/RS_nz
    for x in xrange(RS_nx):
        for z in xrange(RS_nz):
            fx = x*d_x
            fz = z*d_z
            t_file.write('%f %f %f\n' % (fx, fz, float(rho[x,z])))
        t_file.write('\n')
    t_file.close()
    pipe = GPD[key]
    if count != None:
        pipe.write('set title "count=%d"\n'%count)
    else:
        pipe.write('set title "%s"\n'%titles[key])
    pipe.write('splot "%s" title "%s"\n'%(tname,surface_name[key]))
    pipe.flush()
    return
save_gnuplot = """
set terminal epslatex color dashed "normal" 9
set output "%s.tex"
replot
set terminal x11 title "%s"
"""
def Save_Plots(save_names=None):
    import shutil, tkFileDialog as TKFD
    for key in GPD.keys():
        pipe = GPD[key]
        if save_names == None:
            save_tex = TKFD.asksaveasfilename(defaultextension='.tex')
            parts = save_tex.split('.')
            assert len(parts) == 2,'len(parts)=%d'%len(parts)
            base = os.path.basename(parts[0])
        else:
            base = save_names[key]
        pipe.write(save_gnuplot%(base,key))
        pipe.flush()
        File = open(base+'.log','w')
        for key in GVD.keys():
            if GVD[key] != None:
                if type(GVD[key]) == type(1): #isint?
                    print >>File, ' %7s = %7d'%(key,GVD[key])
                else:
                    print >>File, ' %7s = %7.2g'%(key,GVD[key])
    return
############# Service routines for buttons #########################
def NewPrior():
    global Prior_OK, GE, GVD
    t0 = time.time()
    GE = geom.Experiment(GVD['Basis_d'],Energy,GVD['Dx'],GVD['Dy'],GVD['Dz'],
                         GVD['Z0'])
    GE.TV_prior(GVD['Basis_s'],GVD['alpha'],GVD['beta'],GVD['RS_nx'],
                GVD['RS_ny'],GVD['RS_nz'])
    Prior_OK = True
    print "New prior in %6.2f seconds"%(
      time.time()-t0),"RSM.shape=",GE.RSM.shape
    return
def DataPrior():
    global Prior_OK, GE, GVD
    t0 = time.time()
    GE = geom.Experiment(GVD['Basis_d'],Energy,GVD['Dx'],GVD['Dy'],
                         GVD['Dz'],GVD['Z0'])
    GE.TV_prior(GVD['Basis_s'],GVD['alpha'],GVD['beta'],GVD['RS_nx'],
                GVD['RS_ny'],GVD['RS_nz'])
    Prior_OK = True
    print "New prior in %6.2f seconds"%(
      time.time()-t0),"RSM.shape=",GE.RSM.shape
    return
def DataMuons():
    global Muons_OK, GE, Radon_T_name, Radon_B_name
    import scipy.io
    if not Prior_OK:
        return PriorHelp()
    Radon_T_name='tmp/Tdata_dx%.1f_s0%.2f_bin%.1f'%(GE.d_x,GE.s_0,GVD['XYbin'])
    Radon_B_name='tmp/Bdata_dx%.1f_s0%.2f_bin%.1f'%(GE.d_x,GE.s_0,GVD['XYbin'])
    def check_files(name):
        File = open(name,'r')
        File = open(name+'head','r')
        # Format: 'Nmax: %d\nN_mu: %d\nt: %d'.  Second item in each line.
        L = [int(line.split()[1]) for line in File.readlines()]
        return (L[0],L[1],float(L[2]))
    try:
        GE.NT_max,GE.T_muons,GE.t = check_files(Radon_T_name)
        GE.NB_max,GE.B_muons,GE.tb = check_files(Radon_B_name)
    except:
        try:
            data_dict = scipy.io.loadmat('tmp/AG_rotated')
        except:
            print 'Failed to read tmp/AG_rotated. Creating it'
            Gamma = [0.0,    0.0,   0.0,   0.00911,  0.0]
            # Gamma describes rotation that J. Andrew Green reported for data
            name_ = '/home/andy/fuse/green/luis_4Dhist_OverheadBricks_pos%d.txt'
            data_dict = {}
            for i in xrange(4):
                Name = name_%(i+1)
                A = geom.ReadGreen(Name)
                key = 'AG%d'%(i+1)
                data_dict[key] = geom.Reframe(Gamma,A) # [x,y,theta,phi,n]
            scipy.io.savemat('tmp/AG_rotated',data_dict)
        # Now data_dict has rotated data for 4 segments
        data = 4*[None]
        for i in xrange(4):
            key = 'AG%d'%(i+1)
            data[i] = data_dict[key] # Used dict to keep loadmat happy
        dataT,dataB = geom.join_data(data,shifts,times,binspace=GVD['XYbin'])
        rays,t = (None,1.0) # To placate writeRadon
        t0 = time.time()
        GE.NT_max,GE.T_muons = GE.writeRadon(rays,t,Radon_T_name,data=dataT)
        GE.NB_max,GE.B_muons = GE.writeRadon(rays,t,Radon_B_name,data=dataB)
        GE.t = GE.tb = 1.0  # Duration information went into the n_k fields
        print 'Generated Radon vectors from data in %.1f seconds'%(
            time.time()-t0)
    Muons_OK = True
    return
def NewMuons():
    global Muons_OK, GE, V_0, Radon_T_name, Radon_B_name, message
    if not Prior_OK:
        return PriorHelp()
    V_0 = GE.gen_V(GVD['Sp_d'],GVD['Sp_r'],
                   numpy.array([GVD['Sp_x'],GVD['Sp_y'],GVD['Sp_z']]))
    splot_v('V_0',GE,v=V_0)
    Radon_T_name='tmp/T_dx%.1f_s0%.2f_N%d_E%d_D%.2f_r%.1f_x%.1f_y%.1f_z%.1f'%(
        GE.d_x,GE.s_0,GVD['N_mu'],E_index,GVD['Sp_d'],GVD['Sp_r'],GVD['Sp_x'],
        GVD['Sp_y'],GVD['Sp_z'])
    Radon_B_name='tmp/B_dx%.1f_s0%.2f_N%d_E%d'%(GE.d_x,GE.s_0,GVD['N_FF'],
                                                E_index)
    location = numpy.array([GVD['Sp_x'],GVD['Sp_y'],GVD['Sp_z']],numpy.float64)
    GP = geom.Physics(Energy,GVD['Sp_d'],GVD['Sp_r'],location,GVD['Dx'],
                      GVD['Dy'],GVD['Dz'])
    def check_files(name):
        File = open(name,'r')
        File = open(name+'head','r')
        # Format: 'Nmax: %d\nN_mu: %d\nt: %d'.  Second item in each line.
        L = [int(line.split()[1]) for line in File.readlines()]
        return (L[0],L[1],float(L[2]))
    t0 = time.time()
    try:
        GE.NT_max,GE.T_muons,GE.t = check_files(Radon_T_name)
    except:
        GP.gen_data(GVD['N_mu'])
        GE.writeRadonT(GP,Radon_T_name)
    try:
        GE.NB_max,GE.B_muons,GE.tb = check_files(Radon_B_name)
    except:
        GP.FF(N_FF)
        GE.writeRadonB(GP,Radon_B_name)
    t = 'Read or generated muons in %.1f seconds'%(time.time()-t0)
    message.set('Muons\nN_mu = %d\nN_FF   = %d'%(GVD['N_mu'],GVD['N_FF']))
    print t
    Muons_OK = True
    return
def Base_test(GP):
    """ Selects trajectories that are regularly spaced and uses
    deterministic function to get number of hits for each trajectory.
    """
    global Muons_OK, GE, V_0, Radon_T_name, Radon_B_name
    if not Prior_OK:
        return PriorHelp()
    V_0 = GE.gen_V(GVD['Sp_d'],GVD['Sp_r'],numpy.array([GVD['Sp_x'],
                                  GVD['Sp_y'],GVD['Sp_z']]))
    splot_v('V_0',GE,v=V_0)
    def check_files(name):
        File = open(name,'r')
        File = open(name+'head','r')
        # Format: 'Nmax: %d\nN_mu: %d\nt: %d'.  Second item in each line.
        L = [int(line.split()[1]) for line in File.readlines()]
        return (L[0],L[1],float(L[2]))
    try:
        GE.NT_max,GE.T_muons,GE.t = check_files(Radon_T_name)
        GE.NB_max,GE.B_muons,GE.tb = check_files(Radon_B_name)
        print 'After reading, NT_max=%d, T_muons=%d, t=%d'%(
            GE.NT_max,GE.T_muons,GE.t)
    except:
        GE.gen_test(N_mu,Radon_T_name,Radon_B_name,GP=GP)
        print 'After gen_test, NT_max=%d, T_muons=%d, t=%d'%(
            GE.NT_max,GE.T_muons,GE.t)
    Muons_OK = True
    return
def Test_Rrho():
    """ Substitute for NewMuons() that is deterministic.  Uses true integrals.
    """
    global GE, Radon_T_name, Radon_B_name
    Radon_T_name='tmp/Trho_dx%.1f_s0%.2f_N%d_E%d_D%.2f_r%.1f_x%.1f_y%.1f_z%.1f'%(
      GE.d_x,GE.s_0,GVD['N_mu'],E_index,GVD['Sp_d'],GVD['Sp_r'],GVD['Sp_x'],
      GVD['Sp_y'],GVD['Sp_z'])
    Radon_B_name='tmp/Brho_dx%.1f_s0%.2f_N%d_E%d'%(GE.d_x,GE.s_0,GVD['N_FF'],
                                                   E_index)
    location = numpy.array([GVD['Sp_x'],GVD['Sp_y'],GVD['Sp_z']],numpy.float64)
    GP = geom.Physics(Energy,Sp_d,Sp_r,location,GVD['Dx'],GVD['Dy'],GVD['Dz'])
    Base_test(GP)
    message.set('Test rho\nN_mu=%d'%N_mu)
    return
def Test_Rv():
    """ Substitute for NewMuons() that uses Radon*v instead of true
    integrals.
    """
    global GE, Radon_T_name, Radon_B_name
    Radon_T_name='tmp/Test_dx%.1f_s0%.2f_N%d_E%d_D%.2f_r%.1f_x%.1f_y%.1f_z%.1f'%(
      GE.d_x,GE.s_0,GVD['N_mu'],E_index,GVD['Sp_d'],GVD['Sp_r'],GVD['Sp_x'],
      GVD['Sp_y'],GVD['Sp_z'])
    Radon_B_name='tmp/Best_dx%.1f_s0%.2f_N%d_E%d'%(GE.d_x,GE.s_0,GVD['N_FF'],
                                                   E_index)
    Base_test(None)
    message.set('Test v\nN_mu=%d'%GVD['N_mu'])
    return
fmin_count = 0 # Counts iterations for callback()
def Estimate():
    global Muons_OK, fmin_count
    if not Muons_OK:
        return MuonHelp()
    fmin_count = 0
    def callback(v):
        global fmin_count
        splot_v('V_hat',GE,v=v,count=fmin_count)
        print 'fmin_step=%d'%fmin_count
        fmin_count += 1
    t0 = time.time()
    V_hat = GE.fminLM(Radon_T_name,Radon_B_name,callback=callback)
    print '%d background muons, %d target muons, len(v)=%d'%(
        GE.B_muons,GE.T_muons,len(V_hat))
    print '%f seconds for fminLM'%(time.time()-t0)
    splot_v('V_hat',GE,v=V_hat)
    return
    print '%10s in %6.3f seconds'%(msg,ts[-1][0]-ts[-2][0])
def set_E(ev=None):
    global E_index, Energy, Energies, Muons_OK, Prior_OK
    E_index = TK_E.get()
    Energy = Energies[E_index]
    Muons_OK = False
    Prior_OK = False
    return

############# Service routines for new slider values #########################
def newBasis_d(ev=None): # Center Spacing
    global GVD, Prior_OK
    GVD['Basis_d'] = BASIS_D.get()
    Prior_OK = False

def newBasis_s(ev=None): # Basis function range
    global Prior_OK
    GVD['Basis_s'] = BASIS_S.get()
    Prior_OK = False

def newAlpha(ev=None): # Regularization strength
    global Prior_OK
    t = ALPHA.get()
    GVD['alpha'] = 10.0**t
    Prior_OK = False

def newBeta(ev=None): # Correlation range
    global Prior_OK
    t = BETA.get()
    GVD['beta'] = 10.0**t
    Prior_OK = False

def newN_mu(ev=None): # Number of muons
    global GVD
    GVD['N_mu'] = N_MU.get()
    
def newN_FF(ev=None): # Number of muons for flat field
    global GVD
    GVD['N_FF'] = N_FF_slider.get()

def newSp_d(ev=None): # Density of spherical perturbation
    global GVD, Muons_OK
    GVD['Sp_d'] = SP_D.get()
    Muons_OK = False

def newSp_r(ev=None): # Radius of perturbation
    global GVD, Muons_OK
    GVD['Sp_r'] = SP_R.get()
    Muons_OK = False

def newSp_x(ev=None): # X location of perturbation
    global GVD, Muons_OK
    GVD['Sp_x'] = SP_X.get()
    Muons_OK = False

def newSp_y(ev=None): # Y location of perturbation
    global GVD, Muons_OK
    GVD['Sp_y'] = SP_Y.get()
    Muons_OK = False

def newSp_z(ev=None): # Z location of perturbation
    global GVD, Muons_OK
    GVD['Sp_z'] = SP_Z.get()
    Muons_OK = False

def newV_y(ev=None): # Vslice, ie, Y of plane for viewing density
    global GVD
    old_iy = min(int(GVD['V_y']*GVD['RS_nx']),GVD['RS_nx']-1)
    GVD['V_y'] = V_Y.get()
    iy = min(int(GVD['V_y']*GVD['RS_nx']),GVD['RS_nx']-1)
    if iy == old_iy:
        return
    for key in GPD.keys():
        try:
            splot_v(key,GE)
        except:
            print 'Failed to plot %s'%key
    return
def newLim(ev=None): # Change zrange for viewing density
    global GVD
    t = LIM.get()
    GVD['Lim_val'] = 10.0**t
    for key in GPD.keys():
        GPD[key].write('set zrange[%.1g:%.1g]\n'%(-GVD['Lim_val'],
                                                   GVD['Lim_val']))
        try:
            splot_v(key,GE)
        except:
            print 'Failed to plot %s'%key
    return
def newXYBin(ev=None): # Data rebin resolution
    global GVD
    GVD['XYbin'] = BINXY.get()
    return

#######################################################################
# The stuff below defines the look and initial settings of the GUI and
# ties service routines to the items in the GUI.
root = Tkinter.Tk()
root.title('List Mode GUI')

def Quit():
    #for name in [Radon_T_name,Radon_B_name]:
    #  os.system('rm '+name)
    root.quit  

# b: A row of buttons
b = Tkinter.Frame(root)
b.pack(side=Tkinter.TOP)
QB = Tkinter.Button(b,text='QUIT', command=root.quit)
QB.pack(side=Tkinter.LEFT)
if data_mode == 'simulation':
    NP_button = Tkinter.Button(b,text='New Prior', command=NewPrior)
    NP_button.pack(side=Tkinter.LEFT)
    newmuons = Tkinter.Button(b,text='Muons', command=NewMuons)
    newmuons.pack(side=Tkinter.LEFT)
    TestInt = Tkinter.Button(b,text='Test v', command=Test_Rv)
    TestInt.pack(side=Tkinter.LEFT)
    Testrho = Tkinter.Button(b,text='Test rho', command=Test_Rrho)
    Testrho.pack(side=Tkinter.LEFT)
else: # data_mode == 'data'
    NP_button = Tkinter.Button(b,text='Prior', command=DataPrior)
    NP_button.pack(side=Tkinter.LEFT)
    TestInt = Tkinter.Button(b,text='Radon Vectors', command=DataMuons)
    TestInt.pack(side=Tkinter.LEFT)
est = Tkinter.Button(b,text='Estimate', command=Estimate)
est.pack(side=Tkinter.LEFT)
SP = Tkinter.Button(b,text='Save Plots', command=Save_Plots)
SP.pack(side=Tkinter.LEFT)
HB = Tkinter.Button(b,text='Help', command=ButtonHelp)
HB.pack(side=Tkinter.LEFT)
MODES = [
    ("Power Law", 0),
    ("Exponential", 1),
    ("F(x)=x", 2)
    ]
TK_E = Tkinter.IntVar()
TK_E.set(E_index) # initialize
for text, mode in MODES:
    rb = Tkinter.Radiobutton(b, text=text,
                    variable=TK_E, value=mode, command=set_E)
    rb.pack(anchor=Tkinter.W)

# Row 1: basis functions, prior, number of muons, and detector bins
R1 = Tkinter.Frame(root)
R1.pack(side=Tkinter.TOP)

BasisBlock= Tkinter.Frame(R1,bd=5,relief=Tkinter.RIDGE)
BasisBlock.pack(side=Tkinter.LEFT)
BB1 = Tkinter.Frame(BasisBlock)
BB1.pack(side=Tkinter.TOP)
Tkinter.Label(BB1,text="Basis").pack(side=Tkinter.LEFT)
BBhelp = Tkinter.Button(BB1,text='HELP', command=BasisHelp)
BBhelp.pack(side=Tkinter.LEFT)
BB2 = Tkinter.Frame(BasisBlock)
BB2.pack(side=Tkinter.TOP)
BASIS_D = Tkinter.Scale(BB2, label='d_0',length=201,resolution=0.5,
                  to=(0.5), from_=(10.0),command=newBasis_d)
BASIS_D.set(GVD['Basis_d'])
BASIS_D.pack(side=Tkinter.LEFT)
BASIS_S = Tkinter.Scale(BB2, label='s_0',length=201,resolution=0.05,
                  to=(.3), from_=(.7),command=newBasis_s)
BASIS_S.set(GVD['Basis_s'])
BASIS_S.pack(side=Tkinter.LEFT)

RegularizationBlock= Tkinter.Frame(R1,bd=5,relief=Tkinter.RIDGE)
RegularizationBlock.pack(side=Tkinter.LEFT)
CB1 = Tkinter.Frame(RegularizationBlock)
CB1.pack(side=Tkinter.TOP)
Tkinter.Label(CB1,text="TV Reg (log_10)").pack(side=Tkinter.LEFT)
CBhelp = Tkinter.Button(CB1,text='HELP', command=TVHelp)
CBhelp.pack(side=Tkinter.LEFT)
CB2 = Tkinter.Frame(RegularizationBlock)
CB2.pack(side=Tkinter.TOP)
ALPHA = Tkinter.Scale(CB2, label='alpha',length=201,resolution=0.2,
                  to=(-7), from_=(2),command=newAlpha)
ALPHA.set(math.log10(GVD['alpha']))
ALPHA.pack(side=Tkinter.LEFT)
BETA = Tkinter.Scale(CB2, label='beta',length=201,resolution=0.2,
                  to=(-4), from_=(7),command=newBeta)
BETA.set(math.log10(GVD['beta']))
BETA.pack(side=Tkinter.LEFT)

if data_mode=='simulation': #Simulate rather than experimental data
    NMUBlock= Tkinter.Frame(R1,bd=5,relief=Tkinter.RIDGE) # Number of muons
    NMUBlock.pack(side=Tkinter.LEFT)
    NMUB1 = Tkinter.Frame(NMUBlock)
    NMUB1.pack(side=Tkinter.TOP)
    Tkinter.Label(NMUB1,text="N muon").pack(side=Tkinter.LEFT)
    NMUBhelp = Tkinter.Button(NMUB1,text='HELP', command=NMUHelp)
    NMUBhelp.pack(side=Tkinter.LEFT)
    NMUB2 = Tkinter.Frame(NMUBlock)
    NMUB2.pack(side=Tkinter.TOP)
    N_MU = Tkinter.Scale(
        NMUB2, label='Nmu', length=201, resolution=10000, to=(10000),
        from_=(2000000),command=newN_mu)
    N_MU.set(GVD['N_mu'])
    N_MU.pack(side=Tkinter.LEFT)
    N_FF_slider = Tkinter.Scale(
        NMUB2, label='NFF', length=201, resolution=10000, to=(10000),
        from_=(2000000),command=newN_FF)
    N_FF_slider.set(GVD['N_FF'])
    N_FF_slider.pack(anchor=Tkinter.W)
    
    # Row 2: Describe spherical density perturbation and view slices
    R2 = Tkinter.Frame(root)
    R2.pack(side=Tkinter.TOP)

    message = Tkinter.StringVar()
    Message = Tkinter.Message(R2,textvariable=message,width=200)
    Message.pack(side=Tkinter.LEFT)

    SphereBlock= Tkinter.Frame(R2,bd=5,relief=Tkinter.RIDGE)
    SphereBlock.pack(side=Tkinter.LEFT)
    SB1 = Tkinter.Frame(SphereBlock)
    SB1.pack(side=Tkinter.TOP)
    Tkinter.Label(SB1,text="Spherical Perturbation").pack(side=Tkinter.LEFT)
    SBhelp = Tkinter.Button(SB1,text='HELP', command=SPHelp)
    SBhelp.pack(side=Tkinter.LEFT)
    SB2 = Tkinter.Frame(SphereBlock)
    SB2.pack(side=Tkinter.TOP)
    SP_D = Tkinter.Scale(SB2, label='D',length=201,resolution=0.05,
                    to=(-1), from_=(1),command=newSp_d)
    SP_D.set(GVD['Sp_d'])
    SP_D.pack(side=Tkinter.LEFT)
    SP_R = Tkinter.Scale(SB2, label='r',length=201,resolution=0.1,
                    to=(0.1), from_=(20),command=newSp_r)
    SP_R.set(GVD['Sp_r'])
    SP_R.pack(side=Tkinter.LEFT)
    SP_X = Tkinter.Scale(SB2, label='x',length=201,resolution=0.2,
                    to=(0), from_=(GVD['Dx']),command=newSp_x)
    SP_X.set(GVD['Sp_x'])
    SP_X.pack(side=Tkinter.LEFT)
    SP_Y = Tkinter.Scale(SB2, label='y',length=201,resolution=0.1,
                    to=(0), from_=(GVD['Dy']),command=newSp_y)
    SP_Y.set(GVD['Sp_y'])
    SP_Y.pack(side=Tkinter.LEFT)
    SP_Z = Tkinter.Scale(SB2, label='z',length=201,resolution=0.1,
                    to=(0), from_=(GVD['Dz']),command=newSp_z)
    SP_Z.set(GVD['Sp_z'])
    SP_Z.pack(side=Tkinter.LEFT)

    VsliceBlock= Tkinter.Frame(R2,bd=5,relief=Tkinter.RIDGE)

else: # Analyze experimental data
    RebinBlock= Tkinter.Frame(R1,bd=5,relief=Tkinter.RIDGE)
    RebinBlock.pack(side=Tkinter.LEFT)
    RBB1 = Tkinter.Frame(RebinBlock)
    RBB1.pack(side=Tkinter.TOP)
  
    Tkinter.Label(RBB1,text="Rebin").pack(side=Tkinter.LEFT)
    RBBhelp = Tkinter.Button(RBB1,text='HELP', command=RebinHelp)
    RBBhelp.pack(side=Tkinter.LEFT)
    RBB2 = Tkinter.Frame(RebinBlock)
    RBB2.pack(side=Tkinter.TOP)
    BINXY = Tkinter.Scale(RBB2, label='resolution',length=201,resolution=0.5,
                    to=(0), from_=(25),command=newXYBin)
    BINXY.set(GVD['XYbin'])
    BINXY.pack(side=Tkinter.LEFT)

    VsliceBlock= Tkinter.Frame(R1,bd=5,relief=Tkinter.RIDGE)
VsliceBlock.pack(side=Tkinter.LEFT)
VSB1 = Tkinter.Frame(VsliceBlock)
VSB1.pack(side=Tkinter.TOP)
Tkinter.Label(VSB1,text="Vslice").pack(side=Tkinter.LEFT)
VSBhelp = Tkinter.Button(VSB1,text='HELP', command=VsliceHelp)
VSBhelp.pack(side=Tkinter.LEFT)
VSB2 = Tkinter.Frame(VsliceBlock)
VSB2.pack(side=Tkinter.TOP)
V_Y = Tkinter.Scale(VSB2, label='y',length=201,resolution=0.01,
                    to=(0), from_=(1),command=newV_y)
V_Y.set(GVD['V_y'])
V_Y.pack(side=Tkinter.LEFT)
LIM = Tkinter.Scale(VSB2, label='Log range',length=201,resolution=0.2,
                    to=(-3), from_=(2),command=newLim)
LIM.set(math.log10(GVD['Lim_val']))
LIM.pack(side=Tkinter.LEFT)

# ============================================================

if len(sys.argv) < 3:  # Start the GUI
    Tkinter.mainloop()
else: # Code to read commands from a file.  Example formats:
    """
    set Basis_d 5.0
    call NewPrior
    call Test_Rv
    call Estimate
    call Save_Plots plot1 plot2
    """
    calls = {'NewPrior':lambda:NewPrior(), 'DataPrior':lambda:DataPrior(),
             'DataMuons':lambda:DataMuons(),'NewMuons':lambda:NewMuons(),
             'Test_Rho':lambda:Test_Rho(),'Test_Rv':lambda:Test_Rv(),
             'Estimate':lambda:Estimate(),'Save_Plots':lambda:Save_Plots()}
    # Next parse the file of commands from the file named in argv[2]
    for line in open(sys.argv[2],'r').readlines():
        parts = line.split()
        if len(parts) == 0 or parts[0] == '#':
            continue
        assert parts[0] == 'call' or parts[0] == 'set',\
            "line=\n%sCan't interpret parts[0]=%s"%(line,parts[0])
        if parts[0] == 'set':
            if parts[1] == 'Energy':
                E_index = int(parts[2])
                Energy = Energies[E_index]
            if parts[1][0] == 'N':
                GVD[parts[1]] = int(parts[2])
            else:
                GVD[parts[1]] = float(parts[2])
        if parts[0] == 'call':
            if parts[1] == 'Save_Plots':
                assert len(parts)-2 == len(keylist),\
                    'len(parts)=%d, len(keylist)=%d'%(len(parts), len(keylist))
                if len(keylist) == 2:
                    D = {'V_0':parts[2],'V_hat':parts[3]}
                if len(keylist) == 1:
                    D = {'V_hat':parts[2]}
                Save_Plots(save_names=D)
            else:
                calls[parts[1]]()


#Local Variables:
#mode:python
#End:
