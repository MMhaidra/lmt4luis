#include <Python.h>
#include "numpy/arrayobject.h"
#include "numpy/arraybase.h"
#define max(A,B) ( (A) > (B) ? (A):(B)) 
#define min(A,B) ( (A) < (B) ? (A):(B))
double exp_cut = -7.0;
/* http://docs.python.org/extending/
   http://docs.python.org/c-api/memory.html
   http://www.scipy.org/Cookbook/C_Extensions
   http://www.scipy.org/Cookbook/C_Extensions/NumPy_arrays Pecora
   /usr/src/scipy/scipy/sparse/linalg/dsolve/_ssuperlumodule.c
   http://docs.python.org/c-api/arg.html#arg-parsing

From http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html#scipy.sparse.csr_matrix

    csr_matrix((data, indices, indptr), [shape=(M, N)]) is the
    standard CSR representation where the column indices for row i are
    stored in indices[indptr[i]:indices[i+1]] and their corresponding
    values are stored in data[indptr[i]:indptr[i+1]]. If the shape
    parameter is not supplied, the matrix dimensions are inferred from
    the index arrays.
*/

void print3i(char* name, int* vec){
  int i;
  fprintf(stderr,"%s=",name);
  for (i=0;i<3;i++)
    fprintf(stderr," %2d ",vec[i]);
  fprintf(stderr,"\n");
}
void print3d(char* name, double* vec){
  int i;
  fprintf(stderr,"%s=",name);
  for (i=0;i<3;i++)
    fprintf(stderr," %6lg ",vec[i]);
  fprintf(stderr,"\n");
}
double d2oca(double *c, double theta, double phi, double *p){
  /* Calculate squared distance between point c and trajectory of
    incoming ray given by (theta,phi,p)
   */
  int i;
  double r,d2,dx[3],u[3];

  for(i=0;i<3;i++)
    dx[i] = c[i] - p[i];
  u[0] = sin(theta)*cos(phi); u[1] = sin(theta)*sin(phi); u[2] = cos(theta);
  for(r=0,i=0;i<3;i++)
    r += dx[i]*u[i];
  for(d2=0,i=0;i<3;i++){
    double t = dx[i]-r*u[i];
    d2 += t*t;
  }
  return d2;
}

int closest_uz(double x_0, double y_0, double theta, double phi, int iz,
	    double var2, double z, double d_x, double d_y, int N_x,
	    int N_y, PyArrayObject *Py_centers_xyz, int *best){
  /* Find the triple of indices for the center that is closest to the
        ray (x_0,y_0,theta,phi) in the iz plane and put it in *best.
        If the result is r and -r^2/var2 < exp_cut, return -1
  */
  double r,rr,p[3],*c,smallest;
  int ixc, iyc,ix,ix_,iy,iy_;
  
  p[0] = x_0; p[1] = y_0; p[2] = 0;
  r = z/cos(theta); //Distance along the ray to the iz plane
  ixc = ((x_0+r*sin(theta)*cos(phi))/d_x);
  iyc = ((y_0+r*sin(theta)*sin(phi))/d_y); // integer coords of intersection
  ix = min(max(0,ixc),N_x-1);
  iy = min(max(0,iyc),N_y-1);
  c = PyArray_GETPTR4(Py_centers_xyz,ix,iy,iz,0);
  best[0] = ix; best[1] = iy; best[2] = iz;
  smallest = d2oca(c,theta,phi,p);
  for(ix_=ixc-1;ix_<ixc+1;ix_++){
    ix = min(max(0,ix_),N_x-1);
    for(iy_=iyc-1;iy_<iyc+1;iy_++){
      iy = min(max(0,iy_),N_y-1);
      c = PyArray_GETPTR4(Py_centers_xyz,ix,iy,iz,0);
      rr = d2oca(c,theta,phi,p);
      if (rr < smallest){
	smallest = rr;
	best[0] = ix; best[1] = iy; best[2] = iz;
      }
    }
  }
  if (-smallest/var2 < exp_cut)
    return -1;
  return 1;
}

static PyObject *
amf_readR(PyObject *self, PyObject *args)
/* Read indices and data for sparse Radon matrix.  Call is
   amf.readR(N,indices,data,name)
 */
{
  PyArrayObject *Py_indices, *Py_data;
  int N;
  char *name;
  /* end of argument variables */
  int i,index,dumb,*ip;
  float datum,*fp;
  FILE *file;
  if (!PyArg_ParseTuple(args, "iO&O&s",
			&N,
			PyArray_Converter, &Py_indices,
			PyArray_Converter, &Py_data,
			&name
			))
    return NULL;
  file  = fopen(name,"r");
  for(i=0;i<N;i++){
    if (fread(&dumb,sizeof(dumb),1,file) <= 0)
      fprintf(stderr,"FixMe fread failed\n");
    if (fread(&index,sizeof(index),1,file) <= 0)
      fprintf(stderr,"FixMe fread failed\n");
    if (fread(&datum,sizeof(datum),1,file) <= 0)
      fprintf(stderr,"FixMe fread failed\n");
    ip = (int *)PyArray_GETPTR1(Py_indices,i);
    *ip = index;
    fp = (float *)PyArray_GETPTR1(Py_data,i);
    *fp = datum;
  }
  fclose(file);
  Py_DECREF(Py_indices);
  Py_DECREF(Py_data);
  return Py_BuildValue("");
}

static PyObject *
amf_uline(PyObject *self, PyObject *args)
/* Calculate Radon vector for ray specified by (p_0,theta,phi)

   Here is the python call: uline(
       Zs,            # Array of z plane elevations, float64
       D_xy_u,        # Array of integers.  Offsets to consider from poca at z
       theta,
       phi,
       d_x,           # Center spacing in x direction
       d_y,           # Center spacing in y direction
       p_0,           # Array [x,y,0] float64.  Point that ray crosses z=0
       centers_xyz,   # Array with shape (N_centers,3), float64
       norm,          # Normalization for basis functions
       var2,          # Denominator for basis function exponents
       N_x,           # N_centers = N_x*N_y*N_z
       N_y,
       N_z,
       i_u,           # Index of ray
       name           # Write results to file thusly named
       )

  This function writes triples (i_u,j,ev) to the file and returns the
  number of triples written.  i_u is the index of the ray or the index
  of u, j is the index of v and ev = Radon[i_u,j].
 */
{
  PyArrayObject *Py_Zs, *Py_D_xy_u, *Py_p_0, *Py_centers_xyz;
  double theta,phi,d_x,d_y,var2,norm;
  int Nx,Ny,Nz,i_u;
  double *p_0, *Zs;
  /* end of argument variables */
  int i,j,k,sum[3],m,n,NXYZ[3],iz,*dxyz;
  double *center, x_0, y_0, z;
  float ev;
  char *name;
  FILE *file;
  if (!PyArg_ParseTuple(args, "O&O&ddddO&O&ddiiiis",
			PyArray_Converter, &Py_Zs,
			PyArray_Converter, &Py_D_xy_u,
			&theta,
			&phi,
			&d_x,
			&d_y,
			PyArray_Converter, &Py_p_0,
			PyArray_Converter, &Py_centers_xyz,
			&norm,
			&var2,
			&Nx,
			&Ny,
			&Nz,
			&i_u,
			&name
			))
    return NULL;
  // FixMe: check types with PyArray TYPE pg 239 of Oliphant
  file = fopen(name,"a");
  Zs     =(double *) PyArray_DATA(Py_Zs);
  p_0     =(double *) PyArray_DATA(Py_p_0);
  x_0 = p_0[0]; y_0 = p_0[1];
  NXYZ[0] = Nx; NXYZ[1] = Ny; NXYZ[2] = Nz;
  m = 0;
  for(iz=0;iz<Nz;iz++){
    int ixyz[3];
    z = Zs[iz];
    if (closest_uz(x_0, y_0, theta, phi, iz, var2, z, d_x, d_y,  Nx, Ny,
		   Py_centers_xyz, ixyz) < 0) continue;
    n = PyArray_DIM (Py_D_xy_u, 0);
    for (k=0;k<n;k++){
      dxyz = PyArray_GETPTR2(Py_D_xy_u,k,0);
      for(j=0,i=0;i<3;i++){
	sum[i] = ixyz[i] + dxyz[i];
	if ((NXYZ[i] <= sum[i]) || (sum[i]<0)){
	  j = -1;
	  break;
	}
      }
      if (j < 0)
	continue;
      center = (double *)PyArray_GETPTR4(Py_centers_xyz,sum[0],sum[1],sum[2],0);
      ev = -d2oca(center, theta, phi, p_0)/var2;
      if (ev < exp_cut){
	continue;
      }
      j = sum[0]*Ny*Nz + sum[1]*Nz + sum[2];
      ev = norm*exp(ev);
      if (fwrite(&i_u,sizeof(i_u),1,file) <=0)
	fprintf(stderr,"FixMe fwrite failed\n");
      if (fwrite(&j,sizeof(j),1,file) <=0)
	fprintf(stderr,"FixMe fwrite failed\n");
      if (fwrite(&ev,sizeof(ev),1,file) <=0)
	fprintf(stderr,"FixMe fwrite failed\n");
      m += 1;
    }
  }
  fclose(file);
  Py_DECREF(Py_Zs);
  Py_DECREF(Py_D_xy_u);
  Py_DECREF(Py_p_0);
  Py_DECREF(Py_centers_xyz);
  return Py_BuildValue("i",m);
}

static PyObject *
amf_rebin(PyObject *self, PyObject *args)
/* Rebin count data.  Call is amf.rebin(N, pars, counts_in, count_out, flag,
                      Dx, Nx,
                      Dy, Ny,
                      Dtheta, Ntheta,
                      bin_pointers)
 */
{
  PyArrayObject *Py_pars, *Py_counts_in, *Py_count_out, *Py_flag,
    *Py_bin_pointers;
  int N, Nx, Ny, Ntheta, Nphi;
  double Dx, Dy, Dtheta;
  /* end of argument variables */
  float *fp;
  int  *ip, i, j, t, i_x, i_y, i_theta, bin, Nx1, Ny1, Ntheta1, index;
  if (!PyArg_ParseTuple(args, "iO&O&O&O&dididiO&",
			&N,
			PyArray_Converter, &Py_pars,
			PyArray_Converter, &Py_counts_in,
			PyArray_Converter, &Py_count_out,
			PyArray_Converter, &Py_flag,
			&Dx,
			&Nx,
			&Dy,
			&Ny,
			&Dtheta,
			&Ntheta,
			PyArray_Converter, &Py_bin_pointers
			))
    return NULL;
  Nx1 = Nx-1;
  Ny1 = Ny-1;
  Ntheta1 = Ntheta-1;
  for(i=1;i<N;i++){
    fp = (float *)PyArray_GETPTR2(Py_pars,i,0);
    j = (*fp/Dx)*Nx;
    i_x = min(j,Nx1);

    fp = (float *)PyArray_GETPTR2(Py_pars,i,1);
    j = (*fp/Dy)*Ny;
    i_y = min(j,Ny1);
    
    fp = (float *)PyArray_GETPTR2(Py_pars,i,2);
    j = (*fp/Dtheta)*Ntheta;
    i_theta = min(j,Ntheta1);

    ip = (int *)PyArray_GETPTR4(Py_bin_pointers,i_x,i_y,i_theta,0);
    Nphi = ip[0];
    index = ip[1];
    fp = (float *)PyArray_GETPTR2(Py_pars,i,3);
    j = (*fp/(6.2831853071795862))*Nphi; 
    bin = index + min(j,Nphi-1);
    
    ip = (int *)PyArray_GETPTR1(Py_counts_in,i);
    j = *ip;
    t = *(int *)PyArray_GETPTR1(Py_count_out,bin);
    if((j<0) || ((t+j) < 0)){
      fprintf(stderr,"Counter overflow j=%d, t=%d, t+j=%d\n",j,t,t+j);
      return NULL;
    }
    *(int *)PyArray_GETPTR1(Py_count_out,bin) = t+j;
    *(char *)PyArray_GETPTR1(Py_flag,bin) += 1;
  }
  return Py_BuildValue("");
}

static PyObject *
amf_read_count(PyObject *self, PyObject *args)
/* Read count data.  Call is amf.read_count(pars,counts,name)
 */
{
  PyArrayObject *Py_pars, *Py_counts;
  char *name;
  /* end of argument variables */
  int i,N,*ip;
  float *fp;
  FILE *file;
  if (!PyArg_ParseTuple(args, "O&O&s",
			PyArray_Converter, &Py_pars,
			PyArray_Converter, &Py_counts,
			&name
			))
    return NULL;
  file  = fopen(name,"r");
  if (fread(&N,sizeof(N),1,file) <= 0)
    fprintf(stderr,"FixMe fread failed\n");
  for(i=0;i<N;i++){
    fp = (float *)PyArray_GETPTR2(Py_pars,i,0);
    if (fread(fp,sizeof(*fp),4,file) <= 0)
      fprintf(stderr,"FixMe fread failed\n");
    ip = (int *)PyArray_GETPTR1(Py_counts,i);
    if (fread(ip,sizeof(*ip),1,file) <= 0)
      fprintf(stderr,"FixMe fread failed\n");
  }
  fclose(file);
  Py_DECREF(Py_pars);
  Py_DECREF(Py_counts);
  return Py_BuildValue("i",N);
}

static PyObject *
amf_write_count(PyObject *self, PyObject *args)
/* Read count data.  Call is amf.read_count(N,pars,counts,name)
 */
{
  PyArrayObject *Py_pars, *Py_counts;
  int N;
  char *name;
  /* end of argument variables */
  int i,*ip;
  float *fp;
  FILE *file;
  if (!PyArg_ParseTuple(args, "iO&O&s",
			&N,
			PyArray_Converter, &Py_pars,
			PyArray_Converter, &Py_counts,
			&name
			))
    return NULL;
  file  = fopen(name,"w");
  if (fwrite(&N,sizeof(N),1,file) <=0)
    fprintf(stderr,"FixMe fwrite failed\n");
  for(i=0;i<N;i++){
    fp = (float *)PyArray_GETPTR2(Py_pars,i,0);
    if (fwrite(fp,sizeof(*fp),4,file) <=0)
      fprintf(stderr,"FixMe fwrite failed\n");
    ip = (int *)PyArray_GETPTR1(Py_counts,i);
    if (fwrite(ip,sizeof(*ip),1,file) <=0)
      fprintf(stderr,"FixMe fwrite failed\n");
  }
  fclose(file);
  Py_DECREF(Py_pars);
  Py_DECREF(Py_counts);
  return Py_BuildValue("");
}

static PyObject * amf_writeRadon(PyObject *self, PyObject *args)
/* Python call: writeRadon(n,indices,values,theta,nk,file) Writes
   n,indices,values,nk to file which is already open */
{
  PyArrayObject *Py_indices, *Py_values;
  PyObject *pyfile;
  int n;
  double theta, nk;
  FILE *file;
  /* end of argument variables */
  int *indices, rv;
  float *values; 

  if (!PyArg_ParseTuple(args, "iO&O&ddO!",
			&n,
			PyArray_Converter, &Py_indices,
			PyArray_Converter, &Py_values,
			&theta,
			&nk,
			&PyFile_Type, &pyfile
			))
    return NULL;
  file = PyFile_AsFile(pyfile);
  indices = (int *)PyArray_GETPTR1(Py_indices,0);
  values = (float *)PyArray_GETPTR1(Py_values,0);
  rv = fwrite(&n,sizeof(n),1,file);
  if (rv <=0)
    fprintf(stderr,"FixMe fwrite n=%d failed. rv=%d\n",n,rv);
  if (fwrite(indices,sizeof(*indices),n,file) <=0){
    PyErr_SetString(PyExc_IOError,"fwrite indices failed\n");
    Py_DECREF(Py_indices);
    Py_DECREF(Py_values);
    return NULL;
  }    
  if (fwrite(values,sizeof(*values),n,file) <=0)
    fprintf(stderr,"FixMe fwrite values failed\n");
  if (fwrite(&theta,sizeof(theta),1,file) <=0)
    fprintf(stderr,"FixMe fwrite theta failed\n");
  if (fwrite(&nk,sizeof(nk),1,file) <=0)
    fprintf(stderr,"FixMe fwrite nk failed\n");
  Py_DECREF(Py_indices);
  Py_DECREF(Py_values);
  return Py_BuildValue("");
}

static PyObject * amf_readRadon(PyObject *self, PyObject *args)
/* Python call: n,theta,nk = readRadon(indices,values,file) Reads
   n,indices,values,nk to file which is already open */
{
  PyArrayObject *Py_indices, *Py_values;
  PyObject *pyfile;
  int n;
  double theta, nk;
  FILE *file;
  /* end of argument variables */
  int *indices, rv;
  float *values; 

  if (!PyArg_ParseTuple(args, "O&O&O!",
			PyArray_Converter, &Py_indices,
			PyArray_Converter, &Py_values,
			&PyFile_Type, &pyfile
			))
    return NULL;
  file = PyFile_AsFile(pyfile);
  indices = (int *)PyArray_GETPTR1(Py_indices,0);
  values = (float *)PyArray_GETPTR1(Py_values,0);
  rv = fread(&n,sizeof(n),1,file);
  if (rv <=0)
    fprintf(stderr,"FixMe fread n=%d failed. rv=%d\n",n,rv);
  if (fread(indices,sizeof(*indices),n,file) <=0)
    fprintf(stderr,"FixMe fread indices failed\n");
  if (fread(values,sizeof(*values),n,file) <=0)
    fprintf(stderr,"FixMe fread values failed\n");
  if (fread(&theta,sizeof(theta),1,file) <=0)
    fprintf(stderr,"FixMe fread theta failed\n");
  if (fread(&nk,sizeof(nk),1,file) <=0)
    fprintf(stderr,"FixMe fread nk failed\n");
  Py_DECREF(Py_indices);
  Py_DECREF(Py_values);
  return Py_BuildValue("idd",n,theta,nk);
}

static PyObject * amf_dot(PyObject *self, PyObject *args)
/* Python call: f = dot(v,n,indices,values) The sparse vector w with n
   nonzero entries is represented by (indices,values) and v is a dense
   vector.  This function returns the dot product of w and v.  */
{
  PyArrayObject *Py_indices, *Py_values, *Py_v;
  int n;
  /* end of argument variables */
  double rv=0.0;
  int *indices, i;
  float *values;
  double *v;

  if (!PyArg_ParseTuple(args, "O&iO&O&",
			PyArray_Converter, &Py_v,
			&n,
			PyArray_Converter, &Py_indices,
			PyArray_Converter, &Py_values
			))
    return NULL;
  indices = (int *)PyArray_GETPTR1(Py_indices,0);
  values = (float *)PyArray_GETPTR1(Py_values,0);
  v = (double *)PyArray_GETPTR1(Py_v,0);
  for (i=0;i<n;i++)
    rv += v[indices[i]]*values[i];
  Py_DECREF(Py_v);
  Py_DECREF(Py_indices);
  Py_DECREF(Py_values);
  return Py_BuildValue("d",rv);
}

static PyObject * amf_accumulate(PyObject *self, PyObject *args)
/* Python call: accumulate(v,c,n,indices,values) The sparse vector
   w with n nonzero entries is represented by (indices,values), v is a
   dense vector, and c is a scalar.  This function performs the
   following: v += c*w.  */
{
  PyArrayObject *Py_indices, *Py_values, *Py_v;
  int n;
  double c;
  /* end of argument variables */
  int *indices, i;
  float *values;
  double *v;

  if (!PyArg_ParseTuple(args, "O&diO&O&",
			PyArray_Converter, &Py_v,
			&c, &n,
			PyArray_Converter, &Py_indices,
			PyArray_Converter, &Py_values
			))
    return NULL;
  indices = (int *)PyArray_GETPTR1(Py_indices,0);
  values = (float *)PyArray_GETPTR1(Py_values,0);
  v = (double *)PyArray_GETPTR1(Py_v,0);
  for (i=0;i<n;i++)
    v[indices[i]] += c*values[i];
  Py_DECREF(Py_v);
  Py_DECREF(Py_indices);
  Py_DECREF(Py_values);
  return Py_BuildValue("");
}

static PyMethodDef AmfMethods[] = {
    {"uline",  amf_uline, METH_VARARGS, "Calculate row of Radon matrix"},
    {"readR",  amf_readR, METH_VARARGS, "Calculate row of Radon matrix"},
    {"write_count",  amf_write_count, METH_VARARGS, "Write data cube to disk"},
    {"read_count",  amf_read_count, METH_VARARGS, "Read data cube from disk"},
    {"rebin",  amf_rebin, METH_VARARGS, "Calculate counts for coarser bins"},
    {"writeRadon",  amf_writeRadon, METH_VARARGS, "Write data for one hit"},
    {"readRadon",  amf_readRadon, METH_VARARGS, "Read data for one hit"},
    {"dot",amf_dot, METH_VARARGS,
  "f = dot(v,n,indices,values)\n Returns product of dense and sparse vectors"},
    {"accumulate",  amf_accumulate, METH_VARARGS,
     "accumulate(v,c,n,indices,values)\n Adds scalar multiple of sparse vector to dense vector"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};
PyMODINIT_FUNC
initamf(void)
{
    (void) Py_InitModule("amf", AmfMethods);
    import_array();
}
