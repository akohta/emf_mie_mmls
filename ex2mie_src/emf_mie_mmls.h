#if !defined EMF_MIE_MMLS_H_
#define EMF_MIE_MMLS_H_

#include "multi_fbeam.h"
#include "rctbess.h"
#include "gsl/gsl_specfunc.h"
#include "my_utils.h"
#include "osu_mksa.h"

// default sphere datafile name
#define fn_mlsphr "mlsphr.txt"

// for iterative operations
#define ito_max 256
#define ito_eps 1.0e-12
#define ito_breakcount 2

// struct definition
typedef struct discrete_data{
  int    nt,np;              // number of sampling points
  double *xt,*wt,*xp,*wp;    // Gauss-Legendre point and weight data 
  double complex *eri,*hri;  // incident field on sphere surface (gauss-legendre node data). Er(r,theta,phi),Hr(r,theta,phi) 
  double complex *ers,*hrs;  // scattered field on sphere surfrace.
  int l_max;                 // limit of Ricatti-Bessel function order 
  double *cab;               // coefficient for Alm,Blm
  double complex *ceM,*chM;  // characteristic matrix for electromagnetic field coefficient
  double complex *ca,*cb;    // coefficient for alm,blm of scattered field
  double complex *ca0,*cb0;  // coefficient for Alm,Blm of layer number 0 field
  double complex **Alm,**Blm;// all field coefficient
  double complex **alm,**blm;// all field coefficient 
}DDT;

typedef struct sphere_data{
  int bsn;              // basic sampling number for Gauss-Legendre quadrature
  int bdv;              // division number of sphere surface ( per M_PI )
  int l_limit;          // l limit of order number ( limit of summation )
  double xs,ys,zs;      // sphere center (xs,ys,zs)
  int n_l;              // number of layers 
  double *a;            // radius of layers   
  double complex *ns;   // refractive index of layers 
  DDT ddt;
}SPD;

typedef struct sphere_objct{
  int n_sphr;   // number of spheres
  SPD *sp;      // sphere data 
  Bobj bm;      // multi_fbeam data
}MSPD;

// ---- emf_mie_mmls.c ----
void read_data_mmls(MSPD *msp);      // seach the sphere datafile (defined fn_mlsphr) and load 
void print_data_mmls(MSPD *msp);     // print loaded data
void print_data_mmls_mksa(MSPD *msp);// print loaded data in MKSA system of units
void setup_mmls(MSPD *msp);          // allocate memory and setup coefficients
void  free_mmls(MSPD *msp);          // free allocated memory
void iterative_ops_mmls(MSPD *msp);  // solve multiple scattering

// ---- emf_mie_mmls_field.c ----
void  incident_EH_mmls(double complex *e,double complex *h,double *r,MSPD *msp); // calculate incident field 
void  internal_EH_mmls(double complex *e,double complex *h,double *r,MSPD *msp); // calculate internal field 
void scattered_EH_mmls(double complex *e,double complex *h,double *r,MSPD *msp); // calculate scattered field
void     total_EH_mmls(double complex *e,double complex *h,double *r,MSPD *msp); // calculate total field
// e[0]=Ex,e[1]=Ey,e[2]=Ez,h[0]=Hx,h[1]=Hy,h[2]=Hz,r[0]=x,r[1]=y,r[2]=z

// ---- emf_mie_mmls_force.c ----
void  force_mmls(int id,double *vf,MSPD *msp);                 // calculate radiation force of the sphere 
void torque_mmls(int id,double *vn,MSPD *msp);                 // calculate radiation torque of the sphere 
void force_torque_mmls(int id,double *vf,double *vn,MSPD *msp); // calculate radiation force and torque of the sphere
// id : sphere id, vf[0]=Fx,vf[1]=Fy,vf[2]=Fz,vn[0]=Nx,vn[1]=Ny,vn[2]=Vz

// ---- emf_mie_mmls_dat.c ----
void write_dat_mmls(char *fn,MSPD *msp); // write datafile
void  read_dat_mmls(char *fn,MSPD *msp); // read datafile and allocate memory

#endif
