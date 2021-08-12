#include "emf_mie_mmls.h"

void incident_EH_mmls(double complex *e,double complex *h,double *r,MSPD *msp)
{
  calc_mfb_EH(e,h,r,&(msp->bm));
}

void internal_EH_mmls(double complex *e,double complex *h,double *r,MSPD *msp)
{
  int layer_id(double r,SPD *sp); //emf_mie_mmls.c. select inside layer when r is on boundary.
  void internal_EH(int id,double complex *e3,double complex *h3,double *x3,SPD *sp,Bobj *bm); // emf_mie_mmls.c
  
  int i,id;
  double x,y,z,rd;
 
  for(i=0;i<msp->n_sphr;i++){
    x=r[0]-msp->sp[i].xs;    y=r[1]-msp->sp[i].ys;    z=r[2]-msp->sp[i].zs;
    rd=sqrt(x*x+y*y+z*z);
    id=layer_id(rd,&(msp->sp[i]));
    if(id!=msp->sp[i].n_l){
      internal_EH(id,e,h,r,&(msp->sp[i]),&(msp->bm));
      return;
    }
  }
  e[0]=0.0;  e[1]=0.0;  e[2]=0.0;
  h[0]=0.0;  h[1]=0.0;  h[2]=0.0; 
}

void scattered_EH_mmls(double complex *e,double complex *h,double *r,MSPD *msp)
{
  void scattered_EH(double complex *e3,double complex *h3,double *x3,SPD *sp,Bobj *bm); // emf_mie_mmls.c
  
  double complex te[3],th[3];
  int i,j;
  double x,y,z,rd;
  
  for(j=0;j<3;j++){
    e[j]=0.0;    h[j]=0.0;
  }
  // internal
  for(i=0;i<msp->n_sphr;i++){
    x=r[0]-msp->sp[i].xs;    y=r[1]-msp->sp[i].ys;    z=r[2]-msp->sp[i].zs;
    rd=sqrt(x*x+y*y+z*z);
    if(rd<msp->sp[i].a[msp->sp[i].n_l-1]) { // select outside layer ( scattered field ) when rd is on boundary
      e[0]=0.0;      e[1]=0.0;      e[2]=0.0;
      h[0]=0.0;      h[1]=0.0;      h[2]=0.0;
      return;
    }
  }
  // external
  for(i=0;i<msp->n_sphr;i++){
    scattered_EH(te,th,r,&(msp->sp[i]),&(msp->bm));
    for(j=0;j<3;j++){
      e[j]+=te[j];    h[j]+=th[j];
    }
  }
}

void total_EH_mmls(double complex *e,double complex *h,double *r,MSPD *msp)
{
  int layer_id(double r,SPD *sp); //emf_mie_mmls.c. select inside layer when r is on boundary.
  void internal_EH(int id,double complex *e3,double complex *h3,double *x3,SPD *sp,Bobj *bm); // emf_mie_mmls.c
  void scattered_EH(double complex *e3,double complex *h3,double *x3,SPD *sp,Bobj *bm); // emf_mie_mmls.c

  double complex te[3],th[3];
  int id,i,j;
  double x,y,z,rd;
  
  for(j=0;j<3;j++){
    e[j]=0.0;    h[j]=0.0;
  }
  // internal  
  for(i=0;i<msp->n_sphr;i++){
    x=r[0]-msp->sp[i].xs;    y=r[1]-msp->sp[i].ys;    z=r[2]-msp->sp[i].zs;
    rd=sqrt(x*x+y*y+z*z);
    id=layer_id(rd,&(msp->sp[i]));
    if(id!=msp->sp[i].n_l){
      internal_EH(id,e,h,r,&(msp->sp[i]),&(msp->bm));
      return;
    }
  }
  // external
  for(i=0;i<msp->n_sphr;i++){
    scattered_EH(te,th,r,&(msp->sp[i]),&(msp->bm));
    for(j=0;j<3;j++){
      e[j]+=te[j];      h[j]+=th[j];
    }
  }
  // incident field  
  calc_mfb_EH(te,th,r,&(msp->bm));
  for(j=0;j<3;j++){
    e[j]+=te[j];    h[j]+=th[j];
  }
}


