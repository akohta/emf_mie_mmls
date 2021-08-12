#include "emf_mie_mmls.h"

void  force_mmls(int id,double *vf,MSPD *msp)
{
  void force_sphr(double *f3,SPD *sp,Bobj *bm);
  
  if(id>=msp->n_sphr){
    printf("Sphere ID error! id=%d. Exit...\n",id);
    exit(1);
  }
  force_sphr(vf,&(msp->sp[id]),&(msp->bm));  
}

void torque_mmls(int id,double *vn,MSPD *msp)
{
  void torque_sphr(double *n3,SPD *sp,Bobj *bm);
  
  if(id>=msp->n_sphr){
    printf("Sphere ID error! id=%d. Exit...\n",id);
    exit(1);
  }
  torque_sphr(vn,&(msp->sp[id]),&(msp->bm));
}

void force_torque_mmls(int id,double *vf,double *vn,MSPD *msp)
{ 
  void force_torque_sphr(double *f3,double *n3,SPD *sp,Bobj *bm);
  
  if(id>=msp->n_sphr){
    printf("Sphere ID error! id=%d. Exit...\n",id);
    exit(1);
  } 
  force_torque_sphr(vf,vn,&(msp->sp[id]),&(msp->bm));
}

///////////////////////////////////////////////////////////////////////
void force_sphr(double *f3,SPD *sp,Bobj *bm)
{
  double complex pickup_Alm(int l,int m,SPD *sp); // emf_mie_mmls.c
  double complex pickup_Blm(int l,int m,SPD *sp); // emf_mie_mmls.c
  double complex pickup_alm(int l,int m,SPD *sp); // emf_mie_mmls.c
  double complex pickup_blm(int l,int m,SPD *sp); // emf_mie_mmls.c
  
  double complex A00,B00,App,Bpp,Apm,Bpm,A0p,B0p,Ap0,Bp0;
  double complex a00,b00,app,bpp,apm,bpm,a0p,b0p,ap0,bp0;
  double complex tfxy=0.0,tfz=0.0;
  double c0p,c0m,c1p,c1m,c2p,c2m;
  double cep=pow(bm->n_0,2)*epsilon0;
  double mu_0=mu0;
  double ke=2.0*M_PI*bm->n_0/bm->lambda_0;
  int lm=sp->l_limit;
  int l,m;
    
  m=0;
  for(l=1;l<lm;l++){
    c0p=sqrt((double)((l+m+2)*(l+m+1))/(double)((2*l+1)*(2*l+3)))*(double)(l*(l+2));
    c0m=sqrt((double)((l-m+2)*(l-m+1))/(double)((2*l+1)*(2*l+3)))*(double)(l*(l+2));
    c1p=sqrt((double)((l+m+1)*(l-m)))*sqrt(cep*mu_0);
    c1m=sqrt((double)((l-m+1)*(l+m)))*sqrt(cep*mu_0);
    c2p=sqrt((double)((l-m+1)*(l+m+1))/(double)((2*l+3)*(2*l+1)))*(double)(l*(l+2));
    c2m=sqrt((double)((l+m+1)*(l-m+1))/(double)((2*l+3)*(2*l+1)))*(double)(l*(l+2));
    A00=pickup_Alm(l,m,sp);          a00=pickup_alm(l,m,sp);
    B00=pickup_Blm(l,m,sp);          b00=pickup_blm(l,m,sp);
    App=pickup_Alm(l+1,m+1,sp);      app=pickup_alm(l+1,m+1,sp);
    Bpp=pickup_Blm(l+1,m+1,sp);      bpp=pickup_blm(l+1,m+1,sp);
    Apm=pickup_Alm(l+1,m-1,sp);      apm=pickup_alm(l+1,m-1,sp);
    Bpm=pickup_Blm(l+1,m-1,sp);      bpm=pickup_blm(l+1,m-1,sp);
    A0p=pickup_Alm(l,m+1,sp);        a0p=pickup_alm(l,m+1,sp);
    B0p=pickup_Blm(l,m+1,sp);        b0p=pickup_blm(l,m+1,sp);
    Ap0=pickup_Alm(l+1,m,sp);        ap0=pickup_alm(l+1,m,sp);
    Bp0=pickup_Blm(l+1,m,sp);        bp0=pickup_blm(l+1,m,sp);
    tfxy+=c0p*(cep*(2.0*a00*conj(app)+a00*conj(App)+A00*conj(app))+mu_0*(2.0*b00*conj(bpp)+b00*conj(Bpp)+B00*conj(bpp)))
      +   c0m*(cep*(2.0*apm*conj(a00)+apm*conj(A00)+Apm*conj(a00))+mu_0*(2.0*bpm*conj(b00)+bpm*conj(B00)+Bpm*conj(b00)))
      -   c1p*(-2.0*a00*conj(b0p)+2.0*b00*conj(a0p)-a00*conj(B0p)+b00*conj(A0p)+B00*conj(a0p)-A00*conj(b0p));
    tfz +=c2p*(cep*(2.0*ap0*conj(a00)+ap0*conj(A00)+Ap0*conj(a00))+mu_0*(2.0*bp0*conj(b00)+bp0*conj(B00)+Bp0*conj(b00)))
      +sqrt(cep*mu_0)*(double)m*(2.0*a00*conj(b00)+a00*conj(B00)+A00*conj(b00));
  }
  for(m=1;m<lm;m++){
    for(l=m;l<lm;l++){
      c0p=sqrt((double)((l+m+2)*(l+m+1))/(double)((2*l+1)*(2*l+3)))*(double)(l*(l+2));
      c0m=sqrt((double)((l-m+2)*(l-m+1))/(double)((2*l+1)*(2*l+3)))*(double)(l*(l+2));
      c1p=sqrt((double)((l+m+1)*(l-m)))*sqrt(cep*mu_0);
      c1m=sqrt((double)((l-m+1)*(l+m)))*sqrt(cep*mu_0);
      c2p=sqrt((double)((l-m+1)*(l+m+1))/(double)((2*l+3)*(2*l+1)))*(double)(l*(l+2));
      c2m=sqrt((double)((l+m+1)*(l-m+1))/(double)((2*l+3)*(2*l+1)))*(double)(l*(l+2));
      A00=pickup_Alm(l,m,sp);          a00=pickup_alm(l,m,sp);
      B00=pickup_Blm(l,m,sp);          b00=pickup_blm(l,m,sp);
      App=pickup_Alm(l+1,m+1,sp);      app=pickup_alm(l+1,m+1,sp);
      Bpp=pickup_Blm(l+1,m+1,sp);      bpp=pickup_blm(l+1,m+1,sp);
      Apm=pickup_Alm(l+1,m-1,sp);      apm=pickup_alm(l+1,m-1,sp);
      Bpm=pickup_Blm(l+1,m-1,sp);      bpm=pickup_blm(l+1,m-1,sp);
      A0p=pickup_Alm(l,m+1,sp);        a0p=pickup_alm(l,m+1,sp);
      B0p=pickup_Blm(l,m+1,sp);        b0p=pickup_blm(l,m+1,sp);
      Ap0=pickup_Alm(l+1,m,sp);        ap0=pickup_alm(l+1,m,sp);
      Bp0=pickup_Blm(l+1,m,sp);        bp0=pickup_blm(l+1,m,sp);
      tfxy+=c0p*(cep*(2.0*a00*conj(app)+a00*conj(App)+A00*conj(app))+mu_0*(2.0*b00*conj(bpp)+b00*conj(Bpp)+B00*conj(bpp)))
        +   c0m*(cep*(2.0*apm*conj(a00)+apm*conj(A00)+Apm*conj(a00))+mu_0*(2.0*bpm*conj(b00)+bpm*conj(B00)+Bpm*conj(b00)))
        -   c1p*(-2.0*a00*conj(b0p)+2.0*b00*conj(a0p)-a00*conj(B0p)+b00*conj(A0p)+B00*conj(a0p)-A00*conj(b0p));
      tfz +=c2p*(cep*(2.0*ap0*conj(a00)+ap0*conj(A00)+Ap0*conj(a00))+mu_0*(2.0*bp0*conj(b00)+bp0*conj(B00)+Bp0*conj(b00)))
        +sqrt(cep*mu_0)*(double)m*(2.0*a00*conj(b00)+a00*conj(B00)+A00*conj(b00));
      
      A00=pickup_Alm(l,-m,sp);          a00=pickup_alm(l,-m,sp);
      B00=pickup_Blm(l,-m,sp);          b00=pickup_blm(l,-m,sp);
      App=pickup_Alm(l+1,-m+1,sp);      app=pickup_alm(l+1,-m+1,sp);
      Bpp=pickup_Blm(l+1,-m+1,sp);      bpp=pickup_blm(l+1,-m+1,sp);
      Apm=pickup_Alm(l+1,-m-1,sp);      apm=pickup_alm(l+1,-m-1,sp);
      Bpm=pickup_Blm(l+1,-m-1,sp);      bpm=pickup_blm(l+1,-m-1,sp);
      A0p=pickup_Alm(l,-m+1,sp);        a0p=pickup_alm(l,-m+1,sp);
      B0p=pickup_Blm(l,-m+1,sp);        b0p=pickup_blm(l,-m+1,sp);
      Ap0=pickup_Alm(l+1,-m,sp);        ap0=pickup_alm(l+1,-m,sp);
      Bp0=pickup_Blm(l+1,-m,sp);        bp0=pickup_blm(l+1,-m,sp);
      tfxy+=c0m*(cep*(2.0*a00*conj(app)+a00*conj(App)+A00*conj(app))+mu_0*(2.0*b00*conj(bpp)+b00*conj(Bpp)+B00*conj(bpp)))
        +   c0p*(cep*(2.0*apm*conj(a00)+apm*conj(A00)+Apm*conj(a00))+mu_0*(2.0*bpm*conj(b00)+bpm*conj(B00)+Bpm*conj(b00)))
        -   c1m*(-2.0*a00*conj(b0p)+2.0*b00*conj(a0p)-a00*conj(B0p)+b00*conj(A0p)+B00*conj(a0p)-A00*conj(b0p));
      tfz +=c2m*(cep*(2.0*ap0*conj(a00)+ap0*conj(A00)+Ap0*conj(a00))+mu_0*(2.0*bp0*conj(b00)+bp0*conj(B00)+Bp0*conj(b00)))
        -sqrt(cep*mu_0)*(double)m*(2.0*a00*conj(b00)+a00*conj(B00)+A00*conj(b00));
    }
  }
  tfxy*=( 0.25*pow(ke,2)*I);
  tfz *=(-0.5 *pow(ke,2));
  f3[0]=creal(tfxy);
  f3[1]=cimag(tfxy);
  f3[2]=cimag(tfz );
}

void torque_sphr(double *n3,SPD *sp,Bobj *bm)
{
  double complex pickup_Alm(int l,int m,SPD *sp); // emf_mie_mmls.c
  double complex pickup_Blm(int l,int m,SPD *sp); // emf_mie_mmls.c
  double complex pickup_alm(int l,int m,SPD *sp); // emf_mie_mmls.c
  double complex pickup_blm(int l,int m,SPD *sp); // emf_mie_mmls.c
  
  double complex A00,B00,A0p,B0p;
  double complex a00,b00,a0p,b0p;
  double tnx=0.0,tny=0.0,tnz=0.0;
  double c0p,c0m;
  double cep=pow(bm->n_0,2)*epsilon0;
  double mu_0=mu0;
  double ke=2.0*M_PI*bm->n_0/bm->lambda_0;
  int lm=sp->l_limit;
  int l,m;
  
  m=0;
  for(l=1;l<lm;l++){
    c0p=sqrt((double)((l-m)*(l+m+1)))*(double)(l*(l+1));
    A00=pickup_Alm(l,m  ,sp);    a00=pickup_alm(l,m,sp);
    B00=pickup_Blm(l,m  ,sp);    b00=pickup_blm(l,m,sp);
    A0p=pickup_Alm(l,m+1,sp);    a0p=pickup_alm(l,m+1,sp);
    B0p=pickup_Blm(l,m+1,sp);    b0p=pickup_blm(l,m+1,sp);
    tnx+= creal(c0p*(cep*a00*conj(a0p)+mu_0*b00*conj(b0p)+0.5*(cep*(a00*conj(A0p)+a0p*conj(A00))+mu_0*(b00*conj(B0p)+b0p*conj(B00)))));
    tny+= cimag(c0p*(cep*a00*conj(a0p)+mu_0*b00*conj(b0p)+0.5*(cep*(a00*conj(A0p)-a0p*conj(A00))+mu_0*(b00*conj(B0p)-b0p*conj(B00)))));
    tnz+= (double)(m*(l+1)*l)*(cep*pow(cabs(a00),2)+mu_0*pow(cabs(b00),2)+creal(cep*a00*conj(A00)+mu_0*b00*conj(B00)));
  }
  for(m=1;m<lm;m++){
    for(l=m;l<lm;l++){
      c0p=sqrt((double)((l-m)*(l+m+1)))*(double)(l*(l+1));
      c0m=sqrt((double)((l+m)*(l-m+1)))*(double)(l*(l+1));
      A00=pickup_Alm(l,m  ,sp);    a00=pickup_alm(l,m,sp);
      B00=pickup_Blm(l,m  ,sp);    b00=pickup_blm(l,m,sp);
      A0p=pickup_Alm(l,m+1,sp);    a0p=pickup_alm(l,m+1,sp);
      B0p=pickup_Blm(l,m+1,sp);    b0p=pickup_blm(l,m+1,sp);
      tnx+= creal(c0p*(cep*a00*conj(a0p)+mu_0*b00*conj(b0p)+0.5*(cep*(a00*conj(A0p)+a0p*conj(A00))+mu_0*(b00*conj(B0p)+b0p*conj(B00)))));
      tny+= cimag(c0p*(cep*a00*conj(a0p)+mu_0*b00*conj(b0p)+0.5*(cep*(a00*conj(A0p)-a0p*conj(A00))+mu_0*(b00*conj(B0p)-b0p*conj(B00)))));
      tnz+= (double)(m*(l+1)*l)*(cep*pow(cabs(a00),2)+mu_0*pow(cabs(b00),2)+creal(cep*a00*conj(A00)+mu_0*b00*conj(B00)));
      
      A00=pickup_Alm(l,-m  ,sp);    a00=pickup_alm(l,-m,sp);
      B00=pickup_Blm(l,-m  ,sp);    b00=pickup_blm(l,-m,sp);
      A0p=pickup_Alm(l,-m+1,sp);    a0p=pickup_alm(l,-m+1,sp);
      B0p=pickup_Blm(l,-m+1,sp);    b0p=pickup_blm(l,-m+1,sp);
      tnx+= creal(c0m*(cep*a00*conj(a0p)+mu_0*b00*conj(b0p)+0.5*(cep*(a00*conj(A0p)+a0p*conj(A00))+mu_0*(b00*conj(B0p)+b0p*conj(B00)))));
      tny+= cimag(c0m*(cep*a00*conj(a0p)+mu_0*b00*conj(b0p)+0.5*(cep*(a00*conj(A0p)-a0p*conj(A00))+mu_0*(b00*conj(B0p)-b0p*conj(B00)))));
      tnz+=-(double)(m*(l+1)*l)*(cep*pow(cabs(a00),2)+mu_0*pow(cabs(b00),2)+creal(cep*a00*conj(A00)+mu_0*b00*conj(B00)));
    }
  }
  n3[0]=(-0.5*ke*tnx);
  n3[1]=(-0.5*ke*tny);
  n3[2]=(-0.5*ke*tnz);
}

void force_torque_sphr(double *f3,double *n3,SPD *sp,Bobj *bm)
{
  double complex pickup_Alm(int l,int m,SPD *sp); // emf_mie_mmls.c
  double complex pickup_Blm(int l,int m,SPD *sp); // emf_mie_mmls.c
  double complex pickup_alm(int l,int m,SPD *sp); // emf_mie_mmls.c
  double complex pickup_blm(int l,int m,SPD *sp); // emf_mie_mmls.c
  
  double complex A00,B00,App,Bpp,Apm,Bpm,A0p,B0p,Ap0,Bp0;
  double complex a00,b00,app,bpp,apm,bpm,a0p,b0p,ap0,bp0;
  double complex tfxy=0.0,tfz=0.0;
  double tnx=0.0,tny=0.0,tnz=0.0;
  double c0p,c0m,c1p,c1m,c2p,c2m,c3p,c3m;
  double cep=pow(bm->n_0,2)*epsilon0;
  double mu_0=mu0;
  double ke=2.0*M_PI*bm->n_0/bm->lambda_0;
  int lm=sp->l_limit;
  int l,m;
    
  m=0;
  for(l=1;l<lm;l++){
    c0p=sqrt((double)((l+m+2)*(l+m+1))/(double)((2*l+1)*(2*l+3)))*(double)(l*(l+2));
    c0m=sqrt((double)((l-m+2)*(l-m+1))/(double)((2*l+1)*(2*l+3)))*(double)(l*(l+2));
    c1p=sqrt((double)((l+m+1)*(l-m)))*sqrt(cep*mu_0);
    c1m=sqrt((double)((l-m+1)*(l+m)))*sqrt(cep*mu_0);
    c2p=sqrt((double)((l-m+1)*(l+m+1))/(double)((2*l+3)*(2*l+1)))*(double)(l*(l+2));
    c2m=sqrt((double)((l+m+1)*(l-m+1))/(double)((2*l+3)*(2*l+1)))*(double)(l*(l+2));
    c3p=sqrt((double)((l-m)*(l+m+1)))*(double)(l*(l+1));
    A00=pickup_Alm(l,m,sp);          a00=pickup_alm(l,m,sp);
    B00=pickup_Blm(l,m,sp);          b00=pickup_blm(l,m,sp);
    App=pickup_Alm(l+1,m+1,sp);      app=pickup_alm(l+1,m+1,sp);
    Bpp=pickup_Blm(l+1,m+1,sp);      bpp=pickup_blm(l+1,m+1,sp);
    Apm=pickup_Alm(l+1,m-1,sp);      apm=pickup_alm(l+1,m-1,sp);
    Bpm=pickup_Blm(l+1,m-1,sp);      bpm=pickup_blm(l+1,m-1,sp);
    A0p=pickup_Alm(l,m+1,sp);        a0p=pickup_alm(l,m+1,sp);
    B0p=pickup_Blm(l,m+1,sp);        b0p=pickup_blm(l,m+1,sp);
    Ap0=pickup_Alm(l+1,m,sp);        ap0=pickup_alm(l+1,m,sp);
    Bp0=pickup_Blm(l+1,m,sp);        bp0=pickup_blm(l+1,m,sp);
    tfxy+= c0p*(cep*(2.0*a00*conj(app)+a00*conj(App)+A00*conj(app))+mu_0*(2.0*b00*conj(bpp)+b00*conj(Bpp)+B00*conj(bpp)))
          +c0m*(cep*(2.0*apm*conj(a00)+apm*conj(A00)+Apm*conj(a00))+mu_0*(2.0*bpm*conj(b00)+bpm*conj(B00)+Bpm*conj(b00)))
          -c1p*(-2.0*a00*conj(b0p)+2.0*b00*conj(a0p)-a00*conj(B0p)+b00*conj(A0p)+B00*conj(a0p)-A00*conj(b0p));
    tfz += c2p*(cep*(2.0*ap0*conj(a00)+ap0*conj(A00)+Ap0*conj(a00))+mu_0*(2.0*bp0*conj(b00)+bp0*conj(B00)+Bp0*conj(b00)))
          +sqrt(cep*mu_0)*(double)m*(2.0*a00*conj(b00)+a00*conj(B00)+A00*conj(b00));
    tnx += creal(c3p*(cep*a00*conj(a0p)+mu_0*b00*conj(b0p)+0.5*(cep*(a00*conj(A0p)+a0p*conj(A00))+mu_0*(b00*conj(B0p)+b0p*conj(B00)))));
    tny += cimag(c3p*(cep*a00*conj(a0p)+mu_0*b00*conj(b0p)+0.5*(cep*(a00*conj(A0p)-a0p*conj(A00))+mu_0*(b00*conj(B0p)-b0p*conj(B00)))));
    tnz += (double)(m*(l+1)*l)*(cep*pow(cabs(a00),2)+mu_0*pow(cabs(b00),2)+creal(cep*a00*conj(A00)+mu_0*b00*conj(B00)));    
  }
  for(m=1;m<lm;m++){
    for(l=m;l<lm;l++){
      c0p=sqrt((double)((l+m+2)*(l+m+1))/(double)((2*l+1)*(2*l+3)))*(double)(l*(l+2));
      c0m=sqrt((double)((l-m+2)*(l-m+1))/(double)((2*l+1)*(2*l+3)))*(double)(l*(l+2));
      c1p=sqrt((double)((l+m+1)*(l-m)))*sqrt(cep*mu_0);
      c1m=sqrt((double)((l-m+1)*(l+m)))*sqrt(cep*mu_0);
      c2p=sqrt((double)((l-m+1)*(l+m+1))/(double)((2*l+3)*(2*l+1)))*(double)(l*(l+2));
      c2m=sqrt((double)((l+m+1)*(l-m+1))/(double)((2*l+3)*(2*l+1)))*(double)(l*(l+2));
      c3p=sqrt((double)((l-m)*(l+m+1)))*(double)(l*(l+1));
      c3m=sqrt((double)((l+m)*(l-m+1)))*(double)(l*(l+1));
      A00=pickup_Alm(l,m,sp);          a00=pickup_alm(l,m,sp);
      B00=pickup_Blm(l,m,sp);          b00=pickup_blm(l,m,sp);
      App=pickup_Alm(l+1,m+1,sp);      app=pickup_alm(l+1,m+1,sp);
      Bpp=pickup_Blm(l+1,m+1,sp);      bpp=pickup_blm(l+1,m+1,sp);
      Apm=pickup_Alm(l+1,m-1,sp);      apm=pickup_alm(l+1,m-1,sp);
      Bpm=pickup_Blm(l+1,m-1,sp);      bpm=pickup_blm(l+1,m-1,sp);
      A0p=pickup_Alm(l,m+1,sp);        a0p=pickup_alm(l,m+1,sp);
      B0p=pickup_Blm(l,m+1,sp);        b0p=pickup_blm(l,m+1,sp);
      Ap0=pickup_Alm(l+1,m,sp);        ap0=pickup_alm(l+1,m,sp);
      Bp0=pickup_Blm(l+1,m,sp);        bp0=pickup_blm(l+1,m,sp);
      tfxy+= c0p*(cep*(2.0*a00*conj(app)+a00*conj(App)+A00*conj(app))+mu_0*(2.0*b00*conj(bpp)+b00*conj(Bpp)+B00*conj(bpp)))
            +c0m*(cep*(2.0*apm*conj(a00)+apm*conj(A00)+Apm*conj(a00))+mu_0*(2.0*bpm*conj(b00)+bpm*conj(B00)+Bpm*conj(b00)))
            -c1p*(-2.0*a00*conj(b0p)+2.0*b00*conj(a0p)-a00*conj(B0p)+b00*conj(A0p)+B00*conj(a0p)-A00*conj(b0p));
      tfz += c2p*(cep*(2.0*ap0*conj(a00)+ap0*conj(A00)+Ap0*conj(a00))+mu_0*(2.0*bp0*conj(b00)+bp0*conj(B00)+Bp0*conj(b00)))
            +sqrt(cep*mu_0)*(double)m*(2.0*a00*conj(b00)+a00*conj(B00)+A00*conj(b00));
      tnx += creal(c3p*(cep*a00*conj(a0p)+mu_0*b00*conj(b0p)+0.5*(cep*(a00*conj(A0p)+a0p*conj(A00))+mu_0*(b00*conj(B0p)+b0p*conj(B00)))));
      tny += cimag(c3p*(cep*a00*conj(a0p)+mu_0*b00*conj(b0p)+0.5*(cep*(a00*conj(A0p)-a0p*conj(A00))+mu_0*(b00*conj(B0p)-b0p*conj(B00)))));
      tnz += (double)(m*(l+1)*l)*(cep*pow(cabs(a00),2)+mu_0*pow(cabs(b00),2)+creal(cep*a00*conj(A00)+mu_0*b00*conj(B00)));

      A00=pickup_Alm(l,-m,sp);          a00=pickup_alm(l,-m,sp);
      B00=pickup_Blm(l,-m,sp);          b00=pickup_blm(l,-m,sp);
      App=pickup_Alm(l+1,-m+1,sp);      app=pickup_alm(l+1,-m+1,sp);
      Bpp=pickup_Blm(l+1,-m+1,sp);      bpp=pickup_blm(l+1,-m+1,sp);
      Apm=pickup_Alm(l+1,-m-1,sp);      apm=pickup_alm(l+1,-m-1,sp);
      Bpm=pickup_Blm(l+1,-m-1,sp);      bpm=pickup_blm(l+1,-m-1,sp);
      A0p=pickup_Alm(l,-m+1,sp);        a0p=pickup_alm(l,-m+1,sp);
      B0p=pickup_Blm(l,-m+1,sp);        b0p=pickup_blm(l,-m+1,sp);
      Ap0=pickup_Alm(l+1,-m,sp);        ap0=pickup_alm(l+1,-m,sp);
      Bp0=pickup_Blm(l+1,-m,sp);        bp0=pickup_blm(l+1,-m,sp);
      tfxy+= c0m*(cep*(2.0*a00*conj(app)+a00*conj(App)+A00*conj(app))+mu_0*(2.0*b00*conj(bpp)+b00*conj(Bpp)+B00*conj(bpp)))
            +c0p*(cep*(2.0*apm*conj(a00)+apm*conj(A00)+Apm*conj(a00))+mu_0*(2.0*bpm*conj(b00)+bpm*conj(B00)+Bpm*conj(b00)))
            -c1m*(-2.0*a00*conj(b0p)+2.0*b00*conj(a0p)-a00*conj(B0p)+b00*conj(A0p)+B00*conj(a0p)-A00*conj(b0p));
      tfz += c2m*(cep*(2.0*ap0*conj(a00)+ap0*conj(A00)+Ap0*conj(a00))+mu_0*(2.0*bp0*conj(b00)+bp0*conj(B00)+Bp0*conj(b00)))
            -sqrt(cep*mu_0)*(double)m*(2.0*a00*conj(b00)+a00*conj(B00)+A00*conj(b00));
      tnx += creal(c3m*(cep*a00*conj(a0p)+mu_0*b00*conj(b0p)+0.5*(cep*(a00*conj(A0p)+a0p*conj(A00))+mu_0*(b00*conj(B0p)+b0p*conj(B00)))));
      tny += cimag(c3m*(cep*a00*conj(a0p)+mu_0*b00*conj(b0p)+0.5*(cep*(a00*conj(A0p)-a0p*conj(A00))+mu_0*(b00*conj(B0p)-b0p*conj(B00)))));
      tnz +=-(double)(m*(l+1)*l)*(cep*pow(cabs(a00),2)+mu_0*pow(cabs(b00),2)+creal(cep*a00*conj(A00)+mu_0*b00*conj(B00)));
    }
  }
  tfxy*=( 0.25*pow(ke,2)*I);
  tfz *=(-0.5 *pow(ke,2));
  f3[0]=creal(tfxy);
  f3[1]=cimag(tfxy);
  f3[2]=cimag(tfz );
  n3[0]=(-0.5*ke*tnx);
  n3[1]=(-0.5*ke*tny);
  n3[2]=(-0.5*ke*tnz);
}
