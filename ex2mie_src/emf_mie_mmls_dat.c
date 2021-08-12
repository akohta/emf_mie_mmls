#include "emf_mie_mmls.h"

void write_dat_mmls(char *fn,MSPD *msp)
{
  FILE *fp;
  int i,j,np,nt,mn,nl;
  
  if((fp=fopen(fn,"wb"))==NULL){    printf("write_dat_ms(), Failed to open the file %s. Exit...\n",fn);    exit(1);  }
  
  fwrite(msp,sizeof(MSPD),1,fp);
  // beam data
  fwrite(msp->bm.bd.ipw,sizeof(Ipw),msp->bm.n_ipw,fp);
  fwrite(msp->bm.bd.fpw,sizeof(Fpw),msp->bm.n_fpw,fp);
  fwrite(msp->bm.bd.lgb,sizeof(LGb),msp->bm.n_lgb,fp);
  fwrite(msp->bm.bd.bsb,sizeof(Bsb),msp->bm.n_bsb,fp);
  fwrite(msp->bm.bd.blg,sizeof(BsLGb),msp->bm.n_blg,fp);
  fwrite(msp->bm.bd.rab,sizeof(RAb),msp->bm.n_rab,fp);
  // sphere data
  fwrite(msp->sp,sizeof(SPD),msp->n_sphr,fp);
  for(i=0;i<msp->n_sphr;i++){
    np=msp->sp[i].ddt.np;
    nt=msp->sp[i].ddt.nt;
    mn=msp->sp[i].l_limit;
    nl=msp->sp[i].n_l;
    fwrite(msp->sp[i].a,sizeof(double),nl,fp);
    fwrite(msp->sp[i].ns,sizeof(double complex),nl,fp);
    fwrite(msp->sp[i].ddt.eri,sizeof(double complex),np*nt,fp);
    fwrite(msp->sp[i].ddt.hri,sizeof(double complex),np*nt,fp);
    fwrite(msp->sp[i].ddt.cab,sizeof(double),mn+1,fp);
    fwrite(msp->sp[i].ddt.ca,sizeof(double complex),mn+1,fp);
    fwrite(msp->sp[i].ddt.cb,sizeof(double complex),mn+1,fp);
    fwrite(msp->sp[i].ddt.ca0,sizeof(double complex),mn+1,fp);
    fwrite(msp->sp[i].ddt.cb0,sizeof(double complex),mn+1,fp);
    fwrite(msp->sp[i].ddt.ceM,sizeof(double complex),2*nl*(mn+1),fp);
    fwrite(msp->sp[i].ddt.chM,sizeof(double complex),2*nl*(mn+1),fp);
    for(j=0;j<=nl;j++) fwrite(msp->sp[i].ddt.Alm[j],sizeof(double complex),mn*(mn+2),fp);
    for(j=0;j<=nl;j++) fwrite(msp->sp[i].ddt.Blm[j],sizeof(double complex),mn*(mn+2),fp);
    for(j=0;j<=nl;j++) fwrite(msp->sp[i].ddt.alm[j],sizeof(double complex),mn*(mn+2),fp);
    for(j=0;j<=nl;j++) fwrite(msp->sp[i].ddt.blm[j],sizeof(double complex),mn*(mn+2),fp);
  }

  fclose(fp);
}

void  read_dat_mmls(char *fn,MSPD *msp)
{
  void *m_alloc2(size_t num,size_t size, char *txt); // emf_mie_ms.c
  void setup_sp(SPD *sp); // emf_mie_ms.c

  FILE *fp;
  int mn,i,j,nt,np,nl;
  if((fp=fopen(fn,"rb"))==NULL){    printf("read_dat_ms(), Failed to open the %s. Exit...\n",fn);    exit(1);  }
  
  fread(msp,sizeof(MSPD),1,fp);
  // beam data
  msp->bm.bd.ipw=(Ipw *)m_alloc2(msp->bm.n_ipw,sizeof(Ipw),"read_dat_ms(),msp->bm.bd.ipw");
  fread(msp->bm.bd.ipw,sizeof(Ipw),msp->bm.n_ipw,fp);
  msp->bm.bd.fpw=(Fpw *)m_alloc2(msp->bm.n_fpw,sizeof(Fpw),"read_dat_ms(),msp->bm.bd.fpw");
  fread(msp->bm.bd.fpw,sizeof(Fpw),msp->bm.n_fpw,fp);
  msp->bm.bd.lgb=(LGb *)m_alloc2(msp->bm.n_lgb,sizeof(LGb),"read_dat_ms(),msp->bm.bd.lgb");
  fread(msp->bm.bd.lgb,sizeof(LGb),msp->bm.n_lgb,fp);
  msp->bm.bd.bsb=(Bsb *)m_alloc2(msp->bm.n_bsb,sizeof(Bsb),"read_dat_ms(),msp->bm.bd.bsb");
  fread(msp->bm.bd.bsb,sizeof(Bsb),msp->bm.n_bsb,fp);
  msp->bm.bd.blg=(BsLGb *)m_alloc2(msp->bm.n_blg,sizeof(BsLGb),"read_dat_ms(),msp->bm.bd.blg");
  fread(msp->bm.bd.blg,sizeof(BsLGb),msp->bm.n_blg,fp);
  msp->bm.bd.rab=(RAb *)m_alloc2(msp->bm.n_rab,sizeof(RAb),"read_dat_ms(),msp->bm.bd.rab");
  fread(msp->bm.bd.rab,sizeof(RAb),msp->bm.n_rab,fp);
  setup_mfb(&(msp->bm));
  // sphere data
  msp->sp=(SPD *)m_alloc2(msp->n_sphr,sizeof(SPD),"read_dat_ms(),msp->sp");
  fread(msp->sp,sizeof(SPD),msp->n_sphr,fp);
  for(i=0;i<msp->n_sphr;i++){
    np=msp->sp[i].ddt.np;
    nt=msp->sp[i].ddt.nt;
    mn=msp->sp[i].l_limit;
    nl=msp->sp[i].n_l;
    msp->sp[i].a=(double *)m_alloc2(nl,sizeof(double),"read_dat_mmls(),msp->sp[i].a");
    msp->sp[i].ns=(double complex *)m_alloc2(nl,sizeof(double complex),"read_dat_mmls(),msp->sp[i].ns");
    setup_sp(&(msp->sp[i]));   
    fread(msp->sp[i].a,sizeof(double),nl,fp);
    fread(msp->sp[i].ns,sizeof(double complex),nl,fp);
    fread(msp->sp[i].ddt.eri,sizeof(double complex),np*nt,fp);
    fread(msp->sp[i].ddt.hri,sizeof(double complex),np*nt,fp);
    fread(msp->sp[i].ddt.cab,sizeof(double),mn+1,fp);
    fread(msp->sp[i].ddt.ca,sizeof(double complex),mn+1,fp);
    fread(msp->sp[i].ddt.cb,sizeof(double complex),mn+1,fp);
    fread(msp->sp[i].ddt.ca0,sizeof(double complex),mn+1,fp);
    fread(msp->sp[i].ddt.cb0,sizeof(double complex),mn+1,fp);
    fread(msp->sp[i].ddt.ceM,sizeof(double complex),2*nl*(mn+1),fp);
    fread(msp->sp[i].ddt.chM,sizeof(double complex),2*nl*(mn+1),fp);
    for(j=0;j<=nl;j++) fread(msp->sp[i].ddt.Alm[j],sizeof(double complex),mn*(mn+2),fp);
    for(j=0;j<=nl;j++) fread(msp->sp[i].ddt.Blm[j],sizeof(double complex),mn*(mn+2),fp);
    for(j=0;j<=nl;j++) fread(msp->sp[i].ddt.alm[j],sizeof(double complex),mn*(mn+2),fp);
    for(j=0;j<=nl;j++) fread(msp->sp[i].ddt.blm[j],sizeof(double complex),mn*(mn+2),fp);
  } 
  fclose(fp);
}
