#include "emf_mie_mmls.h"

void write_dat_mmls(char *fn,MSPD *msp)
{
  FILE *fp;
  int i,j,np,nt,mn,nl;
  
  if((fp=fopen(fn,"wb"))==NULL){    printf("write_dat_ms(), Failed to open the file %s. Exit...\n",fn);    exit(1);  }
  
  if(fwrite(msp,sizeof(MSPD),1,fp)!=1){
    printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the msp. exit...\n");
    exit(1);
  }
  // beam data
  if(fwrite(msp->bm.bd.ipw,sizeof(Ipw),msp->bm.n_ipw,fp)!=msp->bm.n_ipw){
    printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the ipw. exit...\n");
    exit(1);
  }
  if(fwrite(msp->bm.bd.fpw,sizeof(Fpw),msp->bm.n_fpw,fp)!=msp->bm.n_fpw){
    printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the fpw. exit...\n");
    exit(1);
  }
  if(fwrite(msp->bm.bd.lgb,sizeof(LGb),msp->bm.n_lgb,fp)!=msp->bm.n_lgb){
    printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the lgb. exit...\n");
    exit(1);
  }
  if(fwrite(msp->bm.bd.bsb,sizeof(Bsb),msp->bm.n_bsb,fp)!=msp->bm.n_bsb){
    printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the bsb. exit...\n");
    exit(1);
  }
  if(fwrite(msp->bm.bd.blg,sizeof(BsLGb),msp->bm.n_blg,fp)!=msp->bm.n_blg){
    printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the blg. exit...\n");
    exit(1);
  }
  if(fwrite(msp->bm.bd.rab,sizeof(RAb),msp->bm.n_rab,fp)!=msp->bm.n_rab){
    printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the rab. exit...\n");
    exit(1);
  }
  // sphere data
  if(fwrite(msp->sp,sizeof(SPD),msp->n_sphr,fp)!=msp->n_sphr){
    printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the sp. exit...\n");
    exit(1);
  }
  for(i=0;i<msp->n_sphr;i++){
    np=msp->sp[i].ddt.np;
    nt=msp->sp[i].ddt.nt;
    mn=msp->sp[i].l_limit;
    nl=msp->sp[i].n_l;
    if(fwrite(msp->sp[i].a,sizeof(double),nl,fp)!=nl){
      printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the a. exit...\n");
      exit(1);
    }
    if(fwrite(msp->sp[i].ns,sizeof(double complex),nl,fp)!=nl){
      printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the ns. exit...\n");
      exit(1);
    }
    if(fwrite(msp->sp[i].ddt.eri,sizeof(double complex),np*nt,fp)!=np*nt){
      printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the eri. exit...\n");
      exit(1);
    }
    if(fwrite(msp->sp[i].ddt.hri,sizeof(double complex),np*nt,fp)!=np*nt){
      printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the hri. exit...\n");
      exit(1);
    }
    if(fwrite(msp->sp[i].ddt.cab,sizeof(double),mn+1,fp)!=mn+1){
      printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the cab. exit...\n");
      exit(1);
    }
    if(fwrite(msp->sp[i].ddt.ca,sizeof(double complex),mn+1,fp)!=mn+1){
      printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the ca. exit...\n");
      exit(1);
    }
    if(fwrite(msp->sp[i].ddt.cb,sizeof(double complex),mn+1,fp)!=mn+1){
      printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the cb. exit...\n");
      exit(1);
    }
    if(fwrite(msp->sp[i].ddt.ca0,sizeof(double complex),mn+1,fp)!=mn+1){
      printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the ca0. exit...\n");
      exit(1);
    }
    if(fwrite(msp->sp[i].ddt.cb0,sizeof(double complex),mn+1,fp)!=mn+1){
      printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the cb0. exit...\n");
      exit(1);
    }
    if(fwrite(msp->sp[i].ddt.ceM,sizeof(double complex),2*nl*(mn+1),fp)!=2*nl*(mn+1)){
      printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the ceM. exit...\n");
      exit(1);
    }
    if(fwrite(msp->sp[i].ddt.chM,sizeof(double complex),2*nl*(mn+1),fp)!=2*nl*(mn+1)){
      printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the chM. exit...\n");
      exit(1);
    }
    for(j=0;j<=nl;j++){
      if(fwrite(msp->sp[i].ddt.Alm[j],sizeof(double complex),mn*(mn+2),fp)!=mn*(mn+2)){
        printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the Alm[j]. exit...\n");
        exit(1);
      }
    }
    for(j=0;j<=nl;j++){
      if(fwrite(msp->sp[i].ddt.Blm[j],sizeof(double complex),mn*(mn+2),fp)!=mn*(mn+2)){
        printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the Blm[j]. exit...\n");
        exit(1);
      }
    }
    for(j=0;j<=nl;j++){
      if(fwrite(msp->sp[i].ddt.alm[j],sizeof(double complex),mn*(mn+2),fp)!=mn*(mn+2)){
        printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the alm[j]. exit...\n");
        exit(1);
      }
    }
    for(j=0;j<=nl;j++){
      if(fwrite(msp->sp[i].ddt.blm[j],sizeof(double complex),mn*(mn+2),fp)!=mn*(mn+2)){
        printf("emf_mie_mmls_dat.c, write_dat_mmls(), failed to write the blm[j]. exit...\n");
        exit(1);
      }
    }
  }

  fclose(fp);
}

void  read_dat_mmls(char *fn,MSPD *msp)
{
  void setup_sp(SPD *sp); // emf_mie_ms.c

  FILE *fp;
  int mn,i,j,nt,np,nl;
  if((fp=fopen(fn,"rb"))==NULL){    printf("read_dat_ms(), Failed to open the %s. Exit...\n",fn);    exit(1);  }
  
  if(fread(msp,sizeof(MSPD),1,fp)!=1){
    printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the msp. exit...\n");
    exit(1);
  }
  // beam data
  msp->bm.bd.ipw=(Ipw *)m_alloc2(msp->bm.n_ipw,sizeof(Ipw),"read_dat_ms(),msp->bm.bd.ipw");
  if(fread(msp->bm.bd.ipw,sizeof(Ipw),msp->bm.n_ipw,fp)!=msp->bm.n_ipw){
    printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the ipw. exit...\n");
    exit(1);
  }
  msp->bm.bd.fpw=(Fpw *)m_alloc2(msp->bm.n_fpw,sizeof(Fpw),"read_dat_ms(),msp->bm.bd.fpw");
  if(fread(msp->bm.bd.fpw,sizeof(Fpw),msp->bm.n_fpw,fp)!=msp->bm.n_fpw){
    printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the fpw. exit...\n");
    exit(1);
  }
  msp->bm.bd.lgb=(LGb *)m_alloc2(msp->bm.n_lgb,sizeof(LGb),"read_dat_ms(),msp->bm.bd.lgb");
  if(fread(msp->bm.bd.lgb,sizeof(LGb),msp->bm.n_lgb,fp)!=msp->bm.n_lgb){
    printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the lgb. exit...\n");
    exit(1);
  }
  msp->bm.bd.bsb=(Bsb *)m_alloc2(msp->bm.n_bsb,sizeof(Bsb),"read_dat_ms(),msp->bm.bd.bsb");
  if(fread(msp->bm.bd.bsb,sizeof(Bsb),msp->bm.n_bsb,fp)!=msp->bm.n_bsb){
    printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the bsb. exit...\n");
    exit(1);
  }
  msp->bm.bd.blg=(BsLGb *)m_alloc2(msp->bm.n_blg,sizeof(BsLGb),"read_dat_ms(),msp->bm.bd.blg");
  if(fread(msp->bm.bd.blg,sizeof(BsLGb),msp->bm.n_blg,fp)!=msp->bm.n_blg){
    printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the blg. exit...\n");
    exit(1);
  }
  msp->bm.bd.rab=(RAb *)m_alloc2(msp->bm.n_rab,sizeof(RAb),"read_dat_ms(),msp->bm.bd.rab");
  if(fread(msp->bm.bd.rab,sizeof(RAb),msp->bm.n_rab,fp)!=msp->bm.n_rab){
    printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the rab. exit...\n");
    exit(1);
  }
  setup_mfb(&(msp->bm));
  // sphere data
  msp->sp=(SPD *)m_alloc2(msp->n_sphr,sizeof(SPD),"read_dat_ms(),msp->sp");
  if(fread(msp->sp,sizeof(SPD),msp->n_sphr,fp)!=msp->n_sphr){
    printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the sp. exit...\n");
    exit(1);
  }
  for(i=0;i<msp->n_sphr;i++){
    np=msp->sp[i].ddt.np;
    nt=msp->sp[i].ddt.nt;
    mn=msp->sp[i].l_limit;
    nl=msp->sp[i].n_l;
    msp->sp[i].a=(double *)m_alloc2(nl,sizeof(double),"read_dat_mmls(),msp->sp[i].a");
    msp->sp[i].ns=(double complex *)m_alloc2(nl,sizeof(double complex),"read_dat_mmls(),msp->sp[i].ns");
    setup_sp(&(msp->sp[i]));   
    if(fread(msp->sp[i].a,sizeof(double),nl,fp)!=nl){
      printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the a. exit...\n");
      exit(1);
    }
    if(fread(msp->sp[i].ns,sizeof(double complex),nl,fp)!=nl){
      printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the ns. exit...\n");
      exit(1);
    }
    if(fread(msp->sp[i].ddt.eri,sizeof(double complex),np*nt,fp)!=np*nt){
      printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the eri. exit...\n");
      exit(1);
    }
    if(fread(msp->sp[i].ddt.hri,sizeof(double complex),np*nt,fp)!=np*nt){
      printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the hri. exit...\n");
      exit(1);
    }
    if(fread(msp->sp[i].ddt.cab,sizeof(double),mn+1,fp)!=mn+1){
      printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the cab. exit...\n");
      exit(1);
    }
    if(fread(msp->sp[i].ddt.ca,sizeof(double complex),mn+1,fp)!=mn+1){
      printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the ca. exit...\n");
      exit(1);
    }
    if(fread(msp->sp[i].ddt.cb,sizeof(double complex),mn+1,fp)!=mn+1){
      printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the cb. exit...\n");
      exit(1);
    }
    if(fread(msp->sp[i].ddt.ca0,sizeof(double complex),mn+1,fp)!=mn+1){
      printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the ca0. exit...\n");
      exit(1);
    }
    if(fread(msp->sp[i].ddt.cb0,sizeof(double complex),mn+1,fp)!=mn+1){
      printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the cb0. exit...\n");
      exit(1);
    }
    if(fread(msp->sp[i].ddt.ceM,sizeof(double complex),2*nl*(mn+1),fp)!=2*nl*(mn+1)){
      printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the ceM. exit...\n");
      exit(1);
    }
    if(fread(msp->sp[i].ddt.chM,sizeof(double complex),2*nl*(mn+1),fp)!=2*nl*(mn+1)){
      printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the chM. exit...\n");
      exit(1);
    }
    for(j=0;j<=nl;j++){
      if(fread(msp->sp[i].ddt.Alm[j],sizeof(double complex),mn*(mn+2),fp)!=mn*(mn+2)){
        printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the Alm[j]. exit...\n");
        exit(1);
      }
    }
    for(j=0;j<=nl;j++){
      if(fread(msp->sp[i].ddt.Blm[j],sizeof(double complex),mn*(mn+2),fp)!=mn*(mn+2)){
        printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the Blm[j]. exit...\n");
        exit(1);
      }
    }
    for(j=0;j<=nl;j++){
      if(fread(msp->sp[i].ddt.alm[j],sizeof(double complex),mn*(mn+2),fp)!=mn*(mn+2)){
        printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the alm[j]. exit...\n");
        exit(1);
      }
    }
    for(j=0;j<=nl;j++){
      if(fread(msp->sp[i].ddt.blm[j],sizeof(double complex),mn*(mn+2),fp)!=mn*(mn+2)){
        printf("emf_mie_mmls_dat.c, read_dat_mmls(), failed to read the blm[j]. exit...\n");
        exit(1);
      }
    }
  } 
  fclose(fp);
}
