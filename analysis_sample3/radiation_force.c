// calculation example of radiation force. 
// trapping force and trappable point can be analyzed ( single sphere only )
#include "emf_mie_mmls.h"

void set_coef(MSPD *msp);

int main()
{
  MSPD msp;
  FILE *ff;
  double rang,dr,r[3],vf[3];
  int max,i,j;
  
  read_data_mmls(&msp);
  print_data_mmls(&msp);
  setup_mmls(&msp);
  
  max=15;
  rang=1.0*msp.bm.lambda_0;
  dr=rang*2/(double)(max-1);
  
  // x=0 plane 
  if((ff=fopen("force_yz.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(ff,"%s\n","# y z Fx[N] Fy[N] Fz[N]");

  r[0]=0.0;
  for(i=0;i<max;i++){
    r[1]=-rang+(double)i*dr;
    for(j=0;j<max;j++){
      r[2]=-rang+(double)j*dr;
      // set shpere position
      msp.sp[0].xs=r[0];      msp.sp[0].ys=r[1];      msp.sp[0].zs=r[2];
      set_coef(&msp);
      // force 
      force_mmls(0,vf,&msp);
      fprintf(ff,"%g %g %g %g %g\n",r[1],r[2],OSUtoMKSA_Force(vf[0]),OSUtoMKSA_Force(vf[1]),OSUtoMKSA_Force(vf[2]));
     }
    fprintf(ff,"\n");

  }
  fclose(ff);
  printf("x=0 plane finished\n"); fflush(stdout);
  
  
  // y=0 plane 
  if((ff=fopen("force_xz.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(ff,"%s\n","# x z Fx[N] Fy[N] Fz[N]");

  r[1]=0.0;
  for(i=0;i<max;i++){
    r[0]=-rang+(double)i*dr;
    for(j=0;j<max;j++){
      r[2]=-rang+(double)j*dr;
      // set shpere position
      msp.sp[0].xs=r[0];      msp.sp[0].ys=r[1];      msp.sp[0].zs=r[2];
      set_coef(&msp);
      // force 
      force_mmls(0,vf,&msp);
      fprintf(ff,"%g %g %g %g %g\n",r[0],r[2],OSUtoMKSA_Force(vf[0]),OSUtoMKSA_Force(vf[1]),OSUtoMKSA_Force(vf[2]));
     }
    fprintf(ff,"\n");
  }
  fclose(ff);
  printf("y=0 plane finished\n"); fflush(stdout);

  
  // z=0 plane
  if((ff=fopen("force_xy.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(ff,"%s\n","# x y Fx[N] Fy[N] Fz[N]");

  r[2]=0.0;
  for(i=0;i<max;i++){
    r[0]=-rang+(double)i*dr;
    for(j=0;j<max;j++){
      r[1]=-rang+(double)j*dr;
      // set shpere position
      msp.sp[0].xs=r[0];      msp.sp[0].ys=r[1];      msp.sp[0].zs=r[2];
      set_coef(&msp);
      // force 
      force_mmls(0,vf,&msp);
      fprintf(ff,"%15.14g %15.14g %15.14g %15.14g %15.14g\n",r[0],r[1],OSUtoMKSA_Force(vf[0]),OSUtoMKSA_Force(vf[1]),OSUtoMKSA_Force(vf[2]));
    }
    fprintf(ff,"\n");
  }
  fclose(ff);
  printf("analysis of radiation force is finished\n");
  
  free_mmls(&msp);
  return 0;
}

void set_coef(MSPD *msp)
{
  void initialize_eh_r(SPD *sp,Bobj *bm);       // emf_mie_ms.c
  void coefficient(SPD *sp);                    // emf_mie_ms.c
  void characteristic_matrix(SPD *sp,Bobj *bm); // emf_mie_ms.c
  
  initialize_eh_r(&(msp->sp[0]),&(msp->bm));
  characteristic_matrix(&(msp->sp[0]),&(msp->bm)); 
  coefficient(&(msp->sp[0]));
}
  
