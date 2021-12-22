#include "emf_mie_mmls.h"

void read_data_mmls(MSPD *msp)
{
  FILE *fp;
  char buf[256]="";
  double tmpd,tmpd2;
  int s,i,tmpi;
  
  if((fp=fopen(fn_mlsphr,"rt"))==NULL){    printf("Can not open the '%s' file. Exit...\n",fn_mlsphr);    exit(1);  }
  if(fgets(buf,256,fp)==NULL){
    printf("emf_mie_mmls.c, read_data_mmls(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("emf_mie_mmls.c, read_data_mmls(), failed to read the line. exit...\n");
    exit(1);
  }
  
  if(fscanf(fp,"%d\n",&tmpi)!=1){
    printf("emf_mie_mmls.c, read_data_mmls(), failed to read the n_sphr. exit...\n");
    exit(1);
  }
  msp->n_sphr=tmpi;
  if(msp->n_sphr==0) {
    printf("Sphere number is 0. Exit...\n");
    exit(0);
  }

  msp->sp=(SPD *)m_alloc2(tmpi,sizeof(SPD),"read_data_mmls(),msp->sp"); // malloc

  for(s=0;s<msp->n_sphr;s++){
    if(fgets(buf,256,fp)==NULL){
      printf("emf_mie_mmls.c, read_data_mmls(), failed to read the line. exit...\n");
      exit(1);
    }

    if(fscanf(fp,"%d",&tmpi)!=1){
      printf("emf_mie_mmls.c, read_data_mmls(), failed to read the bsn. exit...\n");
      exit(1);
    }
    msp->sp[s].bsn      =tmpi;    
    if(fscanf(fp,"%d",&tmpi)!=1){
      printf("emf_mie_mmls.c, read_data_mmls(), failed to read the bdv. exit...\n");
      exit(1);
    }
    msp->sp[s].bdv      =tmpi;
    if(fscanf(fp,"%d",&tmpi)!=1){
      printf("emf_mie_mmls.c, read_data_mmls(), failed to read the l_limit. exit...\n");
      exit(1);
    }
    msp->sp[s].l_limit  =tmpi; 
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("emf_mie_mmls.c, read_data_mmls(), failed to read the xs. exit...\n");
      exit(1);
    }
    msp->sp[s].xs       =tmpd;
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("emf_mie_mmls.c, read_data_mmls(), failed to read the ys. exit...\n");
      exit(1);
    }
    msp->sp[s].ys       =tmpd; 
    if(fscanf(fp,"%lf\n",&tmpd)!=1){
      printf("emf_mie_mmls.c, read_data_mmls(), failed to read the zs. exit...\n");
      exit(1);
    }
    msp->sp[s].zs       =tmpd;
  
    if(fgets(buf,256,fp)==NULL){
      printf("emf_mie_mmls.c, read_data_mmls(), failed to read the line. exit...\n");
      exit(1);
    }
    if(fscanf(fp,"%d\n",&tmpi)!=1){
      printf("emf_mie_mmls.c, read_data_mmls(), failed to read the n_l. exit...\n");
      exit(1);
    }
    msp->sp[s].n_l      =tmpi;
    if(tmpi<1){
      printf("The number of layers must be positive integer. Exit...\n");
      exit(1);
    }
    else if(tmpi==1){
      printf("The number of layers must be greater than or equal to 2.\n");
      printf("For non-layered sphere, set 2 layers as the same refractive index. Exit...\n");
      exit(1);
    }
    if(fgets(buf,256,fp)==NULL){
      printf("emf_mie_mmls.c, read_data_mmls(), failed to read the line. exit...\n");
      exit(1);
    }

    msp->sp[s].a=(double *)m_alloc2(tmpi,sizeof(double),"read_data_mmls(),msp->sp[s].a"); // malloc
    msp->sp[s].ns=(double complex *)m_alloc2(tmpi,sizeof(double complex),"read_data_mmls(),msp->sp[s].ns"); // malloc  
    for(i=0;i<tmpi;i++){
      if(fscanf(fp,"%lf",&tmpd)!=1){
        printf("emf_mie_mmls.c, read_data_mmls(), failed to read the a. exit...\n");
        exit(1);
      }
      msp->sp[s].a [i]=tmpd;
      if(fscanf(fp,"%lf",&tmpd)!=1){
        printf("emf_mie_mmls.c, read_data_mmls(), failed to read the real(ns). exit...\n");
        exit(1);
      }
      if(fscanf(fp,"%lf\n",&tmpd2)!=1){
        printf("emf_mie_mmls.c, read_data_mmls(), failed to read the imag(ns). exit...\n");
        exit(1);
      }
      msp->sp[s].ns[i]=tmpd+I*tmpd2;  
      if(i>0){
        if(msp->sp[s].a[i-1]>msp->sp[s].a[i]){
          printf("layer data must be defined in order from the inside\n");
          printf("sphere id = %d. a[%d]=%g,a[%d]=%g. Exit...\n",s,i-1,msp->sp[s].a[i-1],i,msp->sp[s].a[i]);
          exit(1);
        }
        if(msp->sp[s].a[i-1]==msp->sp[s].a[i]){
          printf("layer is overlapped.\n");
          printf("sphere id = %d. a[%d]=%g,a[%d]=%g. Exit...\n",s,i-1,msp->sp[s].a[i-1],i,msp->sp[s].a[i]);
          exit(1);
        }
      }
    }
  }
  fclose(fp);
  
  // multi fbeam
  init_mfb(&(msp->bm));       // initialize
  read_data_mfb(&(msp->bm));  // search and read beam datafile
}

void print_data_mmls(MSPD *msp)
{
  int s,i;
  
  print_data_mfb(&(msp->bm)); // print beam data
  
  printf("--------- sphere data ---------\n");
  for(s=0;s<msp->n_sphr;s++){
    printf("sphere id %d\n",s);
    printf("  basic sampling number on sphere surface       : %16d\n",msp->sp[s].bsn);
    printf("  division number of sphere surface     (per PI): %16d\n",msp->sp[s].bdv);
    printf("  limit of order number l                       : %16d\n",msp->sp[s].l_limit);
    printf("  x-coordinate of sphere center                 : %16.15g\n",msp->sp[s].xs);
    printf("  y-coordinate of sphere center                 : %16.15g\n",msp->sp[s].ys);
    printf("  z-coordinate of sphere center                 : %16.15g\n",msp->sp[s].zs);
    for(i=0;i<msp->sp[s].n_l;i++){
      printf("  layer id %d\n",i);
      printf("    radius of layer                             : %16.15g\n",msp->sp[s].a[i]);
      printf("    refractive index of layer                   : %7.6g+%7.6gI\n",creal(msp->sp[s].ns[i]),cimag(msp->sp[s].ns[i]));    }
  }
  printf("\n");
}

void print_data_mmls_mksa(MSPD *msp)
{
  int s,i;
  
  print_data_mfb_mksa(&(msp->bm)); // print beam data
  
  printf("--------- sphere data, MKSA system ---------\n");
  for(s=0;s<msp->n_sphr;s++){
    printf("sphere id %d\n",s);
    printf("  basic sampling number on sphere surface       : %16d\n",msp->sp[s].bsn);
    printf("  division number of sphere surface     (per PI): %16d\n",msp->sp[s].bdv);
    printf("  limit of order number l                       : %16d\n",msp->sp[s].l_limit);
    printf("  x-coordinate of sphere center              [m]: %16.15g\n",OSUtoMKSA_length(msp->sp[s].xs));
    printf("  y-coordinate of sphere center              [m]: %16.15g\n",OSUtoMKSA_length(msp->sp[s].ys));
    printf("  z-coordinate of sphere center              [m]: %16.15g\n",OSUtoMKSA_length(msp->sp[s].zs));
    for(i=0;i<msp->sp[s].n_l;i++){
      printf("  layer id %d\n",i);
      printf("    radius of layer                          [m]: %16.15g\n",OSUtoMKSA_length(msp->sp[s].a[i]));
      printf("    refractive index of layer                   : %7.6g+%7.6gI\n",creal(msp->sp[s].ns[i]),cimag(msp->sp[s].ns[i]));    }
  }
  printf("\n");
}

void setup_mmls(MSPD *msp)
{
  void check_data(MSPD *msp);
  void setup_sp(SPD *sp);
  void setup_cab(SPD *sp,Bobj *bm);
  void initialize_eh_r(SPD *sp,Bobj *bm);
  void characteristic_matrix(SPD *sp,Bobj *bm);
  void coefficient(SPD *sp,Bobj *bm);
   
  int i;
  
  // check sphere data
  check_data(msp);

  // multi_fbeam
  setup_mfb(&(msp->bm));

  // spheres
  for(i=0;i<msp->n_sphr;i++){
    setup_sp(&(msp->sp[i]));
    setup_cab(&(msp->sp[i]),&(msp->bm));
    initialize_eh_r(&(msp->sp[i]),&(msp->bm));
    characteristic_matrix(&(msp->sp[i]),&(msp->bm));   
    coefficient(&(msp->sp[i]),&(msp->bm));
  }

}

void  free_mmls(MSPD *msp)
{
  void free_sp(SPD *sp);
  
  int i;

  // spheres  
  for(i=0;i<msp->n_sphr;i++){
    free_sp(&(msp->sp[i]));
  }
  free(msp->sp);  msp->n_sphr=0;
  
  // multi_fbeam
  free_mfb(&(msp->bm));
}

void iterative_ops_mmls(MSPD *msp)
{
  void coefficient(SPD *sp,Bobj *bm);
  void all_coefficient(SPD *sp,Bobj *bm);
  void field_s_ehr(int src,int obj,MSPD *msp);
  
  int i,j,t,nn,sbc,num,s,*bc;
  double f1,f3[3],f_delta,*f0;
  
  num=msp->n_sphr;
  
  if(num<2){
    all_coefficient(&(msp->sp[0]),&(msp->bm));
    return;
  }
  
  bc=(int *)m_alloc2(num,sizeof(int),"iterative_ops_mmls(),bc");
  f0=(double *)m_alloc2(num,sizeof(double),"iterative_ops_mmls(),f0");

  for(t=0;t<num;t++){
    force_mmls(t,f3,msp);
    f0[t]=f3[0]*f3[0]+f3[1]*f3[1]+f3[2]*f3[2];
    bc[t]=ito_breakcount;
  }
  if(num>1){
    printf("iterative operation start (convergence criterion : cv < %g)\n",ito_eps);
    for(nn=0;nn<ito_max;nn++){
      for(i=0;i<num;i++)
        for(j=0;j<num;j++)
          if(i!=j)      field_s_ehr(i,j,msp);
      for(i=0;i<num;i++) coefficient(&(msp->sp[i]),&(msp->bm));
      printf("%3d, cv : ",nn);
      for(t=0;t<num;t++){
        force_mmls(t,f3,msp);
        f1=f3[0]*f3[0]+f3[1]*f3[1]+f3[2]*f3[2];
        f_delta=fabs(f1/f0[t]-1.0);
        if(f_delta<ito_eps)  bc[t]--;
        printf("%g\t",f_delta);
        f0[t]=f1;
      }
      printf("\n");
      sbc=0;
      for(t=0;t<num;t++) if(bc[t]<=0) sbc++;
      if(sbc==num) break;
    }
    
    if(nn==ito_max){
      printf("The maximum number of iterations has been reached (The result has not converged).\n");
    }
  }
  for(s=0;s<num;s++)    all_coefficient(&(msp->sp[s]),&(msp->bm));
  printf("\n");

  free(bc);  free(f0);
}

void output_node_particles(char *fname,MSPD *msp)
{
  FILE *fp;
  double a,st,ct,sp,cp,x,y,z;
  int s1,oid,i,j;
  char *sd,fo[256]={},tf[200]={};

  s1=strlen(fname);
  if(s1>200){
    printf("emf_mie_mmls.c, output_node_particles(), file name is too long. exit...\n");
    exit(1);
  }
  sprintf(fo,"%s",fname);
  sd=strrchr(fo,'.');
  if(sd!=NULL){
    strncpy(tf,fname,s1-strlen(sd));
    sprintf(fo,"%s.particles",tf);
  }
  
  if((fp=fopen(fo,"wt"))==NULL){    printf("Can not open the %s file.\n",fo);    exit(1);  }
  fprintf(fp,"# x y z object_id\n");
  
  for(oid=0;oid<msp->n_sphr;oid++){
    a=msp->sp[oid].a[msp->sp[oid].n_l-1];
    for(i=0;i<msp->sp[oid].ddt.nt;i++){
      st=sin(msp->sp[oid].ddt.xt[i]);
      ct=cos(msp->sp[oid].ddt.xt[i]);
      for(j=0;j<msp->sp[oid].ddt.np;j++){
        sp=sin(msp->sp[oid].ddt.xp[j]);
        cp=cos(msp->sp[oid].ddt.xp[j]);

        x=a*st*cp+msp->sp[oid].xs;
        y=a*st*sp+msp->sp[oid].ys;
        z=a*ct   +msp->sp[oid].zs;
        fprintf(fp,"%15.14e %15.14e %15.14e %d\n",x,y,z,oid);
      }
    }
  }
  
  fclose(fp);
}

/////////////////////////////////////////////////////////////////////
void check_data(MSPD *msp)
{
  double r,rs;
  int i,j,s;
  
  // position check
  for(i=0;i<msp->n_sphr;i++){
    for(j=i+1;j<msp->n_sphr;j++){
      r =msp->sp[i].a[msp->sp[i].n_l-1]+msp->sp[j].a[msp->sp[j].n_l-1];
      rs=sqrt(pow(msp->sp[j].xs-msp->sp[i].xs,2)+pow(msp->sp[j].ys-msp->sp[i].ys,2)+pow(msp->sp[j].zs-msp->sp[i].zs,2));
      if(rs<r){
        printf("Sphere Position Check Error! sphere id=%d and sphere id=%d is overlaped. Exit...\n",i,j);
        exit(1);
      }
    }
  }
  // sphere data check
  for(s=0;s<msp->n_sphr;s++){
    if(msp->sp[s].n_l<1){
      printf("The number of layers must be positive integer.\n");
      printf("sphere id=%d. number of layers=%d. Exit...\n",s,msp->sp[s].n_l);
      exit(1);
    }
    else if(msp->sp[s].n_l==1){
      printf("The number of layers must be greater than or equal to 2.\n");
      printf("For non-layered sphere, set 2 layers as the same refractive index. Exit...\n");
      exit(1);
    }
    for(i=1;i<msp->sp[s].n_l;i++){
      if(msp->sp[s].a[i-1]>msp->sp[s].a[i]){
        printf("layer data must be defined in order from the inside\n");
        printf("sphere id = %d. a[%d]=%g,a[%d]=%g. Exit...\n",s,i-1,msp->sp[s].a[i-1],i,msp->sp[s].a[i]);
        exit(1);
      }
      if(msp->sp[s].a[i-1]==msp->sp[s].a[i]){
        printf("layer is overlapped.\n");
        printf("sphere id = %d. a[%d]=%g,a[%d]=%g. Exit...\n",s,i-1,msp->sp[s].a[i-1],i,msp->sp[s].a[i]);
        exit(1);
      }
    }
  }
}

void setup_sp(SPD *sp)
{
  void gauleg_dv(double a,double b,double *x,double *w,int bn,int dv);
  
  int i,mn,nl,nt,np;
  
  // gauleg
  nt=  sp->bsn*sp->bdv;
  np=2*sp->bsn*sp->bdv;
  sp->ddt.nt=nt;
  sp->ddt.np=np;
  sp->ddt.xp=(double *)m_alloc2(np,sizeof(double),"setup_sp(),sp->ddt,xp");
  sp->ddt.wp=(double *)m_alloc2(np,sizeof(double),"setup_sp(),sp->ddt.wp");
  sp->ddt.xt=(double *)m_alloc2(nt,sizeof(double),"setup_sp(),sp->ddt.xt");
  sp->ddt.wt=(double *)m_alloc2(nt,sizeof(double),"setup_sp(),sp->ddt.wt");
  gauleg_dv(0.0,    M_PI,sp->ddt.xt,sp->ddt.wt,sp->bsn,  sp->bdv);
  gauleg_dv(0.0,2.0*M_PI,sp->ddt.xp,sp->ddt.wp,sp->bsn,2*sp->bdv);
  
  sp->ddt.eri=(double complex *)m_alloc2(np*nt,sizeof(double complex),"setup_sp(),sp->ddt.eri");
  sp->ddt.ers=(double complex *)m_alloc2(np*nt,sizeof(double complex),"setup_sp(),sp->ddt.ers");
  sp->ddt.hri=(double complex *)m_alloc2(np*nt,sizeof(double complex),"setup_sp(),sp->ddt.hri");
  sp->ddt.hrs=(double complex *)m_alloc2(np*nt,sizeof(double complex),"setup_sp(),sp->ddt.hrs");
  // sphere
  mn=sp->l_limit;
  sp->ddt.l_max=mn;
  nl=sp->n_l;
  sp->ddt.cab=(double *)m_alloc2(mn+1,sizeof(double),"setup_sp(),sp->ddt.cab");
  sp->ddt.ca =(double complex *)m_alloc2(mn+1,sizeof(double complex),"setup_sp(),sp->ddt.ca");
  sp->ddt.cb =(double complex *)m_alloc2(mn+1,sizeof(double complex),"setup_sp(),sp->ddt.cb");
  sp->ddt.ca0=(double complex *)m_alloc2(mn+1,sizeof(double complex),"setup_sp(),sp->ddt.ca0");
  sp->ddt.cb0=(double complex *)m_alloc2(mn+1,sizeof(double complex),"setup_sp(),sp->ddt.cb0");
  sp->ddt.ceM=(double complex *)m_alloc2(2*nl*(mn+1),sizeof(double complex),"setup_sp(),sp->ddt.ceM");
  sp->ddt.chM=(double complex *)m_alloc2(2*nl*(mn+1),sizeof(double complex),"setup_sp(),sp->ddt.chM");
  sp->ddt.Alm =(double complex **)m_alloc2(nl+1,sizeof(double complex *),"setup_sp(),sp->ddt.Alm");
  sp->ddt.Blm =(double complex **)m_alloc2(nl+1,sizeof(double complex *),"setup_sp(),sp->ddt.Blm");
  sp->ddt.alm =(double complex **)m_alloc2(nl+1,sizeof(double complex *),"setup_sp(),sp->ddt.alm");
  sp->ddt.blm =(double complex **)m_alloc2(nl+1,sizeof(double complex *),"setup_sp(),sp->ddt.blm");
  for(i=0;i<=nl;i++){
    sp->ddt.Alm [i] =(double complex *)m_alloc2(mn*(mn+2),sizeof(double complex),"setup_sp(),sp->ddt.Alm[i]");
    sp->ddt.Blm [i] =(double complex *)m_alloc2(mn*(mn+2),sizeof(double complex),"setup_sp(),sp->ddt.Blm[i]");
    sp->ddt.alm [i] =(double complex *)m_alloc2(mn*(mn+2),sizeof(double complex),"setup_sp(),sp->ddt.alm[i]");
    sp->ddt.blm [i] =(double complex *)m_alloc2(mn*(mn+2),sizeof(double complex),"setup_sp(),sp->ddt.blm[i]");
  }
}

void gauleg_dv(double a,double b,double *x,double *w,int bn,int dv)
{
  double xt[bn],wt[bn];
  gauleg(-1.0, 1.0,xt,wt,bn);
  
  double h,dh,x0,x1,cx,cc;
  int d,i,j;
  h=b-a;
  dh=h/(double)dv;
  x1=a;
  j=0;
  for(d=0;d<dv;d++){
    x0=x1;
    x1=x0+dh;
    
    cx=0.5*(x1-x0);
    cc=0.5*(x1+x0);
    for(i=0;i<bn;i++){
      x[j]= cx*xt[i]+cc;
      w[j]= cx*wt[i];
      j++;
    }
  }
}

void free_sp(SPD *sp)
{
  int i;
  
  free(sp->ddt.xp);   free(sp->ddt.wp); 
  free(sp->ddt.xt);   free(sp->ddt.wt);   
  free(sp->ddt.eri);  free(sp->ddt.hri);  
  free(sp->ddt.ers);  free(sp->ddt.hrs);

  free(sp->ddt.cab);
  free(sp->ddt.ca);  free(sp->ddt.cb);
  free(sp->ddt.ca0); free(sp->ddt.cb0);
  free(sp->ddt.ceM); free(sp->ddt.chM);

  for(i=0;i<sp->n_l;i++){
    free(sp->ddt.alm[i]);  free(sp->ddt.blm[i]);
    free(sp->ddt.Alm[i]);  free(sp->ddt.Blm[i]);
  }
  free(sp->ddt.alm);  free(sp->ddt.blm);
  free(sp->ddt.Alm);  free(sp->ddt.Blm);
  
  sp->n_l=0;
  free(sp->a);
  free(sp->ns);
}

void setup_cab(SPD *sp,Bobj *bm)
{
  double x,a2;
  double *psi,*dpsi;
  int mn,nn,i,nl;
  
  mn=sp->l_limit;
  nl=sp->n_l;
  psi =(double *)m_alloc2(mn+1,sizeof(double),"setup_cab(),psi");  
  dpsi=(double *)m_alloc2(mn+1,sizeof(double),"setup_cab(),dpsi");
  x=2.0*M_PI*bm->n_0*sp->a[nl-1]/bm->lambda_0;
  rctjd(mn,x,&nn,psi,dpsi);
  if(nn<sp->ddt.l_max) sp->ddt.l_max=nn;
  a2=pow(sp->a[nl-1],2);
  for(i=1;i<=mn;i++){ 
    sp->ddt.cab[i]=a2/((double)(i*(i+1))*psi[i]);
  }
  
  free(psi);
  free(dpsi);
}

void initialize_eh_r(SPD *sp,Bobj *bm)
{
  double complex e[3],h[3];
  double r,theta,phi,x[3],sin_t,cos_t,sin_p,cos_p;
  int i,j,nt,np;

  nt=sp->ddt.nt;
  np=sp->ddt.np;
  r=sp->a[sp->n_l-1];
  #pragma omp parallel for schedule(dynamic) private(theta,sin_t,cos_t,j,phi,sin_p,cos_p,x,e,h) // omp parallel
  for(i=0;i<nt;i++){
    theta=sp->ddt.xt[i];
    sin_t=sin(theta);    cos_t=cos(theta);
    for(j=0;j<np;j++){
      phi=sp->ddt.xp[j];
      sin_p=sin(phi);      cos_p=cos(phi);
      x[0]=r*sin_t*cos_p+sp->xs;
      x[1]=r*sin_t*sin_p+sp->ys;
      x[2]=r*cos_t      +sp->zs;
      calc_mfb_EH(e,h,x,bm);
      sp->ddt.eri[i*np+j]=e[0]*sin_t*cos_p+e[1]*sin_t*sin_p+e[2]*cos_t;
      sp->ddt.hri[i*np+j]=h[0]*sin_t*cos_p+h[1]*sin_t*sin_p+h[2]*cos_t;
    }
  }
}

void check_l_limit_ms(MSPD *msp)
{
  int i;
  for(i=0;i<msp->n_sphr;i++){
    if(msp->sp[i].l_limit>msp->sp[i].ddt.l_max){
      printf("Overflow and underflow problem of Riccati-Bessel function occurred sphere id %d. Check the data precision.\n",i);
      printf("Available order number is less than %d.\n",msp->sp[i].ddt.l_max); 
    }   
  }
}

void characteristic_matrix(SPD *sp,Bobj *bm)
{
  double complex n0,n1,z0,z1,ce,cm,i_n1,teM[4],thM[4],tt[8],tr[2];
  double complex *psic0,*dpsic0,*psic1,*dpsic1,*chic0,*dchic0,*chic1,*dchic1,*reM,*rhM;
  int i,nn,l,mn,nb;

  mn=sp->l_limit;
  nb=sp->n_l;

  reM=(double complex *)m_alloc2(2*(mn+1),sizeof(double complex),"characteristic_matrix(),reM");
  rhM=(double complex *)m_alloc2(2*(mn+1),sizeof(double complex),"characteristic_matrix(),rhM");

  psic0 =(double complex *)m_alloc2(mn+1,sizeof(double complex),"characteristic_matrix(),psic0");
  dpsic0=(double complex *)m_alloc2(mn+1,sizeof(double complex),"characteristic_matrix(),dpsic0");
  psic1 =(double complex *)m_alloc2(mn+1,sizeof(double complex),"characteristic_matrix(),psic1");
  dpsic1=(double complex *)m_alloc2(mn+1,sizeof(double complex),"characteristic_matrix(),dpsic1");
  chic0 =(double complex *)m_alloc2(mn+1,sizeof(double complex),"characteristic_matrix(),chic0");
  dchic0=(double complex *)m_alloc2(mn+1,sizeof(double complex),"characteristic_matrix(),dchic0");
  chic1 =(double complex *)m_alloc2(mn+1,sizeof(double complex),"characteristic_matrix(),chic1");
  dchic1=(double complex *)m_alloc2(mn+1,sizeof(double complex),"characteristic_matrix(),dchic1");
  
  i=0;
  n0=sp->ns[i+0];
  n1=sp->ns[i+1];
  z0=2.0*M_PI*n0*sp->a[i]/bm->lambda_0;
  z1=2.0*M_PI*n1*sp->a[i]/bm->lambda_0;
  
  rctjc (mn,z0,&nn,psic0,dpsic0);  if(nn<sp->ddt.l_max) sp->ddt.l_max=nn;
  rctjc (mn,z1,&nn,psic1,dpsic1);  if(nn<sp->ddt.l_max) sp->ddt.l_max=nn;
  rctyc (mn,z1,&nn,chic1,dchic1);  if(nn<sp->ddt.l_max) sp->ddt.l_max=nn;
  
  i=0;
  i_n1=1.0/n1;
  ce=n0*i_n1*i_n1;
  cm=i_n1;
  for(l=0;l<=mn;l++){
    tt[0]=dpsic0[l]* chic1[l];    tt[1]= psic0[l]*dchic1[l];
    tt[2]= psic0[l]*dpsic1[l];    tt[3]=dpsic0[l]* psic1[l];
    
    teM[0]=ce*(n1*tt[0]-n0*tt[1]);    teM[1]=0.0;
    teM[2]=ce*(n0*tt[2]-n1*tt[3]);    teM[3]=0.0;

    thM[0]=cm*(n0*tt[0]-n1*tt[1]);    thM[1]=0.0;
    thM[2]=cm*(n1*tt[2]-n0*tt[3]);    thM[3]=0.0;
    
    reM[l*2+0]=teM[0];    reM[l*2+1]=teM[2];
    rhM[l*2+0]=thM[0];    rhM[l*2+1]=thM[2];

    sp->ddt.ceM[2*(mn+1)*i+2*l+0]=teM[0];    sp->ddt.ceM[2*(mn+1)*i+2*l+1]=teM[2];
    sp->ddt.chM[2*(mn+1)*i+2*l+0]=thM[0];    sp->ddt.chM[2*(mn+1)*i+2*l+1]=thM[2];
  }
  for(i=1;i<nb-1;i++){
    n0=sp->ns[i+0];
    n1=sp->ns[i+1];
    z0=2.0*M_PI*n0*sp->a[i]/bm->lambda_0;
    z1=2.0*M_PI*n1*sp->a[i]/bm->lambda_0;
    rctjc (mn,z0,&nn,psic0,dpsic0);  if(nn<sp->ddt.l_max) sp->ddt.l_max=nn;
    rctjc (mn,z1,&nn,psic1,dpsic1);  if(nn<sp->ddt.l_max) sp->ddt.l_max=nn;
    rctyc (mn,z0,&nn,chic0,dchic0);  if(nn<sp->ddt.l_max) sp->ddt.l_max=nn;
    rctyc (mn,z1,&nn,chic1,dchic1);  if(nn<sp->ddt.l_max) sp->ddt.l_max=nn;

    i_n1=1.0/n1;
    cm=i_n1;
    ce=n0*i_n1*i_n1;
    for(l=0;l<=mn;l++){
      tt[0]=dpsic0[l]* chic1[l];      tt[1]= psic0[l]*dchic1[l];
      tt[2]=dchic0[l]* chic1[l];      tt[3]= chic0[l]*dchic1[l];
      tt[4]= psic0[l]*dpsic1[l];      tt[5]=dpsic0[l]* psic1[l];
      tt[6]= chic0[l]*dpsic1[l];      tt[7]=dchic0[l]* psic1[l];
      // teM
      teM[0]=ce*(n1*tt[0]-n0*tt[1]);      teM[1]=ce*(n1*tt[2]-n0*tt[3]);
      teM[2]=ce*(n0*tt[4]-n1*tt[5]);      teM[3]=ce*(n0*tt[6]-n1*tt[7]);
      // thM
      thM[0]=cm*(n0*tt[0]-n1*tt[1]);      thM[1]=cm*(n0*tt[2]-n1*tt[3]);
      thM[2]=cm*(n1*tt[4]-n0*tt[5]);      thM[3]=cm*(n1*tt[6]-n0*tt[7]);
      // reM
      tr[0]=teM[0]*reM[l*2+0]+teM[1]*reM[l*2+1];
      tr[1]=teM[2]*reM[l*2+0]+teM[3]*reM[l*2+1];
      reM[2*l+0]=tr[0];
      reM[2*l+1]=tr[1];
      // rhM
      tr[0]=thM[0]*rhM[l*2+0]+thM[1]*rhM[l*2+1];
      tr[1]=thM[2]*rhM[l*2+0]+thM[3]*rhM[l*2+1];
      rhM[2*l+0]=tr[0];
      rhM[2*l+1]=tr[1];

      sp->ddt.ceM[2*(mn+1)*i+2*l+0]=reM[2*l+0];    sp->ddt.ceM[2*(mn+1)*i+2*l+1]=reM[2*l+1];
      sp->ddt.chM[2*(mn+1)*i+2*l+0]=rhM[2*l+0];    sp->ddt.chM[2*(mn+1)*i+2*l+1]=rhM[2*l+1];
    }
  }
  i=nb-1;
  n0=sp->ns[i];
  n1=bm->n_0;
  z0=2.0*M_PI*n0*sp->a[i]/bm->lambda_0;
  z1=2.0*M_PI*n1*sp->a[i]/bm->lambda_0;
  rctjc (mn,z0,&nn,psic0,dpsic0);  if(nn<sp->ddt.l_max) sp->ddt.l_max=nn;
  rctyc (mn,z0,&nn,chic0,dchic0);  if(nn<sp->ddt.l_max) sp->ddt.l_max=nn;
  rctjc (mn,z1,&nn,psic1,dpsic1);  if(nn<sp->ddt.l_max) sp->ddt.l_max=nn;
  rctyc (mn,z1,&nn,chic1,dchic1);  if(nn<sp->ddt.l_max) sp->ddt.l_max=nn;
  
  i_n1=1.0/n1;
  cm=i_n1*I;
  ce=n0*i_n1*i_n1*I;
  for(l=0;l<=mn;l++){
    tt[0]=I* psic0[l]*dchic1[l];    tt[1]=I*dpsic0[l]* chic1[l];
    tt[2]=I* chic0[l]*dchic1[l];    tt[3]=I*dchic0[l]* chic1[l];
    tt[4]=   psic0[l]*dpsic1[l];    tt[5]=  dpsic0[l]* psic1[l];
    tt[6]=   chic0[l]*dpsic1[l];    tt[7]=  dchic0[l]* psic1[l];
    // eM
    teM[0]=ce*(n0*tt[0]-n1*tt[1]);    teM[1]=ce*(n0*tt[2]-n1*tt[3]);
    teM[2]=ce*(n0*tt[4]-n1*tt[5]);    teM[3]=ce*(n0*tt[6]-n1*tt[7]);
    // hM
    thM[0]=cm*(n1*tt[0]-n0*tt[1]);    thM[1]=cm*(n1*tt[2]-n0*tt[3]);
    thM[2]=cm*(n1*tt[4]-n0*tt[5]);    thM[3]=cm*(n1*tt[6]-n0*tt[7]);
    // reM
    tr[0]=teM[0]*reM[l*2+0]+teM[1]*reM[l*2+1];
    tr[1]=teM[2]*reM[l*2+0]+teM[3]*reM[l*2+1];
    reM[2*l+0]=tr[0];
    reM[2*l+1]=tr[1];
    // rhM                                                                                                                                                                        
    tr[0]=thM[0]*rhM[l*2+0]+thM[1]*rhM[l*2+1];
    tr[1]=thM[2]*rhM[l*2+0]+thM[3]*rhM[l*2+1];
    rhM[2*l+0]=tr[0];
    rhM[2*l+1]=tr[1];

    sp->ddt.ceM[2*(mn+1)*i+2*l+0]=reM[2*l+0];    sp->ddt.ceM[2*(mn+1)*i+2*l+1]=reM[2*l+1];
    sp->ddt.chM[2*(mn+1)*i+2*l+0]=rhM[2*l+0];    sp->ddt.chM[2*(mn+1)*i+2*l+1]=rhM[2*l+1];
  }
  
  for(l=0;l<=mn;l++){
    tr[0]=reM[2*l+0];
    tr[1]=reM[2*l+1];
    reM[2*l+0]=(tr[0]+tr[1])/(tr[0]*tr[0]-tr[1]*tr[1]); // A_lm^{(0)} / A_lm
    reM[2*l+1]=reM[2*l+0]*tr[1];                        // a_lm / A_lm

    tr[0]=rhM[2*l+0];
    tr[1]=rhM[2*l+1];
    rhM[2*l+0]=(tr[0]+tr[1])/(tr[0]*tr[0]-tr[1]*tr[1]); // B_lm^{(0)} / B_lm   
    rhM[2*l+1]=rhM[2*l+0]*tr[1];                        // b_lm / B_lm           
  }
  
  // initialize 0-coefficient
  for(l=1;l<=mn;l++){
    sp->ddt.ca0[l]=reM[2*l+0];    sp->ddt.cb0[l]=rhM[2*l+0];
    sp->ddt.ca [l]=reM[2*l+1];    sp->ddt.cb [l]=rhM[2*l+1];
  }

  free(psic0);  free(dpsic0);
  free(psic1);  free(dpsic1);
  free(chic0);  free(dchic0);
  free(chic1);  free(dchic1);

  free(reM);  free(rhM);
}

void coefficient(SPD *sp,Bobj *bm)
{

  int m,l,i,j,tt,ti,t,lm,nl,np,nt;
  size_t ms,lmax;
  
  lmax=(size_t)sp->l_limit;
  ms=gsl_sf_legendre_array_n(lmax);
  lm=(int)lmax;
  nl=sp->n_l;
  np=sp->ddt.np;
  nt=sp->ddt.nt;
  
  for(ti=0;ti<lm*(lm+2);ti++){
    sp->ddt.Alm[nl][ti]=0.0;     sp->ddt.Blm[nl][ti]=0.0;
  }

  #pragma omp parallel private(i,j,l,m,t,tt,ti) // omp parallel
  {
  double complex Yp,Ym;
  double theta,sin_t,cos_t,flag;    
    
  double *sphPlm=(double *)m_alloc2(ms,sizeof(double),"coefficient(),*sphPlm");
  double complex *e_phim=(double complex *)m_alloc2(np,sizeof(double complex),"coefficient(),*i_phim");
  double complex *tmpAlm=(double complex *)m_alloc2(lm*(lm+2),sizeof(double complex),"coefficient(),*tmpAlm");
  double complex *tmpBlm=(double complex *)m_alloc2(lm*(lm+2),sizeof(double complex),"coefficient(),*tmpBlm");
  
  #pragma omp for schedule(dynamic) // omp parallel for
  for(i=0;i<nt;i++){
    theta=sp->ddt.xt[i];
    sin_t=sin(theta);    cos_t=cos(theta);
    tt=0;
    m=0;
    gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM,lmax,cos_t,-1,sphPlm);
    for(l=1;l<=lm;l++){
      Yp=sphPlm[gsl_sf_legendre_array_index(l,m)];
      for(j=0;j<np;j++){
        tmpAlm[tt]+=(sp->ddt.eri[i*np+j]+sp->ddt.ers[i*np+j])*conj(Yp)*sp->ddt.wp[j];
        tmpBlm[tt]+=(sp->ddt.hri[i*np+j]+sp->ddt.hrs[i*np+j])*conj(Yp)*sp->ddt.wp[j];
      }
      tt++;
    }
    for(m=1;m<=lm;m++){
      for(t=0;t<np;t++) e_phim[t]=cos((double)m*sp->ddt.xp[t])+I*sin((double)m*sp->ddt.xp[t]);
      if(m%2==0) flag=1.0;
      else flag=-1.0;
      for(l=m;l<=lm;l++){
        for(j=0;j<np;j++){
          Yp=sphPlm[gsl_sf_legendre_array_index(l,m)]*e_phim[j];
          Ym=flag*conj(Yp);
          tmpAlm[tt  ]+=(sp->ddt.eri[i*np+j]+sp->ddt.ers[i*np+j])*conj(Yp)*sp->ddt.wp[j];
          tmpBlm[tt  ]+=(sp->ddt.hri[i*np+j]+sp->ddt.hrs[i*np+j])*conj(Yp)*sp->ddt.wp[j];
          tmpAlm[tt+1]+=(sp->ddt.eri[i*np+j]+sp->ddt.ers[i*np+j])*conj(Ym)*sp->ddt.wp[j];
          tmpBlm[tt+1]+=(sp->ddt.hri[i*np+j]+sp->ddt.hrs[i*np+j])*conj(Ym)*sp->ddt.wp[j];
        }
        tt+=2;
      }
    }
    for(ti=0;ti<lm*(lm+2);ti++){
      #pragma omp critical // omp critical
      sp->ddt.Alm[nl][ti]+=tmpAlm[ti]*sin_t*sp->ddt.wt[i];
      #pragma omp critical // omp critical
      sp->ddt.Blm[nl][ti]+=tmpBlm[ti]*sin_t*sp->ddt.wt[i];
      tmpAlm[ti]=0.0;      tmpBlm[ti]=0.0;
    }
  }
  free(sphPlm);
  free(e_phim);
  free(tmpAlm);  free(tmpBlm);
  } // omp parallel end

  tt=0;  m=0;
  for(l=1;l<=lm;l++){
    sp->ddt.Alm[nl][tt]*=sp->ddt.cab[l];
    sp->ddt.Blm[nl][tt]*=sp->ddt.cab[l];
    sp->ddt.alm[nl][tt]=sp->ddt.ca[l]*sp->ddt.Alm[nl][tt];
    sp->ddt.blm[nl][tt]=sp->ddt.cb[l]*sp->ddt.Blm[nl][tt];
    tt++;
  }
  for(m=1;m<=lm;m++){
    for(l=m;l<=lm;l++){
      sp->ddt.Alm[nl][tt  ]*=sp->ddt.cab[l];   
      sp->ddt.alm[nl][tt  ] =sp->ddt.ca[l]*sp->ddt.Alm[nl][tt  ];
      sp->ddt.Alm[nl][tt+1]*=sp->ddt.cab[l];
      sp->ddt.alm[nl][tt+1] =sp->ddt.ca[l]*sp->ddt.Alm[nl][tt+1];
      sp->ddt.Blm[nl][tt  ]*=sp->ddt.cab[l];
      sp->ddt.blm[nl][tt  ] =sp->ddt.cb[l]*sp->ddt.Blm[nl][tt  ];
      sp->ddt.Blm[nl][tt+1]*=sp->ddt.cab[l];
      sp->ddt.blm[nl][tt+1] =sp->ddt.cb[l]*sp->ddt.Blm[nl][tt+1];
      tt+=2;
    }
  }
  
  for(i=0;i<nt;i++){
    for(j=0;j<np;j++){
      sp->ddt.ers[i*np+j]=0.0;
      sp->ddt.hrs[i*np+j]=0.0;
    }
  }
} 

void all_coefficient(SPD *sp,Bobj *bm)
{
  int i,l,m,cc;
  const int nb=sp->n_l;
  const int mn=sp->l_limit;

  // 0-coefficient
  cc=0;
  m=0;
  for(l=1;l<=mn;l++){
    sp->ddt.Alm[ 0][cc+0]=sp->ddt.ca0[l]*sp->ddt.Alm[nb][cc+0];    sp->ddt.alm[0][cc+0]=0.0;
    sp->ddt.Blm[ 0][cc+0]=sp->ddt.cb0[l]*sp->ddt.Blm[nb][cc+0];    sp->ddt.blm[0][cc+0]=0.0;
    cc++;
  }
  for(m=1;m<=mn;m++){
    for(l=m;l<=mn;l++){
      sp->ddt.Alm[ 0][cc+0]=sp->ddt.ca0[l]*sp->ddt.Alm[nb][cc+0];      sp->ddt.alm[0][cc+0]=0.0;
      sp->ddt.Blm[ 0][cc+0]=sp->ddt.cb0[l]*sp->ddt.Blm[nb][cc+0];      sp->ddt.blm[0][cc+0]=0.0;
      sp->ddt.Alm[ 0][cc+1]=sp->ddt.ca0[l]*sp->ddt.Alm[nb][cc+1];      sp->ddt.alm[0][cc+1]=0.0;
      sp->ddt.Blm[ 0][cc+1]=sp->ddt.cb0[l]*sp->ddt.Blm[nb][cc+1];      sp->ddt.blm[0][cc+1]=0.0;
      cc+=2;
    }
  }

  // 1-coefficient
  i=0;
  cc=0;
  m=0;
  for(l=1;l<=mn;l++){
    sp->ddt.Alm[i+1][cc+0]=sp->ddt.ceM[2*(mn+1)*i+l*2+0]*sp->ddt.Alm[0][cc+0];    
    sp->ddt.alm[i+1][cc+0]=sp->ddt.ceM[2*(mn+1)*i+l*2+1]*sp->ddt.Alm[0][cc+0];
    sp->ddt.Blm[i+1][cc+0]=sp->ddt.chM[2*(mn+1)*i+l*2+0]*sp->ddt.Blm[0][cc+0];
    sp->ddt.blm[i+1][cc+0]=sp->ddt.chM[2*(mn+1)*i+l*2+1]*sp->ddt.Blm[0][cc+0];
    cc++;
  }
  for(m=1;m<=mn;m++){
    for(l=m;l<=mn;l++){
      sp->ddt.Alm[i+1][cc+0]=sp->ddt.ceM[2*(mn+1)*i+l*2+0]*sp->ddt.Alm[0][cc+0];
      sp->ddt.alm[i+1][cc+0]=sp->ddt.ceM[2*(mn+1)*i+l*2+1]*sp->ddt.Alm[0][cc+0];
      sp->ddt.Blm[i+1][cc+0]=sp->ddt.chM[2*(mn+1)*i+l*2+0]*sp->ddt.Blm[0][cc+0];
      sp->ddt.blm[i+1][cc+0]=sp->ddt.chM[2*(mn+1)*i+l*2+1]*sp->ddt.Blm[0][cc+0];

      sp->ddt.Alm[i+1][cc+1]=sp->ddt.ceM[2*(mn+1)*i+l*2+0]*sp->ddt.Alm[0][cc+1];
      sp->ddt.alm[i+1][cc+1]=sp->ddt.ceM[2*(mn+1)*i+l*2+1]*sp->ddt.Alm[0][cc+1];
      sp->ddt.Blm[i+1][cc+1]=sp->ddt.chM[2*(mn+1)*i+l*2+0]*sp->ddt.Blm[0][cc+1];
      sp->ddt.blm[i+1][cc+1]=sp->ddt.chM[2*(mn+1)*i+l*2+1]*sp->ddt.Blm[0][cc+1];
      cc+=2;
    }
  }
  // 2 to nb-1 coefficient
  for(i=1;i<nb-1;i++){
    cc=0;
    m=0;
    for(l=1;l<=mn;l++){
      sp->ddt.Alm[i+1][cc+0]=sp->ddt.ceM[2*(mn+1)*i+l*2+0]*sp->ddt.Alm[0][cc+0];
      sp->ddt.alm[i+1][cc+0]=sp->ddt.ceM[2*(mn+1)*i+l*2+1]*sp->ddt.Alm[0][cc+0];
      sp->ddt.Blm[i+1][cc+0]=sp->ddt.chM[2*(mn+1)*i+l*2+0]*sp->ddt.Blm[0][cc+0];
      sp->ddt.blm[i+1][cc+0]=sp->ddt.chM[2*(mn+1)*i+l*2+1]*sp->ddt.Blm[0][cc+0];
      cc++;
    }
    for(m=1;m<=mn;m++){
      for(l=m;l<=mn;l++){
        sp->ddt.Alm[i+1][cc+0]=sp->ddt.ceM[2*(mn+1)*i+l*2+0]*sp->ddt.Alm[0][cc+0];
        sp->ddt.alm[i+1][cc+0]=sp->ddt.ceM[2*(mn+1)*i+l*2+1]*sp->ddt.Alm[0][cc+0];
        sp->ddt.Blm[i+1][cc+0]=sp->ddt.chM[2*(mn+1)*i+l*2+0]*sp->ddt.Blm[0][cc+0];
        sp->ddt.blm[i+1][cc+0]=sp->ddt.chM[2*(mn+1)*i+l*2+1]*sp->ddt.Blm[0][cc+0];

        sp->ddt.Alm[i+1][cc+1]=sp->ddt.ceM[2*(mn+1)*i+l*2+0]*sp->ddt.Alm[0][cc+1];
        sp->ddt.alm[i+1][cc+1]=sp->ddt.ceM[2*(mn+1)*i+l*2+1]*sp->ddt.Alm[0][cc+1];
        sp->ddt.Blm[i+1][cc+1]=sp->ddt.chM[2*(mn+1)*i+l*2+0]*sp->ddt.Blm[0][cc+1];
        sp->ddt.blm[i+1][cc+1]=sp->ddt.chM[2*(mn+1)*i+l*2+1]*sp->ddt.Blm[0][cc+1];
        cc+=2;
      }
    }
  }
}

void field_s_ehr(int src,int obj,MSPD *msp)
{
  void scattered_EH(double complex *e3,double complex *h3,double *x3,SPD *sp,Bobj *bm);
  
  double complex es[3],hs[3];
  double xt,xp,cos_t,sin_t,cos_p,sin_p;
  double x[3];
  double r=msp->sp[obj].a[msp->sp[obj].n_l-1];
  int np=msp->sp[obj].ddt.np;
  int nt=msp->sp[obj].ddt.nt;
  int i,j;

  #pragma omp parallel for schedule(dynamic) private(xt,cos_t,sin_t,j,xp,cos_p,sin_p,x,es,hs)
  for(i=0;i<nt;i++){
    xt=msp->sp[obj].ddt.xt[i];
    cos_t=cos(xt);    sin_t=sin(xt);
    for(j=0;j<np;j++){
      xp=msp->sp[obj].ddt.xp[j];
      cos_p=cos(xp);      sin_p=sin(xp);
      x[0]=r*sin_t*cos_p+msp->sp[obj].xs;
      x[1]=r*sin_t*sin_p+msp->sp[obj].ys;
      x[2]=r*cos_t      +msp->sp[obj].zs;
      scattered_EH(es,hs,x,&(msp->sp[src]),&(msp->bm));
      msp->sp[obj].ddt.ers[i*np+j]+=es[0]*sin_t*cos_p+es[1]*sin_t*sin_p+es[2]*cos_t;
      msp->sp[obj].ddt.hrs[i*np+j]+=hs[0]*sin_t*cos_p+hs[1]*sin_t*sin_p+hs[2]*cos_t;
    }
  }
}

void scattered_EH(double complex *e3,double complex *h3,double *x3,SPD *sp,Bobj *bm)
{
  double complex er,et,ep,hr,ht,hp,Yp,Ym,dYp,dYm,dep,expi;
  double r,rxy,r2,cos_t,sin_t,cos_p,sin_p,ker,ke,flag,i_ne,ne,x,y,z,i_sin_t;
  int l,m,tt,nn,lm,ai,nb;
  size_t ms,lmax;
  
  lmax=(size_t)sp->l_limit;
  ms=gsl_sf_legendre_array_n(lmax);
  lm=(int)lmax;
  nb=sp->n_l;
  
  double *sphPlm =(double *)m_alloc2(ms,sizeof(double),"scattered_EH(),sphPlm");
  double *dsphPlm=(double *)m_alloc2(ms,sizeof(double),"scattered_EH(),dsphPlm");
  double complex *xi =(double complex *)m_alloc2(lm+1,sizeof(double complex),"scattered_EH(),xi");
  double complex *dxi=(double complex *)m_alloc2(lm+1,sizeof(double complex),"scattered_EH(),dxi");

  x=x3[0]-sp->xs;  y=x3[1]-sp->ys;  z=x3[2]-sp->zs;
  r2=x*x+y*y+z*z;    r=sqrt(r2);
  rxy=sqrt(x*x+y*y);
  if(rxy==0.0){ // x==0,y==0
    x=z*0.7e-7;
    y=z*0.7e-7;
    r2=x*x+y*y+z*z;    r=sqrt(r2);
    rxy=sqrt(x*x+y*y); 
  }
  cos_t=z/r;    sin_t=rxy/r;   i_sin_t=r/rxy;
  cos_p=x/rxy;  sin_p=y/rxy;
  ke =2.0*M_PI*bm->n_0/bm->lambda_0;
  ker=ke*r;
  ne=bm->n_0;
  i_ne=1.0/(bm->n_0);
  
  rcth1d(lm,ker,&nn,xi,dxi);  if(nn<sp->ddt.l_max) sp->ddt.l_max=nn;
  dep=cos_p+I*sin_p; expi=1.0;
  er=0.0;  et=0.0;  ep=0.0;
  hr=0.0;  ht=0.0;  hp=0.0;
  tt=0;  m=0;
  gsl_sf_legendre_deriv_alt_array_e(GSL_SF_LEGENDRE_SPHARM,lmax,cos_t,-1,sphPlm,dsphPlm);
  for(l=1;l<=lm;l++){
    ai=gsl_sf_legendre_array_index(l,m);
    Yp = sphPlm[ai];
    dYp=dsphPlm[ai];
    er+=(double)(l*(l+1))*sp->ddt.alm[nb][tt]*xi[l]*Yp;
    hr+=(double)(l*(l+1))*sp->ddt.blm[nb][tt]*xi[l]*Yp;
    et+=sp->ddt.alm[nb][tt]*dxi[l]*dYp-(double)m*i_ne*sp->ddt.blm[nb][tt]*xi[l]*Yp*i_sin_t;
    ht+=sp->ddt.blm[nb][tt]*dxi[l]*dYp+(double)m*  ne*sp->ddt.alm[nb][tt]*xi[l]*Yp*i_sin_t;
    ep+=(double)m*sp->ddt.alm[nb][tt]*dxi[l]*Yp*i_sin_t-i_ne*sp->ddt.blm[nb][tt]*xi[l]*dYp;
    hp+=(double)m*sp->ddt.blm[nb][tt]*dxi[l]*Yp*i_sin_t+  ne*sp->ddt.alm[nb][tt]*xi[l]*dYp;
    tt++;
  }
  for(m=1;m<=lm;m++){
    expi*=dep;
    if(m%2==0)flag= 1.0;
    else      flag=-1.0;
    for(l=m;l<=lm;l++){
      ai=gsl_sf_legendre_array_index(l,m);
      Yp = sphPlm[ai]*expi;
      dYp=dsphPlm[ai]*expi;
      Ym =flag*conj( Yp);
      dYm=flag*conj(dYp);
      er+=(double)(l*(l+1))*xi[l]*(sp->ddt.alm[nb][tt  ]*Yp+sp->ddt.alm[nb][tt+1]*Ym);
      et+=dxi[l]*(sp->ddt.alm[nb][tt  ]*dYp+sp->ddt.alm[nb][tt+1]*dYm)
          -(double)m*i_ne*xi[l]*i_sin_t*(sp->ddt.blm[nb][tt  ]*Yp-sp->ddt.blm[nb][tt+1]*Ym);
      ep+=(double)m*dxi[l]*i_sin_t*(sp->ddt.alm[nb][tt  ]*Yp-sp->ddt.alm[nb][tt+1]*Ym)
          -i_ne*xi[l]*(sp->ddt.blm[nb][tt  ]*dYp+sp->ddt.blm[nb][tt+1]*dYm);
      hr+=(double)(l*(l+1))*xi[l]*(sp->ddt.blm[nb][tt  ]*Yp+sp->ddt.blm[nb][tt+1]*Ym); 
      ht+=dxi[l]*(sp->ddt.blm[nb][tt  ]*dYp+sp->ddt.blm[nb][tt+1]*dYm)
          +(double)m*ne*xi[l]*i_sin_t*(sp->ddt.alm[nb][tt  ]*Yp-sp->ddt.alm[nb][tt+1]*Ym);
      hp+=(double)m*dxi[l]*i_sin_t*(sp->ddt.blm[nb][tt  ]*Yp-sp->ddt.blm[nb][tt+1]*Ym)
          +ne*xi[l]*(sp->ddt.alm[nb][tt  ]*dYp+sp->ddt.alm[nb][tt+1]*dYm);
      tt+=2;
    }
  }
  er/=r2;
  et*=ke/r;
  ep*=I*ke/r;
  hr/=r2;
  ht*=ke/r;
  hp*=I*ke/r;

  e3[0]=er*sin_t*cos_p+et*cos_t*cos_p-ep*sin_p;
  e3[1]=er*sin_t*sin_p+et*cos_t*sin_p+ep*cos_p;
  e3[2]=er*cos_t-et*sin_t;
  h3[0]=hr*sin_t*cos_p+ht*cos_t*cos_p-hp*sin_p;
  h3[1]=hr*sin_t*sin_p+ht*cos_t*sin_p+hp*cos_p;
  h3[2]=hr*cos_t-ht*sin_t;

  free(sphPlm);  free(dsphPlm);
  free(xi);      free(dxi);
}

void internal_EH(int id,double complex *e3,double complex *h3,double *x3,SPD *sp,Bobj *bm)
{
  void internal_EH_r0(double complex *e,double complex *h,SPD *sp,Bobj *bm);
  
  double complex er,et,ep,hr,ht,hp,Yp,Ym,dYp,dYm,dep,expi,ke,ker,i_ne,ne;
  double r,rxy,r2,cos_t,sin_t,cos_p,sin_p,flag,x,y,z,i_sin_t;
  int l,m,tt,nn,lm,ai;
  size_t ms,lmax;
  
  lmax=(size_t)sp->ddt.l_max;
  ms=gsl_sf_legendre_array_n(lmax);
  lm=(int)lmax;
  
  double *sphPlm =(double *)m_alloc2(ms,sizeof(double),"internal_EH(),*sphPlm");
  double *dsphPlm=(double *)m_alloc2(ms,sizeof(double),"internal_EH(),*dsphPlm");
  double complex *psi =(double complex *)m_alloc2(lm+1,sizeof(double complex),"internal_EH(),*psi");
  double complex *dpsi=(double complex *)m_alloc2(lm+1,sizeof(double complex),"internal_EH(),*dpsi");
  double complex *chi =(double complex *)m_alloc2(lm+1,sizeof(double complex),"internal_EH(),*chi");
  double complex *dchi=(double complex *)m_alloc2(lm+1,sizeof(double complex),"internal_EH(),*dchi");
  
  x=x3[0]-sp->xs;  y=x3[1]-sp->ys;  z=x3[2]-sp->zs;
  r2=x*x+y*y+z*z;    r=sqrt(r2);
  rxy=sqrt(x*x+y*y);
  if(rxy==0.0){ // x==0,y==0
	  if(z==0.0){
      internal_EH_r0(e3,h3,sp,bm);
      return;
    } 
    x=z*0.7e-7;
    y=z*0.7e-7;
    r2=x*x+y*y+z*z;    r=sqrt(r2);
    rxy=sqrt(x*x+y*y); 
  }
  cos_t=z/r;    sin_t=rxy/r;   i_sin_t=r/rxy;
  cos_p=x/rxy;  sin_p=y/rxy;
  ke =2.0*M_PI*sp->ns[id]/bm->lambda_0;
  ker=ke*r;
  ne=sp->ns[id];
  i_ne=1.0/(sp->ns[id]);

  rctjc(lm,ker,&nn,psi,dpsi);  if(nn<sp->ddt.l_max) sp->ddt.l_max=nn;
  rctyc(lm,ker,&nn,chi,dchi);  if(nn<sp->ddt.l_max) sp->ddt.l_max=nn;
  dep=cos_p+I*sin_p; expi=1.0;
  er=0.0;  et=0.0;  ep=0.0;
  hr=0.0;  ht=0.0;  hp=0.0;
  tt=0;  m=0;
  gsl_sf_legendre_deriv_alt_array_e(GSL_SF_LEGENDRE_SPHARM,lmax,cos_t,-1,sphPlm,dsphPlm);
  for(l=1;l<=lm;l++){
    ai=gsl_sf_legendre_array_index(l,m);
    Yp = sphPlm[ai];
    dYp=dsphPlm[ai];
    er+=(double)(l*(l+1))*(sp->ddt.Alm[id][tt]*psi[l]+sp->ddt.alm[id][tt]*chi[l])*Yp;
    et+=(sp->ddt.Alm[id][tt]*dpsi[l]+sp->ddt.alm[id][tt]*dchi[l])*dYp
        -(double)m*i_ne*(sp->ddt.Blm[id][tt]*psi[l]+sp->ddt.blm[id][tt]*chi[l])*Yp*i_sin_t;
    ep+=(double)m*(sp->ddt.Alm[id][tt]*dpsi[l]+sp->ddt.alm[id][tt]*dchi[l])*Yp*i_sin_t
        -i_ne*(sp->ddt.Blm[id][tt]*psi[l]+sp->ddt.blm[id][tt]*chi[l])*dYp;
    hr+=(double)(l*(l+1))*(sp->ddt.Blm[id][tt]*psi[l]+sp->ddt.blm[id][tt]*chi[l])*Yp;
    ht+=(sp->ddt.Blm[id][tt]*dpsi[l]+sp->ddt.blm[id][tt]*dchi[l])*dYp
        +(double)m*ne*(sp->ddt.Alm[id][tt]*psi[l]+sp->ddt.alm[id][tt]*chi[l])*Yp*i_sin_t;
    hp+=(double)m*(sp->ddt.Blm[id][tt]*dpsi[l]+sp->ddt.blm[id][tt]*dchi[l])*Yp*i_sin_t
        +ne*(sp->ddt.Alm[id][tt]*psi[l]+sp->ddt.alm[id][tt]*chi[l])*dYp;
    tt++;
  }
  for(m=1;m<=lm;m++){
    expi*=dep;
    if(m%2==0)flag= 1.0;
    else      flag=-1.0;
    for(l=m;l<=lm;l++){
      ai=gsl_sf_legendre_array_index(l,m);
      Yp = sphPlm[ai]*expi;
      dYp=dsphPlm[ai]*expi;
      Ym =flag*conj( Yp);
      dYm=flag*conj(dYp);
      er+=(double)(l*(l+1))*(psi[l]*(sp->ddt.Alm[id][tt  ]*Yp+sp->ddt.Alm[id][tt+1]*Ym)
          +chi[l]*(sp->ddt.alm[id][tt]*Yp+sp->ddt.alm[id][tt+1]*Ym)); 
      et+=dpsi[l]*(sp->ddt.Alm[id][tt  ]*dYp+sp->ddt.Alm[id][tt+1]*dYm)
          +dchi[l]*(sp->ddt.alm[id][tt]*dYp+sp->ddt.alm[id][tt+1]*dYm)
          -(double)m*i_ne*i_sin_t*(psi[l]*(sp->ddt.Blm[id][tt  ]*Yp-sp->ddt.Blm[id][tt+1]*Ym)
          +chi[l]*(sp->ddt.blm[id][tt]*Yp-sp->ddt.blm[id][tt+1]*Ym));
      ep+=(double)m*i_sin_t*(dpsi[l]*(sp->ddt.Alm[id][tt  ]*Yp-sp->ddt.Alm[id][tt+1]*Ym)
          +dchi[l]*(sp->ddt.alm[id][tt]*Yp-sp->ddt.alm[id][tt+1]*Ym))
          -i_ne*(psi[l]*(sp->ddt.Blm[id][tt  ]*dYp+sp->ddt.Blm[id][tt+1]*dYm)
          +chi[l]*(sp->ddt.blm[id][tt]*dYp+sp->ddt.blm[id][tt+1]*dYm));
      hr+=(double)(l*(l+1))*(psi[l]*(sp->ddt.Blm[id][tt  ]*Yp+sp->ddt.Blm[id][tt+1]*Ym)
          +chi[l]*(sp->ddt.blm[id][tt]*Yp+sp->ddt.blm[id][tt+1]*Ym));
      ht+=dpsi[l]*(sp->ddt.Blm[id][tt  ]*dYp+sp->ddt.Blm[id][tt+1]*dYm)
          +dchi[l]*(sp->ddt.blm[id][tt]*dYp+sp->ddt.blm[id][tt+1]*dYm)
          +(double)m*ne*i_sin_t*(psi[l]*(sp->ddt.Alm[id][tt  ]*Yp-sp->ddt.Alm[id][tt+1]*Ym)
          +chi[l]*(sp->ddt.alm[id][tt]*Yp-sp->ddt.alm[id][tt+1]*Ym));
      hp+=(double)m*i_sin_t*(dpsi[l]*(sp->ddt.Blm[id][tt  ]*Yp-sp->ddt.Blm[id][tt+1]*Ym)
          +dchi[l]*(sp->ddt.blm[id][tt]*Yp-sp->ddt.blm[id][tt+1]*Ym))
          +ne*(psi[l]*(sp->ddt.Alm[id][tt  ]*dYp+sp->ddt.Alm[id][tt+1]*dYm)
          +chi[l]*(sp->ddt.alm[id][tt]*dYp+sp->ddt.alm[id][tt+1]*dYm));
      tt+=2;
    }
  }
  er/=r2;
  et*=ke/r;
  ep*=I*ke/r;
  hr/=r2;
  ht*=ke/r;
  hp*=I*ke/r;

  e3[0]=er*sin_t*cos_p+et*cos_t*cos_p-ep*sin_p;
  e3[1]=er*sin_t*sin_p+et*cos_t*sin_p+ep*cos_p;
  e3[2]=er*cos_t-et*sin_t;
  h3[0]=hr*sin_t*cos_p+ht*cos_t*cos_p-hp*sin_p;
  h3[1]=hr*sin_t*sin_p+ht*cos_t*sin_p+hp*cos_p;
  h3[2]=hr*cos_t-ht*sin_t;

  free(sphPlm);  free(dsphPlm);
  free(psi);      free(dpsi);
  free(chi);      free(dchi);
}

// verification function for sphere center field 
void internal_EH_r0(double complex *e,double complex *h,SPD *sp,Bobj *bm)
{
  double complex pickup_Alm_nl(int nl,int l,int m,SPD *sp);
  double complex pickup_Blm_nl(int nl,int l,int m,SPD *sp);

  double complex ke,er,et,ep,hr,ht,hp,A1m1,A1p1,A10,B1m1,B1p1,B10,Y10,Y1p1,Y1m1,dY10,dY1p1,dY1m1;  
  double cos_t,sin_t,i_sin_t,cos_p,sin_p;
  
  ke =2.0*M_PI*sp->ns[0]/bm->lambda_0;
  
  // assume r=(x,0,0) then x to 0
  cos_t=0.0;    sin_t=1.0;   i_sin_t=1.0;
  cos_p=1.0;    sin_p=0.0;

  A10 =pickup_Alm_nl(0,1,0,sp);
  A1p1=pickup_Alm_nl(0,1,1,sp);
  A1m1=pickup_Alm_nl(0,1,-1,sp);
  B10 =pickup_Blm_nl(0,1,0,sp);
  B1p1=pickup_Blm_nl(0,1,1,sp);
  B1m1=pickup_Blm_nl(0,1,-1,sp);

  Y10 = 0.5*sqrt(3.0/M_PI)*cos_t;
  Y1p1=-0.5*sqrt(3.0/(2.0*M_PI))*sin_t*(cos_p+I*sin_p);
  Y1m1= 0.5*sqrt(3.0/(2.0*M_PI))*sin_t/(cos_p+I*sin_p);
  
  dY10 =-0.5*sqrt(3.0/M_PI)*sin_t;
  dY1p1=-0.5*sqrt(3.0/(2.0*M_PI))*cos_t*(cos_p+I*sin_p);
  dY1m1= 0.5*sqrt(3.0/(2.0*M_PI))*cos_t/(cos_p+I*sin_p);
  
  er=2.0*ke*ke/(3.0)*(A10*Y10 +A1p1*Y1p1 +A1m1*Y1m1);
  et=2.0*ke*ke/(3.0)*(A10*dY10+A1p1*dY1p1+A1m1*dY1m1);
  ep=2.0*I*ke*ke*i_sin_t/(3.0)*(A1p1*Y1p1-A1m1*Y1m1);
  hr=2.0*ke*ke/(3.0)*(B10*Y10 +B1p1*Y1p1 +B1m1*Y1m1);
  ht=2.0*ke*ke/(3.0)*(B10*dY10+B1p1*dY1p1+B1m1*dY1m1);
  hp=2.0*I*ke*ke*i_sin_t/(3.0)*(B1p1*Y1p1-B1m1*Y1m1);
  
  e[0]=er*sin_t*cos_p+et*cos_t*cos_p-ep*sin_p;
  e[1]=er*sin_t*sin_p+et*cos_t*sin_p+ep*cos_p;
  e[2]=er*cos_t-et*sin_t;
  h[0]=hr*sin_t*cos_p+ht*cos_t*cos_p-hp*sin_p;
  h[1]=hr*sin_t*sin_p+ht*cos_t*sin_p+hp*cos_p;
  h[2]=hr*cos_t-ht*sin_t;
}

double complex pickup_Alm(int l,int m,SPD *sp)
{
  int nb=sp->n_l;
  int lm=sp->l_limit;
  int am=abs(m);
  int ll;
  
  if(l<am)  return 0.0;
  else {
    if(m==0)
      return sp->ddt.Alm[nb][l-1];
    else {
      ll=2*lm*am-am*am+am-2-lm+2*l;
      if (m>0) return sp->ddt.Alm[nb][ll+0];
      else     return sp->ddt.Alm[nb][ll+1];
    }
  }
}

double complex pickup_Blm(int l,int m,SPD *sp)
{
  int nb=sp->n_l;
  int lm=sp->l_limit;
  int am=abs(m);
  int ll;
  
  if(l<am)  return 0.0;
  else {
    if(m==0)
      return sp->ddt.Blm[nb][l-1];
    else {
      ll=2*lm*am-am*am+am-2-lm+2*l;
      if (m>0) return sp->ddt.Blm[nb][ll+0];
      else     return sp->ddt.Blm[nb][ll+1];
    }
  }
}

double complex pickup_alm(int l,int m,SPD *sp)
{
  int nb=sp->n_l;
  int lm=sp->l_limit;
  int am=abs(m);
  int ll;
  
  if(l<am)  return 0.0;
  else {
    if(m==0)
      return sp->ddt.alm[nb][l-1];
    else {
      ll=2*lm*am-am*am+am-2-lm+2*l;
      if (m>0) return sp->ddt.alm[nb][ll+0];
      else     return sp->ddt.alm[nb][ll+1];
    }
  }
}

double complex pickup_blm(int l,int m,SPD *sp)
{
  int nb=sp->n_l;
  int lm=sp->l_limit;
  int am=abs(m);
  int ll;
  
  if(l<am)  return 0.0;
  else {
    if(m==0)
      return sp->ddt.blm[nb][l-1];
    else {
      ll=2*lm*am-am*am+am-2-lm+2*l;
      if (m>0) return sp->ddt.blm[nb][ll+0];
      else     return sp->ddt.blm[nb][ll+1];
    }
  }
}

double complex pickup_Alm_nl(int nl,int l,int m,SPD *sp)
{
  int lm=sp->l_limit;
  int am=abs(m);
  int ll;
  
  if(l<am)  return 0.0;
  else {
    if(m==0)
      return sp->ddt.Alm[nl][l-1];
    else {
      ll=2*lm*am-am*am+am-2-lm+2*l;
      if (m>0) return sp->ddt.Alm[nl][ll+0];
      else     return sp->ddt.Alm[nl][ll+1];
    }
  }
}

double complex pickup_Blm_nl(int nl,int l,int m,SPD *sp)
{
  int lm=sp->l_limit;
  int am=abs(m);
  int ll;
  
  if(l<am)  return 0.0;
  else {
    if(m==0)
      return sp->ddt.Blm[nl][l-1];
    else {
      ll=2*lm*am-am*am+am-2-lm+2*l;
      if (m>0) return sp->ddt.Blm[nl][ll+0];
      else     return sp->ddt.Blm[nl][ll+1];
    }
  }
}


int layer_id(double r,SPD *sp)
{
  int i;
  for(i=0;i<sp->n_l;i++){
    if(r<=sp->a[i]){  // select inside layer when r is on boundary.
      return i;
    }
  }
  return sp->n_l;
}
