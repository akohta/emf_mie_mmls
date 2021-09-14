#include "emf_mie_mmls.h"

int main(int argc, char *argv[])
{
  if(argc!=2){
    printf("This program needs command line argument.\n");
    printf("Usage : %s 'output datafile name'\n",argv[0]);
    exit(0);
  }
  
  MSPD msp;
  
  read_data_mmls(&msp);               // seach the sphere datafile and load 
  print_data_mmls(&msp);              // print loaded data
  setup_mmls(&msp);                   // allocate memory and setup coefficients
  output_node_particles(argv[1],&msp);// output point cloud data of the nodes
  iterative_ops_mmls(&msp);           // solve multiple scattering
  write_dat_mmls(argv[1],&msp);       // write datafile
  free_mmls(&msp);                    // free allocated memory
  
  return 0;
}
