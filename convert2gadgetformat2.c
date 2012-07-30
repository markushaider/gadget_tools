#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "data_ops.h"
#include "snapshot_io.h"

int main(int argc, char** argv) {

  char *infile;
  infile=argv[1];
  char *outfile;
  outfile=argv[2];

    
  if(argc!=3) {
    printf("Usage: ./convert2gadgetformat2 <infile> <outfile>\n");
    return 0;
  }

  /* load snapshot file */
  io_header header;
  particle_data * P;
  int nTot;
  nTot = load_snapshot(infile, &header, &P);


 
  /* output information about the simulation */
  printf(">> Output from Header \n");
  printf("Boxsize of Simulation: %g\n", header.BoxSize);
  printf("Number of gas particles: %i\n", header.npart[0]);
  printf("Number of high-res dm particles: %i\n", header.npart[1]);
  printf("Number of low-res dm particles: %i\n", header.npart[5]);
  printf("Number of star particles: %i\n", header.npart[4]);

  /* write reduced snapshot file */
  write_snapshot2(outfile, &header, P);

  return 0;
}
