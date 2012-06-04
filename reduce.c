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
    printf("Usage: ./reduce <infile> <outfile>\n");
    return 0;
  }

  /* load snapshot file */
  io_header header1;
  particle_data * P1;
  int nTot;
  nTot = load_snapshot(infile, &header1, &P1);

  /* output information about the simulation */
  printf("Boxsize: %g\n", header1.BoxSize);
  printf("Number of gas particles: %i\n", header1.npart[0]);
  printf("Number of high-res dm particles: %i\n", header1.npart[1]);
  printf("Number of low-res dm particles: %i\n", header1.npart[5]);
  printf("Number of star particles: %i\n", header1.npart[4]);

  /* use only every xth particle */
  int reduction = 7;
  io_header header2;
  particle_data * P2;
  P2 = (particle_data *)calloc(nTot/reduction, sizeof(particle_data));
  
  int i;
  for (i=0;i<nTot/reduction;i++) {
      P2[i] = P1[i*reduction];
  }

  header2 = header1;
  header2.npart[0]/=reduction;
  header2.npart[1]/=reduction;
  header2.npart[4]/=reduction;
  header2.npart[5]/=reduction;
  header2.npartTotal[0]/=reduction;
  header2.npartTotal[1]/=reduction;
  header2.npartTotal[4]/=reduction;
  header2.npartTotal[5]/=reduction;

  /* write reduced snapshot file */
  write_snapshot(outfile, &header2, P2);

  return 0;
}
