#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "data_ops.h"
#include "snapshot_io.h"

int main() {


  io_header header;
  particle_data * P;


  int nTot = load_snapshot("../snapshot_349", &header, &P);

  /* output information about the simulation */
  printf("Boxsize: %g\n", header.BoxSize);
  printf("Number of gas particles: %i\n", header.npart[0]);
  printf("Number of high-res dm particles: %i\n", header.npart[1]);
  printf("Number of low-res dm particles: %i\n", header.npart[5]);
  printf("Number of star particles: %i\n", header.npart[4]);

  float box[6];
  box[0] = 0;
  box[1] = 0;
  box[2] = 0;
  box[3] = 1;
  box[4] = 1;
  box[5] = 1;
  writeIfritParticleFile("./349ifrit", P, nTot , box); //box: boundary box xl,yl,zl, xh,yh,zh

  return 0;
}
