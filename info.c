#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "data_ops.h"
#include "snapshot_io.h"

int main(int argc, char** argv) {

  char *snapfile;
  snapfile=argv[1];

  if(argc!=2) {
    printf("Usage: ./info <snapshot>\n");
    return 0;
  }
  /* load snapshots */
  io_header header;
  particle_data * P;

  load_snapshot(snapfile , &header, &P);

  /* output information about the simulation */
  printf("\n Information about the snapshot: \n");
  printf("Boxsize: %g\n", header.BoxSize);
  printf("Number of gas particles: %i\n", header.npart[0]);
  printf("Number of high-res dm particles: %i\n", header.npart[1]);
  printf("Number of low-res dm particles: %i\n", header.npart[5]);
  printf("Number of star particles: %i\n", header.npart[4]);

  /* find the coordinate values */
  int i,j;
  int counter=0;
  float xMin, xMax, yMin, yMax, zMin, zMax;
  xMin = P[0].Pos[0];
  xMax = P[0].Pos[0];
  yMin = P[0].Pos[1];
  yMax = P[0].Pos[1];
  zMin = P[0].Pos[2];
  zMax = P[0].Pos[2];
  printf(" Values:\n");
  /* printf(" xMin: %g xMax: %f \n", xMin, xMax); */
  /* printf(" yMin: %g yMax: %f \n", yMin, yMax); */
  /* printf(" zMin: %g zMax: %f \n", zMin, zMax); */
  for(i=1; i< header.npart[0]+header.npart[1]; i++) {
    if(P[i].Pos[0] < xMin)
      xMin = P[i].Pos[0];
    if(P[i].Pos[0] > xMax)
      xMax = P[i].Pos[0];
    if(P[i].Pos[1] < yMin)
      yMin = P[i].Pos[1];
    if(P[i].Pos[1] > yMax)
      yMax = P[i].Pos[1];
    if(P[i].Pos[2] < zMin)
      zMin = P[i].Pos[2];
    if(P[i].Pos[2] > zMax)
      zMax = P[i].Pos[2];
   }
 
  printf(" xMin: %g xMax: %f \n", xMin, xMax);
  printf(" yMin: %g yMax: %f \n", yMin, yMax);
  printf(" zMin: %g zMax: %f \n", zMin, zMax);

  return 0;
}
