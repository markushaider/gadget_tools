#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "data_ops.h"
#include "snapshot_io.h"

int main(int argc, char** argv) {

  char *snap_init;
  snap_init=argv[1];
  char *snap_final;
  snap_final=argv[2];
  float bboxX = atof(argv[3]);
  float bboxY = atof(argv[4]);
  float bboxZ = atof(argv[5]);
  float size = atof(argv[6]);

  if(argc!=7) {
    printf("Usage: ./getBBox <snapshot0> <snapshot_final> <x> <y> <z> <extent>\n");
    printf("<x> <y> <z> refer to the coordinates of the center of the desired volume (in kpc)\n");
    printf("extent is the size of the desired volume (in kpc)\n");
    return 0;
  }
  printf(" You specified\n");
  printf(" <x>: %g\n",bboxX);
  printf(" <y>: %g\n",bboxY);
  printf(" <z>: %g\n",bboxZ);
  printf(" extent: %g\n",size);
  printf("\n loading snapshot files...\n");
  /* load snapshots */
  io_header header_ic;
  particle_data * P_ic;
  io_header header;
  particle_data * P;
  load_snapshot(snap_init , &header_ic, &P_ic);
  load_snapshot(snap_final, &header, &P);

  /* output information about the simulation */
  printf("\n Information about the later snapshot: \n");
  printf("Boxsize: %g\n", header.BoxSize);
  printf("Number of gas particles: %i\n", header.npart[0]);
  printf("Number of high-res dm particles: %i\n", header.npart[1]);
  printf("Number of low-res dm particles: %i\n", header.npart[5]);
  printf("Number of star particles: %i\n", header.npart[4]);

  /* find the bounding box at the beginning */
  int i,j;
  int counter=0;
  float iBoxX1, iBoxX2, iBoxY1, iBoxY2, iBoxZ1, iBoxZ2;
  for(i=0; i< header.npart[0]+header.npart[1]; i++) {

    /* check whether particle is within the definded bounding box */
    if (P[i].Pos[0] > bboxX-size/2 && P[i].Pos[0] < bboxX+size/2 &&
	P[i].Pos[1] > bboxY-size/2 && P[i].Pos[1] < bboxY+size/2 &&
	P[i].Pos[2] > bboxZ-size/2 && P[i].Pos[2] < bboxZ+size/2) {
      /* if particle is inside, then get the position of the particle in the ic */
      for(j=0; j< header.npart[0]+header.npart[1]; j++) {
	if(P_ic[j].Id == P[i].Id) {
	  if (counter==0) {
	    iBoxX1 = P_ic[j].Pos[0];
	    iBoxX2 = P_ic[j].Pos[0];
	    iBoxY1 = P_ic[j].Pos[1];
	    iBoxY2 = P_ic[j].Pos[1];
	    iBoxZ1 = P_ic[j].Pos[2];
	    iBoxZ2 = P_ic[j].Pos[2];
	  } else {
	    if (P_ic[j].Pos[0] < iBoxX1)
	      iBoxX1 = P_ic[j].Pos[0];
	    if (P_ic[j].Pos[0] > iBoxX2)
	      iBoxX2 = P_ic[j].Pos[0];
	    if (P_ic[j].Pos[1] < iBoxY1)
	      iBoxY1 = P_ic[j].Pos[1];
	    if (P_ic[j].Pos[1] > iBoxY2)
	      iBoxY2 = P_ic[j].Pos[1];
	    if (P_ic[j].Pos[2] < iBoxZ1)
	      iBoxZ1 = P_ic[j].Pos[2];
	    if (P_ic[j].Pos[2] > iBoxZ2)
	      iBoxZ2 = P_ic[j].Pos[2];
	  }
	  break;
	}
      }
      counter++;
    }
  }

  printf(" Values for the bounding box in the initial conditions:\n");
  printf(" x1: %g x2: %g\n", iBoxX1, iBoxX2);
  printf(" y1: %g y2: %g\n", iBoxY1, iBoxY2);
  printf(" z1: %g z2: %g\n", iBoxZ1, iBoxZ2);
  printf(" x center: %g\n", iBoxX1+(iBoxX2-iBoxX1)/2);
  printf(" y center: %g\n", iBoxY1+(iBoxY2-iBoxY1)/2);
  printf(" z center: %g\n", iBoxZ1+(iBoxZ2-iBoxZ1)/2);
  printf(" x range: %g\n", iBoxX2-iBoxX1);
  printf(" y range: %g\n", iBoxY2-iBoxY1);
  printf(" z range: %g\n", iBoxZ2-iBoxZ1);
  printf(" total number of particles in bbox: %i\n",counter);

  return 0;
}
