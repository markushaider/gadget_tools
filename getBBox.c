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
  float cX = atof(argv[3]);
  float cY = atof(argv[4]);
  float cZ = atof(argv[5]);
  float lRefined = atof(argv[6]);

  if(argc!=7) {
    printf("Usage: ./getBBox <snapshot0> <snapshot_final> <x> <y> <z> <extent>\n");
    printf("<x> <y> <z> refer to the coordinates of the center of the desired volume (in kpc)\n");
    printf("extent is the size of the desired volume (in kpc)\n");
    return 0;
  }
  printf(" You specified\n");
  printf(" <x>: %g\n",cX);
  printf(" <y>: %g\n",cY);
  printf(" <z>: %g\n",cZ);
  printf(" extent: %g\n",lRefined);
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

  printf("\n Information about the initial conditions: \n");
  printf("Boxsize: %g\n", header_ic.BoxSize);
  printf("Number of gas particles: %i\n", header_ic.npart[0]);
  printf("Number of high-res dm particles: %i\n", header_ic.npart[1]);
  printf("Number of low-res dm particles: %i\n", header_ic.npart[5]);
  printf("Number of star particles: %i\n", header_ic.npart[4]);

  /* find the bounding box at the beginning */
  int i,j;

  float xMin, xMax, yMin, yMax, zMin, zMax;
  float halfLength = lRefined/2.;
  float halfBoxSize = header.BoxSize/2;


  /* shift the coordinates around center */
  printf("Shift coordinates around center\n");
  for(i=0; i< header.npart[0]+header.npart[1]; i++) {
    if(P[i].Pos[0]-cX < -halfBoxSize) {
      P[i].Pos[0] = header.BoxSize+P[i].Pos[0]-cX;
    }
    else if(P[i].Pos[0]-cX > halfBoxSize) {
      P[i].Pos[0] = -header.BoxSize+P[i].Pos[0]-cX;
    }
    else {
      P[i].Pos[0] -= cX;
    }
    if(P[i].Pos[1]-cY < -halfBoxSize) {
      P[i].Pos[1] = header.BoxSize+P[i].Pos[1]-cY;
    }
    else if(P[i].Pos[1]-cY > halfBoxSize) {
      P[i].Pos[1] = -header.BoxSize+P[i].Pos[1]-cY;
    }
    else {
      P[i].Pos[1] -= cY;
    }
    if(P[i].Pos[2]-cZ < -halfBoxSize) {
      P[i].Pos[2] = header.BoxSize+P[i].Pos[2]-cZ;
    }
    else if(P[i].Pos[2]-cZ > halfBoxSize) {
      P[i].Pos[2] = -header.BoxSize+P[i].Pos[2]-cZ;
    }
    else {
      P[i].Pos[2] -= cZ;
    }
  }
  printf("Done shifting coordinates\n");

  printf("Shift coordinates of IC around center\n");
  for(i=0; i< header_ic.npart[0]+header_ic.npart[1]; i++) {
    if(P_ic[i].Pos[0]-cX < -halfBoxSize) {
      P_ic[i].Pos[0] = header_ic.BoxSize+P_ic[i].Pos[0]-cX;
    }
    else if(P_ic[i].Pos[0]-cX > halfBoxSize) {
      P_ic[i].Pos[0] = -header_ic.BoxSize+P_ic[i].Pos[0]-cX;
    }
    else {
      P_ic[i].Pos[0] -= cX;
    }
    if(P_ic[i].Pos[1]-cY < -halfBoxSize) {
      P_ic[i].Pos[1] = header_ic.BoxSize+P_ic[i].Pos[1]-cY;
    }
    else if(P_ic[i].Pos[1]-cY > halfBoxSize) {
      P_ic[i].Pos[1] = -header_ic.BoxSize+P_ic[i].Pos[1]-cY;
    }
    else {
      P_ic[i].Pos[1] -= cY;
    }
    if(P_ic[i].Pos[2]-cZ < -halfBoxSize) {
      P_ic[i].Pos[2] = header_ic.BoxSize+P_ic[i].Pos[2]-cZ;
    }
    else if(P_ic[i].Pos[2]-cZ > halfBoxSize) {
      P_ic[i].Pos[2] = -header_ic.BoxSize+P_ic[i].Pos[2]-cZ;
    }
    else {
      P_ic[i].Pos[2] -= cZ;
    }
  }
  printf("Done shifting coordinates of IC\n");

  int counter=0;
  int found=0;
  int notfound=0;
  for(i=0; i< header.npart[0]+header.npart[1]; i++) {
    /* check whether particle is within the definded bounding box */
    if (P[i].Pos[0] < halfLength && P[i].Pos[0] > -halfLength &&
	P[i].Pos[1] < halfLength && P[i].Pos[1] > -halfLength &&
	P[i].Pos[2] < halfLength && P[i].Pos[2] > -halfLength) {
      /* if particle is inside, then get the position of the particle in the ic */
      for(j=0; j< header_ic.npart[0]+header_ic.npart[1]; j++) {

	if(P_ic[j].Id == P[i].Id) {
	  if(counter==0) {
	    xMin = P_ic[j].Pos[0];
	    xMax = P_ic[j].Pos[0];
	    yMin = P_ic[j].Pos[1];
	    yMax = P_ic[j].Pos[1];
	    zMin = P_ic[j].Pos[2];
	    zMax = P_ic[j].Pos[2];

	    found++;
	    break;
	  }
	  found++;
	  if (P_ic[j].Pos[0] < xMin)
	    xMin = P_ic[j].Pos[0];
	  if (P_ic[j].Pos[0] > xMax)
	    xMax = P_ic[j].Pos[0];
	  if (P_ic[j].Pos[1] < yMin)
	    yMin = P_ic[j].Pos[1];
	  if (P_ic[j].Pos[1] > yMax)
	    yMax = P_ic[j].Pos[1];
	  if (P_ic[j].Pos[2] < zMin)
	    zMin = P_ic[j].Pos[2];
	  if (P_ic[j].Pos[2] > zMax)
	    zMax = P_ic[j].Pos[2];
	  break;
	}
      }
      if(j==header_ic.npart[0]+header_ic.npart[1]) {
	notfound++;
      }
      counter++;
    }
  }
  if(notfound!=0) {
    printf(" Problem! There were particles which have not been found in the ICs\n");
  }
  printf(" total number of particles in bbox: %i\n",counter);
  printf(" Bounding Box coordinates in shifted coordinates:\n");
  printf(" xMin: %g xMax: %g\n", xMin, xMax);
  printf(" yMin: %g yMax: %g\n", yMin, yMax);
  printf(" zMin: %g zMax: %g\n", zMin, zMax);
  printf(" x-range: %g\n", xMax-xMin);
  printf(" y-range: %g\n", yMax-yMin);
  printf(" z-range: %g\n", zMax-zMin);
  printf(" Bounding Box coordinates in real coordinates:\n");
  printf(" xMin: %g xMax: %g\n", xMin+cX, xMax+cX);
  printf(" yMin: %g yMax: %g\n", yMin+cY, yMax+cY);
  printf(" zMin: %g zMax: %g\n", zMin+cZ, zMax+cZ);
  printf(" center of the region in the ics:\n");
  printf(" x-center: %g\n", xMin+cX + (xMax-xMin)/2);
  printf(" y-center: %g\n", yMin+cY + (yMax-yMin)/2);
  printf(" z-center: %g\n", zMin+cZ + (zMax-zMin)/2);
  printf(" specified center:\n");
  printf(" x-center at z=0: %g\n", cX);
  printf(" y-center at z=0: %g\n", cY);
  printf(" z-center at z=0: %g\n", cZ);
  printf(" corrected bounding box:\n");
  xMin = xMin+cX;
  xMax = xMax+cX;
  yMin = yMin+cY;
  yMax = yMax+cY;
  zMin = zMin+cZ;
  zMax = zMax+cZ;
  if(xMax<cX+halfLength)
    xMax=cX+halfLength;
  if(xMin>cX-halfLength)
    xMin=cX-halfLength;
  if(yMax<cY+halfLength)
    yMax=cY+halfLength;
  if(yMin>cY-halfLength)
    yMin=cY-halfLength;
  if(zMax<cZ+halfLength)
    zMax=cZ+halfLength;
  if(zMin>cZ-halfLength)
    zMin=cZ-halfLength;
  printf(" xMin: %g xMax: %g\n", xMin, xMax);
  printf(" yMin: %g yMax: %g\n", yMin, yMax);
  printf(" zMin: %g zMax: %g\n", zMin, zMax);


  return 0;
}
