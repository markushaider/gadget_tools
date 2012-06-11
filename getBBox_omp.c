#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <omp.h>

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

  /* find the bounding box at the beginning */
  int i,j;
  int counter=0;
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

  /* set first values */
  xMin = P[0].Pos[0];
  xMax = P[0].Pos[0];
  yMin = P[0].Pos[1];
  yMax = P[0].Pos[1];
  zMin = P[0].Pos[2];
  zMax = P[0].Pos[2];
  
  float xMinR, xMaxR, yMinR, yMaxR, zMinR, zMaxR;
  xMinR = xMin;
  xMaxR = xMax;
  yMinR = yMin;
  yMaxR = yMax;
  zMaxR = zMax;
  zMinR = zMin;
  int th_id;
  printf("Finding particle in initial condition\n");
  #pragma omp parallel private(xMin,xMax,yMin,yMax,zMin,zMax)
  {
    //th_id = omp_get_thread_num();
    //printf("Hello World from thread %d\n", th_id);
    #pragma omp for private(i) reduction(+:counter)
    for(i=1; i< header.npart[0]+header.npart[1]; i++) {
      /* check whether particle is within the definded bounding box */
      if (P[i].Pos[0] < halfLength && P[i].Pos[0] > -halfLength &&
	  P[i].Pos[1] < halfLength && P[i].Pos[1] > -halfLength &&
	  P[i].Pos[2] < halfLength && P[i].Pos[2] > -halfLength) {
	/* if particle is inside, then get the position of the particle in the ic */
	for(j=0; j< header.npart[0]+header.npart[1]; j++) {

	  if(P_ic[j].Id == P[i].Id) {
	    if (P[i].Pos[0] < xMin)
	      xMin = P[i].Pos[0];
	    if (P[i].Pos[0] > xMax)
	      xMax = P[i].Pos[0];
	    if (P[i].Pos[1] < yMin)
	      yMin = P[i].Pos[1];
	    if (P[i].Pos[1] > yMax)
	      yMax = P[i].Pos[1];
	    if (P[i].Pos[2] < zMin)
	      zMin = P[i].Pos[2];
	    if (P[i].Pos[2] > zMax)
	      zMax = P[i].Pos[2];
	    break;
	  }
	}
	counter++;
      }
    }

    if (xMin<xMinR) {
      #pragma critical
      {
        if ( xMin < xMinR) xMinR = xMin;
      }
    }
    if (xMax>xMaxR) {
      #pragma critical
      {
        if ( xMax > xMaxR) xMaxR = xMax;
      }
    }
    if (yMin<yMinR) {
      #pragma critical
      {
        if ( yMin < yMinR) yMinR = yMin;
      }
    }
    if (yMax>yMaxR) {
      #pragma critical
      {
        if ( yMax > yMaxR) yMaxR = yMax;
      }
    }
    if (zMin<zMinR) {
      #pragma critical
      {
        if ( zMin < zMinR) zMinR = zMin;
      }
    }
    if (zMax>xMaxR) {
      #pragma critical
      {
        if ( zMax > zMaxR) zMaxR = zMax;
      }
    }



  }

  printf(" Values for the distance to center:\n");
  printf(" x1: %g x2: %g\n", xMinR, xMaxR);
  printf(" y1: %g y2: %g\n", yMinR, yMaxR);
  printf(" z1: %g z2: %g\n", zMinR, zMaxR);
  printf(" x range: %g\n", xMaxR-xMinR);
  printf(" y range: %g\n", yMaxR-yMinR);
  printf(" z range: %g\n", zMaxR-zMinR);
  printf(" total number of particles in bbox: %i\n",counter);

  return 0;
}
