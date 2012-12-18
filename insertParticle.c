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
    printf("Usage: ./cutout <infile> <outfile> <x coordinate of clump [Mpc]> <y coordinate> <z coordinate>\n");
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

  float Box; 
  Box = (float)header.BoxSize;
 
  /* insert 1 particle at x=3.5 Mpc */
  P[nTot].Pos[0] = Box/2.0+3500.0;
  P[nTot].Pos[1] = 0.0;
  P[nTot].Pos[2] = 0.0;
  P[nTot].Vel[0] = 0.0;
  P[nTot].Vel[1] = 0.0;
  P[nTot].Vel[2] = 0.0;
  P[nTot].Type=5;


  int counterType [6] = {0,0,0,0,0,0};
  /* get the numbers of the particle types */
  int i;
  for (i=0;i<=nTot;i++) {
    if(P[i].Type==0) {
      counterType[0]+=1;
    }
    if(P[i].Type==1) {
      counterType[1]+=1;
    }
    if(P[i].Type==2) {
      counterType[2]+=1;
    }
    if(P[i].Type==3) {
      counterType[3]+=1;
    }
    if(P[i].Type==4) {
      counterType[4]+=1;
    }
    if(P[i].Type==5) {
      counterType[5]+=1;
    }
  }

  header.npart[0]=counterType[0];
  header.npartTotal[0]=counterType[0];
  header.npart[1]=counterType[1];
  header.npartTotal[1]=counterType[1];
  header.npart[4]=counterType[4];
  header.npartTotal[4]=counterType[4];
  header.npart[5]=counterType[5];
  header.npartTotal[5]=counterType[5];
  header.mass[5]=1E2;

  /* output information about the simulation */
  printf(">> Output Cutout Box\n");
  printf("Boxsize: %g\n", header.BoxSize);
  printf("Number of gas particles: %i\n", header.npart[0]);
  printf("Number of high-res dm particles: %i\n", header.npart[1]);
  printf("Number of low-res dm particles: %i\n", header.npart[5]);
  printf("Number of star particles: %i\n", header.npart[4]);
 
  /* write reduced snapshot file */
  write_snapshot(outfile, &header, P);

  return 0;
}
