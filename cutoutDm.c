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
  float gasPMass=atof(argv[3]);
  if(argc!=4) {
    printf("Usage: ./cutoutDm <infile> <outfile> <gas particle mass>\n");
    return 0;
  }
  printf("Got %g for gas particle mass\n",gasPMass);


  /* load snapshot file */
  io_header header;
  particle_data * P;
  int nTot;
  nTot = load_snapshot(infile, &header, &P);

  /* output information about the simulation */
  printf("--------- Input File --------- \n");
  printf("Boxsize: %g\n", header.BoxSize);
  printf("Number of gas particles: %i\n", header.npart[0]);
  printf("Number of high-res dm particles: %i\n", header.npart[1]);
  printf("Number of low-res dm particles: %i\n", header.npart[5]);
  printf("Number of star particles: %i\n", header.npart[4]);

  int i;
  /* throw everything away except for the fine dark matter particles */
  particle_data * P_out;
  P_out = (particle_data *)calloc(header.npart[1],sizeof(particle_data));
  if(!P_out) {
      fprintf(stderr,"failed to allocate memory.\n");
      return -1;
  }
  /* use only the fine dm particles */
  for (i=0;i<header.npart[1];i++) {
    P_out[i] = P[i+header.npart[0]];
  }

  /* correct the masses */
  printf("Changing DM mass from %g to %g (+%g)\n",header.mass[1],header.mass[1]+gasPMass,gasPMass);
  header.mass[1]+=gasPMass; /* add gas mass to dm mass in order to comply with rockstar */
  header.mass[0]=0;
  header.mass[4]=0;
  header.mass[5]=0;


  /* shift the coordinates */
  /* not sure wheather this makes sense */
  float extent=20000;
  float center = header.BoxSize/2;
  for (i=0;i<nTot;i++) {
      P_out[i].Pos[0]-= center-extent/2;
      P_out[i].Pos[1]-= center-extent/2;
      P_out[i].Pos[2]-= center-extent/2;
  }

  header.npart[0]=0;
  header.npartTotal[0]=0;
  header.npart[4]=0;
  header.npartTotal[4]=0;
  header.npart[5]=0;
  header.npartTotal[5]=0;
  
  /* write reduced snapshot file */
  write_snapshot(outfile, &header, P_out);

  /* output information about the simulation */
  printf("--------- Output File --------- \n");
  printf("Boxsize: %g\n", header.BoxSize);
  printf("Number of gas particles: %i\n", header.npart[0]);
  printf("Number of high-res dm particles: %i\n", header.npart[1]);
  printf("Number of low-res dm particles: %i\n", header.npart[5]);
  printf("Number of star particles: %i\n", header.npart[4]);


  return 0;
}
