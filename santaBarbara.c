//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "snapshot_io.h"
#include "data_ops.h"


/*inline*/ float dist(float *a, float *b)
{
	float d[3] = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
	return sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
}


void check_omega(particle_data *P, int NumPart, const float BoxSize)
{
	
  //const float BoxSize = 100000;
  const float Hubble = 0.1;
  const float G = 43007.1;
  const float Omega0 = 0.276;
  
  double mass = 0, masstot = 0, omega;
  int i;

  for(i = 0; i < NumPart; i++)
    masstot += P[i].Mass;

  //MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  omega =
    masstot / (BoxSize * BoxSize * BoxSize) / (3 * Hubble * Hubble / (8 * M_PI * G));

  //if(fabs(omega - Omega0) > 1.0e-3)
    {
     // if(ThisTask == 0)
	{
	  printf("\n\nI've found something odd!\n");
	  printf
	    ("The mass content accounts only for Omega=%g,\nbut you specified Omega=%g in the parameterfile.\n",
	     omega, Omega0);
	  printf("\nI better stop.\n");

	  fflush(stdout);
	}
      //endrun(1);
    }
}


int main(int argc, char **argv)
{
	//load file
	particle_data *P, *p, *pi;
	io_header header;
	int tot;
	
	int i, j, k, format;

	const float a = 0.015625;
	//const float a = 1.0;
	const float sq_a = sqrt(a);
	const float h = 0.5;
	const float h_inv = 1.0/h;
	const float h_inv_sq = h_inv*h_inv;
	
	//int j, start, end, type, tot;
		
	if(argc < 4)
	{
		fprintf(stderr, "usage: %s <IC file> <num part.> <out file>\n", argv[0]);
		exit(-1);
	}
	
	tot = atoi(argv[2]);

	
	//format x[Mpc], vx[km/s], y, vy, z, vz, mass[M_s], tag //single precision
	FILE *f = fopen(argv[1], "r");
	P = (particle_data*)calloc(tot, sizeof(particle_data));
	for(j = 0; j < tot; ++j)
	{
		float data[8];
				
		fread(data, sizeof(float), 8, f);
		
		for(k = 0; k < 3; ++k)
		{
			P[j].Pos[k] = data[k*2]*1000 * h;
			P[j].Vel[k] = data[k*2 + 1] * sq_a;
		}
		P[j].Mass = data[6]*1e-10 * h;
		//P[j].Id = (int)data[7];
		P[j].Id = j;
		P[j].Type = 1;
		
		if(j < 10)
			printf("Id: %d %g %g %g %g %g %g %g\n", P[j].Id, P[j].Mass, P[j].Pos[0], P[j].Pos[1], P[j].Pos[2], P[j].Vel[0], P[j].Vel[1], P[j].Vel[2]);
	}
	fclose(f);
	
	//test
	check_omega(P, tot, 64000*h);
	
	//set header to 0
	memset(&header, 0, sizeof(header));
	
	generateHeader(&header, P, tot, 1);
	write_snapshot(argv[3], &header, P);
	
	return 0;
}

