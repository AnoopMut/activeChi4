#ifndef GLOBAL_H
	#define GLOBAL_H
	#include <stdio.h>
	#include <math.h>
	#include <time.h>
	#include <stdlib.h>
	#include <assert.h>
	#include <unistd.h>
	#include <stdbool.h>
	#include "mpi.h"
	#define N					50000								//No of spheres
	#define NA					(int)(.800*(double)Ns)			//No of type A spheres
	#define NB					(int)(.200*(double)Ns)			//No of type B spheres
	#define pi      		    3.14159265358979		
	#define TWO_PI				2.0*pi
	#define Sigma				1.0
	#define LENGTH				1.0								//scaled domain size
	#define HALF_LENGTH			0.5								//half of it
	#define time_step			0.005
   	#define half_time_step      ( time_step*0.5 )
	#define DIVIDE_TO			2.0								//frac of maximal log interval
	#define RHO					1.2
	typedef struct{
		double x , y,z;
		double x0, y0,z0; 											//time 0 coordinates
		int    type,type0;
	}SPHERE;
	SPHERE *sphere;
    double DENSITY;
	int num_task,task_id,job_no;
 //   double L ;
	double T,P;
	int serial, origin , num_of_frames;
	int *frame_index_value;
	double *Allx;
	double *Ally;
	double *Allz;
	int *Alltype;
	double *Alloriginx;
	double *Alloriginy;
	double *Alloriginz;
	int *Allorigintype;
	double *all_data_x;
	double *all_data_y;
	double *all_data_z;
	int *all_data_type;
	double L;
	int node,Nprocs;												//present block
	int NPP;
	int tot_origin;
    double *avg_msd[3];
    double *avg_relax[3];
    int ActTag[N],ActNum;
	double invL,tau_alpha;
	char main_dir[128];
	char sub_dir[128];
	char buff[1024];
	int nq;
	int NT;
	int q[12*23+100];
	int *Vecx, *Vecy,*Vecz;
#endif
