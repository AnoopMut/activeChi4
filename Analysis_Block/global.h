#ifndef GLOBAL_H
	#define GLOBAL_H
	#include <stdio.h>
	#include <math.h>
	#include <time.h>
	#include <stdlib.h>
	#include <assert.h>
	#include <unistd.h>
	#include <stdbool.h>
	#define Ns				50000								//No of spheres
	#define pi      		    3.14159265358979		
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
	SPHERE sphere[Ns];
	double pe,ke,rot_e,virial;
    double DENSITY,DIAMETER;
	int num_task,task_id,job_no;
    double L,invL ;
	double T,P;
	int serial, origin , num_of_frames;
	int *frame_index_value;
	double *all_data_x;
	double *all_data_y;
	double *all_data_z;
	int *all_data_type;
	double *avg_msd[3];
	double *avg_relax[7];
	int N,nb;
	int tot_origin;
	double tau_alpha;
	char main_dir[128];
	char sub_dir[128];
	char buff[1024];

#endif
