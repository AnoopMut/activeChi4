#include "global.h"
#include "subroutine.h"


void mean_square_displacement(double factor, int is_avg_req,int flag){
	FILE *msd_file;
	char msd_file_name[128];
	int k,i,current,other,origin_number,offset,maximal_interval;
	int current_frame,other_frame,time_index;
	int num_of_origins,num_of_intervals;
	int *counters;
	double *msd[3], *delta_t;
	double sum,dx,dy,dz,r2,r4,r3,sumSq,sumCube;			//saperate sum for rod and sphere
	double SUM,SUMSQ,SUMCUBE;	
	sphere =	calloc(NPP,sizeof(SPHERE)); 
	//count number of origins
	current = 0; origin_number = 0;
	do{
		while (current<num_of_frames && frame_index_value[current] != origin_number*origin)
			current++;
		if (current<num_of_frames)
			origin_number++;
	}while(current<num_of_frames);
	double TotPart=(double)N;
	if (flag==1) TotPart=0.8*N;

	
	maximal_interval = origin_number*origin/DIVIDE_TO;
	num_of_intervals = 10+(int)( log((double)maximal_interval)/log(factor) )+10; 
	if(node==0){
		msd[0] = (double *)malloc(sizeof(double)*num_of_intervals);
		msd[1] = (double *)malloc(sizeof(double)*num_of_intervals);
		msd[2] = (double *)malloc(sizeof(double)*num_of_intervals);
		counters = (int *)malloc(sizeof(int)*num_of_intervals);
		delta_t = (double *)malloc(sizeof(double)*num_of_intervals);
		for (k=0 ; k<num_of_intervals ; k++){
			msd[0][k] = 0.0;
			msd[1][k] = 0.0;
			counters[k] = 0;
			delta_t[k] = 0.0;
		}

	}	
	current = 0; origin_number = 0;

	MPI_Barrier(MPI_COMM_WORLD);
	//printf("HERE\n");
	do{
		current_frame = frame_index_value[current];
		LoadFramePP(current,0);				//0 to load data for t=0 slot
		other = current;
		offset = 1; //the smallest interval
		time_index = 1;
		while (other < num_of_frames && offset < maximal_interval){
			other_frame = current_frame + offset; //this is the frame at later time 
			//look for that frame
			while (other < num_of_frames && frame_index_value[other] != other_frame)
				other++;
			if (other < num_of_frames){
				//at this point other and current have the indices of the frames.
				LoadFramePP(other,1);			//1 to load data for normal slot at time t 
				MPI_Barrier(MPI_COMM_WORLD);
				sum=0.0;sumSq=0.0;sumCube=0.0;
				for (i=0 ; i<NPP ; i++){
					if(sphere[i].type!=1 || flag==0){
						dx = sphere[i].x- sphere[i].x0;
						dy = sphere[i].y- sphere[i].y0;
						dz = sphere[i].z- sphere[i].z0;
						r2=L*L*(dx*dx+dy*dy+dz*dz);	
						r3= r2*r2;
						r4=L*L*L*L*(dx*dx*dx*dx + dy*dy*dy*dy + dz*dz*dz*dz);
						sum += r2;
						sumCube+=r3;
						sumSq+=r4;
					}
				}
				//printf("%lf\n",L);
                SUM=0.0;
                SUMCUBE=0.0;
                SUMSQ=0.0;
                MPI_Reduce(&sum, &SUM,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Reduce(&sumSq, &SUMSQ,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Reduce(&sumCube, &SUMCUBE,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);
				if(node==0){
					msd[0][time_index] += SUM/TotPart;
					msd[1][time_index] += SUMSQ/TotPart;
					msd[2][time_index] += SUMCUBE/TotPart;
					counters[time_index]++;
					delta_t[time_index] = time_step*(double)(other_frame - current_frame);
				}
				//advance offset and timeIndex
				if ( (int)(offset*factor) == offset )
					offset++;
				else
					offset = (int)(((double)offset)*factor);
				time_index++;
			}
		}
		//go to next origin
		origin_number++;
		while (current<num_of_frames && frame_index_value[current] != origin_number*origin)
			current++;
		//printf("origin_number is %d\n",origin_number);
	} while (current < num_of_frames);
    MPI_Barrier(MPI_COMM_WORLD);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if(node==0){	
		sprintf(msd_file_name,"%s/msd.dat",sub_dir);
		if(flag==1) sprintf(msd_file_name,"%s/msdA.dat",sub_dir);
		msd_file = fopen(msd_file_name,"wb");
		fprintf(msd_file,"0.0\t0.0\t0.0\t0.0\n");
		for (k=1;k<num_of_intervals;k++){
			if ( counters[k] )
				fprintf(msd_file,"%.14lf\t%.14lf\t%.14lf\t%.14lf\n",delta_t[k],msd[0][k]/(double)counters[k],msd[1][k]/(double)counters[k],msd[2][k]/(double)counters[k]);
		}
		fclose(msd_file);
		free(msd[0]);free(msd[1]);free(msd[2]);free(delta_t); free(counters);
	}
	free(sphere);
	return;
}
