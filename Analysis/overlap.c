#include "global.h"
#include "subroutine.h"

void overlap_correlation_function(double factor, double PARAMETER, int is_avg_req, int flag ){
	double parameter=PARAMETER*PARAMETER;
	char c_file_name[128];
	int k,i,current,other,origin_number,offset,maximal_interval;
	int current_frame,other_frame,time_index,max_time_index;
	int num_of_origins,num_of_intervals;
	int *counters;
	double *c, *delta_t, *X4;
	double sum,dx,dy,dz,r2;
	sphere =	calloc(NPP,sizeof(SPHERE)); 
	num_of_intervals=num_of_frames;
	double SUM;
	current = 0; origin_number = 0;
	do{
		while (current<num_of_frames && frame_index_value[current] != origin_number*origin+1)
			current++;
		//printf("originNumber is %d\n",originNumber);
		if (current<num_of_frames)
			origin_number++;
	}while(current<num_of_frames);

	double TotPart=(double)N;
	if (flag==1) TotPart=0.8*N;
	
	maximal_interval = origin_number*origin/DIVIDE_TO;
	num_of_intervals = 10+(int)( log((double)maximal_interval)/log(factor) ); //just to make sure
	if(node==0){
		c = (double *)malloc(sizeof(double)*num_of_intervals);
		X4 = (double *)malloc(sizeof(double)*num_of_intervals);
		counters = (int *)malloc(sizeof(int)*num_of_intervals);
		delta_t = (double *)malloc(sizeof(double)*num_of_intervals);
		for (k=0;k<num_of_intervals;k++){
			c[k] = 0.0;
			X4[k] = 0.0;
			counters[k] = 0;
			delta_t[k] = 0.0;
		}
	}
	current = 0; origin_number = 0;
	max_time_index = 0;

	do{
        current_frame = frame_index_value[current];
		LoadFramePP(current,0);
		//printf("loaded frame number %d which is %d\n",current,node);
		other = current;
		offset = 1; //the smallest interval
		time_index = 1;
		while (other < num_of_frames && offset < maximal_interval){
			other_frame = current_frame + offset; //this is the frame to look for
			//look for that frame
			while (other < num_of_frames && frame_index_value[other] != other_frame)
				other++;
			if (other < num_of_frames){
				LoadFramePP(other,1);
    			MPI_Barrier(MPI_COMM_WORLD);
				sum=0.0;
				for (i=0 ; i<NPP ; i++){
					if(sphere[i].type!=1 || flag==0){
						dx = sphere[i].x- sphere[i].x0;
						dy = sphere[i].y- sphere[i].y0;
						dz = sphere[i].z- sphere[i].z0;
						// mess due to periodic boundary conditions 
    	        	    if ( dz >= HALF_LENGTH )
    	        	        dz -= LENGTH;
    	        	    else if ( dz < -HALF_LENGTH )
    	        	        dz += LENGTH;
    	        	    if ( dy >= HALF_LENGTH )
    	        	        dy -= LENGTH;
    	        	    else if ( dy < -HALF_LENGTH )
    	        	        dy += LENGTH;
    	        	    if ( dx >= HALF_LENGTH )
    	        	        dx -= LENGTH;
    	        	    else if ( dx < -HALF_LENGTH )
    	        	        dx += LENGTH;
						// end of mess 
						r2 = L*L*( dx*dx + dy*dy + dz*dz);
						if (r2 < parameter)
							sum += 1.0;
					}
				}
				sum = sum/TotPart;
 				SUM=0.0;
		  		MPI_Reduce(&sum, &SUM,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   				MPI_Barrier(MPI_COMM_WORLD);
				if(node==0) {
					c[time_index] += SUM;
					X4[time_index] += SUM*SUM;
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
		origin_number++;
		while (current<num_of_frames && frame_index_value[current] != origin_number*origin)
			current++;
	} while (current<num_of_frames);
	/* -----------------------------------------------------------------*/
	if(node==0){
		FILE *c_file;
		sprintf(c_file_name,"%s/relax_a_%.2lf.dat",sub_dir,PARAMETER);
		if(flag==1) sprintf(c_file_name,"%s/relaxAA_a_%.2lf.dat",sub_dir,PARAMETER);
		c_file = fopen(c_file_name,"wb");
		fprintf(c_file,"0.0\t1.0\t0.0\n");
		for (k=1;k<num_of_intervals;k++){
			if ( counters[k] )
				fprintf(c_file,"%.14lf\t%.14lf\t%.14lf\n",delta_t[k],
									c[k]/(double)counters[k],
									X4[k]/(double)counters[k] - c[k]*c[k]/((double)(counters[k]*counters[k])));
		}
		fclose(c_file);
		free(c); free(X4);
		free(delta_t); free(counters);
	printf("Overlap Calculated\n");
	}

	free(sphere);
	return;
}
