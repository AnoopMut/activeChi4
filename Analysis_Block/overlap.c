#include "global.h"
#include "subroutine.h"

void overlap_correlation_function(double factor, double parameter, int is_avg_req ,double CellCut, int scale){
	FILE *relax_file;
	char relax_file_name[128];
	int k,i,current,other,origin_number,offset,maximal_interval;
	int current_frame,other_frame,time_index,max_time_index;
	int num_of_origins,num_of_intervals;
	int **counters, *countersT;
	double **c, **delta_t, **X4, *X4T, *cT;
	double dx,dy,dz,r2;
//	sphere =	calloc(Ns,sizeof(SPHERE)); 
	//count number of origins
	//printf("HERE\n");	
	current = 0; origin_number = 0;
	do{
		while (current<num_of_frames && frame_index_value[current] != origin_number*origin+1)
			current++;
		//printf("originNumber is %d\n",originNumber);
		if (current<num_of_frames)
			origin_number++;
	}while(current<num_of_frames);

	maximal_interval = origin_number*origin/DIVIDE_TO;
	num_of_intervals = 10+(int)( log((double)maximal_interval)/log(factor) ); //just to make sure

	if(is_avg_req==1){ 
		for(i=0;i<2;i++){
			avg_relax[i]= (double *)malloc(sizeof(double)*num_of_intervals);
			for(k=0 ; k<num_of_intervals;k++) {
				avg_relax[i][k] = 0.0;
			}
		}
	}

	int Ncellx=(int)(L/CellCut);
	int Ncell=Ncellx*Ncellx*Ncellx;
	int Ncellxy=Ncellx*Ncellx;
	double InvCellSize=(double)Ncellx;
	int NumInCell[Ncell];
	int Cell[Ns];
	maximal_interval = origin_number*origin/DIVIDE_TO;
	num_of_intervals = 10+(int)( log((double)maximal_interval)/log(factor) ); 
	c = (double **)malloc(sizeof(double *)*Ncell);
	X4 = (double **)malloc(sizeof(double *)*Ncell);
	counters = (int **)malloc(sizeof(int *)*Ncell);
	delta_t = (double **)malloc(sizeof(double *)*Ncell);

	cT = (double *)malloc(sizeof(double)*num_of_intervals);
	X4T = (double *)malloc(sizeof(double)*num_of_intervals);
	countersT = (int *)malloc(sizeof(int)*num_of_intervals);

	for(int i=0;i<Ncell;i++){
		c[i] = (double *)malloc(sizeof(double)*num_of_intervals);
		X4[i] = (double *)malloc(sizeof(double)*num_of_intervals);
		counters[i] = (int *)malloc(sizeof(int)*num_of_intervals);
		delta_t[i] = (double *)malloc(sizeof(double)*num_of_intervals);
	}
	for(int i=0;i<Ncell;i++){
		for (k=0 ; k<num_of_intervals ; k++){
			c[i][k] = 0.0;
			X4[i][k] = 0.0;
			counters[i][k] = 0;
			delta_t[i][k] = 0.0;
		}
	}	
	for (k=0 ; k<num_of_intervals ; k++){
		cT[k] = 0.0;
		X4T[k] = 0.0;
		countersT[k] = 0;
	}
	double sum[Ncell];
	current = 0; origin_number = 0;
	max_time_index = 0;
	do{
		current_frame = frame_index_value[current];
		load_frame(current,0);				//0 to load data for t=0 slot
		other = current;
		offset = 1; //the smallest interval
		time_index = 1;
		//Lets Distribute them to cells
		for (int i=0;i<Ncell;i++)
			NumInCell[i]=0;
		for(int i=0; i<Ns; i++){
			while (sphere[i].x0 >= LENGTH) sphere[i].x0 -=LENGTH;
			while (sphere[i].x0 < 0.0 )    sphere[i].x0 +=LENGTH;
			while (sphere[i].y0 >= LENGTH) sphere[i].y0 -=LENGTH;
			while (sphere[i].y0 < 0.0 )    sphere[i].y0 +=LENGTH;
			while (sphere[i].z0 >= LENGTH) sphere[i].z0 -=LENGTH;
			while (sphere[i].z0 < 0.0 )    sphere[i].z0 +=LENGTH;            
		}

		
		for(int i=0; i<Ns; i++){
        	int cell_index = floor(sphere[i].x0*InvCellSize) + floor(sphere[i].y0*InvCellSize)*Ncellx +
                        floor(sphere[i].z0*InvCellSize)*Ncellxy;
			Cell[i]=cell_index;
            NumInCell[cell_index]++;
        }
		while (other < num_of_frames && offset < maximal_interval){
			other_frame = current_frame + offset; //this is the frame to look for
			//look for that frame
			while (other < num_of_frames && frame_index_value[other] != other_frame)
				other++;
			if (other < num_of_frames){
				//at this point other and current have the indices of the frames.
				//printf("comparing %d with %d	(indices %d with %d)\n",frameIndexArray[current],frameIndexArray[other],current,other);
				load_frame(other,1);
				for(int i=0; i<Ns; i++){
					while (sphere[i].x >= LENGTH) sphere[i].x -=LENGTH;
					while (sphere[i].x < 0.0 )    sphere[i].x +=LENGTH;
					while (sphere[i].y >= LENGTH) sphere[i].y -=LENGTH;
					while (sphere[i].y < 0.0 )    sphere[i].y +=LENGTH;
					while (sphere[i].z >= LENGTH) sphere[i].z -=LENGTH;
					while (sphere[i].z < 0.0 )    sphere[i].z +=LENGTH;            
				}	

                for(i=0;i<Ncell;i++){
                    sum[i]=0.0;
                }
				for (i=0 ; i<Ns ; i++){
					dx = sphere[i].x - sphere[i].x0;
					dy = sphere[i].y - sphere[i].y0;
					dz = sphere[i].z - sphere[i].z0;
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
						sum[Cell[i]] += 1.0;
				}
                for(i=0;i<Ncell;i++){
                    c[i][time_index] += (sum[i])/(double)NumInCell[i];
                    X4[i][time_index] += (sum[i]/(double)NumInCell[i])*(sum[i]/(double)NumInCell[i]);
                    counters[i][time_index]++;
                    delta_t[i][time_index] = time_step*(double)(other_frame - current_frame);
                    cT[time_index] += (sum[i])/(double)NumInCell[i];
				//	X4T[time_index] += (sum[i]/(double)NumInCell[i])*(sum[i]/(double)NumInCell[i]);
                    countersT[time_index]++;
                }
				if ( (int)(offset*factor) == offset )
					offset++;
				else
					offset = (int)(((double)offset)*factor);
				time_index++;
				if(max_time_index < time_index)
					max_time_index = time_index;
			}
		}
		//go to next origin
		origin_number++;
		while (current<num_of_frames && frame_index_value[current] != origin_number*origin)
			current++;
		//printf("originNumber is %d\n",originNumber);
	} while (current<num_of_frames);
	
	/* -----------------------------------------------------------------*/

//	Averaging over blocks
	double BlkAvgc[num_of_intervals];
	double BlkAvgX4[num_of_intervals];
	int BlkAvgCounter[num_of_intervals];

	for(int j=1;j<num_of_intervals;j++){
		BlkAvgc[j]=0.0;
		BlkAvgX4[j]=0.0;
		BlkAvgCounter[j]=0.0;
		for(i=0;i<Ncell;i++){
			if(counters[i][j]){
				BlkAvgc[j]+=c[i][j]/(double)counters[i][j];
				BlkAvgX4[j]+=(X4[i][j]/(double)counters[i][j] - c[i][j]*c[i][j]/((double)(counters[i][j]*counters[i][j])))*NumInCell[i];
				BlkAvgCounter[j]++;
			X4T[j]+=X4[i][j];
			}
		}
	}
	double avgNumInCell=0.0;
	for(i=0;i<Ncell;i++){
//		sprintf(relax_file_name,"%s/BlockAnalysis/Relax/%d/RelaxBlockNum%d.dat",sub_dir,scale,i);
//		relax_file = fopen(relax_file_name,"wb");
//		for(int j=1;j<num_of_intervals;j++){
//			if(counters[i][j]){
//				fprintf(relax_file,"%.14lf\t%.14lf\t%.14lf\n",delta_t[0][j],c[i][j]/(double)counters[i][j],(X4[i][j]/(double)counters[i][j] - c[i][j]*c[i][j]/((double)(counters[i][j]*counters[i][j])))*NumInCell[i]);			
//			}

//		}
//		fclose(relax_file);
		avgNumInCell+=(double)NumInCell[i];
	}
	avgNumInCell/=(double)Ncell;
	printf("%lf\n",avgNumInCell);
	sprintf(relax_file_name,"%s/BlockAnalysis/Relax/RelaxLb%.3lf.dat",sub_dir,L/InvCellSize);
	relax_file = fopen(relax_file_name,"wb");
//	fprintf(relax_file,"0.0\t1.0\t0.0\n");
	sprintf(relax_file_name,"%s/BlockAnalysis/Relax/RelaxTotLb%.3lf.dat",sub_dir,L/InvCellSize);
	FILE *c_file;
	c_file = fopen(relax_file_name,"wb");
	for (k=1;k<num_of_intervals;k++){
		if ( BlkAvgCounter[k] )
			fprintf(relax_file,"%.14lf\t%.14lf\t%.14lf\n",delta_t[0][k],BlkAvgc[k]/(double)BlkAvgCounter[k],BlkAvgX4[k]/(double)BlkAvgCounter[k]);
		if ( countersT[k] )
			fprintf(c_file,"%.14lf\t%.14lf\t%.14lf\n",delta_t[0][k],
                                    cT[k]/(double)countersT[k], 
                                    avgNumInCell*(X4T[k]/(double)countersT[k]-cT[k]*cT[k]/((double)countersT[k]*countersT[k]) ));
	}
	fclose(relax_file);
	fclose(c_file);
/*	for (k=1;k<num_of_intervals;k++){
		if ( counters[k] ){
			avg_relax[1][k] +=BlkAvgMSD[k]/(double)BlkAvgCounter[k]; 
			avg_relax[2][k] +=BlkAvgMSDSq[k]/(double)BlkAvgCounter[k];
			avg_msd[0][k] =delta_t[0][k];
		}
	}	
	*/for(i=0;i<Ncell;i++)
	{	free(c[i]);free(X4[i]);free(delta_t[i]);free(counters[i]);}
	free(c);free(X4);free(delta_t);free(counters);
//	free(sphere);
	return;
}
