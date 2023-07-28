#include "global.h"
#include "subroutine.h"
void mean_square_displacement(double factor, int is_avg_req, double CellCut){
	FILE *msd_file;
	char msd_file_name[128];
	int k,i,current,other,origin_number,offset,maximal_interval;
	int current_frame,other_frame,time_index;
	int num_of_origins,num_of_intervals;
	int **counters,*CounterBinder;
	double **msd,*msdSqCell,*msdCell, **delta_t;
	double dx,dy,dz,r2;			//saperate sum for rod and sphere
//	sphere =(SPHERE *)malloc(sizeof(SPHERE)*Ns); 
	//count number of origins
	current = 0; origin_number = 0;
	do{
		while (current<num_of_frames && frame_index_value[current] != origin_number*origin+1)
			current++;
		if (current<num_of_frames)
			origin_number++;
	}while(current<num_of_frames);
/*	if(is_avg_req==1){ 
		avg_msd[0]= (double *)malloc(sizeof(double)*num_of_intervals);
		avg_msd[1]= (double *)malloc(sizeof(double)*num_of_intervals);
		avg_msd[2]= (double *)malloc(sizeof(double)*num_of_intervals);
		for(k=0 ; k<num_of_intervals;k++) {
			avg_msd[0][k] =0.0;
			avg_msd[1][k] =0.0;
			avg_msd[2][k] =0.0;
		}
	}
*/
	//Cell Division
	int Ncellx=(int)(L/CellCut);
	int Ncell=Ncellx*Ncellx*Ncellx;
	int Ncellxy=Ncellx*Ncellx;
	printf("# of cells equals %d Length of Block %lf\n",Ncell,L/Ncellx);
	double InvCellSize=(double)Ncellx;
	int NumInCell[Ncell];
	double sum[Ncell];
	double COMx0[Ncell],COMy0[Ncell],COMz0[Ncell];
	double COMx[Ncell],COMy[Ncell],COMz[Ncell];
	int Cell[Ns];
	maximal_interval = origin_number*origin/DIVIDE_TO;
	num_of_intervals = 10+(int)( log((double)maximal_interval)/log(factor) ); 
	msd = (double **)malloc(sizeof(double *)*Ncell);
//	msdCell = (double **)malloc(sizeof(double *)*Ncell);
//	msdSqCell = (double **)malloc(sizeof(double *)*Ncell);
	counters = (int **)malloc(sizeof(int *)*Ncell);
	delta_t = (double **)malloc(sizeof(double *)*Ncell);
	for(int i=0;i<Ncell;i++){
		msd[i] = (double *)malloc(sizeof(double)*num_of_intervals);
//		msdCell[i] = (double *)malloc(sizeof(double)*num_of_intervals);
//		msdSqCell[i] = (double *)malloc(sizeof(double)*num_of_intervals);
		counters[i] = (int *)malloc(sizeof(int)*num_of_intervals);
		delta_t[i] = (double *)malloc(sizeof(double)*num_of_intervals);
	}
		msdCell= (double *)malloc(sizeof(double)*num_of_intervals);
		msdSqCell = (double *)malloc(sizeof(double)*num_of_intervals);
		 CounterBinder = (int *)malloc(sizeof(int)*num_of_intervals);
	for(int i=0;i<Ncell;i++){
		for (k=0 ; k<num_of_intervals ; k++){
			msd[i][k] = 0.0;
			msdCell[k] = 0.0;
			msdSqCell[k] = 0.0;
			counters[i][k] = 0;
			delta_t[i][k] = 0.0;
			CounterBinder[k]=0;
		}
	}	
	current = 0; origin_number = 0;
	do{
		current_frame = frame_index_value[current];
		load_frame(current,0);				//0 to load data for t=0 slot
		other = current;
		offset = 1; //the smallest interval
		time_index = 1;
		//Lets Distribute them to cells
		for (int i=0;i<Ncell;i++){
			NumInCell[i]=0;
			COMx0[i]=0.0;
			COMy0[i]=0.0;
			COMz0[i]=0.0;
		}
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
			if(cell_index>Ncell) printf("Fuck\n%lf\t%lf\t%lf\n",sphere[i].x0,sphere[i].y0,sphere[i].z0);
			Cell[i]=cell_index;
			NumInCell[cell_index]++;
		}	

//		load_frame(current,0);				//0 to load data for t=0 slot
		for(int i=0; i<Ns; i++){
			COMx0[Cell[i]]+=sphere[i].x0;
			COMy0[Cell[i]]+=sphere[i].y0;
			COMz0[Cell[i]]+=sphere[i].z0;
		}
		while (other < num_of_frames && offset < maximal_interval){
			other_frame = current_frame + offset; //this is the frame at later time 
			//look for that frame
			while (other < num_of_frames && frame_index_value[other] != other_frame)
				other++;
			if (other < num_of_frames){
				//at this point other and current have the indices of the frames.
				load_frame(other,1);			//1 to load data for normal slot at time t 
				
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
					//sumSq[i]=0.0;
					COMx[i]=0.0;COMy[i]=0.0;COMz[i]=0.0;
				}
				for (i=0 ; i<Ns ; i++){
					dx = sphere[i].x- sphere[i].x0;
					dy = sphere[i].y- sphere[i].y0;
					dz = sphere[i].z- sphere[i].z0;
					
					COMx[Cell[i]]+=sphere[i].x;
					COMy[Cell[i]]+=sphere[i].y;
					COMz[Cell[i]]+=sphere[i].z;
                    if ( dz >= HALF_LENGTH )
						COMz[Cell[i]]-=LENGTH;
                    else if ( dz < -HALF_LENGTH )
						COMz[Cell[i]]+=LENGTH;
                    if ( dy >= HALF_LENGTH )
						COMy[Cell[i]]-=LENGTH;
                    else if ( dy < -HALF_LENGTH )
						COMy[Cell[i]]+=LENGTH;
                    if ( dx >= HALF_LENGTH )
						COMx[Cell[i]]-=LENGTH;
                    else if ( dx < -HALF_LENGTH )
						COMx[Cell[i]]+=LENGTH;


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
					
					r2=L*L*(dx*dx+dy*dy+dz*dz);	
					sum[Cell[i]] += r2;

					//if(r2>.0)printf("%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",r2,i,sphere[i].x,sphere[i].x0,sphere[i].y,sphere[i].y0,sphere[i].z,sphere[i].z0);
				}
				for(i=0;i<Ncell;i++){
					msd[i][time_index] += (sum[i])/(double)NumInCell[i];
					dx=(COMx[i]-COMx0[i])/NumInCell[i];
					dy=(COMy[i]-COMy0[i])/NumInCell[i];
					dz=(COMz[i]-COMz0[i])/NumInCell[i];
					msdCell[time_index] += L*L*(dx*dx+dy*dy+dz*dz);
					msdSqCell[time_index] += L*L*L*L*(dx*dx*dx*dx+dy*dy*dy*dy+dz*dz*dz*dz);
					counters[i][time_index]++;
					CounterBinder[time_index]++;
					delta_t[i][time_index] = time_step*(double)(other_frame - current_frame);
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

//	Averaging over blocks
	double BlkAvgMSD[num_of_intervals];
	int BlkAvgCounter[num_of_intervals];

	for(int j=1;j<num_of_intervals;j++){
		BlkAvgMSD[j]=0.0;
		//BlkAvgMSDCell[j]=0.0;
	//	BlkAvgMSDSqCell[j]=0.0;
		BlkAvgCounter[j]=0.0;
		for(i=0;i<Ncell;i++){
			if(counters[i][j]){
				BlkAvgMSD[j]+=msd[i][j]/(double)counters[i][j];
	//			BlkAvgMSDCell[j]+=msdCell[i][j]/(double)counters[i][j];
	//			BlkAvgMSDSqCell[j]+=msdSqCell[i][j]/(double)counters[i][j];
				BlkAvgCounter[j]++;
			}
		}
	}
	sprintf(msd_file_name,"%s/BlockAnalysis/MSD/MSDLb%.3lf.dat",sub_dir,L/InvCellSize);
	msd_file = fopen(msd_file_name,"wb");
	fprintf(msd_file,"0.0\t0.0\t0.0\t0.0\n");
	for (k=1;k<num_of_intervals;k++){
		if ( BlkAvgCounter[k] )
			fprintf(msd_file,"%.14lf\t%.14lf\t%.14lf\t%.14lf\n",delta_t[1][k],BlkAvgMSD[k]/(double)BlkAvgCounter[k],msdCell[k]/(double)CounterBinder[k],msdSqCell[k]/(double)CounterBinder[k]);
			//fprintf(msd_file,"%.14lf\t%.14lf\t%.14lf\n",delta_t[0][k],msd[0][k]/(double)counters[0][k],msdSq[0][k]/(double)counters[0][k]);
	}
	fclose(msd_file);

/*	for (k=1;k<num_of_intervals;k++){
		if ( BlkAvgCounter[k] ){
			avg_msd[1][k] +=BlkAvgMSD[k]/(double)BlkAvgCounter[k]; 
			avg_msd[2][k] +=BlkAvgMSDSq[k]/(double)BlkAvgCounter[k];
			avg_msd[0][k] =delta_t[0][k];
		}
	}	
*/
	for(i=0;i<Ncell;i++)
	{	free(msd[i]);
	//	free(msdCell[i]);free(msdSqCell[i]);
		free(delta_t[i]);free(counters[i]);}
	free(msd);free(msdCell);free(msdSqCell);free(delta_t);free(counters);
//	free(sphere);
	return;
}
