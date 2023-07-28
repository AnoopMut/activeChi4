#include "global.h"
#include "subroutine.h"
double  sigma[3][3] = {{1.00,0.80,1.0},{0.80,0.88,1.0},{1.0,1.0,1.0}};
void overlap_cage(double factor, double PARAMETER, int is_avg_req ){
	double parameter=PARAMETER*PARAMETER;
	char c_file_name[128];
	int k,i,current,other,origin_number,offset,maximal_interval;
	int current_frame,other_frame,time_index,max_time_index;
	int num_of_origins,num_of_intervals;
	int *counters;
	double *c, *delta_t, *X4;
	double *DX,*dxpp,*DY,*dypp,*DZ,*dzpp;
	double sum,dx,dy,dz,r2;
	sphere =	calloc(N,sizeof(SPHERE)); 
    DX     =   calloc(N,sizeof(double));
    dxpp       =   calloc(NPP,sizeof(double));
    DY     =   calloc(N,sizeof(double));
    dypp       =   calloc(NPP,sizeof(double));
    DZ     =   calloc(N,sizeof(double));
    dzpp       =   calloc(NPP,sizeof(double));
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


	maps();	
	maximal_interval = origin_number*origin/DIVIDE_TO;
	num_of_intervals = 10+(int)( log((double)maximal_interval)/log(factor) ); //just to make sure
	if(node==0){
		c = (double *)malloc(sizeof(double)*num_of_intervals);
		X4 = (double *)malloc(sizeof(double)*num_of_intervals);
		counters = (int *)malloc(sizeof(int)*num_of_intervals);
		delta_t = (double *)malloc(sizeof(double)*num_of_intervals);
		if(is_avg_req==1){ 
			for(i=0;i<3;i++){
				avg_relax[i]= (double *)malloc(sizeof(double)*num_of_intervals);
				for(k=0 ; k<num_of_intervals;k++) {
					avg_relax[i][k] = 0.0;
				}
			}
		}
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
        double r2;
        LoadOrigin(origin_number);                          //Loading the origin frames to all procs

		GetList();
    	MPI_Barrier(MPI_COMM_WORLD);
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
					dxpp[i] = L*dx;
					dypp[i] = L*dy;
					dzpp[i] = L*dz;

					
//					if (r2 < parameter)
//						sum += 1.0;
				}
        		GatherVec(DX,dxpp);
        		GatherVec(DY,dypp);
        		GatherVec(DZ,dzpp);
				for(i=0;i<NPP;i++){
					dx=0.0;
					dy=0.0;
					dz=0.0;
					for(int j=0;j<NeibNum[i];j++){
						int part=Neib[i][j];
						dx+=DX[part];	
						dy+=DY[part];	
						dz+=DZ[part];	
					}
					dx = dx/(double)NeibNum[i];
					dy = dy/(double)NeibNum[i];
					dz = dz/(double)NeibNum[i];
					double x=sphere[i].x*L;
					double y=sphere[i].y*L;
					double z=sphere[i].z*L;
					r2=(DX[i]*DX[i] + DY[i]*DY[i] + DZ[i]*DZ[i] );
					//r2=(x*x + y*y + z*z );
					r2+=(dx*dx+dy*dy+dz*dz);
					r2-=2.00*(DX[i]*dx+DY[i]*dy+DZ[i]*dz);
				//	r2-=2.00*(x*dx+y*dy+z*dz);
					if(r2<parameter) sum+=1.0;
				}
 				SUM=0.0;
		  		MPI_Reduce(&sum, &SUM,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   				MPI_Barrier(MPI_COMM_WORLD);
				SUM = SUM/(double)N;

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
		sprintf(c_file_name,"%s/relaxMPICR1_a_%.2lf.dat",sub_dir,PARAMETER);
		c_file = fopen(c_file_name,"wb");
		fprintf(c_file,"0.0\t1.0\t0.0\n");
		for (k=1;k<num_of_intervals;k++){
			if ( counters[k] )
				fprintf(c_file,"%.14lf\t%.14lf\t%.14lf\n",delta_t[k],
									c[k]/(double)counters[k],
									X4[k]/(double)counters[k] - c[k]*c[k]/((double)(counters[k]*counters[k])));
		}
		fclose(c_file);
		
		avg_relax[0][0] = 0.0; avg_relax[1][0] = 0.0; avg_relax[2][0] = 0.0;
		for (k=1;k<num_of_intervals;k++){
			if ( counters[k] ){
				avg_relax[1][k] += c[k]/(double)counters[k];
				avg_relax[2][k] += X4[k]/(double)counters[k] - c[k]*c[k]/((double)(counters[k]*counters[k]));
				avg_relax[0][k] = delta_t[k];
			}
		}
		free(c); free(X4);
		free(delta_t); free(counters);
	printf("Overlap Calculated\n");
	}

	free(sphere);
	return;
}
void GetList(){
	n_cellx = L/Cell_cut ;		
	n_celly = n_cellx ; n_cellz =n_cellx;
	N_cell	= n_cellx*n_celly*n_cellz;
	double cell_size = LENGTH/(double)n_cellx;
	int cell_index , cell_index2 , k  ;
//	int In_cell[N_cell+1][200];					   //Ncell+1 becauase cell numbering starts from one 
	int tot_part_in_cell[N_cell+1];
	for(int i=0; i<=N_cell; i++){
		for (int j=0; j<max_part_in_cell; j++) In_cell[i][j]=-1;
			tot_part_in_cell[i]=0;
	}
	for(int i=0; i<N ; i++){
		for(int j=0;j<400;j++) Neib[i][j]=-1;
		NeibNum[i]=0;
		int part=i;
		while (sphere[part].x0 >= LENGTH) sphere[part].x0 -=LENGTH;
		while (sphere[part].x0 < 0.0 )    sphere[part].x0 +=LENGTH;
		while (sphere[part].y0 >= LENGTH) sphere[part].y0 -=LENGTH;
		while (sphere[part].y0 < 0.0 )    sphere[part].y0 +=LENGTH;
		while (sphere[part].z0 >= LENGTH) sphere[part].z0 -=LENGTH;
		while (sphere[part].z0 < 0.0 )    sphere[part].z0 +=LENGTH;            
	}

	/***********************distributing particles to respective cells**********************/
	for(int i=0; i<N; i++){
		cell_index = 1+ floor(sphere[i].x0/cell_size) + floor(sphere[i].y0/cell_size)*n_cellx +
						floor(sphere[i].z0/cell_size)*n_cellx*n_celly ;
		In_cell[cell_index][tot_part_in_cell[cell_index]]=i;
		tot_part_in_cell[cell_index] ++;
		cell_index_each_part[i]=cell_index;
	}

   //***********************Generating cell-list for all particles**************************//
	for(int i=0; i<N; i++){
		cell_index = cell_index_each_part[i];
		for(int j=0; j<tot_part_in_cell[cell_index]; j++){				  	//in its own cell
		        k = In_cell[cell_index][j];

			if(k>i){
				bool is_neib=len_check(sphere[i],sphere[k]);
				if(is_neib){
					Neib[i][NeibNum[i]] = k;
					NeibNum[i] ++; 
					Neib[k][NeibNum[k]] = i;
					NeibNum[k] ++; 
				}
			}
		}
    	MPI_Barrier(MPI_COMM_WORLD);
		for(int cell_neibs = 0; cell_neibs < half_neibs; cell_neibs++){
			cell_index2 = MAP[cell_neibs][cell_index];
			for(int j=0; j<tot_part_in_cell[cell_index2]; j++){				//in neighbouring cells
				k = In_cell[cell_index2][j];
				bool is_neib=len_check(sphere[i],sphere[k]);
				if(is_neib){
					Neib[i][NeibNum[i]] = k;
					NeibNum[i]++;
					Neib[k][NeibNum[k]] = i;
					NeibNum[k]++;
				} 

			}
		}
	}
return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool len_check(SPHERE s1,SPHERE s2){
	bool is_neib=false;
	double dx = s1.x0-s2.x0;
	double dy = s1.y0-s2.y0;
	double dz = s1.z0-s2.z0;
	if( dx < -HALF_LENGTH) dx += LENGTH;
	if( dx > HALF_LENGTH)  dx -= LENGTH;
	if( dy < -HALF_LENGTH) dy += LENGTH;
	if( dy > HALF_LENGTH)  dy -= LENGTH;
	if( dz < -HALF_LENGTH) dz += LENGTH;
	if( dz > HALF_LENGTH)  dz -= LENGTH;
	double dr2 = (dx*dx + dy*dy + dz*dz)*L*L;
//	double sig = sigma[s1.type0][s2.type0];
	//if( dr2 < sig*sig*1.3*1.3 ) is_neib = true;
	if( dr2 < 1.3*1.3 ) is_neib = true;
	return is_neib;	
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 			Get Cell index for location Ix Iy Iz
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int		getcell_I(int Ix ,int Iy , int Iz){
	int cellno = 1+ ((Ix-1+n_cellx*3 )%n_cellx) + ((Iy-1+n_celly*3)%n_celly)*n_cellx
										+((Iz-1+n_cellz*3)%n_cellz)*n_cellx*n_celly;
		return cellno;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 			Routine to lay the map of neighbouring cells
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void maps(){
	n_cellx = L/Cell_cut ;
	n_celly = n_cellx ; n_cellz =n_cellx;
	N_cell	= n_cellx*n_celly*n_cellz;
	for(int i=0;i < half_neibs;i++) {								//rem. to free MAP at end of mains and call it just once at begining if not
		MAP[i]=malloc(sizeof(int)*(N_cell+10));	
		if(MAP[i] ==NULL){
			printf("\nMEMORY NOT AVAILABLE\n");
			exit(0);
		}
	}
	for(int Ix=1 ;Ix<=n_cellx;Ix++){
		for(int Iy=1 ; Iy<=n_celly ; Iy++){
			for(int Iz=1 ; Iz<=n_cellz ; Iz++){	
				int IMAP=getcell_I(Ix,Iy,Iz);
				MAP[0][IMAP] = getcell_I(Ix+1, Iy, Iz);				//east
				MAP[1][IMAP] = getcell_I(Ix+1, Iy+1, Iz);			//north-east			
				MAP[2][IMAP] = getcell_I(Ix, Iy+1, Iz);				//north
				MAP[3][IMAP] = getcell_I(Ix-1 , Iy+1, Iz);			//north-west
				MAP[4][IMAP] = getcell_I(Ix+1, Iy, Iz-1);			//east-bottom
				MAP[5][IMAP] = getcell_I(Ix+1, Iy+1, Iz-1);			//north-east-bottom			
				MAP[6][IMAP] = getcell_I(Ix, Iy+1, Iz-1);			//north-bottom
				MAP[7][IMAP] = getcell_I(Ix-1 , Iy+1, Iz-1);		//north-west-bottom
				MAP[8][IMAP] = getcell_I(Ix+1, Iy, Iz+1);			//east-top
				MAP[9][IMAP] = getcell_I(Ix+1, Iy+1, Iz+1);			//north-east-top			
				MAP[10][IMAP] = getcell_I(Ix, Iy+1, Iz+1);			//north-top
				MAP[11][IMAP] = getcell_I(Ix-1 , Iy+1, Iz+1);		//north-west-top
				MAP[12][IMAP] = getcell_I(Ix, Iy, Iz+1);			//top
			}
		}
	}
	return;
}

