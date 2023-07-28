#include "global.h"
#include "subroutine.h"
void VanHove(double factor, int is_avg_req, double CellCut, double Time){
	FILE *msd_file;
	char msd_file_name[128];
	int k,i,current,other,origin_number,offset,maximal_interval;
	int current_frame,other_frame,time_index;
	int num_of_origins,num_of_intervals;
	double *delta_t;
	int FrameNum;
	double dx,dy,dz,r2;			//saperate sum for rod and sphere
	current = 0; origin_number = 0;
	int countvan=0;

	double *vanHove;
	vanHove= (double *)malloc(sizeof(double)*3*Ns*200);



	do{
		while (current<num_of_frames && frame_index_value[current] != origin_number*origin+1)
			current++;
		if (current<num_of_frames)
			origin_number++;
	}while(current<num_of_frames);
	maximal_interval = origin_number*origin/DIVIDE_TO;
	num_of_intervals = 10+(int)( log((double)maximal_interval)/log(factor) ); 
	delta_t = (double *)malloc(sizeof(double)*num_of_intervals);

	time_index=0;
	while (other < num_of_frames && offset < maximal_interval){
		other_frame = current_frame + offset; //this is the frame to look for
		//look for that frame
		while (other < num_of_frames && frame_index_value[other] != other_frame)
			other++;
		if (other < num_of_frames){
			delta_t[time_index] = time_step*(double)(other_frame - current_frame);
			time_index++;
			if ( (int)(offset*factor) == offset )
				offset++;
			else
				offset = (int)(((double)offset)*factor);
		}
	}
	double inter=1000.0;
	for(int i=1;i<num_of_intervals;i++){
		if((fabs(delta_t[i]-Time)<inter)){
			FrameNum=(delta_t[i]*200.0);	
			inter=fabs(delta_t[i]-Time);
		}
	}

	//Cell Division
	int Ncellx=(int)(L/CellCut);
	int Ncell=Ncellx*Ncellx*Ncellx;
	int Ncellxy=Ncellx*Ncellx;
	printf("# of cells equals %d Length of Block %lf\n",Ncell,L/Ncellx);
	double InvCellSize=(double)Ncellx;
	int NumInCell[Ncell];
	double COMx0[Ncell],COMy0[Ncell],COMz0[Ncell];
	double COMx[Ncell],COMy[Ncell],COMz[Ncell];
	int Cell[Ns];

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

		for(int i=0; i<Ns; i++){
			COMx0[Cell[i]]+=sphere[i].x0;
			COMy0[Cell[i]]+=sphere[i].y0;
			COMz0[Cell[i]]+=sphere[i].z0;
		}
//		while (other < num_of_frames && offset < maximal_interval){

			other_frame = current_frame+FrameNum;
			other=current;
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
				}
				for(i=0;i<Ncell;i++){
					vanHove[countvan]=(COMx[i]-COMx0[i])/NumInCell[i];
					vanHove[countvan+1]=(COMy[i]-COMy0[i])/NumInCell[i];
					vanHove[countvan+2]=(COMz[i]-COMz0[i])/NumInCell[i];
					countvan+=3;
				}
//			}
		}
		//go to next origin
		origin_number++;
		while (current<num_of_frames && frame_index_value[current] != origin_number*origin)
			current++;
		//printf("origin_number is %d\n",origin_number);
	} while (current < num_of_frames);

	sprintf(msd_file_name,"%s/BlockAnalysis/vanHove%.3lf.dat",sub_dir,L/InvCellSize);
	msd_file = fopen(msd_file_name,"wb");
	for (k=0;k<countvan;k++){
			fprintf(msd_file,"%.14lf\n",vanHove[k]*L);
	}
	fclose(msd_file);
	free(delta_t);
	free(vanHove);
	return;
}
