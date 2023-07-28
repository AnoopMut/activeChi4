#include "global.h"
#include "subroutine.h"
void guu(double factor, double Time){
	FILE *guu_file;
	char guu_file_name[128];
	int k,i,current,other,origin_number,offset,maximal_interval;
	int current_frame,other_frame,time_index;
	int num_of_origins,num_of_intervals;
	double *delta_t;
	double dx,dy,dz;	
	double MSDt,MDt,MSDtPP,MDtPP,bin;	
	MSDtPP=0.0;
	MDtPP=0.0;
	MSDt=0.0;
	MDt=0.0;
    bin = 0.005;
	double   MaxBin = (int)(0.5*L/bin);
	int FrameNum[6];
//------------------------------------------------------------------------
	/*Global variables needed ony at node 0*/
	sphere 	=	calloc(N,sizeof(SPHERE)); 
    double *hist;
    double *hist1;
    double *hist2;
	if(node==0){
	    hist = (double *)malloc(sizeof(double)*MaxBin);
	    hist1 = (double *)malloc(sizeof(double)*(double)MaxBin);
	    hist2 = (double *)malloc(sizeof(double)*(double)MaxBin);
	}
	origin_number=200;
	maximal_interval = origin_number*origin/DIVIDE_TO;
	num_of_intervals = 10+(int)( log((double)maximal_interval)/log(factor) ); //just to make sure
	delta_t = (double *)malloc(sizeof(double)*num_of_intervals);
//------------------------------------------------------------------------
	/*Per Processor (PP) variables needed ony at node 0*/
    double *histPP;
    double *hist1PP;
    double *hist2PP;
	double *dxt,*dyt,*dzt,*drt;
	double *dxtpp,*dytpp,*dztpp,*drtpp;
    histPP = (double *)malloc(sizeof(double)*MaxBin);
    hist1PP = (double *)malloc(sizeof(double)*(double)MaxBin);
    hist2PP = (double *)malloc(sizeof(double)*(double)MaxBin);
	dxt		=	calloc(N,sizeof(double)); 
	dyt		=	calloc(N,sizeof(double)); 
	dzt		=	calloc(N,sizeof(double)); 
	drt		=	calloc(N,sizeof(double)); 
	dxtpp		=	calloc(NPP,sizeof(double)); 
	dytpp		=	calloc(NPP,sizeof(double)); 
	dztpp		=	calloc(NPP,sizeof(double)); 
	drtpp		=	calloc(NPP,sizeof(double)); 
//------------------------------------------------------------------------
	/*Initialising to zero*/
	for(int i=0;i<MaxBin;i++){
		if(node==0){
			hist[i]=0.0;
			hist1[i]=0.0;
			hist2[i]=0.0;
		}
		histPP[i]=0.0;
		hist1PP[i]=0.0;
		hist2PP[i]=0.0;
	}	
//------------------------------------------------------------------------
//			Getting the times to look to calculate guu
//------------------------------------------------------------------------
	current = 0; origin_number = 0;
    current_frame = frame_index_value[current];
	other = current;
	offset = 1; //the smallest interval
	time_index = 1;
	printf("%d\n",num_of_frames);
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
	MPI_Barrier(MPI_COMM_WORLD);
	double inter1,inter2,inter3,inter4,inter5;
	inter1=1000.0;
	inter2=1000.0;
	inter3=1000.0;
	inter4=1000.0;
	inter5=1000.0;
	for(int i=1;i<time_index;i++){
		if((fabs(delta_t[i]-Time/3.0)<inter1)){
			FrameNum[1]=(delta_t[i]*200.0);
			inter1=fabs(delta_t[i]-Time/3.0);
			}
		if((fabs(delta_t[i]-Time/2.0)<inter2)){
			FrameNum[2]=(delta_t[i]*200.0);
			inter2=fabs(delta_t[i]-Time/2.0);
			}
		if((fabs(delta_t[i]-Time)<inter3)){
			FrameNum[3]=(delta_t[i]*200.0);	
			inter3=fabs(delta_t[i]-Time);
		}
		if((fabs(delta_t[i]-Time*3.0)<inter4)){
			FrameNum[4]=(delta_t[i]*200.0);
			inter4=fabs(delta_t[i]-Time*3.0);
		}
		if((fabs(delta_t[i]-Time*10.0)<inter5)){
			FrameNum[5]=(delta_t[i]*200.0);
			inter5=fabs(delta_t[i]-Time*10.0);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
//------------------------------------------------------------------------

//------------------------------------------------------------------------
	for(int frameindx=3;frameindx<4;frameindx++){
//------------------------------------------------------------------------
	/*Lets get the time displacement U(dt)   (dxt,dyt,dzt,drt)*/
		printf("%d\t%lf\n",FrameNum[frameindx],FrameNum[frameindx]/200.0);
		MSDtPP=0.0;
		MDtPP=0.0;
		MSDt=0.0;
		MDt=0.0;
		for(int i=0;i<MaxBin;i++){
			if(node==0){
				hist[i]=0.0;
				hist1[i]=0.0;
				hist2[i]=0.0;
			}
			histPP[i]=0.0;
			hist1PP[i]=0.0;
			hist2PP[i]=0.0;
		}	
		current = 0; origin_number = 0;
		do{
			current_frame = frame_index_value[current];
			LoadFramePP(current,0);					//Load the PP frame for origin
			other_frame = current_frame+FrameNum[frameindx];
			other=current;
			//look for that frame
			while (other < num_of_frames && frame_index_value[other] != other_frame)
				other++;
			if(other<num_of_frames){
			LoadFramePP(other,1);					//Load the PP frame at required time 
	//------------------------------------------------------------------------
			/*Individual processor would get the time dispacements for assigned particles*/
			for (i=0 ; i<NPP ; i++){
				if(ActTag[i]==1){
					int ii=NPP*node+i;
					dxtpp[i] = (sphere[i].x- sphere[i].x0);
					dytpp[i] = (sphere[i].y- sphere[i].y0);
					dztpp[i] = (sphere[i].z- sphere[i].z0);
	        	    if ( dztpp[i] >= HALF_LENGTH )
	        	        dztpp[i] -= LENGTH;
	        	    else if ( dztpp[i] < -HALF_LENGTH )
	        	        dztpp[i] += LENGTH;
	        	    if ( dytpp[i] >= HALF_LENGTH )
	        	        dytpp[ii] -= LENGTH;
	        	    else if ( dytpp[i] < -HALF_LENGTH )
	        	        dytpp[i] += LENGTH;
	        	    if ( dxtpp[i] >= HALF_LENGTH )
	        	        dxtpp[i] -= LENGTH;
	        	    else if ( dxtpp[i] < -HALF_LENGTH )
	        	        dxtpp[i] += LENGTH;
					double r2=(dxtpp[i]*dxtpp[i]+dytpp[i]*dytpp[i]+dztpp[i]*dztpp[i]);
					double r=L*sqrt(r2);
					drtpp[i]=r;	
					MSDtPP+=L*L*r2;									//MSD at time t
					MDtPP+=r;										//Mean displacement at time t
					dxtpp[i]=L*dxtpp[i];	
					dytpp[i]=L*dytpp[i];	
					dztpp[i]=L*dztpp[i];
	/*				dxtpp[i]=L*dxtpp[i]/r;	
					dytpp[i]=L*dytpp[i]/r;	
					dztpp[i]=L*dztpp[i]/r;*/
				}
			}
	//------------------------------------------------------------------------
			/*Lets get time displacement of each particle to each processor (required for next step)*/
			GatherVec(dxt,dxtpp);
			GatherVec(dyt,dytpp);
			GatherVec(dzt,dztpp);
			GatherVec(drt,drtpp);
	//------------------------------------------------------------------------
			/*Now lets get the spatial variation*/
			double r2;
			LoadOrigin(origin_number);							//Loading the origin frames to all procs
	
			for(int ii=0;ii<NPP;ii++){
				int i=Nprocs*ii+node;							//Distributing the particles sequentially to procs (like part 0 to node 0, part 1 to node 1 ....)
				if(i<N-1){										// In this way procs would get the work evenly  
					for(int j=i+1;j<N;j++){
						if(ActTag[j]==1&&ActTag[i]==1){
							dx = (sphere[j].x0- sphere[i].x0);
							dy = (sphere[j].y0- sphere[i].y0);
							dz = (sphere[j].z0- sphere[i].z0);
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
							int index=(int)(sqrt(r2)/bin);
							if(index<MaxBin){
								histPP[index]+=2.0;					// the normal g(r)
								hist1PP[index]+= 2.0*(drt[i]*drt[j]);	// No orientation info
								hist2PP[index]+= 2.0*(dxt[i]*dxt[j]+dyt[i]*dyt[j]+dzt[i]*dzt[j]); // with orientation info
							}
						}
					}
				}
			}
	
			origin_number++;									// the new origin
			while (current<num_of_frames && frame_index_value[current] != origin_number*origin)
				current++;
//			if(node==0)	printf("origin_number is %d\n",origin_number);
		}
		else break;
		} while (current < num_of_frames);
	
	//------------------------------------------------------------------------
		/*Lets combine the work done by all processors*/
		MPI_Reduce(histPP, hist, MaxBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Barrier(MPI_COMM_WORLD);
	
		MPI_Reduce(hist1PP, hist1, MaxBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Barrier(MPI_COMM_WORLD);
	
		MPI_Reduce(hist2PP, hist2, MaxBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Barrier(MPI_COMM_WORLD);
	  
	  	MPI_Reduce(&MSDtPP, &MSDt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	
	    MPI_Reduce(&MDtPP, &MDt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	
		if(node==0){
			MSDt=MSDt/((double)N*(double)origin_number);
			MDt=MDt/((double)N*(double)origin_number);
			printf("HERE\n")	;
			double norm=(4.0*pi*(RHO)*bin*(double)N*(double)origin_number);
			double norm1=(4.0*pi*RHO*bin*(double)N*(double)origin_number*(double)(MDt*MDt));
			double norm2=(4.0*pi*RHO*bin*(double)N*(double)origin_number*(double)(MSDt));
			if(frameindx==1)	
			sprintf(guu_file_name,"%s/guuTauAlphaBy3.dat",sub_dir);
			if(frameindx==2)	
			sprintf(guu_file_name,"%s/guuTauAlphaBy2.dat",sub_dir);
			if(frameindx==3)	
			sprintf(guu_file_name,"%s/guuTauAlphaAct.dat",sub_dir);
			if(frameindx==4)	
			sprintf(guu_file_name,"%s/guu3TauAlpha.dat",sub_dir);
			if(frameindx==5)	
			sprintf(guu_file_name,"%s/guu10TauAlpha.dat",sub_dir);
			guu_file = fopen(guu_file_name,"wb");
			for(int i=1;i<MaxBin;i++){
				double r=(double)(i+1)*bin;
			//	if(hist[i]!=0.0)
				fprintf(guu_file,"%0.10g\t%.10g\t%.10g\t%.10g\t%.10g\n",r,((double)hist[i]/(norm*r*r)),((double)hist1[i]/(norm1*r*r)),hist2[i]/(norm1*r*r),
				hist2[i]/(norm2*r*r));
	//			fprintf(guu_file,"%0.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\n",r,((double)hist[i]/(norm*r*r)),((double)hist2[i]/(norm*r*r)),hist1[i]/(norm1*r*r),
	//			((((double)hist1[i])/((double)hist[i])*norm/norm1)-1.0), ((double)hist2[i])/(norm2*r*r)+((((double)hist1[i])/((double)hist[i])*norm/norm1)-1.0));
			}
			fclose(guu_file);
		}
	    MPI_Barrier(MPI_COMM_WORLD);
	}
//------------------------------------------------------------------------
	/*Lets free some space*/
	free(dxtpp);
	free(dytpp);	
	free(dztpp);
	free(dxt);
	free(dyt);	
	free(dzt);
	free(sphere);
//------------------------------------------------------------------------
	/*Lets free some more space*/
	if(node==0){
		free(hist);
		free(hist1);
		free(hist2);
	}
	free(histPP);
	free(hist1PP);
	free(hist2PP);
	return;
}
/******************************************************************************************************************/
void GatherVec(double *RecvVec,double *SendVec){
	int DispVec[Nprocs];
    int NumVec[Nprocs];
    for(int i=0;i<Nprocs;i++) NumVec[i]=NPP;
	if(node==0){
		DispVec[0]=0;
		for(int j=0;j<NPP;j++)
			RecvVec[j]=SendVec[j];
		for (int i=1;i<Nprocs;i++)
			DispVec[i]=DispVec[i-1]+NPP;
	    MPI_Gatherv( MPI_IN_PLACE, NPP,MPI_DOUBLE,RecvVec,NumVec,DispVec,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
	else MPI_Gatherv( SendVec, NumVec[node],MPI_DOUBLE,RecvVec,NumVec,DispVec,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(RecvVec,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
	return;
}

