#include "global.h"
#include "subroutine.h"

void Fsqt(double factor,double PARAMETER,int flag){
	double parameter=PARAMETER*PARAMETER;							//Sice we always mave MSD
	int k,i,current,other,origin_number,offset,maximal_interval;	//time variables
	int max_time_index,current_frame,other_frame,time_index;		//time variables
	int num_of_origins,num_of_intervals;							//time variables
	double qx,qy,qz,dx,dy,dz,r2;									//Displacement variables	
	int N_Q;														//Total num of Q vecs
	Vecx=(int *) malloc(sizeof(int)*12*23+100);
	Vecy=(int *) malloc(sizeof(int)*12*23+100);							// Allocating q vec
	Vecz=(int *) malloc(sizeof(int)*12*23+100);
	N_Q=GetQ();														//Get the Q's
	double *Q;
	double **SinSum,**CosSum,**CosSqSum,**SinSqSum;					//Double for time data
	double *SINQ,*COSQ;
	double *delta_t;
	int *counters;
	double TotPart=(double)N;
	if (flag==1) TotPart=0.8*TotPart;
	num_of_intervals=num_of_frames;
	current = 0; origin_number = 0;
	do{
		while (current<num_of_frames && frame_index_value[current] != origin_number*origin+1)
			current++;
		//printf("originNumber is %d\n",originNumber);
		if (current<num_of_frames)
			origin_number++;
	}while(current<num_of_frames);
	maximal_interval = origin_number*origin/DIVIDE_TO;
	if(node==0){
		SinSum = (double **)malloc(sizeof(double)*num_of_intervals);
		CosSum = (double **)malloc(sizeof(double)*num_of_intervals);
		CosSqSum = (double **)malloc(sizeof(double)*num_of_intervals);
		SinSqSum = (double **)malloc(sizeof(double)*num_of_intervals);
	}
	if(node==0){
		for(int k=0;k<num_of_intervals;k++){
			SinSum[k] = (double *)malloc(sizeof(double)*N_Q+20);
			CosSum[k] = (double *)malloc(sizeof(double)*N_Q+20);
			CosSqSum[k] = (double *)malloc(sizeof(double)*N_Q+20);
			SinSqSum[k] = (double *)malloc(sizeof(double)*N_Q+20);
		}
		SINQ = (double *)malloc(sizeof(double)*N_Q+20);
		COSQ = (double *)malloc(sizeof(double)*N_Q+20);
		counters = (int *)malloc(sizeof(int)*num_of_intervals);
		delta_t = (double *)malloc(sizeof(double)*num_of_intervals);
	}
    double *cosq;
    double *sinq;
	int  *Degen;
	Q = (double *)malloc(sizeof(double)*N_Q+20);
	cosq = (double *)malloc(sizeof(double)*N_Q+20);
	sinq = (double *)malloc(sizeof(double)*N_Q+20);
	Degen = (int *)malloc(sizeof(int)*N_Q+20);
	

	sphere =	calloc(NPP,sizeof(SPHERE)); 
	//count number of origins
	current = 0; origin_number = 0;
	if(node==0){
		for(int k=0;k<num_of_intervals;k++){
			for(int n=0;n<N_Q;n++){
	    		SinSum[k][n]=0.0;
	    		CosSum[k][n]=0.0;
	    		SinSqSum[k][n]=0.0;
	    		CosSqSum[k][n]=0.0;
				counters[k] = 0;
				delta_t[k] = 0.0;
			}
	    }
	}
	for(int n=0;n<N_Q;n++){
		if((Vecy[n]!=Vecz[n]&&Vecy[n]!=Vecx[n]&&Vecz[n]!=Vecx[n])) Degen[n]=6;
		else if(Vecy[n]==Vecz[n]&&Vecz[n]==Vecx[n])	Degen[n]=1;
		else Degen[n]=3;
	}

//	int FrameNum=(int)(Time*200.00);
	current = 0; origin_number = 0;
	max_time_index = 0;

	do{
        current_frame = frame_index_value[current];
		LoadFramePP(current,0);										// Load the origin frame
		other = current;
		offset = 1; 
		time_index = 1;
		while (other < num_of_frames && offset < maximal_interval){
			other_frame = current_frame + offset; //this is the frame to look for
			//look for that frame
			while (other < num_of_frames && frame_index_value[other] != other_frame)
				other++;
			if (other < num_of_frames){
				for(int n=0;n<N_Q;n++){														//Loop over all q's
					cosq[n]=0.0;
					sinq[n]=0.0;

				}
				LoadFramePP(other,1);														// Load the Del_t frame
    			MPI_Barrier(MPI_COMM_WORLD);
				for (i=0 ; i<NPP ; i++){
					if(sphere[i].type!=1 || flag==0){
						if(sphere[i].x0>=1.0) sphere[i].x0-=1.0;
						if(sphere[i].y0>=1.0) sphere[i].y0-=1.0;
						if(sphere[i].z0>=1.0) sphere[i].z0-=1.0;
						if(sphere[i].x0<0.0) sphere[i].x0+=1.0;
						if(sphere[i].y0<.0) sphere[i].y0+=1.0;
						if(sphere[i].z0<.0) sphere[i].z0+=1.0;			//PBC
						if(sphere[i].x>=1.0) sphere[i].x-=1.0;
						if(sphere[i].y>=1.0) sphere[i].y-=1.0;
						if(sphere[i].z>=1.0) sphere[i].z-=1.0;
						if(sphere[i].x<0.0) sphere[i].x+=1.0;
						if(sphere[i].y<.0) sphere[i].y+=1.0;
						if(sphere[i].z<.0) sphere[i].z+=1.0;
					}
				}
				for (i=0 ; i<NPP ; i++){
					if(sphere[i].type!=1 || flag==0){
						dx = (sphere[i].x- sphere[i].x0);
						dy = (sphere[i].y- sphere[i].y0);
						dz = (sphere[i].z- sphere[i].z0);
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

				//	if(r2<parameter){
						for(int n=0;n<N_Q;n++){														//Loop over all q's
							if(Degen[n]==1){
								qx = TWO_PI*invL*(double)(Vecx[n]);
								qy = TWO_PI*invL*(double)(Vecy[n]);
								qz = TWO_PI*invL*(double)(Vecz[n]);
								Q[n] = qx*qx + qy*qy + qz*qz;
								cosq[n]=cosq[n]+((double)(cos((qx*L*dx) + (qy*L*dy) + (qz*L*dz))));
								sinq[n]=sinq[n]+((double)(sin((qx*L*dx) + (qy*L*dy) + (qz*L*dz))));
							}	
							else if(Degen[n]==3){	
								double cosine=0.0;
								double sine=0.0;
								qx = TWO_PI*invL*(double)(Vecx[n]);
								qy = TWO_PI*invL*(double)(Vecy[n]);
								qz = TWO_PI*invL*(double)(Vecz[n]);
								Q[n] = qx*qx + qy*qy + qz*qz;
								cosine+=((double)(cos((qx*L*dx) + (qy*L*dy) + (qz*L*dz))));
								sine+=	((double)(sin((qx*L*dx) + (qy*L*dy) + (qz*L*dz))));
							
								qx = TWO_PI*invL*(double)(Vecy[n]);
								qy = TWO_PI*invL*(double)(Vecz[n]);
								qz = TWO_PI*invL*(double)(Vecx[n]);
								cosine+=((double)(cos((qx*L*dx) + (qy*L*dy) + (qz*L*dz))));
								sine+=	((double)(sin((qx*L*dx) + (qy*L*dy) + (qz*L*dz))));

								qx = TWO_PI*invL*(double)(Vecz[n]);
								qy = TWO_PI*invL*(double)(Vecx[n]);
								qz = TWO_PI*invL*(double)(Vecy[n]);
								
								cosine+=((double)(cos((qx*L*dx) + (qy*L*dy) + (qz*L*dz))));
								sine+=	((double)(sin((qx*L*dx) + (qy*L*dy) + (qz*L*dz))));
								cosq[n]+=(cosine)/3.0;
								sinq[n]+=(sine)/3.0;
							}
							else{
								double cosine=0.0;
								double sine=0.0;
								qx = TWO_PI*invL*(double)(Vecx[n]);
								qy = TWO_PI*invL*(double)(Vecy[n]);
								qz = TWO_PI*invL*(double)(Vecz[n]);
								Q[n] = qx*qx + qy*qy + qz*qz;
								cosine+=((double)(cos((qx*L*dx) + (qy*L*dy) + (qz*L*dz))));
								sine+=	((double)(sin((qx*L*dx) + (qy*L*dy) + (qz*L*dz))));

								qx = TWO_PI*invL*(double)(Vecx[n]);
								qy = TWO_PI*invL*(double)(Vecz[n]);
								qz = TWO_PI*invL*(double)(Vecy[n]);
								Q[n] = qx*qx + qy*qy + qz*qz;
								cosine+=((double)(cos((qx*L*dx) + (qy*L*dy) + (qz*L*dz))));
								sine+=	((double)(sin((qx*L*dx) + (qy*L*dy) + (qz*L*dz))));
								
								qx = TWO_PI*invL*(double)(Vecy[n]);
								qy = TWO_PI*invL*(double)(Vecx[n]);
								qz = TWO_PI*invL*(double)(Vecz[n]);
								Q[n] = qx*qx + qy*qy + qz*qz;
								cosine+=((double)(cos((qx*L*dx) + (qy*L*dy) + (qz*L*dz))));
								sine+=	((double)(sin((qx*L*dx) + (qy*L*dy) + (qz*L*dz))));
								
								qx = TWO_PI*invL*(double)(Vecy[n]);
								qy = TWO_PI*invL*(double)(Vecz[n]);
								qz = TWO_PI*invL*(double)(Vecx[n]);
								Q[n] = qx*qx + qy*qy + qz*qz;
								cosine+=((double)(cos((qx*L*dx) + (qy*L*dy) + (qz*L*dz))));
								sine+=	((double)(sin((qx*L*dx) + (qy*L*dy) + (qz*L*dz))));

								qx = TWO_PI*invL*(double)(Vecz[n]);
								qy = TWO_PI*invL*(double)(Vecx[n]);
								qz = TWO_PI*invL*(double)(Vecy[n]);
								Q[n] = qx*qx + qy*qy + qz*qz;
								cosine+=((double)(cos((qx*L*dx) + (qy*L*dy) + (qz*L*dz))));
								sine+=	((double)(sin((qx*L*dx) + (qy*L*dy) + (qz*L*dz))));
								
								qx = TWO_PI*invL*(double)(Vecz[n]);
								qy = TWO_PI*invL*(double)(Vecy[n]);
								qz = TWO_PI*invL*(double)(Vecx[n]);
								Q[n] = qx*qx + qy*qy + qz*qz;
								cosine+=((double)(cos((qx*L*dx) + (qy*L*dy) + (qz*L*dz))));
								sine+=	((double)(sin((qx*L*dx) + (qy*L*dy) + (qz*L*dz))));

								cosq[n]+=(cosine)/6.0;
								sinq[n]+=(sine)/6.0;
							}	
						}
					}
				}
				if(node==0){
					for(int n=0;n<N_Q;n++){														//Loop over all q's
						COSQ[n]=0.0;
						SINQ[n]=0.0;
					}
				}
				MPI_Reduce(cosq, COSQ,N_Q, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Reduce(sinq, SINQ,N_Q, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);
				if(node==0){
					counters[time_index]++;
					delta_t[time_index] = time_step*(double)(other_frame - current_frame);
			        for(int n=0;n<N_Q;n++){
			        	SinSum[time_index][n]+=SINQ[n];
			        	CosSum[time_index][n]+=COSQ[n];
			        	SinSqSum[time_index][n]+=SINQ[n]*SINQ[n];
			        	CosSqSum[time_index][n]+=COSQ[n]*COSQ[n];
					}
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
	if(node==0) printf("%d\n",origin_number);
	}while (current<num_of_frames);
	if(node==0){
		FILE *S4Qt_file;
		char S4QT_file_name[128];
		sprintf(S4QT_file_name,"%s/Fsqt_a_%.2lf.dat",sub_dir,PARAMETER);
		if(flag==1) sprintf(S4QT_file_name,"%s/FsqtAA_a_%.2lf.dat",sub_dir,PARAMETER);
		S4Qt_file = fopen(S4QT_file_name,"wb");
		for(int n=0;n<N_Q;n++){
	    	for(int k=0;k<num_of_frames;k++){
				if(counters[k]!=0){
					CosSum[k][n]=CosSum[k][n]/(double)(counters[k]);
    			    CosSqSum[k][n]=CosSqSum[k][n]/(double)(counters[k]);
    			    CosSqSum[k][n]=(double)(CosSqSum[k][n]-CosSum[k][n]*CosSum[k][n])/TotPart;
    			    CosSum[k][n]=(double)CosSum[k][n]/TotPart;
			        SinSum[k][n]=SinSum[k][n]/(double)(counters[k]);
   				    SinSqSum[k][n]=SinSqSum[k][n]/(double)(counters[k]);
    			   	SinSqSum[k][n]=(double)(SinSqSum[k][n]-SinSum[k][n]*SinSum[k][n])/TotPart;
    		    	SinSum[k][n]=(double)SinSum[k][n]/TotPart;
    		    	fprintf(S4Qt_file,"%0.10g\t%0.10g\t%.10g\t%.10g\t%.10g\n",delta_t[k],Q[n],CosSum[k][n],SinSum[k][n],(SinSqSum[k][n]+CosSqSum[k][n]));
				}
			}
    		fprintf(S4Qt_file,"\n\n");
		}
	    for(int k=0;k<num_of_frames;k++){
			free(CosSum[k]);
			free(CosSqSum[k]);
			free(SinSum[k]);
			free(SinSqSum[k]);
		}
		free(CosSum);
		free(CosSqSum);
		free(SinSum);
		free(SinSqSum);
		free(COSQ);
		free(SINQ);
		fclose(S4Qt_file);
	}
	free(sphere);
	free(cosq);
	free(sinq);
	free(Degen);
	free(Vecx);
	free(Vecy);
	free(Vecz);
	MPI_Barrier(MPI_COMM_WORLD);

	return;
}
