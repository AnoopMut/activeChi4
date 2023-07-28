#include "global.h"
#include "subroutine.h"
int main(int n,char **inputStrings){
	double factor = 1.4000000;
//	double PARAMETER=0.18;
	double PARAMETER=0.3;
	double F0,tauAlpha;
	MPI_Init(&n,&inputStrings);
	MPI_Comm_size(MPI_COMM_WORLD,&Nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&node);
	PARAMETER=PARAMETER*PARAMETER;
	sscanf(inputStrings[1],"%d",&serial);
	sscanf(inputStrings[2],"%lf",&T);
	sscanf(inputStrings[3],"%lf",&tauAlpha);
	L = pow((double)N/RHO,1.0/3.0);
	invL=1.0/L;
//	T = 0.900;
	int jobno;
	int start_ensemble, end_ensemble;
	int tot_ensemble;
	jobno=serial;
	//###### $$ Reading system parameters form file $$ ###########
	FILE *para_file;
	para_file=fopen("para_file","r");
	if(para_file==NULL){
		fprintf(stderr, "\nSorry man PARA file doesn't exit\n");
		exit(0);
	}
	int dummy;
	double dummyd;
	fgets(buff,1024,(FILE*)para_file ); //first line reads discription
	fgets(buff,1024,(FILE*)para_file ); 
	sscanf(buff,"%lf\t%lf\t%lf\t%lf\t%d\t%lf",&dummyd,&dummy, &dummyd,&dummyd,&dummy,&dummy);
	fgets(buff,1024,(FILE*)para_file ); //third line reads discription
	fgets(buff,1024,(FILE*)para_file ); 
	sscanf(buff,"%d\t%d",&start_ensemble,&end_ensemble);
	fclose(para_file);
	tot_ensemble = end_ensemble-start_ensemble+1;
	sprintf(main_dir, "../%.3lf",T);
	int is_avg_req =1;
	nq=15;
	NT=23*12+100;
	for(jobno=serial; jobno <= serial; jobno++){
		sprintf(sub_dir, "%s/data_00%d",main_dir,jobno);
		if(node==0) read_data_file();
		BCastConstants();
		if(node==0)printf("Data Read\n");
		GetTheOriginConfigs();									// Get the origin frames to each processor (for guu)
		SplitTheSystem();										// Split the system between processors according to particle #
		if(node==0) printf("Splitting Done\n");
		
//		mean_square_displacement(factor,0,1);
		mean_square_displacement(factor,0,0);
//		overlap_correlation_function(factor,0.3,0,1);
		overlap_correlation_function(factor,0.3,0,0);
//		overlap_correlation_function(factor,0.18,0);
//		overlap_cage(factor, 0.3,  is_avg_req );
//		overlap_correlation_function(factor,0.18,is_avg_req );
		S4Qt(factor,0.3);
//		S4Qt(factor,0.18);
//		Fsqt(factor,0.3,0);
//		Fsqt(factor,0.3,1);

//		guu(factor,tauAlpha);										
		if(node==0) printf("guu done\n");
//		guu(factor,1096.720);										
	//	if(is_avg_req ==1) is_avg_req=0;
		free(Allx);
		free(Ally);
		free(Allz);
		free(Alltype);
		free(Alloriginx);
		free(Alloriginy);
		free(Alloriginz);
	}
//	if(node==0)	print_ensemble_avg(tot_ensemble);
    MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
/************************************************************************************************/
void BCastConstants(){
    MPI_Bcast(&num_of_frames, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
	if(node!=0) frame_index_value = (int *)malloc((num_of_frames+10)*sizeof(int));
    MPI_Bcast(frame_index_value, num_of_frames, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&origin, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&ActTag, N, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
	return;
}
/************************************************************************************************/
void SplitTheSystem(){
//-----------------------------------------------------------------------------------------------
	/*Lets Split the system evenly to all processors*/
	NPP=N/Nprocs;														//Num of particles per proc.
	if(node==0)	printf("#of particles Per Procs %d\n",NPP);
	double *Buffz;	
	double *Buffx;	
	double *Buffy;	
	int *Bufftype;	
    MPI_Status sStatus,rStatus;
    MPI_Request rRequest,sRequest;
	Allx=(double *)malloc(sizeof(double)*(num_of_frames)*NPP);	
	Ally=(double *)malloc(sizeof(double)*(num_of_frames)*NPP);	
	Allz=(double *)malloc(sizeof(double)*(num_of_frames)*NPP);	
	Alltype=(int *)malloc(sizeof(int)*(num_of_frames)*NPP);	
//-----------------------------------------------------------------------------------------------
	/*Data nened for node=0 */
	if(node==0){
		Buffx=(double *)malloc(sizeof(double)*(num_of_frames)*NPP);	
		Buffy=(double *)malloc(sizeof(double)*(num_of_frames)*NPP);	
		Buffz=(double *)malloc(sizeof(double)*(num_of_frames)*NPP);	
		Bufftype=(int *)malloc(sizeof(int)*(num_of_frames)*NPP);	
		if(node==0){
			for(int i=0;i<num_of_frames;i++){
				for(int j=0;j<NPP;j++){
					Allx[i*NPP+j]=all_data_x[i*N+j];	
					Ally[i*NPP+j]=all_data_y[i*N+j];	
					Allz[i*NPP+j]=all_data_z[i*N+j];	
					Alltype[i*NPP+j]=all_data_type[i*N+j];	
				}
			}
		}
	}
//-----------------------------------------------------------------------------------------------
	/*Spliting and Sending the splitted data to diffferent nodes from root */
	for(int task=1;task<Nprocs;task++){
		if(node==0){
			for(int i=0;i<num_of_frames;i++){
				for(int j=0;j<NPP;j++){
					Buffx[i*NPP+j]=all_data_x[i*N+NPP*task+j];	
					Buffy[i*NPP+j]=all_data_y[i*N+NPP*task+j];	
					Buffz[i*NPP+j]=all_data_z[i*N+NPP*task+j];	
					Bufftype[i*NPP+j]=all_data_type[i*N+NPP*task+j];	
				}
			}
		}
    	MPI_Barrier(MPI_COMM_WORLD);
		int No=NPP*num_of_frames;
	    if(node==0) MPI_Send(Buffx,No,MPI_DOUBLE,task,0,MPI_COMM_WORLD);
        if(node==task) MPI_Recv(Allx,No,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&rStatus);
	  	MPI_Barrier(MPI_COMM_WORLD);
	    if(node==0) MPI_Send(Buffy,No,MPI_DOUBLE,task,0,MPI_COMM_WORLD);
        if(node==task) MPI_Recv(Ally,No,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&rStatus);
    	MPI_Barrier(MPI_COMM_WORLD);
	    if(node==0) MPI_Send(Buffz,No,MPI_DOUBLE,task,0,MPI_COMM_WORLD);
        if(node==task) MPI_Recv(Allz,No,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&rStatus);
    	MPI_Barrier(MPI_COMM_WORLD);
	    if(node==0) MPI_Send(Bufftype,No,MPI_INT,task,0,MPI_COMM_WORLD);
        if(node==task) MPI_Recv(Alltype,No,MPI_INT,0,0,MPI_COMM_WORLD,&rStatus);
    	MPI_Barrier(MPI_COMM_WORLD);
    }
//-----------------------------------------------------------------------------------------------
	/*Freeing the buffers and the all_data*/
	if(node==0){
		free(Buffx);
		free(Buffy);
		free(Buffz);
		free(Bufftype);
		free(all_data_x);
		free(all_data_y);
		free(all_data_z);
		free(all_data_type);
	}
//-----------------------------------------------------------------------------------------------
return;
}
/************************************************************************************************/
void GetTheOriginConfigs(){
//-----------------------------------------------------------------------------------------------
	/*Lets get the origin frames*/
	int current = 0; int origin_number = 1;
	int originFrame[200];
	originFrame[0]=0;
	if(node==0){
		do{
			while (current<num_of_frames && frame_index_value[current] != origin_number*origin)
				current++;
			originFrame[origin_number]=current;
			if (current<num_of_frames)
				origin_number++;
		}while(current<num_of_frames);
		printf("TotOrigin=%d\n",origin_number);
	}
    MPI_Bcast(&origin_number, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&originFrame, 200, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
//-----------------------------------------------------------------------------------------------
	/*Lets save them*/
	Alloriginx=(double *)malloc(sizeof(double)*(origin_number+1)*N);	
	Alloriginy=(double *)malloc(sizeof(double)*(origin_number+1)*N);	
	Alloriginz=(double *)malloc(sizeof(double)*(origin_number+1)*N);	
	Allorigintype=(int *)malloc(sizeof(int)*(origin_number+1)*N);	

	if(node==0){
		for(int i=0;i<origin_number;i++){
			for(int j=0;j<N;j++){
				Alloriginx[i*N+j]=all_data_x[originFrame[i]*N+j];	
				Alloriginy[i*N+j]=all_data_y[originFrame[i]*N+j];	
				Alloriginz[i*N+j]=all_data_z[originFrame[i]*N+j];	
				Allorigintype[i*N+j]=all_data_type[originFrame[i]*N+j];	

			}
		}
	}
    MPI_Barrier(MPI_COMM_WORLD);
//-----------------------------------------------------------------------------------------------
	/*Lets Bcast them*/
	MPI_Bcast(Alloriginx,origin_number*N,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(Alloriginy,origin_number*N,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(Alloriginz,origin_number*N,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(Allorigintype,origin_number*N,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
//-----------------------------------------------------------------------------------------------
	if(node==0)	printf("Origin frames loaded to all processors\n");
return;
}


