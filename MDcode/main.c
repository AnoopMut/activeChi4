#include "global.h"
#include "subroutine.h"
FILE *fid4;
unsigned long long int *nextStepIndex;
char energy_file_name[128];
char data_file_name[128];
char len_file_name[128];
char stuff_file_name[128];
char eq_file_coo_name[128];
char eq_file_angle_name[128];
char ini_file_coo_name[128];
char ini_file_angle_name[128];
int frame_step;
double eqDuration,prodDuration; 

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void shellSort(int n, unsigned long long int *a){ 
	/*Sorts an array a[1..n] into ascending numerical order by Shell's method (diminishing increment
	  sort). n is input; a is replaced on output by its sorted rearrangement.*/
	int i,j,inc;
	unsigned long long int v;
	inc=1; //Determine the starting increment.
	do {
	inc *= 3;
	inc++;
	} while (inc <= n);
	do { //Loop over the partial sorts.
	inc /= 3;
	for (i=inc+1;i<=n;i++) { //Outer loop of straight insertion.
		v=a[i];
		j=i;
		while (a[j-inc] > v) { //Inner loop of straight insertion.
		a[j]=a[j-inc];
		j -= inc; 
		if (j <= inc) break;
		}
		a[j]=v;
	}
	} while (inc > 1);
	return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//	Array for Logarithmic storage 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


void	readSavingArray(){
        char nextStepIndexFileName[128];
        FILE *outFile;
		int num_of_frames;
		char buff[512];
		sprintf(nextStepIndexFileName,"%s/time_index_file.dat",sub_dir);
		outFile = fopen(nextStepIndexFileName,"rb");
		fgets(buff,1024,(FILE*)outFile);
		sscanf(buff,"total timesteps equals %d", &num_of_frames);
		for(int i = 0; i<num_of_frames; i++){
			fgets(buff,1024,(FILE*)outFile);
			sscanf(buff,"%llu", &nextStepIndex[i]);
		}
return;
}
void prepareSavingArray(unsigned long long int runLength, int numOfOrigins, double factor){
    unsigned long long int current;
	int i;
        char nextStepIndexFileName[128];
        FILE *outFile;
	if(node==0){
        int i,j;
        unsigned long long int k,offset,linearInterval,maximalInterval,index;

        linearInterval = runLength/(unsigned long long int)numOfOrigins;
        //maximalInterval = runLength/DIVIDE_TO;
        maximalInterval = runLength;
        if(node == 0){
                printf("no of origins = %d\n",numOfOrigins);
                printf("linearInterval = %lld\n",linearInterval);
                printf("maximalInterval = %lld\n",maximalInterval);
        }
        current = 0;
        for (k=0;k<numOfOrigins;k++){
                nextStepIndex[current] = (unsigned long long int )k*linearInterval;
                current++;
                offset = 1; //the smallest interval
                while (offset < maximalInterval){
                        index = k*linearInterval + offset;
                        if (index<runLength){
                                nextStepIndex[current] = index;
                                current++;
                        }
                        if ( (unsigned long long int)(offset*factor) == offset )
                                offset++;
                        else
                                offset = (unsigned long long int)((double)offset*factor);
                }
        }
        shellSort(current, nextStepIndex);
        j=0; i=0;
        while (j<current){
                while ( nextStepIndex[j] == nextStepIndex[i] )
                        j++;
                i++;
                nextStepIndex[i] = nextStepIndex[j];
        }
        current = i;
	}
 	MPI_Bcast(&current, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(nextStepIndex, current, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if(node==0){
		sprintf(nextStepIndexFileName,"%s/time_index_file.dat",sub_dir);
		outFile = fopen(nextStepIndexFileName,"wb");
		fprintf(outFile, "total timesteps equals %d\n",current);
		for (i=0;i<current;i++)
			fprintf(outFile,"%llu\n",nextStepIndex[i]);
			fclose(outFile);
	}
    return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void equilibrate(){
	unsigned long int steps = eqDuration/time_step;
	unsigned long int i;
	int start; 														//Starting Index
	char restart_file_name[256];
	sprintf(restart_file_name,"%s/restart_file.dat",sub_dir);
    char ActDirectionFileName[128];
    sprintf(ActDirectionFileName,"%s/ActiveDirection.dat",sub_dir);
    char ActParticleFileName[128];
    sprintf(ActParticleFileName,"%s/ActiveFileIndex.dat",sub_dir);
	InitialiseNoseHoover();
//	InitialiseNoseHooverNPT();


	if(wanaa_restart_eq==0){
		start=0;
		for(int i=0;i<N;i++){
			sphere[i].x_u=sphere[i].x;
			sphere[i].y_u=sphere[i].y;
			sphere[i].z_u=sphere[i].z;
		}
		fixDrift();
	//	update_cell_list();
	//	calculate_forces();
	}
	if(wanaa_restart_eq==1){
		int iter, last_saved_index;
		if(node==0){
			read_restart_state(restart_file_name, &iter, &last_saved_index);
			printf("Reading Done\n");
//			InititializeRandomVel();
//			printf("Randomizing velocity Done int_T =%lf\n",inst_T);
		}
 	   	MPI_Bcast(&iter, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
 	       	readActiveDirection(ActDirectionFileName);
 	       	readActiveParticles(ActParticleFileName);
		if(node==0)	printf("Reading Active Done\n");
		initialise_restart_state();
		start = iter +1;
		fixDrift();
		InitialiseNoseHoover();
		if(node==1) printf("Sytem would now get restart equilibriating from time %.5lf\n", start*time_step );
	}
	MPI_Barrier(MPI_COMM_WORLD);


	time_t t0,t1; 
	t0 = time(0);
	for(i=0;i<1000;i++){
		AdvanceTimeNVTNoseHoover();
//		update_cell_list();
	}
	t1 = time(0);
	if(node==0) printf("Time for 10000 Steps= %ld\n",(t1-t0));
    print_stuff(stuff_file_name);
	for(i=start;i<steps;i++){
		AdvanceTimeNVTNoseHoover();
	  if(i%frame_step == 0){
            print_stuff(stuff_file_name);
        }
		if(i%20000 ==0){ 
			write_restart_state(restart_file_name, i, 0);
			saveActiveDirection(ActDirectionFileName);
		}
		if(i%100000 == 0) fixDrift();
        if(!(i%TAU_per)){															//shuffle directions after persistance time
            ShuffleActDirec();
        }	
	}
	//for(int i=0; i<half_neibs ; i++)
	//	free(MAP[i]); 
	/*********SAVE AFTER EQUILIBRATION....************/
//	print_eq(eq_file_coo_name, eq_file_angle_name);
	ExchangeParticlesInBlocks();
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		sphere[part].x_u=sphere[part].x;
		sphere[part].y_u=sphere[part].y;
		sphere[part].z_u=sphere[part].z;
	}
	CommunicateCoordinates(0);
return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void run(double duration){
	unsigned  long int steps = (int)(duration/time_step);
	unsigned  long int i;
	int start, index;									//starting iteration no of run 
	char restart_file_name[256];
	sprintf(restart_file_name,"%s/restart_file.dat",sub_dir);
    char ActDirectionFileName[128];
    sprintf(ActDirectionFileName,"%s/ActiveDirectionRun.dat",sub_dir);
    char ActParticleFileName[128];
    sprintf(ActParticleFileName,"%s/ActiveFileIndex.dat",sub_dir);
	InitialiseNoseHoover();
	if(wanaa_restart_run){
		int iter, last_saved_index;
		if(node==0){
			read_restart_state(restart_file_name, &iter, &last_saved_index);
	//	chop_off_data_files(last_saved_index);
			printf("Reading Done\n");
		}
        readActiveDirection(ActDirectionFileName);
        readActiveParticles(ActParticleFileName);
			printf("Reading Active Done\n");
		initialise_restart_state();
		start = iter +1;
		index = last_saved_index +1;
		start=0;
		index=0;
		InitialiseNoseHoover();
		if(node==0) printf("Sytem would now get restarted from time %.5lf\n", start*time_step );
	}
	else{
		index=0;
		fixDrift();
		start= 0;
		InitialiseNoseHoover();
	}
	if(node==0) printf("Starting run and generating hell lot of data\n");
	for(i=start;i<steps;i++){
        if(!(i%TAU_per)){															//shuffle directions after persistance time
            ShuffleActDirec();
        }	
		AdvanceTimeNVTNoseHoover();
	//	AdvanceTimeNPTNoseHoover();
	//	AdvanceTimeNVTGaussOSM();
	//	advance_time_NVEFrog();
		unsigned long int temp = nextStepIndex[index];
		if((unsigned long int)i==temp){
			print_data(data_file_name);
		//	print_col(len_file_name,L);	
			index++;
		}
		if(i%frame_step == 0){
			print_stuff(stuff_file_name);
		}
		if(i%100000 == 0) fixDrift();
		if(i%100000 ==0) {
			write_restart_state(restart_file_name, i, index-1);
			saveActiveDirection(ActDirectionFileName);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(node==0) printf("Here it ends\n");
return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void write_restart_state(char *outFileName, int iter, int last_saved_index){
	GatherCoordAndMomentaToRoot();
	MPI_Barrier(MPI_COMM_WORLD);
	if(node==0){
		int i;
		FILE *file;
		file = fopen(outFileName,"wb");
		fprintf(file,"%.14lf\t%d\t%d\n",L,iter,last_saved_index); 
		for (i=0; i<N; i++){
			fprintf(file,"%.14lf\t%.14lf\t%.14lf\t%.14lf\t%.14lf\t%.14lf\t%d\n",sphere[i].x_u,sphere[i].y_u,sphere[i].z_u,
										sphere[i].px,sphere[i].py,sphere[i].pz,sphere[i].type);
		}
		fclose(file);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void read_restart_state(char *restart_file_name, int* iter, int* last_saved_index){
	int i ;
	FILE *file;
	int dummy_iter,dummy_index;
	file = fopen(restart_file_name, "rb");
	if(file==NULL) {printf("File Not found %s",restart_file_name);exit(0);}
	fscanf(file,"%lf\t%d\t%d",&L,&dummy_iter,&dummy_index);
	for (i=0; i<N; i++){
		fscanf(file,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n",&sphere[i].x_u,&sphere[i].y_u,&sphere[i].z_u,
									&sphere[i].px,&sphere[i].py,&sphere[i].pz,&sphere[i].type);
	}
	fclose(file);
	*iter= dummy_iter;
	*last_saved_index= dummy_index;
	return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void initialise_restart_state(){
	if(node==0){
		for(int i=0;i<N;i++){
			sphere[i].x=sphere[i].x_u;
			sphere[i].y=sphere[i].y_u;
			sphere[i].z=sphere[i].z_u;
		}
		applypbc();
	}
	initializeGsl() ;
    GetDataTypeSLAB();
    GetDataTypeSLABF();
    AllocateArrays();
    BcastInitState();
    inv_L=1.00/L;
	for(int i=0;i<N;i++){
        sphere[i].x_u=sphere[i].x;
        sphere[i].y_u=sphere[i].y;
        sphere[i].z_u=sphere[i].z;
    }
    SetUpParallel();
    SetUpCellList();
    CommunicateCoordinates(0);                          //Lets get the slabs for first time 
	calculate_temp();
	update_cell_list();
	calculate_forces();
	return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void chop_off_data_files(int last_saved_index){
	FILE *data_file ;
	FILE *len_file ;
	FILE *temp_file;
	int intDummy , no_of_frames, typ;
	double dummy ,x ,y,z,phix,phiy,phiz;
	int k=0;
	char moveDataFile[256];
	char temp_file_name[256];
	sprintf(temp_file_name, "%s/temp.dat",sub_dir);
	data_file = fopen(data_file_name, "rb");
	while ( fscanf(data_file,"%lf\t%lf\t%lf\t%d",&(dummy),&(dummy),&(dummy),&(intDummy)) != EOF )
					k++;
	no_of_frames=(int)k/(double)(N);
	if(no_of_frames < last_saved_index ){
		printf("data file does'nt have data till restart point%d\t%d\n",no_of_frames,last_saved_index);
		exit(0);
	}
	rewind(data_file);
	temp_file = fopen(temp_file_name,"wb");
	for(k=0;k<=last_saved_index;k++){
		for(int i=0;i<N;i++){
			fscanf(data_file,"%lf\t%lf\t%lf\t%d",&x,&y,&z,&typ);
			fprintf(temp_file,"%.5lf\t%.5lf\t%.5lf\t%d\n",x,y,z,typ);
		}
	}
	fclose(data_file);
	fclose(temp_file);
	sprintf(moveDataFile,"mv %s %s",temp_file_name,data_file_name);
	system(moveDataFile);

	len_file=fopen(len_file_name,"rb");
	k=0;
	while ( fscanf(len_file,"%lf",&(dummy)) != EOF )
		k++;
	no_of_frames =k;
	if(no_of_frames < last_saved_index ){
		fprintf(stderr,"length file does'nt have data till restart point\n");
		exit(0);
	}

	rewind(len_file);
	temp_file = fopen(temp_file_name,"wb");
	for(k=0;k<=last_saved_index;k++){
		fscanf(len_file,"%lf",&x);
		fprintf(temp_file,"%.14lf\n",x);
	}
	fclose(len_file);
	fclose(temp_file);
	
	sprintf(moveDataFile,"mv %s %s",temp_file_name,len_file_name);
		system(moveDataFile);

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(int argc, char *argv[]){
	int serial;
	if(argc==1){
		fprintf(stderr,"\n\nmayday.. requires an argument\n\n");
		exit(0);
	}
	if(argc==2){
		fprintf(stderr, "\nGive argument 0 or 1 to tell weather to start equilibrium fresh or restart\n" );
	}
	if(argc==3){
		fprintf(stderr, "\nGive argument 0 or 1 to tell weather to start run fresh or restart\n" );
	}
	sscanf(argv[1],"%d",&serial);
	sscanf(argv[2],"%d",&wanaa_restart_eq);	
	sscanf(argv[3],"%d",&wanaa_restart_run);	
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&Nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&node);
	//###### $$ Reading system parameters form file $$ ###########
	FILE *para_file;
	para_file=fopen("para_file","r");
	if(para_file==NULL){
		fprintf(stderr, "\nSorry man PARA file doesn't exit\n");
		exit(0);
	}
	char dummy[1028];
	fgets(dummy,1028,(FILE*)para_file ); //first line reads discription
	fgets(dummy,1028,(FILE*)para_file );
	sscanf(dummy,"%lf\t%lf\t%d",&eqDuration,&prodDuration,&frame_step);
	fgets(dummy,1028,(FILE*)para_file );//second line also reads discription
	fgets(dummy,1028,(FILE*)para_file ); 
	sscanf(dummy,"%lf\t%lf\t%lf",&DENSITY,&T,&P);
	fgets(dummy,1028,(FILE*)para_file ); 
	fgets(dummy,1028,(FILE*)para_file ); 
	sscanf(dummy,"%lf",&f0);
	fclose(para_file);

	if(node==0){
		printf("\n\n -----------------------------------------------------------\n\n");
		printf("STARTING MOLECULAR DYNAMICS WITH %d SPHERES \n\n",N);
	}
	jobno = serial;

	sprintf(main_dir, "%.3lf",T);
	sprintf(sub_dir, "%s/data_00%d",main_dir,jobno);

	xold = calloc(N,sizeof(double)); 
	yold = calloc(N,sizeof(double)); 
	zold = calloc(N,sizeof(double)); 
	for(int i=0;i<N;i++){
		xold[i]=0.0;
		yold[i]=0.0;
		zold[i]=0.0;
	}
	//different file names 
	sprintf(data_file_name,		"%s/data_coor.out",sub_dir);
	sprintf(energy_file_name,	"%s/energy.dat",sub_dir);
	sprintf(stuff_file_name,	"%s/stuff.dat",sub_dir);
	sprintf(len_file_name,		"%s/length.dat",sub_dir);
	sprintf(eq_file_coo_name,	"%s/eq_coor.dat",sub_dir);
	sprintf(ini_file_coo_name,	"%s/ini_coor.dat",sub_dir);

	if(wanaa_restart_run==0 & wanaa_restart_eq==0){
		//##### $$ MAKE ALL THE DIRECTORIES $$ ##########
		char mkdir[128];
		sprintf(mkdir, "mkdir -p %s", main_dir);
		if (node==0) system(mkdir);
		char remove[128];
		sprintf(remove, "rm -rf %s",sub_dir);	 //removing the preexisting same DIRECTORIES
		if(node==0) system(remove);
		MPI_Barrier(MPI_COMM_WORLD);
		sprintf(mkdir, "mkdir -p %s", sub_dir);
		if(node==0) system(mkdir);
		//!!!!!! all directories for hell lot of data are just made
		L=pow((double) N/DENSITY,0.33333333);					//starting length of box
		inv_L = 1.0/L;
	    int for_ran1 = -(int)(time(0)+jobno);
  	  	srand( jobno + 100 - for_ran1);           
		//print_eq(ini_file_coo_name,ini_file_angle_name);
		initializeSystem();
		if(node == 0){ 
			printf( "\nEquilibriation started and will go for %d steps rest in peace\n",(int)(eqDuration/time_step));
		}
	}
	if(wanaa_restart_run==0) equilibrate();
	if(node == 0){ 
		printf("\nsystem is now equilibriated\n");
	}
	unsigned long long int runLength = (prodDuration/time_step);
	double factor = 1.4000000000000;							//
	int numOfOrigins = 200;
	long int nextStepIndexSize = numOfOrigins*(10+(int)( log((double)(runLength/DIVIDE_TO))/log(factor) ) );
	nextStepIndex = (unsigned long long int*)malloc(sizeof(unsigned long long int)*nextStepIndexSize);
//	readSavingArray();
	prepareSavingArray(runLength, numOfOrigins, factor);
	run(prodDuration);
	MPI_Barrier(MPI_COMM_WORLD);
//		for(int i=0; i<half_neibs ; i++)
//		free(MAP[i]);
	MPI_Finalize();
	return 0;
}
