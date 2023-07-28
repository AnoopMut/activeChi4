#include "global.h"
#include "subroutine.h"

//************************initializing gsl random variable***********************************//
void initializeGsl(){
    rng_ptr = gsl_rng_alloc (gsl_rng_mt19937); // allocate the rng
    unsigned long int seed;
    seed = (unsigned long)time(NULL)+jobno+rand()+node;
    gsl_rng_set (rng_ptr, seed);
    return;
}
//************************Getting normal random number using gsl*****************************//
double normal(double mean, double sigma){
    double ranNum=gsl_ran_gaussian_ziggurat(rng_ptr,sigma);
    ranNum += mean;
return ranNum;
}
//************************Getting uniform random number using gsl*****************************//
double uniform(double a, double b){
	return gsl_ran_flat(rng_ptr,a,b);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//********Generate xal initial config. and rondom momentaaas with zero total drift***********//
double drand(){   /* uniform distribution, (0..1] */
    double  ran=(rand()+1.0)/(RAND_MAX+1.0);
    return ran ;
}

double random_normal(){
 /* normal distribution, centered on 0, std dev 1 */
         return sqrt(-2*log(drand())) * cos(2*pi*drand());
}

void GetFccGrid(){
	int troot_n=1;
	while(troot_n*troot_n*troot_n<N)
		troot_n +=1;													//getting maximum no which is perfect cube and less than N

	double dr = (LENGTH)/(double)(troot_n);
	double X[N]	;
	double Y[N] ;
	double Z[N] ;
	bool typ[N] ;
	for(int i=0 ; i<N ; i++){
		X[i] = dr*(double)((i%(troot_n*troot_n))%troot_n);
		Y[i] = dr*(double)(i%(troot_n*troot_n))/(double)(troot_n);		//X,Y,Z are now points of almost square lattice with dimensions L cross L
		Z[i] = dr*(double)(i/(double)(troot_n*troot_n));
		typ[i]=true;
	}
	int  k=0;
	do{
		int i =floor((float)rand()/(float)(RAND_MAX)*N);
		if(typ[i] == true){
			sphere[k].x = X[i];
			sphere[k].y = Y[i];
			sphere[k].z = Z[i];
			k++;
			typ[i]=false;
		}

	}while(k<N);
	if(k!=N) fprintf(stderr,"Initialisation Failed \n couldn't place all of the spheres on lattice");
	// making 80:20 3D Kob Anderson glass forming mixture.
    for (int i=0; i<N; i++) sphere[i].type = 0;               //Initialisig all of them to type A
    k=0;
    do{
		int ran =floor((float)rand()/(float)(RAND_MAX)*N);
        if(sphere[ran].type == 0){                              //assigning randomy 20% sphere as B
            sphere[ran].type = 1;                               //even one can do that sequentially 
            k++;                                                //since already they are randomly placed
        }
    }while(k < (N-NA));
	//Initialise momentaaas.........
	double px = 0.0;
	double py = 0.0;												//to fix net drift
	double pz = 0.0e0;
	for (int i=0; i<N;i++){
		sphere[i].px = ((double)rand()/(double)RAND_MAX- 0.5 );
		sphere[i].py = ((double)rand()/(double)RAND_MAX- 0.5 );
		sphere[i].pz = ((double)rand()/(double)RAND_MAX- 0.5 );
		px += sphere[i].px;
		py += sphere[i].py;
		pz += sphere[i].pz;
	}
	//fixing drift
	px = px/(double)N;
	py = py/(double)N;
	py = pz/(double)N;
	for (int i=0; i<N ;i++ ) {
		sphere[i].px -= px;
		sphere[i].py -= py;
		sphere[i].pz -= pz;
	}
	char file1[128];
	char file2[128];
	sprintf(file1,"%s/ini.dat",sub_dir);
	print_eq(file1);
 
	//Hope this fix is good...................
	if(node==0) fprintf(stderr, "\n\nMomentum and position initialisation done\n\n");
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void fixDrift(){
	double Px,Py,Pz;
	double Pxblk,Pyblk,Pzblk;
	double drift=0.0,driftblk=0.0;
	int i;
	Pxblk = 0.0; Pyblk = 0.0; Pzblk = 0.0;
	Px = 0.0; Py = 0.0; Pz = 0.0;
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		Pxblk += sphere[part].px ; Pyblk += sphere[part].py; Pzblk += sphere[part].pz;
	}
    MPI_Reduce(&Pxblk, &Px, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Px, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&Pyblk, &Py, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Py, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&Pzblk, &Pz, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Pz, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	Px = Px/(double)N;
	Py = Py/(double)N;
	Pz = Pz/(double)N;
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		sphere[part].px -= Px;
		sphere[part].py -= Py;
		sphere[part].pz -= Pz;
		drift += sphere[part].px;
	}
	MPI_Reduce(&driftblk, &drift, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if(node==0)	printf("Drift After FixDrift=%lf\n",drift);
	return;
}		
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void adjustInitialTemp(double final_temp){
	calculate_temp();
	//scaling factor
	double xi_T = sqrt(final_temp/inst_T);
	for(int i=0; i<N; i++){
		sphere[i].px = sphere[i].px *xi_T;
		sphere[i].py = sphere[i].py *xi_T;
		sphere[i].pz = sphere[i].pz *xi_T;
	}    
	calculate_temp();
	return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void calculate_temp(){
	double keBlk=0.0;
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		keBlk += ( sphere[part].px*sphere[part].px + sphere[part].py*sphere[part].py + sphere[part].pz*sphere[part].pz);        
	}
	keBlk    = keBlk*0.5;
    MPI_Reduce(&keBlk, &ke, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ke, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	inst_T = 2.0*(ke)/(double)DOF;
	return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void breakxal(double rate, double final_temp){
	adjustInitialTemp(final_temp * 1.50);
	T=inst_T;
	update_cell_list();
	calculate_forces();	
	if(node==0) printf("\n High Temp start to break xal at T =  %lf; PE = %lf; it would start cool now :) \n",T,pe);
	while ( T > final_temp ){
		calculate_temp(); 
		adjustInitialTemp(T);
		advance_time_NVEFrog();
		T -= time_step*rate;
	}
	T=final_temp;
	if(node==0)  printf("\n Temprature after cooling %lf \n",T);
	return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void coldstart(double rate, double final_temp){
	adjustInitialTemp(final_temp*0.5);
	T=inst_T;
	CommunicateCoordinates(0);							//Lets get the slabs for first time 
	update_cell_list();
	calculate_forces();	
	if(node==0) printf("\n Lets start it slow  at temprature %lf and PE %lf and heat it so that nothing blows up\n",T,pe);
	while ( T < final_temp ){
		calculate_temp(); 
		adjustInitialTemp(T);
		advance_time_NVEFrog();
		T += time_step*rate;	
		//MPI_Reduce(&peblk, &pe, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		//if(node==0)printf("%lf\t%lf\t%lf\n",pe,ke,pe+ke);
	}
	if(node==0)  printf("\nTemprature after heating %lf \n",T);
	T=final_temp;
	return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void initializeSystem(){
	initializeGsl() ;
	if(node==0) GetFccGrid();
	MPI_Barrier(MPI_COMM_WORLD);
	GetDataTypeSLAB();
	GetDataTypeSLABF();
	AllocateArrays();
	if (node==0) selectRandomActiveParticles();
	BcastInitState();
	for(int i=0;i<N;i++){
		sphere[i].x_u=sphere[i].x;
		sphere[i].y_u=sphere[i].y;
		sphere[i].z_u=sphere[i].z;
	}
	SetUpParallel();
	//printf("F\n\n");
	SetUpCellList();
	CommunicateCoordinates(0);							//Lets get the slabs for first time 
	applypbc1();
	coldstart(.1,T);
	breakxal(.1,T);
	adjustInitialTemp(T);

	return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void swap(int *a, int *b){
	int temp;
	temp = *a;
	*a=*b;
	*b=temp;
return;
}

void InititializeRandomVel(){
	double px = 0.0;
	double py = 0.0;												//to fix net drift
	double pz = 0.0e0;
	for (int i=0; i<N;i++){
		sphere[i].px = ((double)rand()/(double)RAND_MAX- 0.5 );
		sphere[i].py = ((double)rand()/(double)RAND_MAX- 0.5 );
		sphere[i].pz = ((double)rand()/(double)RAND_MAX- 0.5 );
		px += sphere[i].px;
		py += sphere[i].py;
		pz += sphere[i].pz;
	}
	//fixing drift
	px = px/(double)N;
	py = py/(double)N;
	pz = pz/(double)N;
	for (int i=0; i<N ;i++ ) {
		sphere[i].px -= px;
		sphere[i].py -= py;
		sphere[i].pz -= pz;
	}
	double keBlk=0.0;
	for(int part=0; part<N; part++){
		keBlk += ( sphere[part].px*sphere[part].px + sphere[part].py*sphere[part].py + sphere[part].pz*sphere[part].pz);        
	}
	keBlk    = keBlk*0.5;
	inst_T = 2.0*(keBlk)/(double)DOF;
	double xi_T = sqrt(T/inst_T);
	for(int i=0; i<N; i++){
		sphere[i].px = sphere[i].px *xi_T;
		sphere[i].py = sphere[i].py *xi_T;
		sphere[i].pz = sphere[i].pz *xi_T;
	}    
	keBlk=0.0;
	for(int part=0; part<N; part++){
		keBlk += ( sphere[part].px*sphere[part].px + sphere[part].py*sphere[part].py + sphere[part].pz*sphere[part].pz);        
	}
	keBlk    = keBlk*0.5;
	inst_T = 2.0*(keBlk)/(double)DOF;
return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void initializeOverDampedBD(){
	initializeGsl();
	//printf("%lf\n",normal(0.0,1.0));
	fxRandomOld = (double *)malloc(sizeof(double)*N);	
	fyRandomOld = (double *)malloc(sizeof(double)*N);	
	fzRandomOld = (double *)malloc(sizeof(double)*N);	
	for(int i=0;i<N;i++)
		fxRandomOld[i]=normal(0.0,1.0);	
	SqrtVariance=sqrt(2.0*T*T/(time_step*D0));
	InvD0=1.0/D0;
	D0ByT=D0/T;
	TByD0=T/D0;
	expFactor = exp(-time_step*TByD0);
	return;
}
