#include "global.h"
#include "subroutine.h"
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double print_data(char infile1[]){
    GatherCoordinatesToRoot();
	if(node==0){
		FILE *fid1;
		fid1 = fopen( infile1, "ab" );  //a for appending b for binary
		for(int i=0 ;i<N; i++){
			fprintf(fid1, "%.5lf\t%.5lf\t%.5lf\t%d\n", sphere[i].x_u,sphere[i].y_u,sphere[i].z_u,sphere[i].type);
		}
		fclose(fid1);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	return 0.01;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double print_eq(char infile1[]){
	FILE *fid1 ;
	fid1 = fopen(infile1, "wb");
	for (int i =0; i<N; i++){
		fprintf(fid1, "%.14lf\t%.14lf\t%.14lf\t%d\n",sphere[i].x,sphere[i].y,sphere[i].z,sphere[i].type );
	}
	fclose(fid1);
	return 0.01;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Printing data for movie in gnuplot  
double print_movie(){
	FILE *fid1;
	char mov_file_name[128];
	sprintf(mov_file_name,"%s/formovie.dat",sub_dir);
	fid1 = fopen(mov_file_name, "ab");
	for (int i = 0; i < N; ++i){
		fprintf(fid1, "%.14lf\t%.14lf\t%.14lf\t%d\n",sphere[i].x,sphere[i].y,sphere[i].z,sphere[i].type );
	}
	fprintf(fid1, "\n\n" );	// two empty lines for indexing
	fclose(fid1);
	return 0.01;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double print_2col(char infile[], double dummy1, double dummy2){
	int i;
	FILE *fid;
	fid = fopen(infile,"ab");   
	fprintf(fid,"%g\t%g\n",dummy1,dummy2);
	fclose(fid);    
	return 0.01;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double print_pe(){
	int i;
	FILE *fid;
	fid = fopen("pe.dat","ab");   
	fprintf(fid,"%g\n",pe);
	fclose(fid);    
	return 0.01;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double print_stuff(char infile[]){
    int i;
    FILE *fid;
	double drift=0.0, driftblk=0.0;
	keblk=0.0;
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
  		driftblk+=sphere[part].px;
		keblk += ( sphere[part].px*sphere[part].px + 
					sphere[part].py*sphere[part].py + 
					sphere[part].pz*sphere[part].pz );
	}
	MPI_Reduce(&virialblk, &virial, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&peblk, &pe, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&driftblk, &drift, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&keblk, &ke, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	ke    = ke * 0.5;
	inst_T=2.0*ke/DOF;
	Heff=0.0;	
	pressure = DENSITY*inst_T + virial/(L*L*L);
	if(node==0){	
		Heff=ke+pe+eta[0]*T*DOF;																//For Node Hoover NVT
		for(int i=0;i<M;i++) Heff=Heff+p_eta[i]*p_eta[i]/(2.0*q[i]);
		for(int i=1;i<M;i++) Heff=Heff+T*eta[i];
	
//		Heff=ke+pe+eta[0]*T*DOF+PV+p_eps*p_eps/(2.0*w_eps);  //Heff is the conserved quantity   // For Nose Hoover NPT
//		for(int i=0;i<M;i++) Heff=Heff+p_eta[i]*p_eta[i]/(2.0*q[i])+p_eta_baro[i]*p_eta_baro[i]/(2.0*q_baro[i])+T*eta_baro[i];
//		for(int i=1;i<M;i++) Heff=Heff+T*eta[i];
		fid = fopen(infile,"ab");
	    fprintf(fid,"%g\t%g\t%g\t%g\t%.15lf\t%g\n",pe/(double)N,ke/(double)N,inst_T,pressure,drift,Heff);
	    fclose(fid);
	}
    return 0.01;

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void print_col(char infile[], double dummy){
	if(node==0){
		FILE *fid;
	    fid = fopen(infile,"ab");
	    fprintf(fid,"%g\n",dummy);
	    fclose(fid);
	}
	MPI_Barrier(MPI_COMM_WORLD);
return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
