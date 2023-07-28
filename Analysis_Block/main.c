#include "global.h"
#include "subroutine.h"


int main(int n,char **inputStrings){
	double factor = 1.40000000;
//	double PARAMETER=0.18;
	double PARAMETER=0.3;
	double F0;
	PARAMETER=PARAMETER*PARAMETER;
	sscanf(inputStrings[1],"%d",&serial);
	sscanf(inputStrings[2],"%lf",&T);
	sscanf(inputStrings[3],"%lf",&F0);
	L = pow((double)Ns/RHO,1.0/3.00);
//	T = 0.900;
	invL=1.0/L;
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
	sscanf(buff,"%lf\t%lf\t%lf\t%lf\t%d\t%lf",&dummyd,&dummyd,&dummyd,&dummyd,&dummy,&dummyd);
	fgets(buff,1024,(FILE*)para_file ); //third line reads discription
	fgets(buff,1024,(FILE*)para_file ); 
	sscanf(buff,"%d\t%d",&start_ensemble,&end_ensemble);
	fclose(para_file);
	tot_ensemble = end_ensemble-start_ensemble+1;
	sprintf(main_dir, "../%.3lf",T);
	printf("%s\n",main_dir);
	int is_avg_req =1;

//	printf("are we here?\n");
	for(jobno=start_ensemble; jobno <= end_ensemble; jobno++){
		sprintf(sub_dir, "%s/data_00%d",main_dir,jobno);
		double scale = 3.0;
		double cellCutOff = L/scale;
//	printf("are we here?\n");
		read_data_file();
		char mkdir[1024];
		sprintf(mkdir, "mkdir -p %s/BlockAnalysis/MSD",sub_dir);
		system(mkdir);


		sprintf(mkdir, "mkdir -p %s/BlockAnalysis/Relax",sub_dir);
		system(mkdir);
//		for(int i=3;i<22;i++){
//			sprintf(mkdir, "mkdir -p %s/BlockAnalysis/Relax/%d",sub_dir,i);
//			system(mkdir);
//		}
		while(scale<20.0){
			mean_square_displacement(factor,is_avg_req,cellCutOff);
			overlap_correlation_function(factor, PARAMETER,is_avg_req,cellCutOff,(int)scale);
			VanHove(factor, 0 , cellCutOff, 1500.0);
			scale = scale + 1.0;
			cellCutOff = L/scale;

		}
		if(is_avg_req ==1) is_avg_req=0;
		free(all_data_x);
		free(all_data_y);
		free(all_data_z);
		free(all_data_type);
		free(frame_index_value);
	}
//	print_ensemble_avg(tot_ensemble);
	return 0;
}

