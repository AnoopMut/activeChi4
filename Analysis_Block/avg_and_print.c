#include "global.h"
#include "subroutine.h"

void print_ensemble_avg(int tot_ensemble){
	char avg_msd_file_name[128], avg_relax_file_name[128];
	char avg_msd_theta_file_name[128], avg_msd_phi_file_name[128];
	FILE* avg_msd_file;
	FILE* avg_msd_theta_file;
	FILE* avg_msd_phi_file;
	FILE* avg_relax_file;
	int k ;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//							msd 
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	sprintf(avg_msd_file_name,"%s/avg_msd.dat",main_dir);
	avg_msd_file = fopen(avg_msd_file_name,"wb");
	fprintf(avg_msd_file,"0.0\t0.0\n");
	k=1;
	while(avg_msd[1][k] != 0.0){
		fprintf(avg_msd_file,"%.14lf\t%.14lf\n",avg_msd[0][k],avg_msd[1][k]/(double)tot_ensemble);
		k++;
	}
	fclose(avg_msd_file);
	free(avg_msd[0]);free(avg_msd[1]);
	free(avg_msd[2]);free(avg_msd[3]);



	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//							relax
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	sprintf(avg_relax_file_name,"%s/avg_relax_a_0.18.dat",main_dir);
	avg_relax_file = fopen(avg_relax_file_name,"wb");
	fprintf(avg_relax_file,"0.0\t1.0\t0.0\n");
	k=1;
	while(avg_relax[1][k] != 0.0){
		fprintf(avg_relax_file,"%.10g\t%.10g\t%.10g\n",avg_relax[0][k],
							avg_relax[1][k]/(double)tot_ensemble,
							avg_relax[2][k]/(double)tot_ensemble);
		k++;
	}
	fclose(avg_relax_file);
	free(avg_relax[0]);free(avg_relax[1]); free(avg_relax[2]);free(avg_relax[3]);
	free(avg_relax[4]);free(avg_relax[5]); free(avg_relax[6]);
	return;
}


