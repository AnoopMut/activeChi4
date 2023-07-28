#include "global.h"
#include "subroutine.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~reading hell lot of data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void read_data_file(){
	int index;
	char frame_index_file_name[128];
	char data_file_name[128], orient_file_name[128];
	char length_file_name[128];
	FILE *frame_index_file;
	FILE *data_file, *orient_file;
	FILE *length_file;									//for NPT simulations
	double dummy_double;
	int dummy_int;
	char dummyD[512],dummyC[512];
	sprintf(frame_index_file_name, "%s/time_index_file.dat",sub_dir);
	frame_index_file = fopen(frame_index_file_name,"rb");
	int f=0;
	if(frame_index_file == NULL){
		fprintf(stderr,"Frameindex file does'nt exist %s\n",frame_index_file_name);
		exit(0);
	}
	fgets(buff,1024,(FILE*)frame_index_file);
	while ( fscanf(frame_index_file,"%d",&dummy_int) != EOF ){
		f++;}
	num_of_frames=f;
//	num_of_frames=517;
	printf("%d\n",f);
	fclose(frame_index_file);
	frame_index_value = (int *)malloc((num_of_frames+10)*sizeof(int));
	if(frame_index_value == NULL){ fprintf(stderr, "MEMORY NOT AVAILABLE!\n"); exit(0);}
	frame_index_file = fopen(frame_index_file_name, "rb");
	fgets(buff,1024,(FILE*)frame_index_file);
	for(int i = 0; i<num_of_frames; i++){
		fgets(buff,1024,(FILE*)frame_index_file);
		sscanf(buff,"%d", &frame_index_value[i]);
	}
	fclose(frame_index_file);

//	num_of_frames=7005;
	/*getting the origin difference (which is linear)*/

	index = 0;
	while ( frame_index_value[index]+1 == frame_index_value[index+1] )
		index++;
	int dummy = index;
	int consecutive = 0;
	while ( consecutive < dummy ){
		if ( frame_index_value[index]+1 == frame_index_value[index+1] )
			consecutive++;
		else
			consecutive = 0;
		index++;
	}
	origin = frame_index_value[index - dummy];


	
	all_data_x = (double *)malloc(sizeof(double)*(num_of_frames+10)*Ns);
	if(all_data_x == NULL){ fprintf(stderr, "MEMORY NOT AVAILABLE@\n"); exit(0);}
	
	all_data_y = (double *)malloc(sizeof(double)*(num_of_frames+10)*Ns);
	if(all_data_y == NULL){ fprintf(stderr, "MEMORY NOT AVAILABLE#\n"); exit(0);}
	
	all_data_z = (double *)malloc(sizeof(double)*(num_of_frames+10)*Ns);
	if(all_data_z == NULL){ fprintf(stderr, "MEMORY NOT AVAILABLE#\n"); exit(0);}
	
	all_data_type = (int *)malloc(sizeof(int)*(num_of_frames+10)*Ns);
	if(all_data_type == NULL){ fprintf(stderr, "MEMORY NOT AVAILABLE$\n"); exit(0);}
	
	sprintf(data_file_name,"%s/data_coor.out",sub_dir);
	data_file = fopen(data_file_name,"rb");	
	if(data_file == NULL){
		fprintf(stderr, "DATA FILE DOES NOT EXIST %s\n",data_file_name);
		exit(0);
	}
	int i = 0;
        int k = 0;
        for(int j=0;j<num_of_frames;j++){
                for(i=0;i<Ns;i++){
                        fscanf(data_file,"%lf %lf %lf %d",&(all_data_x[k]),&(all_data_y[k]),&(all_data_z[k]),&(all_data_type[k]));
                        k++;
                }
        }
    fclose(data_file);
	printf("%d\t%d\n",k,num_of_frames);
	printf("FUCfdddK\n");
	fprintf(stdout,"all data read!!!\n");
	//count number of origins
	int current = 0;tot_origin = 0;
	do{
		while (current<num_of_frames && frame_index_value[current] != tot_origin*origin)
			current++;
		if (current<num_of_frames)
			tot_origin++;
	}while(current<num_of_frames);
	return;

}
void load_frame(int frame_to_load, int is_t_0){
	int i,j;
	if(is_t_0 == 0){
		for (j=0,i=frame_to_load*Ns; i < (frame_to_load+1)*Ns; i++){
			sphere[j].x0 = all_data_x[i];
			sphere[j].y0 = all_data_y[i];
			sphere[j].z0 = all_data_z[i];
			sphere[j].type0 = all_data_type[i];
			j++;

		}
	}
	else{
		for (j=0,i=frame_to_load*Ns; i < (frame_to_load+1)*Ns; i++){
			sphere[j].x = all_data_x[i];
			sphere[j].y = all_data_y[i];
			sphere[j].z = all_data_z[i];
			sphere[j].type = all_data_type[i];
			j++;
		}
	}
	return;	
}
