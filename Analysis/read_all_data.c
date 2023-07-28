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
	sprintf(frame_index_file_name, "%s/time_index_file.dat",sub_dir);
	frame_index_file = fopen(frame_index_file_name,"rb");
	if(frame_index_file == NULL){
		fprintf(stderr,"Frameindex file does'nt exist %s\n",frame_index_file_name);
		exit(0);
	}
	fgets(buff,1024,(FILE*)frame_index_file);
	sscanf(buff,"total timesteps equals %d", &num_of_frames);
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
	int Num;
//	num_of_frames=1000;
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
	//printf("origin gap? %d\n", origin);
for(int i=0;i<N;i++)
        ActTag[i]=0;
    sprintf(frame_index_file_name, "%s/ActiveFileIndex.dat",sub_dir);
    frame_index_file = fopen(frame_index_file_name, "rb");
    ActNum=0;
    while ( fscanf(frame_index_file,"%d",&Num) != EOF ){
        ActTag[Num]=1;
        ActNum++;
    }
    fclose(frame_index_file);
	printf("%d\n",ActNum);
	sprintf(data_file_name,"%s/data_coor.out",sub_dir);
	sprintf(orient_file_name,"%s/data_orient.out",sub_dir);
	sprintf(length_file_name,"%s/length.dat",sub_dir);
	
	all_data_x = (double *)malloc(sizeof(double)*(num_of_frames+600)*N);
	if(all_data_x == NULL){ fprintf(stderr, "MEMORY NOT AVAILABLE@\n"); exit(0);}
	
	all_data_y = (double *)malloc(sizeof(double)*(num_of_frames+600)*N);
	if(all_data_y == NULL){ fprintf(stderr, "MEMORY NOT AVAILABLE#\n"); exit(0);}
	
	all_data_z = (double *)malloc(sizeof(double)*(num_of_frames+600)*N);
	if(all_data_z == NULL){ fprintf(stderr, "MEMORY NOT AVAILABLE#\n"); exit(0);}
	
	all_data_type = (int *)malloc(sizeof(int)*(num_of_frames+600)*N);
	if(all_data_type == NULL){ fprintf(stderr, "MEMORY NOT AVAILABLE$\n"); exit(0);}
//	L = (double *)malloc(sizeof(double)*(num_of_frames+600));
	

	data_file = fopen(data_file_name,"rb");	
	if(data_file == NULL){
		fprintf(stderr, "DATA FILE DOES NOT EXIST %s\n",data_file_name);
		exit(0);
	}
	int i = 0;
	while ( fscanf(data_file,"%lf\t%lf\t%lf\t%d",&all_data_x[i],&all_data_y[i],&all_data_z[i],&all_data_type[i]) != EOF ){
		i++;
		if((i/N)>num_of_frames){
			printf("# of frames larger than expected");
			break;
		}
	}
	fclose(data_file);
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


		for (j=0,i=frame_to_load*N; i < (frame_to_load+1)*N; i++){
			sphere[j].x0 = all_data_x[i];
			sphere[j].y0 = all_data_y[i];
			sphere[j].z0 = all_data_z[i];
			sphere[j].type0 = all_data_type[i];
			j++;

		}
	}
	else{
		for (j=0,i=frame_to_load*N; i < (frame_to_load+1)*N; i++){
			sphere[j].x = all_data_x[i];
			sphere[j].y = all_data_y[i];
			sphere[j].z = all_data_z[i];
			sphere[j].type = all_data_type[i];
			j++;
		}
	}
	return;	
}

void LoadFramePP(int frame_to_load, int is_t_0){
	int i,j;
	if(is_t_0 == 0){
		for (j=0,i=frame_to_load*NPP; i < (frame_to_load+1)*NPP; i++){
			sphere[j].x0 = Allx[i];
			sphere[j].y0 = Ally[i];
			sphere[j].z0 = Allz[i];
			sphere[j].type0 = Alltype[i];
			j++;
		}
	}
	else{
		for (j=0,i=frame_to_load*NPP; i < (frame_to_load+1)*NPP; i++){
			sphere[j].x = Allx[i];
			sphere[j].y = Ally[i];
			sphere[j].z = Allz[i];
			sphere[j].type = Alltype[i];
			j++;
		}
	}
	return;	
}
void LoadOrigin(int frame_to_load){
	int i,j;
	for (j=0,i=frame_to_load*N; i < (frame_to_load+1)*N; i++){
		sphere[j].x0 = Alloriginx[i];
		sphere[j].y0 = Alloriginy[i];
		sphere[j].z0 = Alloriginz[i];
		sphere[j].type0 = Allorigintype[i];
		j++;
	}
	for(int part=0;part<N;part++){
		while (sphere[part].x0 >= LENGTH) sphere[part].x0 -=LENGTH;
		while (sphere[part].x0 < 0.0 )    sphere[part].x0 +=LENGTH;
		while (sphere[part].y0 >= LENGTH) sphere[part].y0 -=LENGTH;
		while (sphere[part].y0 < 0.0 )    sphere[part].y0 +=LENGTH;
		while (sphere[part].z0 >= LENGTH) sphere[part].z0 -=LENGTH;
		while (sphere[part].z0 < 0.0 )    sphere[part].z0 +=LENGTH;            
	}
	return;	
}



