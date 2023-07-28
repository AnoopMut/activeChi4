#include "global.h"
#include "subroutine.h"
//#################################################################################################
void selectRandomActiveParticles(){
    int i,temp,k;
    int indices[N];
    FILE *activeFile;
    char activeFileName[128];
    sprintf(activeFileName,"%s/ActiveFileIndex.dat",sub_dir);
    activeFile = fopen(activeFileName,"wb");

    for(i=0;i<N;i++)
        indices[i] = i;

    //permute indices:
    for (i=0; i<N; i++){
        k =floor((float)rand()/(float)(RAND_MAX)*N);//choose random index to switch with
        temp = indices[i];
        indices[i] = indices[k];
        indices[k] = temp;
        sphere[i].actTag = -1;
    }

    for(i=0; i<Nact; i++){
        activeParticleIndex[i] = indices[i];
        sphere[activeParticleIndex[i]].actTag = i;
        fprintf(activeFile,"%d\n",activeParticleIndex[i]);
    }
    fclose(activeFile);
    for (i=0;i<Nact;i++){
        if (i<Nact/2){
            kx[i] = 1;
            ky[i] = 1;
            kz[i] = 1;
        }
        else{
            kx[i] = -1;
            ky[i] = -1;
            kz[i] = -1;
        }
    }
    return;
}

void ShuffleActDirec(){
/*	if(node==0){
	//	printf("Shufflling \n");
		gsl_ran_shuffle (rng_ptr, kx, Nact, sizeof (int));
		gsl_ran_shuffle (rng_ptr, ky, Nact, sizeof (int));
		gsl_ran_shuffle (rng_ptr, kz, Nact, sizeof (int));

	}
*/	
	if(node==0){
    int i,j,temp,N1;

    	//permutate indices
		for (i=0;i<Nact;i++){
			j = (int)((double)Nact*uniform(0.0,1.0));
			temp = kx[i];
			kx[i] = kx[j];
			kx[j] = temp;
		}
    	for (i=0;i<Nact;i++){
			j = (int)((double)Nact*uniform(0.0,1.0));
			temp = ky[i];
			ky[i] = ky[j];
			ky[j] = temp;
        }
    	for (i=0;i<Nact;i++){
			j = (int)((double)Nact*uniform(0.0,1.0));
			temp = kz[i];
			kz[i] = kz[j];
			kz[j] = temp;
        }
	}
	MPI_Barrier(MPI_COMM_WORLD);// better to wait for master node
	MPI_Bcast(kx, Nact, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(ky, Nact, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(kz, Nact, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	return;
}

void saveActiveDirection(char *outFileName){
	if(node==0){
    int i;
    FILE *file;
    file = fopen(outFileName,"wb");
    for (i=0; i<Nact; i++)
        fprintf(file,"%d\t%d\t%d\n",kx[i],ky[i],kz[i]);

    fclose(file);
	}
	MPI_Barrier(MPI_COMM_WORLD);
    return;
}

void readActiveDirection(char *inFileName){
	if(node==0){
    	int i;
    	double dummy;
    	FILE *file;
    	file = fopen(inFileName,"rb");
    	for (i=0;i<Nact;i++){
   	    	fscanf(file, "%d",&(kx[i]));
   	    	fscanf(file, "%d",&(ky[i]));
			fscanf(file, "%d",&(kz[i]));
	    }
	    fclose(file);
	}
	MPI_Barrier(MPI_COMM_WORLD);// better to wait for master node
	MPI_Bcast(kx, Nact, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(ky, Nact, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(kz, Nact, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
    return;
}
void readActiveParticles(char *activeFile){
	if(node==0){
    	int i;
    	FILE *file;
    	file = fopen(activeFile,"rb");
    	for(i=0; i<Nact; i++)
    	    fscanf(file, "%d",&(activeParticleIndex[i]));
    	fclose(file);
		for(i=0;i<N;i++)
			sphere[i].actTag=-1;
    	for(i=0; i<Nact; i++){
        	sphere[activeParticleIndex[i]].actTag = i;
    }

	}
	MPI_Barrier(MPI_COMM_WORLD);// better to wait for master node
	MPI_Bcast(&activeParticleIndex, Nact, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
    return;
}
