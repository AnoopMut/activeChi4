#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#define	N					50000
#define	DOUBLE_N			50000.0
#define	RHO					1.20

/* stuff that we don't touch usually */
#define	TIME_STEP				0.005
#define	DELTA(a,b)				((a)==(b)?1.0:0.0)
#define	MAX_NEBZ				64
#define	MAX_CELLS				(N/2)
#define	X_COMP					0
#define Y_COMP 					1
#define	Z_COMP					2
#define	CUTOFF_SQRD				(2.50*2.50)
#define	DIVIDE_TO				2
#define	PARAMETER				0.0324
/* constant constants... :P*******************************************************/
#define	LENGTH					1.0
#define	HALF_LENGTH				0.5

#define	TWO_PI					6.28318530717958
#define	UNI					((double)rand()/((double)RAND_MAX + 1.0))

/*globals*/
double L; 			/*	length of <square> box 	*/
double invL;			/* inverse length of box	*/
double V;			/*	volume			*/
double T;			/*	temperature			*/
double t;			/*	time				*/
double u;			/*	potential energy 		*/
double kinetic;		/*	kinetic energy 		*/
double sxz;
double xi_T;

double rx[N];		/*	x component of position	*/
double ry[N];		/*	y component of position	*/
double rz[N];
double fx[N];		/*	x component of position	*/
double fy[N];		/*	y component of position	*/
double fz[N];
int type[N]; 		/*	for binary system	*/

int serial,origin,numOfFrames;
int *frameIndexArray;
double *allDataX;
double *allDataY;
double *allDataZ;
double F0;
int *allDataType;
char main_dir[128],sub_dir[128];
//notice this initializes the strain increment to the default value... just so its something, but this can be changed.

		char buff[1024];


void readDataFiles(){

        int k,dummy,index,consecutive,j,i;
        char dummyC[10];
        double dummyD;
        char frame_index_file_name[128];
        char dataFileName[128];
        FILE *frame_index_file, *dataFile;

        //sprintf(frameIndexFileName,"nextIndex_%d_%.2f.dat",N,T);
		sprintf(frame_index_file_name, "%s/time_index_file.dat",main_dir);
		frame_index_file = fopen(frame_index_file_name,"rb");
		if(frame_index_file == NULL){
			fprintf(stderr,"Frameindex file does'nt exist %s\n",frame_index_file_name);
			exit(0);
		}
		fgets(buff,1024,(FILE*)frame_index_file);
		sscanf(buff,"total timesteps equals %d", &numOfFrames);
		fclose(frame_index_file);
        frameIndexArray = (int *)malloc(sizeof(int)*numOfFrames);
		frame_index_file = fopen(frame_index_file_name, "rb");
		fgets(buff,1024,(FILE*)frame_index_file);
		for( i = 0; i<numOfFrames; i++){
			fgets(buff,1024,(FILE*)frame_index_file);
			sscanf(buff,"%d", &frameIndexArray[i]);
		}
		fclose(frame_index_file);
        index = 0;
        while ( frameIndexArray[index]+1 == frameIndexArray[index+1] )
                index++;
        dummy = index;
        consecutive = 0;
        while ( consecutive < dummy ){
                if ( frameIndexArray[index]+1 == frameIndexArray[index+1] )
                        consecutive++;
                else
                        consecutive = 0;
                index++;
        }
        printf("origin? %d and No of Frames %d\n", frameIndexArray[index - dummy],numOfFrames);
        origin = frameIndexArray[index - dummy];


        //sprintf(dataFileName,"data_%d_%.2f_%.3d.dat",N,T,serial);
        //sprintf(dataFileName,"dump.liq");
		sprintf(dataFileName,"%s/data_coor.out",sub_dir);
        dataFile = fopen(dataFileName,"rb");
        allDataX = (double *)malloc(sizeof(double)*numOfFrames*N);
        allDataY = (double *)malloc(sizeof(double)*numOfFrames*N);
        allDataZ = (double *)malloc(sizeof(double)*numOfFrames*N);
        allDataType = (int *)malloc(sizeof(int)*numOfFrames*N);
		i = 0;
		while ( fscanf(dataFile,"%lf\t%lf\t%lf\t%d",&allDataX[i],&allDataY[i],&allDataZ[i],&allDataType[i]) != EOF ){
		i++;
		}
		fclose(dataFile);
		printf("FUCfdddK\n");
        return;
}



void loadFrame(int frameToLoad, double *rx, double *ry, double *rz, int *type){
	int i,j;
	
	for (j=0,i=frameToLoad*N; i < (frameToLoad+1)*N; i++){
		//printf("%g\n",allDataX[i]);
		rx[j] = allDataX[i];
		ry[j] = allDataY[i];
		rz[j] = allDataZ[i];
		type[j] = allDataType[i];
		j++;
	}
	return;
}


void applyPBC(double *rx, double *ry, double *rz ){
        int i;

        for(i=0;i<N;i++){
                if (rx[i] >= LENGTH){
                        while(rx[i] >= LENGTH)
                                rx[i] = rx[i] - LENGTH;
                }
                else if (rx[i] < 0.0){
                        while(rx[i] < 0.0)
                                rx[i] = rx[i] + LENGTH;
                }

                if (ry[i] >= LENGTH){
                        while(ry[i] >= LENGTH)
                                ry[i] = ry[i] - LENGTH;
                }
                else if (ry[i] < 0.0){
                        while(ry[i] < 0.0)
                                ry[i] = ry[i] + LENGTH;
                }

                if (rz[i] >= LENGTH){
                        while(rz[i] >= LENGTH)
                                rz[i] = rz[i] - LENGTH;
                }
                else if (rz[i] < 0.0){
                        while(rz[i] < 0.0)
                                rz[i] = rz[i] + LENGTH;
                }

        }
        return;
}
void overlapCorrelationFunction(double factor, double parameter, double cellCutOff ){
	//FILE *cFile;
	//char cFileName[164];
	int k,i,current,other,originNumber,offset,maximalInterval;
	int currentFrame,otherFrame,timeIndex,maxTimeIndex;
	int numOfOrigins,numOfIntervals;
	double *rx_0, *ry_0, *rz_0, *rx_t, *ry_t, *rz_t;
	int *type_0, *type_t;
	int **counters;
	double **c, **X4;
	double **delta_t;
	double dx,dy,dz,r2;
	double invCellSize;
	int m, m2,m3, index1;
	int *numInCell, *tag;
	double *sum;

	//count number of origins
	
	rx_0 = (double *)malloc(sizeof(double)*N);
	ry_0 = (double *)malloc(sizeof(double)*N);
	rz_0 = (double *)malloc(sizeof(double)*N);
	rx_t = (double *)malloc(sizeof(double)*N);
	ry_t = (double *)malloc(sizeof(double)*N);
	rz_t = (double *)malloc(sizeof(double)*N);
	type_0 = (int *)malloc(sizeof(int)*N);
	type_t = (int *)malloc(sizeof(int)*N);
	tag = (int *)malloc(sizeof(int)*N);
	
	current = 0; originNumber = 0;
	do{
		while (current<numOfFrames && frameIndexArray[current] != originNumber*origin+1)
			current++;
		//printf("originNumber is %d\n",originNumber);
		if (current<numOfFrames)
			originNumber++;
	}while(current<numOfFrames);
	
	maximalInterval = originNumber*origin/DIVIDE_TO;
	numOfIntervals = 10+(int)( log((double)maximalInterval)/log(factor) ); //just to make sure
//	printf("%d\n",numOfIntervals);
	// lets do the cell division
	m =(int)(L/cellCutOff);
        m3 = m*m*m;             //the length of cells and numInCell...
//	printf("no of cells = %d\n",m3);
        m2 = m*m;
        invCellSize = (double)m;        // for reduced coordinates...
	printf("relax_%.3f\n",L/invCellSize);	
	c = (double **)malloc(sizeof(double *)*m3);
	X4 = (double **)malloc(sizeof(double *)*m3);
	counters = (int **)malloc(sizeof(int *)*m3);
	delta_t = (double **)malloc(sizeof(double *)*m3);
	for(i=0;i<m3;i++){
		c[i] = (double *)malloc(sizeof(double)*numOfIntervals);
		X4[i] = (double *)malloc(sizeof(double)*numOfIntervals);
		counters[i] = (int *)malloc(sizeof(int)*numOfIntervals);
		delta_t[i] = (double *)malloc(sizeof(double)*numOfIntervals);
	}
	sum = (double *)malloc(sizeof(double)*m3);
	numInCell = (int *)malloc(sizeof(int)*m3);

	//printf("allocation is done\n");	

	for (k=0;k<numOfIntervals;k++){
		for(i=0;i<m3;i++){
			c[i][k] = 0.0;
			X4[i][k] = 0.0;
			counters[i][k] = 0;
			delta_t[i][k] = 0.0;
		}
	}
	//printf("initialization is done\n");	
	
	current = 0; originNumber = 0;
	maxTimeIndex = 0;
	do{
		//printf("current = %d\n",current);
		currentFrame = frameIndexArray[current];
		loadFrame(current,rx_0,ry_0,rz_0,type_0);
		applyPBC(rx_0, ry_0, rz_0);
                //set number in each cell to zero
                for (i=0; i<m3; i++)
                        numInCell[i] = 0;

                //put the particles in respective cells : sorting
                for (i=0;i<N;i++){
                        index1 = m2*(int)(rz_0[i]*invCellSize)+m*(int)(ry_0[i]*invCellSize)+(int)(rx_0[i]*invCellSize);
			if(index1 > m3-1)printf("index1 = %d\n",index1);
			tag[i] = index1;
                        numInCell[index1]++;
                }
		//printf("loaded frame number %d which is %d\n",current,frameIndexArray[current]);
		
		other = current;
		offset = 1; //the smallest interval
		timeIndex = 1;
		while (other < numOfFrames && offset < maximalInterval){
			otherFrame = currentFrame + offset; //this is the frame to look for
//			printf("otherFrame is %d\n",otherFrame);
			
			//look for that frame
			while (other < numOfFrames && frameIndexArray[other] != otherFrame)
				other++;
			//printf("other is %d\n",other);
			if (other < numOfFrames){
				//at this point other and current have the indices of the frames.
				//printf("comparing %d with %d  (indices %d with %d)\n",frameIndexArray[current],frameIndexArray[other],current,other);
				loadFrame(other,rx_t,ry_t,rz_t,type_t);
				applyPBC(rx_t, ry_t, rz_t);
				for(i=0;i<m3;i++)
					sum[i] = 0.0;
				for (i=0;i<N;i++){
					dx = rx_t[i] - rx_0[i];
					dy = ry_t[i] - ry_0[i];
					dz = rz_t[i] - rz_0[i];
					// mess due to periodic boundary conditions 
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
					// end of mess 
					r2 = L*L*( dx*dx + dy*dy + dz*dz );
					if (r2 < parameter)
						sum[tag[i]] += 1.0;
				}
				for(i=0;i<m3;i++){
					sum[i] = sum[i]/numInCell[i];
					c[i][timeIndex] += sum[i];
					X4[i][timeIndex] += sum[i]*sum[i];
					counters[i][timeIndex]++;
					//printf("I am here block = %d and counter = %d\n ",i,counters[i][timeIndex]);
					delta_t[i][timeIndex] = TIME_STEP*(double)(otherFrame - currentFrame);
				}
				//advance offset and timeIndex
				if ( (int)(offset*factor) == offset )
					offset++;
				else
					offset = (int)(((double)offset)*factor);
				timeIndex++;
				if(maxTimeIndex < timeIndex)
					maxTimeIndex = timeIndex;
			}
		}
		//go to next origin
		originNumber++;
		while (current<numOfFrames && frameIndexArray[current] != originNumber*origin){
			//printf("current %d and frame = %d and origin*originNo = %d\n",current, frameIndexArray[current], originNumber*origin);
			current++;
		}
		//printf("originNumber is %d and current = %d \n",originNumber,current);
	} while (current<numOfFrames);
	
	//printf("completed the calculations, need to write the data\n");
	for(i=0;i<m3;i++){
		char cFileName[164];
//		printf("relax_%.3f\n",L/invCellSize);
		sprintf(cFileName,"%s/BlockRelax/relax_%.3f_%.4f_%.6d_%.3d.dat",sub_dir,L/invCellSize,T,i,serial);
		FILE *cFile = fopen(cFileName,"w");
		fprintf(cFile,"0.0\t1.0\t0.0\t%d\n",numInCell[i]);
		for (k=1;k<numOfIntervals;k++){
			if ( counters[i][k] )
				fprintf(cFile,"%.10g\t%.10g\t%.10g\t%d\n",delta_t[i][k],c[i][k]/(double)counters[i][k],
				(X4[i][k]/(double)counters[i][k] - c[i][k]*c[i][k]/((double)(counters[i][k]*counters[i][k])))*numInCell[i], numInCell[i]);
		}
		//printf("I am here %d\n",i);
		fclose(cFile);
	}
	double avgc[numOfIntervals];
	double avgX4[numOfIntervals];
	int avgCounter[numOfIntervals];
//Averaging
	for (k=1;k<numOfIntervals;k++){
		avgCounter[k]=0;
		avgc[k]=0.0;
		avgX4[k]=0.0;
		for(i=0;i<m3;i++){
			if ( counters[i][k] ){
				avgc[k]+=c[i][k]/(double)counters[i][k];
				avgX4[k]+=(X4[i][k]/(double)counters[i][k] - c[i][k]*c[i][k]/((double)(counters[i][k]*counters[i][k])))*numInCell[i];
				avgCounter[k]++;
			}
		}
	}
		char cFileName[164];
		sprintf(cFileName,"%s/BlockRelax/Avgrelax_%.3f_%.4f.dat",sub_dir,L/invCellSize,T);
		FILE *cFile = fopen(cFileName,"w");
		fprintf(cFile,"0.0\t1.0\t0.0\n");
		for (k=1;k<numOfIntervals;k++){
			if ( avgCounter[k] )
				fprintf(cFile,"%.10g\t%.10g\t%.10g\n",delta_t[0][k],avgc[k]/(double)avgCounter[k],avgX4[k]/(double)avgCounter[k]);
		}
		//printf("I am here %d\n",i);
		fclose(cFile);
	
	for(i=0;i<m3;i++){
		free(c[i]); free(X4[i]); free(counters[i]);
		free(delta_t[i]);
	}
	free(c); free(X4); free(counters);
	free(delta_t);
	free(rx_0); free(rx_t);
	free(ry_0); free(ry_t);
	free(rz_0); free(rz_t);
	free(type_0); free(type_t);
	free(sum); free(tag); free(numInCell);

	return;
}


void overlapCorrelationFunction1(double factor, double parameter, double cellCutOff ){
	//FILE *cFile;
	//char cFileName[164];
	int k,i,current,other,originNumber,offset,maximalInterval;
	int currentFrame,otherFrame,timeIndex,maxTimeIndex;
	int numOfOrigins,numOfIntervals;
	double *rx_0, *ry_0, *rz_0, *rx_t, *ry_t, *rz_t;
	int *type_0, *type_t;
	int *counters;
	double *c, *X4;
	double *delta_t;
	double dx,dy,dz,r2;
	double invCellSize;
	int m, m2,m3, index1;
	int numInCell;
	int *tag;
	double sum, rCellCutOff;

	//count number of origins
	
	rx_0 = (double *)malloc(sizeof(double)*N);
	ry_0 = (double *)malloc(sizeof(double)*N);
	rz_0 = (double *)malloc(sizeof(double)*N);
	rx_t = (double *)malloc(sizeof(double)*N);
	ry_t = (double *)malloc(sizeof(double)*N);
	rz_t = (double *)malloc(sizeof(double)*N);
	type_0 = (int *)malloc(sizeof(int)*N);
	type_t = (int *)malloc(sizeof(int)*N);
	tag = (int *)malloc(sizeof(int)*N);
	
	current = 0; originNumber = 0;
	do{
		while (current<numOfFrames && frameIndexArray[current] != originNumber*origin+1)
			current++;
		//printf("originNumber is %d\n",originNumber);
		if (current<numOfFrames)
			originNumber++;
	}while(current<numOfFrames);
	
	maximalInterval = originNumber*origin/DIVIDE_TO;
	numOfIntervals = 10+(int)( log((double)maximalInterval)/log(factor) ); //just to make sure

	// lets do the cell division
	m = 1;
        m3 = m*m*m;             //the length of cells and numInCell...
	printf("no of cells = %d\n",m3);
        m2 = m*m;
	
	c = (double *)malloc(sizeof(double)*numOfIntervals);
	X4 = (double *)malloc(sizeof(double)*numOfIntervals);
	counters = (int *)malloc(sizeof(int)*numOfIntervals);
	delta_t = (double *)malloc(sizeof(double)*numOfIntervals);

	//printf("allocation is done\n");	

	for (k=0;k<numOfIntervals;k++){
		c[k] = 0.0;
		X4[k] = 0.0;
		counters[k] = 0;
		delta_t[k] = 0.0;
	}
	//printf("initialization is done\n");	
	
	current = 0; originNumber = 0;
	maxTimeIndex = 0;
	do{
		//printf("current = %d\n",current);
		currentFrame = frameIndexArray[current];
		loadFrame(current,rx_0,ry_0,rz_0,type_0);
		applyPBC(rx_0, ry_0, rz_0);
                //put the particles in respective cells : sorting
		rCellCutOff = cellCutOff/L;
                numInCell= 0;
                for (i=0;i<N;i++){
			tag[i] = 0;
			if(rx_0[i] > 1.0/L && rx_0[i] < rCellCutOff + 1.0/L){
				if(ry_0[i] > 1.0/L && ry_0[i] < rCellCutOff + 1.0/L){
					if(rz_0[i] > 1.0/L && rz_0[i] < rCellCutOff + 1.0/L){
						tag[i] = 1;
                        			numInCell++;
					}
				}
			}
                }
		printf("numInCell = %d\n",numInCell);
		//printf("loaded frame number %d which is %d\n",current,frameIndexArray[current]);
		
		other = current;
		offset = 1; //the smallest interval
		timeIndex = 1;
		while (other < numOfFrames && offset < maximalInterval){
			otherFrame = currentFrame + offset; //this is the frame to look for
			//printf("otherFrame is %d\n",otherFrame);
			
			//look for that frame
			while (other < numOfFrames && frameIndexArray[other] != otherFrame)
				other++;
			//printf("other is %d\n",other);
			if (other < numOfFrames){
				//at this point other and current have the indices of the frames.
				//printf("comparing %d with %d  (indices %d with %d)\n",frameIndexArray[current],frameIndexArray[other],current,other);
				loadFrame(other,rx_t,ry_t,rz_t,type_t);
				applyPBC(rx_t, ry_t, rz_t);
				sum = 0.0;
				for (i=0;i<N;i++){
					dx = rx_t[i] - rx_0[i];
					dy = ry_t[i] - ry_0[i];
					dz = rz_t[i] - rz_0[i];
					// mess due to periodic boundary conditions 
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
					// end of mess 
					r2 = L*L*( dx*dx + dy*dy + dz*dz );
					if (r2 < parameter && tag[i] == 1)
						sum += 1.0;
				}
				sum = sum/numInCell;
				c[timeIndex] += sum;
				X4[timeIndex] += sum*sum;
				counters[timeIndex]++;
				//printf("I am here block = %d and counter = %d\n ",i,counters[i][timeIndex]);
				delta_t[timeIndex] = TIME_STEP*(double)(otherFrame - currentFrame);
				//advance offset and timeIndex
				if ( (int)(offset*factor) == offset )
					offset++;
				else
					offset = (int)(((double)offset)*factor);
				timeIndex++;
				if(maxTimeIndex < timeIndex)
					maxTimeIndex = timeIndex;
			}
		}
		//go to next origin
		originNumber++;
		while (current<numOfFrames && frameIndexArray[current] != originNumber*origin){
			//printf("current %d and frame = %d and origin*originNo = %d\n",current, frameIndexArray[current], originNumber*origin);
			current++;
		}
		//printf("originNumber is %d and current = %d \n",originNumber,current);
	} while (current<numOfFrames);
	
	printf("completed the calculations, need to write the data\n");
	char cFileName[164];
	sprintf(cFileName,"%s/relax_%.3f_%.4f_%.3d_%.3d.dat",sub_dir,cellCutOff,T,0,serial);
	FILE *cFile = fopen(cFileName,"w");
	fprintf(cFile,"0.0\t1.0\t0.0\t%d\n",numInCell);
	for (k=1;k<numOfIntervals;k++){
		if ( counters[k] )
			fprintf(cFile,"%.10g\t%.10g\t%.10g\t%d\n",delta_t[k],c[k]/(double)counters[k],
			(X4[k]/(double)counters[k] - c[k]*c[k]/((double)(counters[k]*counters[k])))*numInCell, numInCell);
	}
	fclose(cFile);
	
	free(c); free(X4); free(counters);
	free(delta_t);
	free(rx_0); free(rx_t);
	free(ry_0); free(ry_t);
	free(rz_0); free(rz_t);
	free(type_0); free(type_t);
	free(tag);

	return;
}

int main(int n,char **inputStrings){
	double factor = 1.20000000;
	double cellCutOff,scale;


	FILE *para_file;
	para_file=fopen("para_file","r");
	if(para_file==NULL){
		fprintf(stderr, "\nSorry man PARA file doesn't exit\n");
		exit(0);
	}
	

	double dummyd;
	int dummy;
	fgets(buff,1024,(FILE*)para_file ); //first line reads discription
	fgets(buff,1024,(FILE*)para_file ); 
	sscanf(buff,"%lf\t%lf\t%lf\t%lf\t%d\t%lf",&dummyd,&T, &dummyd,&dummyd,&dummy,&F0);
	fgets(buff,1024,(FILE*)para_file ); //third line reads discription
	fgets(buff,1024,(FILE*)para_file ); 
	sscanf(buff,"%d\t%d",&dummy,&dummy);
	fclose(para_file);
	sprintf(main_dir, "../data_N_%d_T_%.2lf_F0_%.1lf",N,T,F0);

	sscanf(inputStrings[1],"%d",&serial);
	sprintf(sub_dir, "%s/data_00%d",main_dir,serial);
	char mkdir[1024];
	sprintf(mkdir, "mkdir %s/BlockRelax",sub_dir);
	system(mkdir);
	L = pow(DOUBLE_N/RHO,1.0/3.0);
//	T = 0.470; 
        invL = 1.0/L;
        printf("box size is %.3f and inverse is %.10f\n",L,invL);
	readDataFiles();
	printf("Read done\n");
	//meanSquareDisplacement(factor);


/// Eta mane nicherta off korte pari but check before off
/*	scale = 0.60;
	cellCutOff = L*scale;
	while(cellCutOff + 2.0 < L){
		overlapCorrelationFunction1(factor,PARAMETER,cellCutOff);
		scale = scale + 0.050;
		cellCutOff = L*scale;
	}*/


///Eta thakbe
	scale = 2.0;
	cellCutOff = L/scale;
	while(cellCutOff > 2.0){
		overlapCorrelationFunction(factor,PARAMETER,cellCutOff);
		scale = scale + 1.0;
		cellCutOff = L/scale;
//                printf("relax_%.3f",L/invCellSize);
	}
	printf("are we here?\n");
	
	
	free(allDataX);
	free(allDataY);
	free(allDataType);
	free(frameIndexArray);
	return 0;
}


