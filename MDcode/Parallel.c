#include "global.h"
#include "subroutine.h"
// Setting up the number of procs in each direction such that npx>npy>npz and npx*npy*npz=Nprocs
//Also Nprocs should be a perfect power of 2
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
void SetUpParallel(){
	int pow2=0;									//largest pow of 2
	while(pow(2,pow2)<Nprocs)
		pow2++;
	npx=0;
	npy=0;
	npz=0;
	int i=0;
	while(i<pow2){
		npx++;	
		i++;
		if(i<pow2){
			npy++;
			i++;
			if(i<pow2){
				npz++;
				i++;
			}
		}	
	}
	npx=(int)pow(2,npx);
	npy=(int)pow(2,npy);
	npz=(int)pow(2,npz);
	if(node==0) printf("Following is the proocessor grid" 
					"for %d Nprocs\n Nprocsx=%d\tNprocsy=%d\tNprocsz=%d\n",Nprocs,npx,npy,npz);
	//Dividing the space to procs
	Lppx=L/(double)npx;									//Length Allocated per processor
	Lppy=L/(double)npy;
	Lppz=L/(double)npz;
	Me[0]=node%npx;										//My block vector
	Me[1]=(node/npx)%npy;
	Me[2]=(node/npy)/npx;
	MyBoundS[0]=(double)Me[0]*Lppx/L;					//Starting x-boundary of my block
	MyBoundE[0]=((double)Me[0]+1.0)*Lppx/L;				//Ending x-boundary of my block
	MyBoundS[1]=(double)Me[1]*Lppy/L;					//Starting y-boundary of my block
	MyBoundE[1]=((double)Me[1]+1.0)*Lppy/L;				//Ending y-boundary of my block
	MyBoundS[2]=(double)Me[2]*Lppz/L;					//Starting z-boundary of my block
	MyBoundE[2]=((double)Me[2]+1.0)*Lppz/L;				//Ending z-boundary of my block
	for(int i=0;i<Nprocs;i++){
		NumInBlock[i]=0;
	}
	for(i=0;i<N;i++){
		int whichBlock=(int)(sphere[i].x*L/Lppx)+npx*(int)(sphere[i].y*L/Lppy)+
						npx*npy*(int)(sphere[i].z*L/Lppz);
		InBlock[whichBlock][NumInBlock[whichBlock]]=i;
		NumInBlock[whichBlock]++;
		Block[i]=whichBlock;
		InMe[i]=true;									//tag for particle to be in Extd block
		if(whichBlock!=node){
			InMe[i]=false;								//limit it to only particle in My-block 
			sphere[i].x=-2.0;
			sphere[i].y=-2.0;
			sphere[i].z=-2.0;
		}
	}
	int TotPartInAllBlocks=0;	
	for(int i=0;i<Nprocs;i++) TotPartInAllBlocks+=NumInBlock[i];
	if(TotPartInAllBlocks!=N){
		fprintf(stderr, "Tot. Particles in All Blocks don't sum upto N,%d\t%d ",N,TotPartInAllBlocks);
	}
	//printf("In Node %d there are %d # of Part out of %d Tot Part\n",node,NumInBlock[node],TotPartInAllBlocks);
//	printf("nx\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",npx,npy,npz,node,Me[0],Me[1],Me[2]);
return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
void SetUpCellList(){
	Ncellx=(int)Lppx/Cell_cut;									//# of cells along x in a block
	Ncelly=(int)Lppy/Cell_cut;									//# of cells along y in a block
	Ncellz=(int)Lppz/Cell_cut;									//# of cells along z in a block
	CellSizex=Lppx/(double)Ncellx;								//cell size along x in system of unit length
	CellSizey=Lppy/(double)Ncelly;								//cell size along y	in system of unit length
	CellSizez=Lppz/(double)Ncellz;								//cell size along z	in system of unit length
	InvCellSizex=L/CellSizex;	
	InvCellSizey=L/CellSizey;									//In Reduced units	
	InvCellSizez=L/CellSizez;	
	
	Ncellx=Ncellx*npx;											//Tot # of cells along x in whole domain
	Ncelly=Ncelly*npy;											//Tot # of cells along y in whole domain
	Ncellz=Ncellz*npz;											//Tot # of cells along z in whole domain
	
	NCell=Ncellx*Ncelly*Ncellz;									//Tot # of cells in whole domain
	NCellOld=NCell;
	Slabdx=Cell_cut/L;											//Tranfer slab which have the dimensions of a unit cell;	
	Slabdy=Cell_cut/L;	
	Slabdz=Cell_cut/L;	
    for(int i=0;i<Nprocs;i++)
        NumCellInBlk[i]=0;
    for(int z=0;z<Ncellz;z++){
        for(int y=0;y<Ncelly;y++){
            for(int x=0;x<Ncellx;x++){
                int cellindx= x + Ncellx*y +z*Ncellx*Ncelly;
                int blkindx= (int)(x*CellSizex/Lppx)+npx*(int)(y*CellSizey/Lppy)+npx*npy*(int)(z*CellSizez/Lppz);
                CellInBlk[blkindx][NumCellInBlk[blkindx]]=cellindx;
                NumCellInBlk[blkindx]++;
            }
        }
    }
//	printf("%d\t%d\t%d\t%d\t%d\n",NumCellInBlk[node],node,Ncellx,Ncelly,Ncellz);
	maps();
return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
//				Routine to transfer the particles to new processor
//						Needed whenever we create a cell list
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
void ExchangeParticlesInBlocks(){
	MPI_Request sRequest, rRequest;
	MPI_Status sStatus, rStatus;
	MigratedPartNum=0;
	int BlockLocal[N];
	int RecvCounter=0;
	for(int i=0;i<N;i++){
		if(Block[i]!=node)
			InMe[i]=false;										//Lets reset the extd block particles to spatially present part 
	}															//Would now use this tag to identify
	for(int i=0;i<N;i++){
		BlockLocal[i]=0;
		Block[i]=0;												//initailise to zero, to Bcast using reduce
	}	

//***********************************************
//Loop to exchange particles along different Dim 
//Scheme: Transfer sequentially x,y, and z
//		  direction particles including the
//		  recieved ones from preceding dimensions 
//***********************************************
	for(int Dim=0;Dim<3;Dim++){
		int npd;												//# of procs in d (x,y,z) dimension
		if(Dim==0) npd=npx;
		else if(Dim==1) npd=npy;
		else npd=npz;
		int MinusId=Me[Dim]-1;
		if(MinusId==-1) MinusId=npd-1;
		int PlusId=Me[Dim]+1;
		if(PlusId==npd) PlusId=0;
		CountPlus=0;
		CountMinus=0;
		GetSlabsForExchange(Dim);
		int SendTo; 
		int RecvFrom;
		if(Dim==0) {
			SendTo = PlusId + Me[1]*npx + Me[2]*npx*npy;			//send right slab 
			RecvFrom = MinusId + Me[1]*npx + Me[2]*npx*npy;
		}
		else if(Dim==1) {
			SendTo = Me[0] + PlusId*npx + Me[2]*npx*npy;			//send top slab 
			RecvFrom = Me[0] + MinusId*npx + Me[2]*npx*npy;
		}
		else{
			SendTo = Me[0] + Me[1]*npx + PlusId*npx*npy;			//send ToMe slab 
			RecvFrom = Me[0] + Me[1]*npx + MinusId*npx*npy;
		}
		if(SendTo!=node){
			MPI_Isend(PlusSlabF,CountPlus,MPI_DOUBLE,SendTo,0,MPI_COMM_WORLD,&sRequest);
			MPI_Irecv(BuffSlabF,MaxSlabPart,MPI_DOUBLE,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&rRequest);
		    MPI_Wait(&sRequest ,&sStatus );
			MPI_Wait(&rRequest ,&rStatus );
			MPI_Get_count(&rStatus ,MPI_DOUBLE,&RecvCounter);
			if(RecvCounter>0){
				for(int i=0;i<RecvCounter;i=i+7){
					int part=(int)BuffSlabF[i];
					sphere[part].x_u=BuffSlabF[i+1];
					sphere[part].y_u=BuffSlabF[i+2];
					sphere[part].z_u=BuffSlabF[i+3];
					sphere[part].x=BuffSlabF[i+1];								//We would be applying PBC 
					sphere[part].y=BuffSlabF[i+2];
					sphere[part].z=BuffSlabF[i+3];
					sphere[part].px=BuffSlabF[i+4];
					sphere[part].py=BuffSlabF[i+5];
					sphere[part].pz=BuffSlabF[i+6];
					InMe[part]=true;
					MigratedPart[MigratedPartNum]=part;
					MigratedPartNum++;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}	

//***********************************************
// Sending to the opposite direction
//***********************************************
		swap(&RecvFrom,&SendTo);													//Roles are now reversed
		if(SendTo!=node){
			MPI_Isend(MinusSlabF,CountMinus,MPI_DOUBLE,SendTo,0,MPI_COMM_WORLD,&sRequest);
			MPI_Irecv(BuffSlabF,MaxSlabPart,MPI_DOUBLE,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&rRequest);
	    	MPI_Wait(&sRequest ,&sStatus );
			MPI_Wait(&rRequest ,&rStatus );
			MPI_Get_count(&rStatus ,MPI_DOUBLE,&RecvCounter);
			if(RecvCounter>0){
				for(int i=0;i<RecvCounter;i=i+7){
					int part=(int)BuffSlabF[i];
					sphere[part].x_u=BuffSlabF[i+1];
					sphere[part].y_u=BuffSlabF[i+2];
					sphere[part].z_u=BuffSlabF[i+3];
					sphere[part].x=BuffSlabF[i+1];								//We would be applying PBC 
					sphere[part].y=BuffSlabF[i+2];
					sphere[part].z=BuffSlabF[i+3];
					sphere[part].px=BuffSlabF[i+4];
					sphere[part].py=BuffSlabF[i+5];
					sphere[part].pz=BuffSlabF[i+6];
					InMe[part]=true;
					MigratedPart[MigratedPartNum]=part;
					MigratedPartNum++;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
//***********************************************
//Updating the particles in InBlock 
//***********************************************
		int NewNum=0;
		for(int i=0;i<NumInBlock[node];i++){
				int part=InBlock[node][i];
				if(InMe[part]){
						InBlock[node][NewNum]=part;
						NewNum++;
				}
		}
		for(int i=0;i<MigratedPartNum;i++){
				int part=MigratedPart[i];
				if(InMe[part]){
						InBlock[node][NewNum]=part;
						NewNum++;
				}
		}
		NumInBlock[node]=NewNum;
		MigratedPartNum=0;
		applypbc1();
	}
	
// Update Block over all processes and check if we missed some sone to kick out
	for(int j=0;j<NumInBlock[node];j++){
		int i=InBlock[node][j];
		int whichBlock=(int)(sphere[i].x*L/Lppx)+npx*(int)(sphere[i].y*L/Lppy)+
						npx*npy*(int)(sphere[i].z*L/Lppz);
		BlockLocal[i]=whichBlock;
		if(whichBlock!=node) fprintf(stderr,"We missed particle # %d to kick out of node %d\n",i,node);
	}
    MPI_Reduce(BlockLocal, Block, N, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(Block, N, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
void GetSlabsForExchange( int Dim){
	int npd;
	if(Dim==0) npd=npx;
	else if(Dim==1) npd=npy;
	else if(Dim==2) npd=npz;
	int MinusId=Me[Dim]-1;
	if(MinusId==-1) MinusId=npd-1;
	int PlusId=Me[Dim]+1;
	if(PlusId==npd) PlusId=0;
	CountPlus=0;
	CountMinus=0;
	for(int i=0;i<NumInBlock[node];i++){
		int part=InBlock[node][i];
		int NewBlk;
		if(Dim==0)  NewBlk=(int)(sphere[part].x*L/Lppx);
		else if(Dim==1) NewBlk=(int)(sphere[part].y*L/Lppy);
		else NewBlk=(int)(sphere[part].z*L/Lppz);
		if(NewBlk!=Me[Dim]){
			if(PlusId!=MinusId){
				if(NewBlk==PlusId){
					PlusSlabF[CountPlus]=(double)part;
					PlusSlabF[CountPlus+1]=sphere[part].x_u;
					PlusSlabF[CountPlus+2]=sphere[part].y_u;
					PlusSlabF[CountPlus+3]=sphere[part].z_u;
					PlusSlabF[CountPlus+4]=sphere[part].px;
					PlusSlabF[CountPlus+5]=sphere[part].py;
					PlusSlabF[CountPlus+6]=sphere[part].pz;
					CountPlus+=7;
					InMe[part]=false;
				}
				else if(NewBlk==MinusId){
					MinusSlabF[CountMinus]=(double)part;
					MinusSlabF[CountMinus+1]=sphere[part].x_u;
					MinusSlabF[CountMinus+2]=sphere[part].y_u;
					MinusSlabF[CountMinus+3]=sphere[part].z_u;
					MinusSlabF[CountMinus+4]=sphere[part].px;
					MinusSlabF[CountMinus+5]=sphere[part].py;
					MinusSlabF[CountMinus+6]=sphere[part].pz;
					CountMinus+=7;
					InMe[part]=false;
				}	
			}
			else {
				if(NewBlk==PlusId){
					PlusSlabF[CountPlus]=(double)part;
					PlusSlabF[CountPlus+1]=sphere[part].x_u;
					PlusSlabF[CountPlus+2]=sphere[part].y_u;
					PlusSlabF[CountPlus+3]=sphere[part].z_u;
					PlusSlabF[CountPlus+4]=sphere[part].px;
					PlusSlabF[CountPlus+5]=sphere[part].py;
					PlusSlabF[CountPlus+6]=sphere[part].pz;
					CountPlus+=7;
					InMe[part]=false;
				}
			}
		}		
	}
return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
//						Communicate slabs for extended block
//		switch Refill defines if you make new slabs or just refill the old ones
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
void CommunicateCoordinates(int Refill){
    MPI_Request sRequest, rRequest;
    MPI_Status sStatus, rStatus;
	int SendTo, RecvFrom,RecvCounter;
	OtherPartNum=0;
	if(Refill) Refillx();
	else GetSlabsx();
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			Send Right
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	if(Me[0]+1 == npx) SendTo = 0 + Me[1]*npx + Me[2]*npx*npy;			//send right slab 
	else SendTo  = Me[0]+1 + Me[1]*npx + Me[2]*npx*npy;					//send right slab 
	if(Me[0] == 0) RecvFrom = npx-1 + Me[1]*npx + Me[2]*npx*npy;		//recv left slab 
	else RecvFrom  = Me[0]-1 + Me[1]*npx + Me[2]*npx*npy;				//recv left slab 
	if(SendTo!=node){
		MPI_Isend(RightSlab,CountR,MPI_DOUBLE,SendTo,0,MPI_COMM_WORLD,&sRequest);
		MPI_Irecv(BuffSlab,MaxSlabPart,MPI_DOUBLE,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&rRequest);
	    MPI_Wait(&sRequest,&sStatus );
		MPI_Wait(&rRequest,&rStatus );
		MPI_Get_count(&rStatus ,MPI_DOUBLE,&RecvCounter);
		for(int i=0;i<RecvCounter;i=i+4){
			int part=(int)BuffSlab[i];
			sphere[part].x=BuffSlab[i+1];
			sphere[part].y=BuffSlab[i+2];
			sphere[part].z=BuffSlab[i+3];
			InMe[part]=true;
			OtherPart[OtherPartNum]=part;
			OtherPartNum++;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			Send Left
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	swap(&RecvFrom,&SendTo);													//Roles are now reversed
	if(SendTo!=node){
		MPI_Isend(LeftSlab,CountL,MPI_DOUBLE,SendTo,0,MPI_COMM_WORLD,&sRequest );
		MPI_Irecv(BuffSlab,MaxSlabPart,MPI_DOUBLE,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&rRequest );
	    MPI_Wait(&sRequest ,&sStatus );
		MPI_Wait(&rRequest ,&rStatus );
		MPI_Get_count(&rStatus ,MPI_DOUBLE,&RecvCounter);
		for(int i=0;i<RecvCounter;i=i+4){
			int part=(int)BuffSlab[i];
			sphere[part].x=BuffSlab[i+1];
			sphere[part].y=BuffSlab[i+2];
			sphere[part].z=BuffSlab[i+3];
			InMe[part]=true;
			OtherPart[OtherPartNum]=part;
			OtherPartNum++;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if(Refill) Refilly();
	else GetSlabsy();
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			Send Top
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	if(Me[1]+1 == npy) SendTo = Me[0] + 0*npx + Me[2]*npx*npy;				//send top slab 
	else	SendTo  = Me[0] + (Me[1]+1)*npx + Me[2]*npx*npy;					//send top slab 
	if(Me[1] == 0) RecvFrom = Me[0] + (npy-1)*npx + Me[2]*npx*npy;			//recv bottom slab 
	else RecvFrom  = Me[0] + (Me[1]-1)*npx + Me[2]*npx*npy;					//recv bottom slab 
	if(SendTo!=node){
		MPI_Isend(TopSlab,CountT,MPI_DOUBLE,SendTo,0,MPI_COMM_WORLD,&sRequest );
		MPI_Irecv(BuffSlab,MaxSlabPart,MPI_DOUBLE,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&rRequest );
	    MPI_Wait(&sRequest ,&sStatus );
		MPI_Wait(&rRequest ,&rStatus );
		MPI_Get_count(&rStatus ,MPI_DOUBLE,&RecvCounter);
		for(int i=0;i<RecvCounter;i=i+4){
			int part=(int)BuffSlab[i];
			sphere[part].x=BuffSlab[i+1];
			sphere[part].y=BuffSlab[i+2];
			sphere[part].z=BuffSlab[i+3];
			InMe[part]=true;
			OtherPart[OtherPartNum]=part;
			OtherPartNum++;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			Send Bottom
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	swap(&RecvFrom,&SendTo);													//Roles are now reversed
	if(SendTo!=node){
		MPI_Isend(BottomSlab,CountB,MPI_DOUBLE,SendTo,0,MPI_COMM_WORLD,&sRequest );
		MPI_Irecv(BuffSlab,MaxSlabPart,MPI_DOUBLE,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&rRequest );
	    MPI_Wait(&sRequest ,&sStatus );
		MPI_Wait(&rRequest ,&rStatus );
		MPI_Get_count(&rStatus ,MPI_DOUBLE,&RecvCounter);
	
		
		for(int i=0;i<RecvCounter;i=i+4){
			int part=(int)BuffSlab[i];
			sphere[part].x=BuffSlab[i+1];
			sphere[part].y=BuffSlab[i+2];
			sphere[part].z=BuffSlab[i+3];
			InMe[part]=true;
			OtherPart[OtherPartNum]=part;
			OtherPartNum++;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(Refill) Refillz();
	else GetSlabsz();
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			Send AwayMe
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	if(Me[2]+1 == npz) SendTo = Me[0] + Me[1]*npx + 0*npx*npy;				//send AwayMe slab 
	else SendTo  = Me[0] + Me[1]*npx + (Me[2]+1)*npx*npy;					//send AwayMe slab 
	if(Me[2] == 0) RecvFrom = Me[0] + Me[1]*npx + (npz-1)*npx*npy;			//recv ToMe slab 
	else RecvFrom  = Me[0] + Me[1]*npx + (Me[2]-1)*npx*npy;					//recv ToMe slab 
	MPI_Isend(AwayMeSlab,CountAwayMe,MPI_DOUBLE,SendTo,0,MPI_COMM_WORLD,&sRequest );
	MPI_Irecv(BuffSlab,MaxSlabPart,MPI_DOUBLE,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&rRequest );
    MPI_Wait(&sRequest ,&sStatus );
	MPI_Wait(&rRequest ,&rStatus );
	MPI_Get_count(&rStatus ,MPI_DOUBLE,&RecvCounter);
		for(int i=0;i<RecvCounter;i=i+4){
			int part=(int)BuffSlab[i];
			sphere[part].x=BuffSlab[i+1];
			sphere[part].y=BuffSlab[i+2];
			sphere[part].z=BuffSlab[i+3];
			InMe[part]=true;
			OtherPart[OtherPartNum]=part;
			OtherPartNum++;
		}
	MPI_Barrier(MPI_COMM_WORLD);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			Send ToMe
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	swap(&RecvFrom,&SendTo);													//Roles are now reversed
	MPI_Isend(ToMeSlab,CountToMe,MPI_DOUBLE,SendTo,0,MPI_COMM_WORLD,&sRequest );
	MPI_Irecv(BuffSlab,MaxSlabPart,MPI_DOUBLE,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&rRequest );
    MPI_Wait(&sRequest ,&sStatus );
	MPI_Wait(&rRequest ,&rStatus );
	MPI_Get_count(&rStatus ,MPI_DOUBLE,&RecvCounter);
		for(int i=0;i<RecvCounter;i=i+4){
			int part=(int)BuffSlab[i];
			sphere[part].x=BuffSlab[i+1];
			sphere[part].y=BuffSlab[i+2];
			sphere[part].z=BuffSlab[i+3];
			InMe[part]=true;
			OtherPart[OtherPartNum]=part;
			OtherPartNum++;
		}
	MPI_Barrier(MPI_COMM_WORLD);
return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
void GetSlabsx(){
	CountL=0;
	CountR=0;
	for(int i=0;i<NumInBlock[node];i++){
		int part=InBlock[node][i];
		if(sphere[part].x>=MyBoundS[0] && sphere[part].x<=MyBoundS[0]+Slabdx){
			LeftSlab[CountL]=(double)part;
			LeftSlab[CountL+1]=sphere[part].x;
			LeftSlab[CountL+2]=sphere[part].y;
			LeftSlab[CountL+3]=sphere[part].z;
			CountL+=4;
		}
		else if(sphere[part].x>=MyBoundE[0]-Slabdx && sphere[part].x<=MyBoundE[0]){
			RightSlab[CountR]=(double)part;
			RightSlab[CountR+1]=sphere[part].x;
			RightSlab[CountR+2]=sphere[part].y;
			RightSlab[CountR+3]=sphere[part].z;
			CountR+=4;
		}
	}
return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
void GetSlabsy(){
	CountT=0;
	CountB=0;
	for(int i=0;i<NumInBlock[node];i++){
		int part=InBlock[node][i];
		if(sphere[part].y>=MyBoundS[1] && sphere[part].y<=MyBoundS[1]+Slabdy){
			BottomSlab[CountB]=(double)part;
			BottomSlab[CountB+1]=sphere[part].x;
			BottomSlab[CountB+2]=sphere[part].y;
			BottomSlab[CountB+3]=sphere[part].z;
			CountB+=4;
		}
		else if(sphere[part].y>=MyBoundE[1]-Slabdy && sphere[part].y<=MyBoundE[1]){
			TopSlab[CountT]=(double)part;
			TopSlab[CountT+1]=sphere[part].x;
			TopSlab[CountT+2]=sphere[part].y;
			TopSlab[CountT+3]=sphere[part].z;
			CountT+=4;
		}
	}

	for(int i=0;i<OtherPartNum;i++){
		int part=OtherPart[i];
		if(sphere[part].y>=MyBoundS[1] && sphere[part].y<=MyBoundS[1]+Slabdy){
			BottomSlab[CountB]=(double)part;
			BottomSlab[CountB+1]=sphere[part].x;
			BottomSlab[CountB+2]=sphere[part].y;
			BottomSlab[CountB+3]=sphere[part].z;
			CountB+=4;
		}
		else if(sphere[part].y>=MyBoundE[1]-Slabdy && sphere[part].y<=MyBoundE[1]){
			TopSlab[CountT]=(double)part;
			TopSlab[CountT+1]=sphere[part].x;
			TopSlab[CountT+2]=sphere[part].y;
			TopSlab[CountT+3]=sphere[part].z;
			CountT+=4;
		}
	}
return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
void GetSlabsz(){
	CountToMe=0;
	CountAwayMe=0;
	for(int i=0;i<NumInBlock[node];i++){
		int part=InBlock[node][i];
		if(sphere[part].z>=MyBoundS[2] && sphere[part].z<=MyBoundS[2]+Slabdz){
			ToMeSlab[CountToMe]=(double)part;
			ToMeSlab[CountToMe+1]=sphere[part].x;
			ToMeSlab[CountToMe+2]=sphere[part].y;
			ToMeSlab[CountToMe+3]=sphere[part].z;
			CountToMe+=4;
		}
		else if(sphere[part].z>=MyBoundE[2]-Slabdz && sphere[part].z<=MyBoundE[2]){
			AwayMeSlab[CountAwayMe]=(double)part;
			AwayMeSlab[CountAwayMe+1]=sphere[part].x;
			AwayMeSlab[CountAwayMe+2]=sphere[part].y;
			AwayMeSlab[CountAwayMe+3]=sphere[part].z;
			CountAwayMe+=4;
		}
	}
	for(int i=0;i<OtherPartNum;i++){
		int part=OtherPart[i];
		if(sphere[part].z>=MyBoundS[2] && sphere[part].z<=MyBoundS[2]+Slabdz){
			ToMeSlab[CountToMe]=(double)part;
			ToMeSlab[CountToMe+1]=sphere[part].x;
			ToMeSlab[CountToMe+2]=sphere[part].y;
			ToMeSlab[CountToMe+3]=sphere[part].z;
			CountToMe+=4;
		}
		else if(sphere[part].z>=MyBoundE[2]-Slabdz && sphere[part].z<=MyBoundE[2]){
			AwayMeSlab[CountAwayMe]=(double)part;
			AwayMeSlab[CountAwayMe+1]=sphere[part].x;
			AwayMeSlab[CountAwayMe+2]=sphere[part].y;
			AwayMeSlab[CountAwayMe+3]=sphere[part].z;
			CountAwayMe+=4;
		}
	}
return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
void Refillx(){
	for(int i=0;i<CountR;i=i+4){
		int part=(int)RightSlab[i];
		RightSlab[i+1]=sphere[part].x;
		RightSlab[i+2]=sphere[part].y;
		RightSlab[i+3]=sphere[part].z;

	}

	for(int i=0;i<CountL;i=i+4){
		int part=(int)LeftSlab[i];
		LeftSlab[i+1]=sphere[part].x;
		LeftSlab[i+2]=sphere[part].y;
		LeftSlab[i+3]=sphere[part].z;
	}
return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
void Refilly(){
	for(int i=0;i<CountT;i=i+4){
		int part=(int)TopSlab[i];
		TopSlab[i+1]=sphere[part].x;
		TopSlab[i+2]=sphere[part].y;
		TopSlab[i+3]=sphere[part].z;
	}
	for(int i=0;i<CountB;i=i+4){
		int part=(int)BottomSlab[i];
		BottomSlab[i+1]=sphere[part].x;
		BottomSlab[i+2]=sphere[part].y;
		BottomSlab[i+3]=sphere[part].z;
	}
return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
void Refillz(){
	for(int i=0;i<CountToMe;i=i+4){
		int part=(int)ToMeSlab[i];
		ToMeSlab[i+1]=sphere[part].x;
		ToMeSlab[i+2]=sphere[part].y;
		ToMeSlab[i+3]=sphere[part].z;
	}
	for(int i=0;i<CountAwayMe;i=i+4){
		int part=(int)AwayMeSlab[i];
		AwayMeSlab[i+1]=sphere[part].x;
		AwayMeSlab[i+2]=sphere[part].y;
		AwayMeSlab[i+3]=sphere[part].z;
	}
return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
void BcastInitState(){
	MPI_Bcast(&L, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	SLABF *Buff;
	Buff=calloc(N,sizeof(SLABF));
	if(node==0){
		for(int i=0;i<N;i++){
			Buff[i].part=sphere[i].type;
			Buff[i].x=sphere[i].x;
			Buff[i].y=sphere[i].y;
			Buff[i].z=sphere[i].z;
			Buff[i].px=sphere[i].px;
			Buff[i].py=sphere[i].py;
			Buff[i].pz=sphere[i].pz;
		}
	}
	int ierr=MPI_Bcast(Buff,N,MPI_SLABF,0,MPI_COMM_WORLD);
	for(int i=0;i<N;i++){
		sphere[i].type=Buff[i].part;
		sphere[i].x=Buff[i].x;
		sphere[i].y=Buff[i].y;
		sphere[i].z=Buff[i].z;
		sphere[i].px=Buff[i].px;
		sphere[i].py=Buff[i].py;
		sphere[i].pz=Buff[i].pz;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	free(Buff);
    MPI_Bcast(activeParticleIndex, Nact, MPI_INT, 0, MPI_COMM_WORLD);
	int Tag[N];
	if(node==0){
		for(int i=0;i<N;i++)
			Tag[i]=sphere[i].actTag;	
	}
	MPI_Bcast(Tag, N, MPI_INT, 0, MPI_COMM_WORLD);
	for(int i=0;i<N;i++)
		sphere[i].actTag=Tag[i];	
return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
void AllocateArrays(){
	MaxSlabPart=10000;
	//MaxSlabPart=(int)(L/(double)npz)*(L/(double)npz)*DENSITY*Cell_cut+10000;
//	printf("%d\n",MaxSlabPart);
	TopSlab = (double *)malloc(4*MaxSlabPart*sizeof(double));
	BottomSlab = (double *)malloc(4*MaxSlabPart*sizeof(double));
	LeftSlab = (double *)malloc(4*MaxSlabPart*sizeof(double));
	RightSlab = (double *)malloc(4*MaxSlabPart*sizeof(double));
	ToMeSlab = (double *)malloc(4*MaxSlabPart*sizeof(double));
	AwayMeSlab = (double *)malloc(4*MaxSlabPart*sizeof(double));
	BuffSlab = (double *)malloc(4*MaxSlabPart*sizeof(double));
	
	PlusSlabF = (double *)malloc(7*MaxSlabPart*sizeof(double));
	MinusSlabF = (double *)malloc(7*MaxSlabPart*sizeof(double));
	BuffSlabF =(double *) malloc(7*MaxSlabPart*sizeof(double));

	OtherPart = (int *)malloc(N*sizeof(int));
	MigratedPart =(int *) malloc(1000*sizeof(int));

return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
void GatherCoordinatesToRoot(){
	MPI_Barrier(MPI_COMM_WORLD);
	if(node==0)
		MPI_Gather(MPI_IN_PLACE, 1, MPI_INT, NumInBlock, 1, MPI_INT, 0, MPI_COMM_WORLD);
	else
		MPI_Gather(&NumInBlock[node], 1, MPI_INT, NumInBlock, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast( NumInBlock, Nprocs, MPI_INT, 0, MPI_COMM_WORLD);

	SLAB *BuffSend;
	if(node!=0)BuffSend =(SLAB*)calloc(NumInBlock[node],sizeof(SLAB)); 
	else BuffSend=(SLAB*)calloc(N,sizeof(SLAB));	
	if(BuffSend==NULL) {fprintf(stderr,"RAM not available\nExiting Now \n"); exit(0);}
	for (int i=0;i<NumInBlock[node];i++){
		int part=InBlock[node][i];
		BuffSend[i].part=part;
		BuffSend[i].x=sphere[part].x_u;
		BuffSend[i].y=sphere[part].y_u;
		BuffSend[i].z=sphere[part].z_u;
	}
	int DispVec[Nprocs];
	if(node==0){
		DispVec[0]=0;
		for (int i=1;i<Nprocs;i++)
			DispVec[i]=DispVec[i-1]+NumInBlock[i-1];
		MPI_Gatherv( MPI_IN_PLACE, NumInBlock[node],MPI_SLAB,BuffSend,NumInBlock,DispVec,MPI_SLAB,0,MPI_COMM_WORLD);	
	}	
	else MPI_Gatherv( BuffSend, NumInBlock[node],MPI_SLAB,BuffSend,NumInBlock,DispVec,MPI_SLAB,0,MPI_COMM_WORLD);	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(node==0){
		for (int i=0;i<N;i++){
			int part=BuffSend[i].part;
			sphere[part].x_u=BuffSend[i].x;
			sphere[part].y_u=BuffSend[i].y;
			sphere[part].z_u=BuffSend[i].z;
			sphere[part].x=BuffSend[i].x;
			sphere[part].y=BuffSend[i].y;
			sphere[part].z=BuffSend[i].z;
		}
	applypbc();
	}
	free(BuffSend);

return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
void GatherCoordAndMomentaToRoot(){
	if(node==0)
		MPI_Gather(MPI_IN_PLACE, 1, MPI_INT, NumInBlock, 1, MPI_INT, 0, MPI_COMM_WORLD);
	else
		MPI_Gather(&NumInBlock[node], 1, MPI_INT, NumInBlock, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast( NumInBlock, Nprocs, MPI_INT, 0, MPI_COMM_WORLD);

	SLABF *BuffSend;
	if(node!=0)BuffSend =calloc(NumInBlock[node],sizeof(SLABF)); 
	else BuffSend=calloc(N,sizeof(SLABF));	


	for (int i=0;i<NumInBlock[node];i++){
		int part=InBlock[node][i];
		BuffSend[i].part=part;
		BuffSend[i].x=sphere[part].x_u;
		BuffSend[i].y=sphere[part].y_u;
		BuffSend[i].z=sphere[part].z_u;
		BuffSend[i].px=sphere[part].px;
		BuffSend[i].py=sphere[part].py;
		BuffSend[i].pz=sphere[part].pz;
	}
	int DispVec[Nprocs];
	if(node==0){
		DispVec[0]=0;
		for (int i=1;i<Nprocs;i++)
			DispVec[i]=DispVec[i-1]+NumInBlock[i-1];
		MPI_Gatherv( MPI_IN_PLACE, NumInBlock[node],MPI_SLABF,BuffSend,NumInBlock,DispVec,MPI_SLABF,0,MPI_COMM_WORLD);	
	}	
	else MPI_Gatherv( BuffSend, NumInBlock[node],MPI_SLABF,BuffSend,NumInBlock,DispVec,MPI_SLABF,0,MPI_COMM_WORLD);	
	MPI_Barrier(MPI_COMM_WORLD);
	if(node==0){
		for (int i=0;i<N;i++){
			int part=BuffSend[i].part;
			sphere[part].x=BuffSend[i].x;
			sphere[part].y=BuffSend[i].y;
			sphere[part].z=BuffSend[i].z;
			sphere[part].x_u=BuffSend[i].x;
			sphere[part].y_u=BuffSend[i].y;
			sphere[part].z_u=BuffSend[i].z;
			sphere[part].px=BuffSend[i].px;
			sphere[part].py=BuffSend[i].py;
			sphere[part].pz=BuffSend[i].pz;
		}
	applypbc();
	}
	free(BuffSend);

return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
void GetDataTypeSLAB(){
    //for struct SLAB
    int nroItems = 4;
    int blockLengths[4] = { 1, 1, 1, 1};
    MPI_Datatype types[4] = { MPI_INT,MPI_DOUBLE, MPI_DOUBLE,MPI_DOUBLE};
    MPI_Aint     offsets[4];
    offsets[0] = offsetof(SLAB, part);
    offsets[1] = offsetof(SLAB, x);
    offsets[2] = offsetof(SLAB, y);
    offsets[3] = offsetof(SLAB, z);
    MPI_Type_create_struct(nroItems, blockLengths, offsets, types, &MPI_SLAB);
    MPI_Type_commit(&MPI_SLAB);
    return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
void GetDataTypeSLABF(){
    //for struct SLABF
    int nroItems = 7;
    int blockLengths[7] = { 1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype types[7] = { MPI_INT,MPI_DOUBLE, MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE, MPI_DOUBLE,MPI_DOUBLE};
    MPI_Aint     offsets[7];
    offsets[0] = offsetof(SLABF, part);
    offsets[1] = offsetof(SLABF, x);
    offsets[2] = offsetof(SLABF, y);
    offsets[3] = offsetof(SLABF, z);
    offsets[4] = offsetof(SLABF, px);
    offsets[5] = offsetof(SLABF, py);
    offsets[6] = offsetof(SLABF, pz);
    MPI_Type_create_struct(nroItems, blockLengths, offsets, types, &MPI_SLABF);
    MPI_Type_commit(&MPI_SLABF);
    return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

