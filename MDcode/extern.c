#include "global.h"
#include "subroutine.h"
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	SPHERE sphere[N];
	int NeibListCount[N];
	int NeibList[N][MaxNeibs];
	int InCellCount[MaxCells];
	int InCell[MaxCells][MaxPartInCell];
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		double pe,ke,virial,peblk,keblk,virialblk;
		int neib_list_counter;							//counter for no of times neib list have'nt be calculated
		int	Ncellx,Ncelly,Ncellz;						//# of cells in a block along x,y,z direction
		int	NCell;										//Total # of cells in whole simulation box
		double CellSizex,CellSizey,CellSizez;			//Cell size in absolute units
		double InvCellSizex,InvCellSizey,InvCellSizez;	//Inv Cell size in relative units
		int	InExtdBlk[N/2],NumInExtdBlk;
	 	int cell_index_each_part[N];
		double maxD,maxmaxD,globalMaxD;
		double inst_T;
		double inv_L;
    	double pressure;
    	double DENSITY;
    	double L ;
	int	MAP[MaxCells][13];
		double T,P;
		char main_dir[128];
	 int NCellOld;
		char sub_dir[128];
		int wanaa_restart_run, wanaa_restart_eq;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//				Nose-Hoover Stuff
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     double eta[3],p_eta[3],q[3];
     double Heff,tau,tau_baro;
     double eta_baro[3],p_eta_baro[3],q_baro[3];
     double p_eps,w_eps,L0,eps,PV;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//				MPI-Parallel Stuff
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	 int jobno;
     gsl_rng *rng_ptr; // pointer to random number generator (rng) 
	 int npx,npy,npz;											//number of processes along x,y,z
	 double Lppx,Lppy,Lppz;									//Length per processor along x,y,z (in actual units)
	 int node,Nprocs;												//present block
	 int Me[3];												//Address of my present block
	 double MyBoundS[3],MyBoundE[3];							//Start & End Bounds of by block (in reduced units)
	 int InBlock[MaxProcs][N]	;								//particle index in block
	 int NumInBlock[MaxProcs];								//# of particles in block
	 int Block[N];											//Block # of a particle
	 bool InMe[N];											//Particle tag if its in block 
	int NumCellInBlk[MaxProcs];								//# of cells (for listing) in a block
	int CellInBlk[MaxProcs][MaxCellInBlk];
	 int MaxSlabPart;											//Max particles in slab 
	 double *TopSlab,*BottomSlab,*LeftSlab;						//Different Slabs
	 double *RightSlab,*ToMeSlab,*AwayMeSlab,*BuffSlab;			//Different Slabs
	 int CountR,CountL,CountT,CountB,CountToMe,CountAwayMe;	//Counter for various slabs
	 double *PlusSlabF, *MinusSlabF,*BuffSlabF;				//Different Slabs contained pos + momenta's
	 int CountPlus,CountMinus;								//counter for no partices to migrace from Plos and Minus slab
	 int *OtherPart,OtherPartNum;								//other particles in block (comming from slab)
	 int MigratedPartNum,*MigratedPart;						//migrated particles
	 double Slabdx,Slabdy,Slabdz;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double vir_actblk,vir_act;
	int activeParticleIndex[Nact];
	int kx[Nact];
	int ky[Nact];
	int kz[Nact];
	double f0;
	double *xold,*yold,*zold;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	double fxRandom[N],fyRandom[N],fzRandom[N];
	double fxold[N],fyold[N],fzold[N];
	double rxold[N],ryold[N],rzold[N];
	double InvD0;
	double D0ByT;
	double TByD0;
	double expFactor;
	double *fxRandomOld; 
	double *fyRandomOld; 
	double *fzRandomOld; 
	double SqrtVariance;
	double SqrtActVariance;
