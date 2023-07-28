#ifndef GLOBAL_H
	#define GLOBAL_H
	#include <stdio.h>
	#include <math.h>
	#include <time.h>
	#include <stdlib.h>
	#include <assert.h>
	#include <unistd.h>
	#include <stdbool.h>
    #include "gsl/gsl_rng.h"
    #include "gsl/gsl_randist.h"
	#include "mpi.h"
//	#include "/apps/openmpi/include/mpi.h"
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#define N					50000							//No of spheres
	#define NA					(int)(.80*(double)N)			//No of type A spheres
	#define NB					(int)(.200*(double)N)			//No of type B spheres
	#define pi					3.14159265358979		
	#define MaxNeibs 			500								//anaticipated no of neibs
	#define MaxPartInCell		600								//anticipates no of particles in cell
	#define MaxCells			(int)(N/15)
	#define MaxCellInBlk		(int)(N/15)
	#define half_neibs	  		13								// for 3D
	#define Sigma				1.00							//
	#define Lcuttoff			( 2.5*Sigma )					//pot, cuttoff length
	#define Lcutoffsqr			( 6.25*Sigma*Sigma )					
	#define LENGTH				1.0								//scaled domain size
	#define HALF_LENGTH			0.5								//hlaf of it
	#define skin				( .35*Sigma )					//skin for cell list
	#define	Cell_cut 			( Lcuttoff+skin )				
	#define Cell_cut_sq			(Cell_cut*Cell_cut)
//	#define Cell_cut_sq			8.1225
	#define time_step			0.005
    #define half_time_step      ( time_step*0.5 )
    #define DOF					(3.0 * (double)N-3.0)			//Degree Of Freedom
	//#define C0					0.161960944271413609385490417480468750						
	//#define C2					-3.88062021229416131973266601562500000E-0002 
	//#define C4					2.48050456866621971130371093750000000E-0003
	#define C0					0.0162665594880					
	#define C2					-0.001949973872640
    #define tau_T               0.2								//Therostat stuff
    #define tau_P               2.0								//Barostat stuff
	#define DIVIDE_TO			2.0								//frac of maximal log interval
	#define	DIM					3	
	#define M					3
	#define MaxProcs			64
	


	#define D0						1.0 
	#define UNI                 	((double)rand()/((double)RAND_MAX + 1.0))
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//													Active Stuff
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#define ACTIVE_DENSITY			0.1000 
	#define Nact					(int)((double)N*ACTIVE_DENSITY) 
//	#define f0						2.500 							//active force factor
	#define TAU_per					200								//persistence time steps, time=Tau_per*time_step
	#define D0						1.0 

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	typedef struct{
		double x , y , z;
		double x_u, y_u, z_u ; 									//unfolded coordinates
		double px , py , pz;
		double Fx , Fy , Fz;
		int    type,actTag;
	}SPHERE;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	typedef struct{
		double x , y , z;
		double px , py , pz;
		int	part;
	}SLABF;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	typedef struct{
		double x , y , z;
		int	part;
	}SLAB;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	extern	SPHERE sphere[N];
	extern	int NeibListCount[N];
	extern 	int NeibList[N][MaxNeibs];
	extern	int InCellCount[MaxCells];
	extern	int InCell[MaxCells][MaxPartInCell];
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    extern gsl_rng *rng_ptr; // pointer to random number generator (rng) 
	extern	double pe,ke,virial,peblk,keblk,virialblk;
	extern	int neib_list_counter;							//counter for no of times neib list have'nt be calculated
	extern	int	Ncellx,Ncelly,Ncellz;						//# of cells in a block along x,y,z direction
	extern	int	NCell;										//Total # of cells in whole simulation box
	extern	double CellSizex,CellSizey,CellSizez;			//Cell size in absolute units
	extern	double InvCellSizex,InvCellSizey,InvCellSizez;	//Inv Cell size in relative units
	extern	int	InExtdBlk[N/2],NumInExtdBlk;
	extern 	int cell_index_each_part[N];
	extern	double maxD,maxmaxD,globalMaxD;
	extern	double inst_T;
	extern	double inv_L;
    extern	double pressure;
    extern	double DENSITY;
    extern	double L ;
	extern	int	MAP[MaxCells][13];
	extern	double T,P;
	extern	char main_dir[128];
	extern int NCellOld;
	extern int jobno;
	extern	char sub_dir[128];
	extern	int wanaa_restart_run, wanaa_restart_eq;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//				Nose-Hoover Stuff
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    extern double eta[3],p_eta[3],q[3];
    extern double Heff,tau,tau_baro;
    extern double eta_baro[3],p_eta_baro[3],q_baro[3];
    extern double p_eps,w_eps,L0,eps,PV;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//				MPI-Parallel Stuff
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	MPI_Datatype MPI_SLAB,MPI_SLABF;
	extern int npx,npy,npz;											//number of processes along x,y,z
	extern double Lppx,Lppy,Lppz;									//Length per processor along x,y,z (in actual units)
	extern int node,Nprocs;												//present block
	extern int Me[3];												//Address of my present block
	extern double MyBoundS[3],MyBoundE[3];							//Start & End Bounds of by block (in reduced units)
	extern int InBlock[MaxProcs][N]	;								//particle index in block
	extern int NumInBlock[MaxProcs];								//# of particles in block
	extern int Block[N];											//Block # of a particle
	extern bool InMe[N];											//Particle tag if its in block 
	extern int NumCellInBlk[MaxProcs];								//# of cells (for listing) in a block
	extern int CellInBlk[MaxProcs][MaxCellInBlk];
	extern int MaxSlabPart;											//Max particles in slab 
	extern double *TopSlab,*BottomSlab,*LeftSlab;						//Different Slabs
	extern double *RightSlab,*ToMeSlab,*AwayMeSlab,*BuffSlab;			//Different Slabs
	extern int CountR,CountL,CountT,CountB,CountToMe,CountAwayMe;	//Counter for various slabs
	extern double *PlusSlabF, *MinusSlabF,*BuffSlabF;				//Different Slabs contained pos + momenta's
	extern int CountPlus,CountMinus;								//counter for no partices to migrace from Plos and Minus slab
	extern int *OtherPart,OtherPartNum;								//other particles in block (comming from slab)
	extern int MigratedPartNum,*MigratedPart;						//migrated particles
	extern double Slabdx,Slabdy,Slabdz;
	extern int exch;
//										Active Variables
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	extern double vir_actblk,vir_act;
	extern int activeParticleIndex[Nact];
	extern int kx[Nact];
	extern int ky[Nact];
	extern int kz[Nact];
	extern double f0;
	extern	double *xold,*yold,*zold;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	extern	double InvD0;
	extern	double D0ByT;
	extern	double TByD0;
	extern	double expFactor;
	extern	double fxRandom[N],fyRandom[N],fzRandom[N];
	extern 	double fxold[N],fyold[N],fzold[N];
	extern 	double rxold[N],ryold[N],rzold[N];
	extern double *fxRandomOld; 
	extern double *fyRandomOld; 
	extern double *fzRandomOld; 
	extern double SqrtVariance;
	extern double SqrtActVariance;
#endif
