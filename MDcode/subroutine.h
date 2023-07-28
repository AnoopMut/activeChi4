#ifndef SUBROUTINE_H
	#define SUBROUTINE_H
	
	//in file cell_list+forces.c	
	int 	getcell_I();
	void	maps();
	void	update_cell_list();
	bool	len_check();
	void	calculate_forces();

	//in file initialise.c
	void initializeGsl();
	int 	randint();
	void 	initialise_system();	
	void 	initializeSystem();
	void 	fixDrift();
	void    adjustInitialTemp();
	void    calculate_temp();
	void	breakxal();
	void	coldstart();
	double random_normal();    
	double drand();
   	double normal();
	double uniform(); 
	void	swap(); 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	void initializeOverDampedBD();
	//in time_evolve.c
	void 	applypbc1();
	void 	applypbc();
	void 	advance_time_NVEVV();
	void 	advance_time_NVEFrog();
	void    advance_time_NVT_Berendsen();
	void    advance_time_NPT();
	void	U1Propagator();
	void	U2Propagator();
	void	U3Propagator();
	void	U4Propagator();
	void	InitialiseNoseHoover();
	void AdvanceTimeNVTNoseHoover();
	void	U1PropagatorNPT();
	void	U2PrimePropagatorNPT();
	void	U2PropagatorNPT();
	void	U3PropagatorNPT();
	void	U4PropagatorNPT();
	void	InitialiseNoseHooverNPT();
	void AdvanceTimeNPTNoseHoover();
	void UpdateBoxForNPT();
	
	//in file print.c
	double print_data();
	double print_eq();
	double print_movie();
	double print_2col();
	double print_pe();
	double print_stuff();
	void print_col();
	
	
    //in file main.c
	void shellSort();
	void prepareSavingArray();
	void equilibrate();
	void run();
	void write_restart_state();
	void read_restart_state();
	void initialise_restart_state();
	void chop_off_data_files();
	int main();
	void AdvanceTimeNVTGaussOSM();

	//in File parallel.c
	void SetUpParallel();
	void SetUpCellList();
	void ExchangeParticlesInBlocks();
	void CommunicateCoordinates( );
	void GetSlabsx();
	void GetSlabsy();
	void GetSlabsz();
	void Refillx();
	void Refilly();
	void Refillz();
	void GetSlabsForExchange();
	void BcastInitState();
	void AllocateArrays();
	void GetDataTypeSLAB();
	void GetDataTypeSLABF();
	void GatherCoordinatesToRoot();
	void GatherCoordAndMomentaToRoot();
	//in file active.c
	void selectRandomActiveParticles()	;
	void ShuffleActDirec();
	void saveActiveDirection();
	void readActiveDirection();
	void readActiveParticles();
#endif
