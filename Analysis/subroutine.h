#ifndef SUBROUTINE_H
	#define SUBROUTINE_H
	void LoadFramePP();

    //in file main.c
	int main();
	void read_data_file();
	void load_frame();
	void mean_square_displacement();
	void mean_square_theta_displacement();
	void overlap_correlation_function();
	void print_ensemble_avg();
	void print_time_cut_things();
	void overlap_cage();
	void S4Qt();
	void guu();
	void GatherVec();
	int	GetQ();
	void Fsqt();
	void SplitTheSystem();
	void LoadOrigin();
	void GetTheOriginConfigs();
	void BCastConstants();
#endif
