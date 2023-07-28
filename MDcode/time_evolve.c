		/*********************************************************************************************************
		**	This bit of code is used to time evolve the system using various integrating schemes i.e.         	**
		**	velocity verlet(NVE); NVT(using two different thermostats Brown-Clarke and Berendsen); and NPT  	**
		**										$$$ROUTINES INCLUDED$$$											**
		**	applypbc1		==> check if new pos. of spheres are in 0-(LENGTH=1)								**
		** 						it contains  inf while loops so do check this if it takes much time				**
		**	advance_time_NVE==> uses velocity verlet to give new pos. of spheres  								**
		**						using the forces on spheres and new vel.										**
		**	advance_time_NVT_Berendsen 		==>	same as above but is NVT with Leap-Frog							**
		**	advance_time_NPT				==>	same as above but is NPT with Leap-Frog							**
		**********************************************************************************************************/
#include "global.h"
#include "subroutine.h"
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void applypbc1(){
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		while (sphere[part].x >= LENGTH) sphere[part].x -=LENGTH;
		while (sphere[part].x < 0.0 )    sphere[part].x +=LENGTH;
		while (sphere[part].y >= LENGTH) sphere[part].y -=LENGTH;
		while (sphere[part].y < 0.0 )    sphere[part].y +=LENGTH;
		while (sphere[part].z >= LENGTH) sphere[part].z -=LENGTH;
		while (sphere[part].z < 0.0 )    sphere[part].z +=LENGTH;            
	}
return;
}
void applypbc(){
	for(int part=0; part<N; part++){
		while (sphere[part].x >= LENGTH) sphere[part].x -=LENGTH;
		while (sphere[part].x < 0.0 )    sphere[part].x +=LENGTH;
		while (sphere[part].y >= LENGTH) sphere[part].y -=LENGTH;
		while (sphere[part].y < 0.0 )    sphere[part].y +=LENGTH;
		while (sphere[part].z >= LENGTH) sphere[part].z -=LENGTH;
		while (sphere[part].z < 0.0 )    sphere[part].z +=LENGTH;            
	}
return;
}
//**********************************Velocity-Verlet********************************************//
void advance_time_NVEVV(){                                  
	double dx,dy,dz, temp;
	maxD=0.0;
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		dx = (sphere[part].px+0.5*sphere[i].Fx*time_step)*time_step; 
		dy = (sphere[part].py+0.5*sphere[i].Fy*time_step)*time_step; 
		dz = (sphere[part].pz+0.5*sphere[i].Fz*time_step)*time_step; 
		sphere[part].x		+= dx*inv_L;
		sphere[part].x_u 	+= dx*inv_L;													//unfolded
		sphere[part].y  	+= dy*inv_L;
		sphere[part].y_u	+= dy*inv_L;													//unfolded
		sphere[part].z  	+= dz*inv_L;
		sphere[part].z_u 	+= dz*inv_L;

		sphere[part].px += sphere[part].Fx*half_time_step;
		sphere[part].py += sphere[part].Fy*half_time_step; 
		sphere[part].pz += sphere[part].Fz*half_time_step; 
	}
	applypbc1();
	maxD=0.0;
	maxmaxD=0.0;
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		dx=sphere[part].x_u-xold[part];
		dy=sphere[part].y_u-yold[part];
		dz=sphere[part].z_u-zold[part];
        temp = dx*dx + dy*dy +dz*dz;
        if(maxD < temp) maxD = temp;
	}
	maxD*=L*L;
	MPI_Allreduce(&maxD, &maxmaxD, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	int refill=1;
	//if(node==0)printf("%d\t%lf\t%lf\n",node,maxmaxD,skin*skin);
	if(3.0*maxmaxD>skin*skin){
		refill=0;
		ExchangeParticlesInBlocks();
		CommunicateCoordinates(refill);
		update_cell_list();
	}
	else CommunicateCoordinates(refill);
	calculate_forces();
	keblk=0.0e0;
//	ke=0.0;
	for(int i=0 ; i<NumInBlock[node] ; i++ ){
		int part=InBlock[node][i];
		sphere[part].px += sphere[part].Fx*half_time_step;
		sphere[part].py += sphere[part].Fy*half_time_step; 
		sphere[part].pz += sphere[part].Fz*half_time_step; 
		keblk += ( sphere[part].px*sphere[part].px + 
					sphere[part].py*sphere[part].py + 
					sphere[part].pz*sphere[part].pz );
	}

//	MPI_Reduce(&keblk, &ke, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//	MPI_Bcast(&ke, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//	ke    = ke * 0.5;
//	inst_T=2.0*ke/DOF;
	return;
}



//********************************************LEAP-FROG*************************************************//
void advance_time_NVEFrog(){                                  
	double dx,dy,dz,temp;
	maxD=0.0;
	//-------------------------------half-kick----------------------------------------
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		sphere[part].px += sphere[part].Fx*half_time_step;
		sphere[part].py += sphere[part].Fy*half_time_step; 
		sphere[part].pz += sphere[part].Fz*half_time_step; 
	}
	//-------------------------------free-drift---------------------------------------
	
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		dx = sphere[part].px*time_step; 
		dy = sphere[part].py*time_step; 
		dz = sphere[part].pz*time_step; 
		sphere[part].x		+= dx*inv_L;
		sphere[part].x_u 	+= dx*inv_L;													//unfolded
		sphere[part].y  	+= dy*inv_L;
		sphere[part].y_u	+= dy*inv_L;													//unfolded
		sphere[part].z  	+= dz*inv_L;
		sphere[part].z_u 	+= dz*inv_L;
	}
	applypbc1();
	maxD=0.0;
	maxmaxD=0.0;
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		dx=sphere[part].x_u-xold[part];
		dy=sphere[part].y_u-yold[part];
		dz=sphere[part].z_u-zold[part];
        temp = dx*dx + dy*dy +dz*dz;
        if(maxD < temp) maxD = temp;
	}
	maxD*=L*L;
	MPI_Allreduce(&maxD, &maxmaxD, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	int refill=1;
	//if(node==0)printf("%d\t%lf\t%lf\n",node,maxmaxD,skin*skin);
	if(3.0*maxmaxD>skin*skin){
		refill=0;
		ExchangeParticlesInBlocks();
		CommunicateCoordinates(refill);
		update_cell_list();
	}
	else CommunicateCoordinates(refill);
	calculate_forces();

	//-------------------------------half-kick----------------------------------------
	keblk=0.0e0;
//	ke=0.0;
	for(int i=0 ; i<NumInBlock[node] ; i++ ){
		int part=InBlock[node][i];
		sphere[part].px += sphere[part].Fx*half_time_step;
		sphere[part].py += sphere[part].Fy*half_time_step; 
		sphere[part].pz += sphere[part].Fz*half_time_step; 
		keblk += ( sphere[part].px*sphere[part].px + 
					sphere[part].py*sphere[part].py + 
					sphere[part].pz*sphere[part].pz );
	}

//	MPI_Reduce(&keblk, &ke, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//	MPI_Bcast(&ke, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//	ke    = ke * 0.5;
//	inst_T=2.0*ke/DOF;
	return;
}



// *****************************Leap-Frog with Berendsen thermostat****************************
void advance_time_NVT_Berendsen(){
	double dx , dy , dz, temp;
	double chi_T;
	maxD=0.00e0;
	//~~~~~~~~~~~~~~~~~~~~~HALF UPDATE VELOCITIES AT TIME t from time t-dt/2~~~~~~~~~~~~~~~~~~~//
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		sphere[part].px += sphere[part].Fx*half_time_step;
		sphere[part].py += sphere[part].Fy*half_time_step; 
		sphere[part].pz += sphere[part].Fz*half_time_step; 
	}
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~get positions using that velocities~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		dx = sphere[part].px*time_step; 
		dy = sphere[part].py*time_step; 
		dz = sphere[part].pz*time_step; 
		sphere[part].x		+= dx*inv_L;
		sphere[part].x_u 	+= dx*inv_L;													//unfolded
		sphere[part].y  	+= dy*inv_L;
		sphere[part].y_u	+= dy*inv_L;													//unfolded
		sphere[part].z  	+= dz*inv_L;
		sphere[part].z_u 	+= dz*inv_L;
		temp = dx*dx + dy*dy + dz*dz;
		if(maxD < temp) maxD = temp;
	}
	maxD=sqrt(maxD);
	applypbc1();
	neib_list_counter ++;
	MPI_Reduce(&maxD, &maxmaxD, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(node==0)
		globalMaxD+=maxmaxD;
	MPI_Bcast(&globalMaxD,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	int refill=1;
	

	CommunicateCoordinates(refill);
	//printf("Comm Done\n")	;
	if(2.0*globalMaxD>skin){
		refill=0;
		//printf("Exchange Started\n")	;
		ExchangeParticlesInBlocks();
//		if(node==0)printf("Exchange Done\n")	;
		CommunicateCoordinates(refill);
		update_cell_list();
	}

	calculate_forces();
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~getting velocities at time t+dt/2~~~~~~~~~~~~~~~~~~~~~~~~~//
	keblk=0.0e0;
//	ke=0.0;
	for(int i=0 ; i<NumInBlock[node] ; i++ ){
		int part=InBlock[node][i];
		sphere[part].px += sphere[part].Fx*half_time_step;
		sphere[part].py += sphere[part].Fy*half_time_step; 
		sphere[part].pz += sphere[part].Fz*half_time_step; 
		keblk += ( sphere[part].px*sphere[part].px + 
					sphere[part].py*sphere[part].py + 
					sphere[part].pz*sphere[part].pz );
	}

	MPI_Reduce(&keblk, &ke, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ke, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	ke    = ke * 0.5;
	//~~~~~~~~~~~~~~~~CALCULATE THE VELOCITY SCALING FACTOR CHI FOR CNST. TEMP~~~~~~~~~~~~~~~//
	inst_T=2.0*ke/DOF;
	chi_T=sqrt(1.00 + (T/inst_T-1.00) * time_step/tau_T);
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~rescaling velocities~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	for(int i=0 ; i<NumInBlock[node] ; i++ ){
		int part=InBlock[node][i];
		sphere[part].px *= chi_T;
		sphere[part].py *= chi_T;
		sphere[part].pz *= chi_T;
	}
	return;
}

// *********************Leap-Frog with Berendsen thermostat and Barostat************************
void advance_time_NPT(){
	double dx,dy,dz,temp;
	double chi_T;
    MPI_Reduce(&virialblk, &virial, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&virial, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	calculate_temp();
	pressure = DENSITY*inst_T + virial/(L*L*L);
	//~~~~~~~~~~~~~~length scaling factor~~~~~~~~~~~~~~~//
	double chi_P = pow((1.0 + time_step*(pressure-P)/tau_P),0.33333333);
	UpdateBoxForNPT();
	L = chi_P*L;
	calculate_forces();
	maxD=0.00e0;
	//~~~~~~~~~~~~~~~~~~~~~HALF UPDATE VELOCITIES AT TIME t from time t-dt/2~~~~~~~~~~~~~~~~~~~//
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		sphere[part].px += sphere[part].Fx*half_time_step;
		sphere[part].py += sphere[part].Fy*half_time_step; 
		sphere[part].pz += sphere[part].Fz*half_time_step; 
	}
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~get positions using that velocities~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		dx = sphere[part].px*time_step; 
		dy = sphere[part].py*time_step; 
		dz = sphere[part].pz*time_step; 
		sphere[part].x		+= dx*inv_L;
		sphere[part].x_u 	+= dx*inv_L;													//unfolded
		sphere[part].y  	+= dy*inv_L;
		sphere[part].y_u	+= dy*inv_L;													//unfolded
		sphere[part].z  	+= dz*inv_L;
		sphere[part].z_u 	+= dz*inv_L;
		temp = dx*dx + dy*dy + dz*dz;
		if(maxD < temp) maxD = temp;
	}
	maxD=sqrt(maxD);
	applypbc1();
	neib_list_counter ++;
	MPI_Reduce(&maxD, &maxmaxD, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(node==0)
		globalMaxD+=maxmaxD;
	MPI_Bcast(&globalMaxD,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	int refill=1;
	

	CommunicateCoordinates(refill);
	//printf("Comm Done\n")	;
	if(2.0*globalMaxD>skin){
		refill=0;
		//printf("Exchange Started\n")	;
		ExchangeParticlesInBlocks();
//		if(node==0)printf("Exchange Done\n")	;
		CommunicateCoordinates(refill);
		update_cell_list();
	}

	calculate_forces();
	ke=0.0e0;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~getting velocities at time t+dt/2~~~~~~~~~~~~~~~~~~~~~~~~~//
	keblk=0.0e0;
//	ke=0.0;
	for(int i=0 ; i<NumInBlock[node] ; i++ ){
		int part=InBlock[node][i];
		sphere[part].px += sphere[part].Fx*half_time_step;
		sphere[part].py += sphere[part].Fy*half_time_step; 
		sphere[part].pz += sphere[part].Fz*half_time_step; 
		keblk += ( sphere[part].px*sphere[part].px + 
					sphere[part].py*sphere[part].py + 
					sphere[part].pz*sphere[part].pz );
	}

	MPI_Reduce(&keblk, &ke, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ke, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	ke    = ke * 0.5;
	//~~~~~~~~~~~~~~~~CALCULATE THE VELOCITY SCALING FACTOR CHI FOR CNST. TEMP~~~~~~~~~~~~~~~//
	inst_T=2.0*ke/DOF;
	chi_T=sqrt(1.00 + (T/inst_T-1.00) * time_step/tau_T);
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~rescaling velocities~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	for(int i=0 ; i<NumInBlock[node] ; i++ ){
		int part=InBlock[node][i];
		sphere[part].px *= chi_T;
		sphere[part].py *= chi_T;
		sphere[part].pz *= chi_T;
	}
	return;
}

// *********************Gaussian Operator splitting NVT ensemble************************
void AdvanceTimeNVTGaussOSM(){
	double dx , dy , dz,temp;
	double tempX1blk = 0.0;
    double tempX2blk = 0.0;
	//maxD=0.0;
	double K0blk=0.e0;
	double tempX1 = 0.0;
    double tempX2 = 0.0;
	//maxD=0.0;
	double K0=0.e0;
	for(int i=0 ; i<NumInBlock[node] ; i++ ){
		int part=InBlock[node][i];
        tempX1blk += sphere[part].px*sphere[part].Fx + sphere[part].py*sphere[part].Fy + sphere[part].pz*sphere[part].Fz;
        tempX2blk += sphere[part].Fx*sphere[part].Fx + sphere[part].Fy*sphere[part].Fy + sphere[part].Fz*sphere[part].Fz;
		K0blk	+= sphere[part].px*sphere[part].px +  sphere[part].py*sphere[part].py + sphere[part].pz*sphere[part].pz ;
	}
    MPI_Reduce(&tempX1blk, &tempX1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tempX1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&tempX2blk, &tempX2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tempX2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&K0blk, &K0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&K0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	double alpha = sqrt(tempX2/K0);
    double zeta0 = -tempX1/K0;
	double beta = exp(-alpha*time_step*0.5);
    double gamma = (zeta0 - alpha)/(zeta0 + alpha);
   	double factor1 = (1.0-gamma)/(beta - gamma/beta);
   	double factor2 = (1.0+gamma-beta-gamma/beta)/((1.0-gamma)*alpha);
	
    //-------------------------------half-kick----------------------------------------
	keblk=0.0e0;
	for(int i=0 ; i<NumInBlock[node] ; i++ ){
		int part=InBlock[node][i];
		sphere[part].px = factor1*(sphere[part].px + factor2*sphere[part].Fx);
		sphere[part].py = factor1*(sphere[part].py + factor2*sphere[part].Fy);
		sphere[part].pz = factor1*(sphere[part].pz + factor2*sphere[part].Fz);
	}
	//-------------------------------free-drift---------------------------------------
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		dx = sphere[part].px*time_step; 
		dy = sphere[part].py*time_step; 
		dz = sphere[part].pz*time_step; 
		sphere[part].x		+= dx*inv_L;
		sphere[part].x_u 	+= dx*inv_L;													//unfolded
		sphere[part].y  	+= dy*inv_L;
		sphere[part].y_u	+= dy*inv_L;													//unfolded
		sphere[part].z  	+= dz*inv_L;
		sphere[part].z_u 	+= dz*inv_L;
	}
	applypbc1();
	maxD=0.0;
	maxmaxD=0.0;
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		dx=sphere[part].x_u-xold[part];
		dy=sphere[part].y_u-yold[part];
		dz=sphere[part].z_u-zold[part];
        temp = dx*dx + dy*dy +dz*dz;
        if(maxD < temp) maxD = temp;
	}
	maxD*=L*L;
	MPI_Allreduce(&maxD, &maxmaxD, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	int refill=1;
	//if(node==0)printf("%d\t%lf\t%lf\n",node,maxmaxD,skin*skin);
	if(3.0*maxmaxD>skin*skin){
		refill=0;
		ExchangeParticlesInBlocks();
		CommunicateCoordinates(refill);
		update_cell_list();
	}
	else CommunicateCoordinates(refill);
	calculate_forces();

	tempX1blk = 0.0;
    tempX2blk = 0.0;
	K0blk=0.0e0;
	for(int i=0 ; i<NumInBlock[node] ; i++ ){
		int part=InBlock[node][i];
        tempX1blk += sphere[part].px*sphere[part].Fx + sphere[part].py*sphere[part].Fy + sphere[part].pz*sphere[part].Fz;
        tempX2blk += sphere[part].Fx*sphere[part].Fx + sphere[part].Fy*sphere[part].Fy + sphere[part].Fz*sphere[part].Fz;
		K0blk	+= sphere[part].px*sphere[part].px +  sphere[part].py*sphere[part].py + sphere[part].pz*sphere[part].pz ;
	}
    MPI_Reduce(&tempX1blk, &tempX1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tempX1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&tempX2blk, &tempX2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tempX2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&K0blk, &K0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&K0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	alpha = sqrt(tempX2/K0);
    zeta0 = -tempX1/K0;
	beta = exp(-alpha*time_step*0.5);
    gamma = (zeta0 - alpha)/(zeta0 + alpha);
    factor1 = (1.0-gamma)/(beta - gamma/beta);
    factor2 = (1.0+gamma-beta-gamma/beta)/((1.0-gamma)*alpha);
    //-------------------------------half-kick----------------------------------------
	keblk=0.0e0;
	for(int i=0 ; i<NumInBlock[node] ; i++ ){
		int part=InBlock[node][i];
		sphere[part].px = factor1*(sphere[part].px + factor2*sphere[part].Fx);
		sphere[part].py = factor1*(sphere[part].py + factor2*sphere[part].Fy);
		sphere[part].pz = factor1*(sphere[part].pz + factor2*sphere[part].Fz);
		keblk += ( sphere[part].px*sphere[part].px + 
					sphere[part].py*sphere[part].py + 
					sphere[part].pz*sphere[part].pz );
	}
	MPI_Reduce(&keblk, &ke, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ke, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	inst_T=2.0*ke/DOF;
	ke = ke * 0.5;
	return;
}


//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//                      M-Chain Nose-Hoover NVT
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
void InitialiseNoseHoover(){
	tau=.20;
	for(int i=0;i<M;i++){
		q[i]=T*tau*tau;
		eta[i]=0.0e0;	
	}
	q[0]=DOF*T*tau*tau;
	if(node==0){
		int for_ran1 = -(int)(time(0)+node);
	  	srand( node + 100 - for_ran1);           
		for(int i=0;i<M;i++){
			double abc=random_normal();
			p_eta[i]=sqrt(T*q[i])*abc;
		}
	}
	MPI_Bcast(&p_eta,M,MPI_DOUBLE,0,MPI_COMM_WORLD);
return;
}

void AdvanceTimeNVTNoseHoover(){
	Heff=0;	
	double dx,dy,dz,dsx,dsy,dsz,temp;
	maxD=0.0;
	U4Propagator(time_step*0.25,M-1,0);
	U3Propagator(time_step*0.50);
	U4Propagator(time_step*0.25,0,M-1);
	
	U2Propagator();
	U1Propagator();
	
	int refill=1;
	//if(node==0)printf("%d\t%lf\t%lf\n",node,maxmaxD,skin*skin);
	if(3.0*maxmaxD>skin*skin){
		refill=0;
		ExchangeParticlesInBlocks();
		CommunicateCoordinates(refill);
		update_cell_list();
	}
	else CommunicateCoordinates(refill);
	calculate_forces();
	U2Propagator();
	U4Propagator(time_step*0.25,M-1,0);
	U3Propagator(time_step*0.50);
	U4Propagator(time_step*0.25,0,M-1);
	ke    = ke * 0.5;
//	Heff=ke+pe+eta[0]*T*DOF;
//	for(int i=0;i<M;i++) Heff=Heff+p_eta[i]*p_eta[i]/(2.0*q[i]);
//	for(int i=1;i<M;i++) Heff=Heff+T*eta[i];
	inst_T=2.0*ke/DOF;
	return;	
}
void U1Propagator(){
	double dx,dy,dz, temp;
	//-------------------------------free-drift---------------------------------------
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		dx = sphere[part].px*time_step; 
		dy = sphere[part].py*time_step; 
		dz = sphere[part].pz*time_step; 
		sphere[part].x		+= dx*inv_L;
		sphere[part].x_u 	+= dx*inv_L;													//unfolded
		sphere[part].y  	+= dy*inv_L;
		sphere[part].y_u	+= dy*inv_L;													//unfolded
		sphere[part].z  	+= dz*inv_L;
		sphere[part].z_u 	+= dz*inv_L;
	}
	applypbc1();
	maxD=0.0;
	maxmaxD=0.0;
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		dx=sphere[part].x_u-xold[part];
		dy=sphere[part].y_u-yold[part];
		dz=sphere[part].z_u-zold[part];
        temp = dx*dx + dy*dy +dz*dz;
        if(maxD < temp) maxD = temp;
	}
	maxD*=L*L;
	MPI_Allreduce(&maxD, &maxmaxD, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	return;
}
void U2Propagator(){
	//-------------------------------half-kick----------------------------------------
	for(int i=0 ; i<NumInBlock[node] ; i++ ){
		int part=InBlock[node][i];
		sphere[part].px += sphere[part].Fx*half_time_step;
		sphere[part].py += sphere[part].Fy*half_time_step; 
		sphere[part].pz += sphere[part].Fz*half_time_step; 
	}
	return;
}
void U3Propagator(double t){
	for(int i=0 ; i<NumInBlock[node] ; i++ ){
		int part=InBlock[node][i];
		sphere[part].px = sphere[part].px*exp(-t*p_eta[0]/q[0]);
		sphere[part].py = sphere[part].py*exp(-t*p_eta[0]/q[0]);
		sphere[part].pz = sphere[part].pz*exp(-t*p_eta[0]/q[0]);
	}	
	for(int i=0;i<M;i++) eta[i]=eta[i]+t*p_eta[i]/q[i];
	return;
}
void U4Propagator(double t,int jStart,int jEnd){
	keblk=0;
	int dj=1;
	double gj=0.0;	
	if(jStart>jEnd) dj = -1;
	int j=jStart;
	do {
		if(j==0){
			for(int i=0 ; i<NumInBlock[node] ; i++ ){
				int part=InBlock[node][i];
				keblk += sphere[part].px*sphere[part].px + sphere[part].py*sphere[part].py + sphere[part].pz*sphere[part].pz;
			}
    		MPI_Reduce(&keblk, &gj, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    		MPI_Bcast(&gj, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			ke=gj;
			gj -=DOF*T;
		}
		else{
			gj = (p_eta[j-1]*p_eta[j-1])/q[j-1]-T;
		}
		if(j==M-1){
			p_eta[j]=p_eta[j]+t*gj;
		}
		else{
			double x=t*p_eta[j+1]/q[j+1];	
			p_eta[j]=p_eta[j]*exp(-x)+t*gj*(1.0-exp(-x))/x;
		}
	j=j+dj;
	}while( j-dj!=jEnd );
	return;	
}

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
				//NoseHoover-NPT  Integrator
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
void InitialiseNoseHooverNPT(){
	tau=.20;
	tau_baro=2.0;	
	L0=L;	
	w_eps=DOF*T*tau_baro*tau_baro;
	eps=0.0;
	for(int i=0;i<M;i++){
		q[i]=T*tau*tau;
		eta[i]=0.0e0;	
	}
	q[0]=DOF*T*tau*tau;
	if(node==0){
		int for_ran1 = -(int)(time(0)+node);
	  	srand( node + 100 - for_ran1);           
		for(int i=0;i<M;i++){
			double abc=random_normal();
			p_eta[i]=sqrt(T*q[i])*abc;
		}
	}
	for(int i=0;i<M;i++){
		q_baro[i]=T*tau_baro*tau_baro;
		eta_baro[i]=0.0e0;	
	}
	if(node==0){	
		for(int i=0;i<M;i++){
			double abc=random_normal();
			p_eta_baro[i]=sqrt(T*q_baro[i])*abc;
		}
	}
	if(node==0){
		double abc=random_normal();
		p_eps=abc*sqrt(T*w_eps);
	}
	MPI_Bcast(&p_eta,M,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&p_eta_baro,M,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&p_eps,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

return;
}

void AdvanceTimeNPTNoseHoover(){
	Heff=0;	
	double dx,dy,dz,temp;
	maxD=0.0;
	U4PropagatorNPT(time_step*0.25,M-1,0);
	U3PropagatorNPT(time_step*0.50);
	U4PropagatorNPT(time_step*0.25,0,M-1);
	
	U2PrimePropagatorNPT();
	U2PropagatorNPT();
	U1PropagatorNPT();
	
	maxD=sqrt(maxD);
	applypbc1();
	neib_list_counter ++;
	MPI_Reduce(&maxD, &maxmaxD, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(node==0)
		globalMaxD+=maxmaxD;
	MPI_Bcast(&globalMaxD,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	int refill=1;
	CommunicateCoordinates(refill);
	if(2.0*globalMaxD>skin){
		refill=0;
		ExchangeParticlesInBlocks();
		CommunicateCoordinates(refill);
		update_cell_list();
	}
	calculate_forces();
	
	U2PropagatorNPT();
	U2PrimePropagatorNPT();

	U4PropagatorNPT(time_step*0.25,M-1,0);
	U3PropagatorNPT(time_step*0.50);
	U4PropagatorNPT(time_step*0.25,0,M-1);
	ke    = ke * 0.5;
	//Heff=ke+pe+eta[0]*T*DOF+PV+p_eps*p_eps/(2.0*w_eps);  //Heff is the conserved quantity
	//for(int i=0;i<M;i++) Heff=Heff+p_eta[i]*p_eta[i]/(2.0*q[i])+p_eta_baro[i]*p_eta_baro[i]/(2.0*q_baro[i])+T*eta_baro[i];
	//for(int i=1;i<M;i++) Heff=Heff+T*eta[i];
	inst_T=2.0*ke/DOF;
	return;	
}
void U1PropagatorNPT(double t){
	double dx,dy,dz, temp;
	//-------------------------------free-drift---------------------------------------
	double x=time_step*p_eps/w_eps;
	double factor=(1.0-exp(-x))/x;
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		dx = sphere[part].px*time_step; 
		dy = sphere[part].py*time_step; 
		dz = sphere[part].pz*time_step; 
		sphere[part].x		+= dx*inv_L;
		sphere[part].x_u 	+= dx*inv_L;													//unfolded
		sphere[part].y  	+= dy*inv_L;
		sphere[part].y_u	+= dy*inv_L;													//unfolded
		sphere[part].z  	+= dz*inv_L;
		sphere[part].z_u 	+= dz*inv_L;
		temp = dx*dx + dy*dy + dz*dz;
		if(maxD < temp) maxD = temp;
	}
	maxD=sqrt(maxD);
	applypbc1();
	eps +=x ;
	L=L0*exp(eps);
	UpdateBoxForNPT();
	return;
}
void U2PropagatorNPT(double t){
	//-------------------------------half-kick----------------------------------------
	double alpha=1.0+3.0/DOF;
	double x=half_time_step*alpha*p_eps/w_eps;
	double factor1=exp(-x);
	double factor2=(1.0-factor1)/x;
	
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		sphere[part].px = sphere[part].px*factor1 + factor2*sphere[part].Fx*half_time_step;
		sphere[part].py = sphere[part].py*factor1 + factor2*sphere[part].Fy*half_time_step;
		sphere[part].pz = sphere[part].pz*factor1 + factor2*sphere[part].Fz*half_time_step;
	}
	return;
	
}
void U2PrimePropagatorNPT(double t){
	double alpha=1.0+3.0/DOF;
	double sum=0.0e0;
	double sum1=0.0e0;
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		sum1 += ( sphere[part].px*sphere[part].px + sphere[part].py*sphere[part].py + sphere[part].pz*sphere[part].pz )  ;
	}
    MPI_Reduce(&virialblk, &virial, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&virial, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&sum1, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	PV = alpha*(sum)/(3.00) + virial;
	p_eps += 3.00*(PV-P*L*L*L)*half_time_step;
	return;
}

void U3PropagatorNPT(double t){
	for(int i=0 ; i<NumInBlock[node] ; i++ ){
		int part=InBlock[node][i];
		sphere[part].px = sphere[part].px*exp(-t*p_eta[0]/q[0]);
		sphere[part].py = sphere[part].py*exp(-t*p_eta[0]/q[0]);
		sphere[part].pz = sphere[part].pz*exp(-t*p_eta[0]/q[0]);
	}	
	for(int i=0;i<M;i++) {
		eta[i]=eta[i]+t*p_eta[i]/q[i];
		eta_baro[i] += t*p_eta_baro[i]/q_baro[i];
	}
	p_eps = p_eps*exp(-t*p_eta_baro[0]/q_baro[0]);
	return;
}
void U4PropagatorNPT(double t,int jStart,int jEnd){
	keblk=0.0;
	int dj=1;
	double gj=0.0;	
	if(jStart>jEnd) dj = -1;
	int j=jStart;
	do {
		if(j==0){
			for(int i=0 ; i<NumInBlock[node] ; i++ ){
				int part=InBlock[node][i];
				keblk += sphere[part].px*sphere[part].px + sphere[part].py*sphere[part].py + sphere[part].pz*sphere[part].pz;
			}
    		MPI_Reduce(&keblk, &gj, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    		MPI_Bcast(&gj, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			ke=gj;
			gj -=DOF*T;
		}
		else{
			gj = (p_eta[j-1]*p_eta[j-1])/q[j-1]-T;
		}
		if(j==M-1){
			p_eta[j]=p_eta[j]+t*gj;
		}
		else{
			double x=t*p_eta[j+1]/q[j+1];	
			p_eta[j]=p_eta[j]*exp(-x)+t*gj*(1.0-exp(-x))/x;
		}
	j=j+dj;
	}while( j-dj!=jEnd );
	j=jStart;
	do {
		if(j==0){
			gj = p_eps*p_eps/w_eps-T;
		}
		else{
			gj = (p_eta_baro[j-1]*p_eta_baro[j-1])/q_baro[j-1]-T;
		}
		if(j==M-1){
			p_eta_baro[j]=p_eta_baro[j]+t*gj;
		}
		else{
			double x=t*p_eta_baro[j+1]/q_baro[j+1];	
			p_eta_baro[j]=p_eta_baro[j]*exp(-x) + t*gj*(1.0-exp(-x))/x;
		}
	j=j+dj;
	}while( j-dj!=jEnd );
	return;	
}

// *********************Update the Simulation box for NPT task************************
void UpdateBoxForNPT(){
	DENSITY = (double)(N)/(L*L*L);
	inv_L = 1.00/L;
    Lppx=L/(double)npx;                                 //Length Allocated per processor
    Lppy=L/(double)npy;
    Lppz=L/(double)npz;
    Ncellx=(int)Lppx/Cell_cut;                                  //# of cells along x in a block
    Ncelly=(int)Lppy/Cell_cut;                                  //# of cells along y in a block
    Ncellz=(int)Lppz/Cell_cut;                                  //# of cells along z in a block
    CellSizex=Lppx/(double)Ncellx;                              //cell size along x in system of unit length
    CellSizey=Lppy/(double)Ncelly;                              //cell size along y in system of unit length
    CellSizez=Lppz/(double)Ncellx;                              //cell size along z in system of unit length
    InvCellSizex=L/CellSizex;
    InvCellSizey=L/CellSizey;                                   //In Reduced units  
    InvCellSizez=L/CellSizez;

    Ncellx=Ncellx*npx;                                          //Tot # of cells along x in whole domain
    Ncelly=Ncelly*npy;                                          //Tot # of cells along y in whole domain
    Ncellz=Ncellz*npz;                                          //Tot # of cells along z in whole domain

    NCell=Ncellx*Ncelly*Ncellz;                                 //Tot # of cells in whole domain
    Slabdx=CellSizex/L;                                         //Tranfer slab which have the dimensions of a unit cell;    
    Slabdy=CellSizey/L;
    Slabdz=CellSizez/L;
  	maps();
	CommunicateCoordinates(1);
return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// *********************Brownian Dynamics with Predictor-Corrector Alg with active forcing ************************
void advance_time_ActiveBrown_PC(){                                  
	double dx,dy,dz,temp;
	double fx,fy,fz;
//	double SqrtVariance=sqrt(2.0*T*T/(time_step*D0));
//	double InvD0=1.0/D0;
//	double D0ByT=D0/T;
//	double TByD0=T/D0;
//	double expFactor = exp(-time_step*TByD0);
	double FxCOMblk = 0.0;
    double FyCOMblk = 0.0;
    double FzCOMblk = 0.0;
	double FxCOM = 0.0;
    double FyCOM = 0.0;
    double FzCOM = 0.0;
	double Netkx,Netky,Netkz;
//	printf("%lf\t%lf\t%lf\n",Netkx,Netky,Netkz);	
	//-------------------------------Getting Random forces----------------------------------------
    for(int i=0; i<NumInBlock[node]; i++){
        int part=InBlock[node][i];
		fxRandom[part]=SqrtVariance*normal(0.0,1.0);
		fyRandom[part]=SqrtVariance*normal(0.0,1.0);
		fzRandom[part]=SqrtVariance*normal(0.0,1.0);
		FxCOMblk+=fxRandom[part];
		FyCOMblk+=fyRandom[part];
		FzCOMblk+=fzRandom[part];
	}
	MPI_Allreduce(&FxCOMblk, &FxCOM, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&FyCOMblk, &FyCOM, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&FzCOMblk, &FzCOM, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	FxCOM /=(double)N;
	FyCOM /=(double)N;
	FzCOM /=(double)N;
    for(int i=0; i<NumInBlock[node]; i++){
        int part=InBlock[node][i];
		fxRandom[part]-=FxCOM;
		fyRandom[part]-=FyCOM;
		fzRandom[part]-=FzCOM;
		rxold[part]=sphere[part].x;
		ryold[part]=sphere[part].y;
		rzold[part]=sphere[part].z;
		fxold[part]=sphere[part].Fx;
		fyold[part]=sphere[part].Fy;
		fzold[part]=sphere[part].Fz;
	}
	maxD=0.0;
	ke=0.0e0;
	int ActPartNum=0;
	//-------------------------------Predictor Step---------------------------------------	
    for(int i=0; i<NumInBlock[node]; i++){
        int part=InBlock[node][i];
		double fx=fxRandom[part]+sphere[part].Fx;
		double fy=fyRandom[part]+sphere[part].Fy;
		double fz=fzRandom[part]+sphere[part].Fz;
		dx = (1.0-expFactor)*sphere[part].px + D0ByT*(TByD0*time_step - (1.0-expFactor))*fx; 
		dy = (1.0-expFactor)*sphere[part].py + D0ByT*(TByD0*time_step - (1.0-expFactor))*fy; 
		dz = (1.0-expFactor)*sphere[part].pz + D0ByT*(TByD0*time_step - (1.0-expFactor))*fz; 
		dx = D0ByT*dx;	
		dy = D0ByT*dy;
		dz = D0ByT*dz;	
		sphere[part].x		+= dx*inv_L;
		sphere[part].y  	+= dy*inv_L;
		sphere[part].z  	+= dz*inv_L;
		temp = dx*dx + dy*dy + dz*dz;
		if(maxD < temp) maxD = temp;
	}
	applypbc1();
    CommunicateCoordinates(1);
	calculate_forces();
	//-------------------------------corrector step---------------------------------------	
    for(int i=0; i<NumInBlock[node]; i++){
        int part=InBlock[node][i];
		sphere[part].x=rxold[part];
		sphere[part].y=ryold[part];
		sphere[part].z=rzold[part];
		double fx=fxRandom[part]+0.5*(sphere[part].Fx+fxold[part]);
		double fy=fyRandom[part]+0.5*(sphere[part].Fy+fyold[part]);
		double fz=fzRandom[part]+0.5*(sphere[part].Fz+fzold[part]);
		dx = (1.0-expFactor)*sphere[part].px + D0ByT*(TByD0*time_step - (1.0-expFactor))*fx; 
		dy = (1.0-expFactor)*sphere[part].py + D0ByT*(TByD0*time_step - (1.0-expFactor))*fy; 
		dz = (1.0-expFactor)*sphere[part].pz + D0ByT*(TByD0*time_step - (1.0-expFactor))*fz; 
		dx = D0ByT*dx;	
		dy = D0ByT*dy;
		dz = D0ByT*dz;	
		sphere[part].x		+= dx*inv_L;
		sphere[part].x_u 	+= dx*inv_L;													//unfolded
		sphere[part].y  	+= dy*inv_L;
		sphere[part].y_u	+= dy*inv_L;													//unfolded
		sphere[part].z  	+= dz*inv_L;
		sphere[part].z_u 	+= dz*inv_L;
		sphere[part].px = expFactor*sphere[part].px + D0ByT*(1.0-expFactor)*fx;
		sphere[part].py = expFactor*sphere[part].py + D0ByT*(1.0-expFactor)*fy;
		sphere[part].pz = expFactor*sphere[part].pz + D0ByT*(1.0-expFactor)*fz;
//		ke += ( sphere[part].px*sphere[part].px + sphere[part].py*sphere[part].py + sphere[part].pz*sphere[part].pz )  ;
	}

	applypbc1();
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		dx=sphere[part].x_u-xold[part];
		dy=sphere[part].y_u-yold[part];
		dz=sphere[part].z_u-zold[part];
        temp = dx*dx + dy*dy +dz*dz;
        if(maxD < temp) maxD = temp;
	}

	maxD*=L*L;
	MPI_Allreduce(&maxD, &maxmaxD, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	int refill=1;

	if(2.50*maxmaxD>skin*skin){
        refill=0;
        ExchangeParticlesInBlocks();
        CommunicateCoordinates(refill);
        update_cell_list();
    }
	else
    	CommunicateCoordinates(refill);
	calculate_forces();
	return;
}
