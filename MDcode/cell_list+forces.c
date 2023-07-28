		/*********************************************************************************************************
		**	This bit of code generates cell list and calculate potential energy, virial, and forces over		**
		**	each particles that includes forces over each sphere, each sphere of rod, and COM force of rod.		**
		**										$$$ROUTINES INCLUDED$$$											**
		**	getcell_I		==> gives back the cell id with input of x,y no of that cell						**
		**	maps			==>	stores the cellno of neighbour cells for in MAPS[neib_no][cellno]. To be called **
		**						onces in code if not NPT otherwise at each step.								**
		**	update_cell_list==>	generates cell list for all spheres including the spheres of rods. To be called **
		**						when tot. displacement of particles exceeds skin of cell list.					**
		**	calculate_forces==>	calculate forces over each particles that includes forces over each sphere,		**
		**						each sphere of rod, COM force of rod and torque on rod.							**
		**																										**
		*********************************************************************************************************/

#include "global.h"
#include "subroutine.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double 	sigma[3][3] = {{1.00,0.80,1.0},{0.80,0.88,1.0},{1.0,1.0,1.0}};
double	epsilon[3][3] = {{1.0,1.5,1.0},{1.5,0.50,1.0},{1.0,1.0,.50}};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 			Get Cell index for location Ix Iy Iz
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int		getcell_I(int Ix ,int Iy , int Iz){
	int cellno = ((Ix+Ncellx*3 )%Ncellx) + ((Iy+Ncelly*3)%Ncelly)*Ncellx
										+((Iz+Ncellz*3)%Ncellz)*Ncellx*Ncelly;
	return cellno;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 			Routine to lay the map of neighbouring cells
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void maps(){
	for(int Ix=0 ;Ix<Ncellx;Ix++){
		for(int Iy=0 ; Iy<Ncelly ; Iy++){
			for(int Iz=0 ; Iz<Ncellz ; Iz++){	
				int IMAP=getcell_I(Ix,Iy,Iz);
				MAP[IMAP][0] = getcell_I(Ix+1, Iy, Iz);				//east
				MAP[IMAP][1] = getcell_I(Ix+1, Iy+1, Iz);			//north-east			
				MAP[IMAP][2] = getcell_I(Ix, Iy+1, Iz);				//north
				MAP[IMAP][3] = getcell_I(Ix-1 , Iy+1, Iz);			//north-west
				MAP[IMAP][4] = getcell_I(Ix+1, Iy, Iz-1);			//east-bottom
				MAP[IMAP][5] = getcell_I(Ix+1, Iy+1, Iz-1);			//north-east-bottom			
				MAP[IMAP][6] = getcell_I(Ix, Iy+1, Iz-1);			//north-bottom
				MAP[IMAP][7] = getcell_I(Ix-1 , Iy+1, Iz-1);		//north-west-bottom
				MAP[IMAP][8] = getcell_I(Ix+1, Iy, Iz+1);			//east-top
				MAP[IMAP][9] = getcell_I(Ix+1, Iy+1, Iz+1);			//north-east-top			
				MAP[IMAP][10] = getcell_I(Ix, Iy+1, Iz+1);			//north-top
				MAP[IMAP][11] = getcell_I(Ix-1 , Iy+1, Iz+1);		//north-west-top
				MAP[IMAP][12] = getcell_I(Ix, Iy, Iz+1);			//top
			}
		}
	}
	return;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//							Lets do the listing
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void update_cell_list(){ 
	int count,i,j,k;
	int part1,part2;
	int cell_index,cell_index2;
	bool is_neib;
	double dx,dy,dz,dr2;
	maxD=0.0;globalMaxD=0.0;										//reseting the MaxD
	for(int i=0;i<NCell;i++)										//Here NCell=#of cell in whole box and not in node
		InCellCount[i]=0;											//Lets initialize
	NumInExtdBlk=0;													//Extd. Block is the block with extra one layer of cell
	for(int i=0; i<N ; i++)
		NeibListCount[i]=0;

	for(int i=0;i<N;i++){
		if(InMe[i]){												//InMe is a tag for extd block
        	cell_index = floor(sphere[i].x*InvCellSizex) + floor(sphere[i].y*InvCellSizey)*Ncellx +
                        floor(sphere[i].z*InvCellSizez)*Ncellx*Ncelly;
			InCell[cell_index][InCellCount[cell_index]]=i;
			InCellCount[cell_index]++;
			InExtdBlk[NumInExtdBlk]=i;
			NumInExtdBlk++;
		}
	}
   //***********************Generating cell-list for all particles in Extd Blocks **************************//
    for(i=0; i<NumCellInBlk[node]; i++){
        cell_index = CellInBlk[node][i];
        count=InCellCount[cell_index];
        for(j=0; j<count-1; j++){               //in its own cell
            part1 = InCell[cell_index][j];
            for(k=j+1;k<count;k++){
                part2 = InCell[cell_index][k];
                dx = sphere[part1].x-sphere[part2].x;
                dy = sphere[part1].y-sphere[part2].y;
                dz = sphere[part1].z-sphere[part2].z;
                if( dx < -HALF_LENGTH) dx += LENGTH;
                if( dx > HALF_LENGTH)  dx -= LENGTH;
                if( dy < -HALF_LENGTH) dy += LENGTH;
                if( dy > HALF_LENGTH)  dy -= LENGTH;
                if( dz < -HALF_LENGTH) dz += LENGTH;
                if( dz > HALF_LENGTH)  dz -= LENGTH;
                dr2 = (dx*dx + dy*dy + dz*dz)*L*L;
                if( dr2 < Cell_cut_sq ){
					NeibList[part1][NeibListCount[part1]] = part2;
					NeibListCount[part1]++;
                }
            }
        }
    } 
	
	for(cell_index=0;cell_index<NCell;cell_index++){
		if(InCellCount[cell_index]>0){
			for(int cell_neibs = 0; cell_neibs < half_neibs; cell_neibs++){
				cell_index2 = MAP[cell_index][cell_neibs];
//				if(Cell[cell_index2].count==0){
//					fprintf(stderr,"Some thing is wrong with slabs \n # in cell %d=0",cell_index2);
//				}
				if(InCellCount[cell_index2]>0){
					for(i=0; i<InCellCount[cell_index]; i++){	
						for(j=0; j<InCellCount[cell_index2]; j++){			//in neighbouring cells
							part1 = InCell[cell_index][i];
							part2 = InCell[cell_index2][j];
							if( Block[part1]==node || Block[part2]==node ){		//atleast one of them should be in me 
								dx = sphere[part1].x-sphere[part2].x;
								dy = sphere[part1].y-sphere[part2].y;
								dz = sphere[part1].z-sphere[part2].z;
								if( dx < -HALF_LENGTH) dx += LENGTH;
								if( dx > HALF_LENGTH)  dx -= LENGTH;
								if( dy < -HALF_LENGTH) dy += LENGTH;
								if( dy > HALF_LENGTH)  dy -= LENGTH;
								if( dz < -HALF_LENGTH) dz += LENGTH;
								if( dz > HALF_LENGTH)  dz -= LENGTH;
								dr2 = (dx*dx + dy*dy + dz*dz)*L*L;
								if( dr2 < Cell_cut_sq ){
									NeibList[part1][NeibListCount[part1]] = part2;
									NeibListCount[part1]++;
								} 
							}
						}
					}
				}
			}				 
		}
	}
	for(int i=0; i<NumInBlock[node]; i++){
		int part=InBlock[node][i];
		xold[part]=sphere[part].x_u;
		yold[part]=sphere[part].y_u;
		zold[part]=sphere[part].z_u;
	}
return;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
	// 				Lets get the forces 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
void calculate_forces(){
	peblk = 0.0 ; virialblk = 0.0; vir_actblk=0.0;												//initialising..........
	int ii,j,tot_neib,A,B,i;
	double factor,dr2,dr,eps,sig,grad;
	double r_over_sig_2 , sig_over_r_2 , sig_over_r_6;
	double dx,dy,dz;
	int part1;
	for(i=0 ; i<NumInExtdBlk ; i++) {
		part1=InExtdBlk[i];
		sphere[part1].Fx=0.0e0 ; sphere[part1].Fy=0.0e0; sphere[part1].Fz=0.0e0;		//initialising...........
	}
	
	for(ii=0; ii<NumInExtdBlk; ii++){
		i=InExtdBlk[ii];
		tot_neib=NeibListCount[i];
		A = sphere[i].type;													  
		for(int k=0; k<tot_neib; k++ ){								  	 	//loop for all neighbours of 'i'
			factor=1.00;
			j= NeibList[i][k];
			if(Block[i]!=node || Block[j]!=node) factor=0.50;
			B = sphere[j].type;
			dx = sphere[i].x - sphere[j].x;
			dy = sphere[i].y - sphere[j].y;
			dz = sphere[i].z - sphere[j].z;
			if( dx < -HALF_LENGTH) dx += LENGTH;
			if( dx > HALF_LENGTH)  dx -= LENGTH;					   		//PBC stuff
			if( dy < -HALF_LENGTH) dy += LENGTH;
			if( dy > HALF_LENGTH)  dy -= LENGTH;
			if( dz < -HALF_LENGTH) dz += LENGTH;
			if( dz > HALF_LENGTH)  dz -= LENGTH;
			dr2 = (dx*dx + dy*dy + dz*dz)*L*L;
			eps = epsilon[A][B];
			sig = sigma[A][B];
			r_over_sig_2 = dr2/(sig*sig);

			
			if(r_over_sig_2 < Lcutoffsqr ){
				sig_over_r_2 = 1.0/r_over_sig_2;
				sig_over_r_6 = sig_over_r_2 * sig_over_r_2 * sig_over_r_2;
				//pe += eps * ( 4.0*(sig_over_r_6*(sig_over_r_6-1.0)) + C0 + (C2+C4*r_over_sig_2)*r_over_sig_2 ); 
				//grad = 2.0 * eps *( 24.0*sig_over_r_6*(sig_over_r_6-0.5)/dr2 - (C2 + 2.0*C4*r_over_sig_2)/(sig*sig) );
				peblk += factor*4.0*eps * ( (sig_over_r_6*(sig_over_r_6-1.0)) + C0 + C2*r_over_sig_2 ); 
				grad = 8.0 * eps *( 6.0*sig_over_r_6*(sig_over_r_6-0.5)/dr2 - C2/(sig*sig) );
				sphere[i].Fx += grad * dx * L;
				sphere[i].Fy += grad * dy * L;
				sphere[i].Fz += grad * dz * L;
				sphere[j].Fx -= grad * dx * L;
				sphere[j].Fy -= grad * dy * L;
				sphere[j].Fz -= grad * dz * L;	
				virialblk += factor*grad*dr2; 
			}

		}
	}
	for(int j=0;j<NumInBlock[node];j++){
		int i=InBlock[node][j];
		if(sphere[i].actTag!=-1){
			int ActIndx=sphere[i].actTag;
		    sphere[i].Fx += f0*(double)kx[ActIndx];
			sphere[i].Fy += f0*(double)ky[ActIndx];
			sphere[i].Fz += f0*(double)kz[ActIndx];
			vir_actblk +=  f0*((double)kx[ActIndx]*sphere[i].x_u 
							+(double)ky[ActIndx]*sphere[i].y_u 
						+(double)kz[ActIndx]*sphere[i].z_u);
		}
	}

	virialblk=virialblk/(3.0e0);
    //need to add the active part of the force
//	MPI_Reduce(&virialblk, &virial, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//	MPI_Reduce(&peblk, &pe, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//	MPI_Reduce(&vir_actblk, &vir_act, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	//fprintf(stderr, "Forces are now calculated\n" );
	return;
}  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool len_check(SPHERE s1,SPHERE s2){
	bool is_neib=false;
	double dx = s1.x-s2.x;
	double dy = s1.y-s2.y;
	double dz = s1.z-s2.z;
	if( dx < -HALF_LENGTH) dx += LENGTH;
	if( dx > HALF_LENGTH)  dx -= LENGTH;
	if( dy < -HALF_LENGTH) dy += LENGTH;
	if( dy > HALF_LENGTH)  dy -= LENGTH;
	if( dz < -HALF_LENGTH) dz += LENGTH;
	if( dz > HALF_LENGTH)  dz -= LENGTH;
	double dr2 = (dx*dx + dy*dy + dz*dz)*L*L;
	double sig = sigma[s1.type][s2.type];
	if( dr2/(sig*sig) < Cell_cut_sq ) is_neib = true;
	return is_neib;	
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
