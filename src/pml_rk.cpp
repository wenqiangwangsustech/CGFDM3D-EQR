/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:pml_rk.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-05
*   Discription:
*
================================================================*/

#include "header.h"
void wave_rk0( float * h_W, float * W, float * t_W, float * m_W, long long num, float DT );
void wave_rk1( float * h_W, float * W, float * t_W, float * m_W, long long num, float DT );
void wave_rk2( float * h_W, float * W, float * t_W, float * m_W, long long num, float DT );
void wave_rk3( float * h_W, float * W, float * t_W, float * m_W, long long num, float DT );


void pmlRk( GRID grid, MPI_BORDER border, int irk, AUX4 Aux4_1, AUX4 Aux4_2, float DT )
{
	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;
	
	int nPML = grid.nPML;

	long long num = 0;
	
	WAVE_RK_FUNC pml_rk[4] = { wave_rk0, wave_rk1, wave_rk2, wave_rk3 };
	long long numx = nPML * ny * nz * WAVESIZE;
	long long numy = nPML * nx * nz * WAVESIZE;
	long long numz = nPML * nx * ny * WAVESIZE;


#ifdef GPU_CUDA
                        
	dim3 thread( 512, 1, 1 );
	dim3 blockX;

	blockX.x = ( numx + thread.x - 1 ) / thread.x;
	blockX.y = 1;
	blockX.z = 1;

	dim3 blockY;
	blockY.x = ( numy + thread.x - 1 ) / thread.x;
	blockY.y = 1;
	blockY.z = 1;

	dim3 blockZ;
	blockZ.x = ( numz + thread.x - 1 ) / thread.x;
	blockZ.y = 1;
	blockZ.z = 1;

	if ( border.isx1 )	pml_rk[irk]<<< blockX, thread >>>( Aux4_1.h_Aux_x.Vx, Aux4_1.Aux_x.Vx, Aux4_1.t_Aux_x.Vx, Aux4_1.m_Aux_x.Vx, numx, DT );
	if ( border.isy1 )  pml_rk[irk]<<< blockY, thread >>>( Aux4_1.h_Aux_y.Vx, Aux4_1.Aux_y.Vx, Aux4_1.t_Aux_y.Vx, Aux4_1.m_Aux_y.Vx, numy, DT );
	if ( border.isz1 )  pml_rk[irk]<<< blockZ, thread >>>( Aux4_1.h_Aux_z.Vx, Aux4_1.Aux_z.Vx, Aux4_1.t_Aux_z.Vx, Aux4_1.m_Aux_z.Vx, numz, DT );
    //                                                                                                                              
	if ( border.isx2 )	pml_rk[irk]<<< blockX, thread >>>( Aux4_2.h_Aux_x.Vx, Aux4_2.Aux_x.Vx, Aux4_2.t_Aux_x.Vx, Aux4_2.m_Aux_x.Vx, numx, DT );
	if ( border.isy2 )  pml_rk[irk]<<< blockY, thread >>>( Aux4_2.h_Aux_y.Vx, Aux4_2.Aux_y.Vx, Aux4_2.t_Aux_y.Vx, Aux4_2.m_Aux_y.Vx, numy, DT );
	//CHECK( cudaDeviceSynchronize( ) );
                                                                                                                                        
#ifndef FREE_SURFACE
	if ( border.isz2 )  pml_rk[irk]<<< blockZ, thread >>>( Aux4_2.h_Aux_z.Vx, Aux4_2.Aux_z.Vx, Aux4_2.t_Aux_z.Vx, Aux4_2.m_Aux_z.Vx, numz, DT );
#endif                                                                                                               

                                                                                                                                   
#else                                                                                                                              
                                                                                                                                   
	
                                                                                                                                   
                                                                                                                                   
	if ( border.isx1 )	pml_rk[irk]( Aux4_1.h_Aux_x.Vx, Aux4_1.Aux_x.Vx, Aux4_1.t_Aux_x.Vx, Aux4_1.m_Aux_x.Vx, numx, DT );
	if ( border.isy1 )  pml_rk[irk]( Aux4_1.h_Aux_y.Vx, Aux4_1.Aux_y.Vx, Aux4_1.t_Aux_y.Vx, Aux4_1.m_Aux_y.Vx, numy, DT );
	if ( border.isz1 )  pml_rk[irk]( Aux4_1.h_Aux_z.Vx, Aux4_1.Aux_z.Vx, Aux4_1.t_Aux_z.Vx, Aux4_1.m_Aux_z.Vx, numz, DT );
    //                                                                                                        
	if ( border.isx2 )	pml_rk[irk]( Aux4_2.h_Aux_x.Vx, Aux4_2.Aux_x.Vx, Aux4_2.t_Aux_x.Vx, Aux4_2.m_Aux_x.Vx, numx, DT );
	if ( border.isy2 )  pml_rk[irk]( Aux4_2.h_Aux_y.Vx, Aux4_2.Aux_y.Vx, Aux4_2.t_Aux_y.Vx, Aux4_2.m_Aux_y.Vx, numy, DT );
#ifndef FREE_SURFACE
	if ( border.isz2 )  pml_rk[irk]( Aux4_2.h_Aux_z.Vx, Aux4_2.Aux_z.Vx, Aux4_2.t_Aux_z.Vx, Aux4_2.m_Aux_z.Vx, numz, DT );
#endif                                                                                                               


#endif

}

