/**************************************************
**** All rights reserved.
**** Author: Wenqiang Wang
**** Department of Earth and Space Sciences, Sustech
**** E-mail: 11849528@mail.sustech.edu.cn
**************************************************/
#include "header.h"


int main( int argc, char ** argv )
{

	PARAMS params;
	GRID grid = { 0 };
		
	MPI_COORD thisMPICoord = { 0 };
	MPI_NEIGHBOR mpiNeighbor = { 0 };
	MPI_Comm comm_cart;


	getParams( &params );
	init_MPI( &argc, &argv, params, &comm_cart, &thisMPICoord, &mpiNeighbor );
	init_grid( params, &grid, thisMPICoord );
	
	createDir( params );

#ifdef GPU_CUDA
	init_gpu( grid.PX, grid.PY, grid.PZ );
#endif


	MPI_Barrier( comm_cart );
	
	run( comm_cart, thisMPICoord, mpiNeighbor, grid, params );
	
	MPI_Barrier( comm_cart );
	MPI_Comm_free( &comm_cart );
	MPI_Finalize( );
	return 0;
}

