/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:init_MPI.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-05
*   Discription:
*
================================================================*/
#include "header.h"

//#include "/public/software/mpich/include/mpi.h"

void finalize_MPI( MPI_Comm * comm_cart  )
{
	MPI_Comm_free( comm_cart );
	MPI_Finalize( );

}

void init_MPI( int *argc, char *** argv, PARAMS params, MPI_Comm * comm_cart,
MPI_COORD * thisMPICoord, MPI_NEIGHBOR * mpiNeighbor )
{
	int PX = params.PX;
	int PY = params.PY;
	int PZ = params.PZ;


	int thisRank, thisMPICoordXYZ[3];
	
	int nDim = 3;
	int mpiDims[3] = { PX, PY, PZ };
	int periods[3] = { 0, 0, 0 };
	int reorder = 0;
	
	MPI_Init( argc, argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &thisRank );
	//MPI_Comm_size( MPI_COMM_WORLD, &nProcs );
	
	MPI_Cart_create( MPI_COMM_WORLD, nDim, mpiDims, periods, reorder, comm_cart );


	MPI_Cart_shift( *comm_cart, 0, 1, &( mpiNeighbor->X1 ), & ( mpiNeighbor->X2 ) );
	MPI_Cart_shift( *comm_cart, 1, 1, &( mpiNeighbor->Y1 ), & ( mpiNeighbor->Y2 ) );
	MPI_Cart_shift( *comm_cart, 2, 1, &( mpiNeighbor->Z1 ), & ( mpiNeighbor->Z2 ) );	

	//if ( mpiNeighbor->X1 < 0 ) mpiNeighbor->X1 = MPI_PROC_NULL;
	//if ( mpiNeighbor->Y1 < 0 ) mpiNeighbor->Y1 = MPI_PROC_NULL;
	//if ( mpiNeighbor->Z1 < 0 ) mpiNeighbor->Z1 = MPI_PROC_NULL;
	//if ( mpiNeighbor->X2 < 0 ) mpiNeighbor->X1 = MPI_PROC_NULL;
	//if ( mpiNeighbor->Y2 < 0 ) mpiNeighbor->Y1 = MPI_PROC_NULL;
	//if ( mpiNeighbor->Z2 < 0 ) mpiNeighbor->Z1 = MPI_PROC_NULL;


	MPI_Cart_coords( *comm_cart, thisRank, 3, thisMPICoordXYZ );

	
	thisMPICoord->X = thisMPICoordXYZ[0]; 
	thisMPICoord->Y = thisMPICoordXYZ[1]; 
	thisMPICoord->Z = thisMPICoordXYZ[2]; 



}

