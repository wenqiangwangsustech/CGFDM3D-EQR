#include "header.h"

void run( MPI_Comm comm_cart, MPI_COORD thisMPICoord, MPI_NEIGHBOR mpiNeighbor, GRID grid, PARAMS params )
{
	int thisRank;
	MPI_Barrier( MPI_COMM_WORLD );
	MPI_Comm_rank( MPI_COMM_WORLD, &thisRank );
	if ( thisRank == 0 )
		printInfo( grid );
	MPI_Barrier( MPI_COMM_WORLD );

	SEND_RECV_DATA sr;
	allocSendRecv( grid, mpiNeighbor, &sr );


	SLICE_DATA sliceData, sliceDataCpu;
	SLICE slice = { 0 };
	locateSlice( params, grid, &slice );
	allocSliceData( grid, slice, &sliceData );
#ifdef GPU_CUDA
	allocSliceData_cpu( grid, slice, &sliceDataCpu );
#endif
	char coordXName[256], coordYName[256], coordZName[256];

	
	//construct coordinate
	COORD coord;
	allocCoord( grid, &coord );
	constructCoord( grid, coord, params );

	{
		sprintf( coordXName, "%s/coordX", params.OUT );
		sprintf( coordYName, "%s/coordY", params.OUT );
		sprintf( coordZName, "%s/coordZ", params.OUT );

		data2D_output_bin( grid, slice, thisMPICoord, coord.x, sliceData, sliceDataCpu, coordXName );
		data2D_output_bin( grid, slice, thisMPICoord, coord.y, sliceData, sliceDataCpu, coordYName );
		data2D_output_bin( grid, slice, thisMPICoord, coord.z, sliceData, sliceDataCpu, coordZName );

		memset( coordXName, 0, 256 );
		memset( coordYName, 0, 256 );
		memset( coordZName, 0, 256 );
	}
	
	//construct medium
	MEDIUM medium;
	allocMedium( grid, &medium );
	constructMedium( grid, medium , params );

	//calculate CFL condition
	calc_CFL( grid, coord, medium, params );

	//translate Vs Vp rho to lambda mu and bouyancy
	vs_vp_rho2lam_mu_bou( medium, grid );

	//solve contravariant matrix and release coordinate memory
	CONTRAVARIANT contravariant;
	float * Jac;
	allocContravariant( grid, &contravariant );
	allocJac( grid, &Jac );

	Mat3x3 _rDZ_DX, _rDZ_DY;
#ifdef FREE_SURFACE
	allocMat3x3( grid, &_rDZ_DX, &_rDZ_DY );
	solveContravariantJac( comm_cart, mpiNeighbor, grid, sr, contravariant, coord, Jac, medium, _rDZ_DX, _rDZ_DY );
#else 
	solveContravariantJac( comm_cart, mpiNeighbor, grid, sr, contravariant, coord, Jac );
#endif
	freeCoord( coord );

		
	MPI_Barrier( MPI_COMM_WORLD );

	WAVE h_W, W, t_W, m_W;
	allocWave( grid, &h_W, &W, &t_W, &m_W );

	propagate( comm_cart, thisMPICoord, mpiNeighbor,
			   grid, params, sr,
			   h_W, W, t_W, m_W,
			   _rDZ_DX, _rDZ_DY,
			   contravariant, Jac, medium,
			   slice, sliceData, sliceDataCpu );


	freeWave( h_W, W, t_W, m_W );

	freeSendRecv( mpiNeighbor, sr );

	freeSliceData( grid, slice, sliceData );
#ifdef GPU_CUDA
	freeSliceData_cpu( grid, slice, sliceDataCpu );
#endif




	//release contravariant Jac medium memory
#ifdef FREE_SURFACE
	freeMat3x3( _rDZ_DX, _rDZ_DY );
#endif
	freeContravariant( contravariant );
	freeJac( Jac );
	freeMedium( medium );
	
	MPI_Barrier( MPI_COMM_WORLD );
	if ( 0 == thisRank )
	{
		printf( "Finish Run function\n"  );
	}
	//system( "sleep 100s"  );

}
