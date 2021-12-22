#include "header.h"

void modelChecking( PARAMS params  )
{
	int thisRank;
	MPI_Barrier( MPI_COMM_WORLD );
	MPI_Comm_rank( MPI_COMM_WORLD, &thisRank );


	if ( params.useMultiSource && params.useSingleSource )
	{
		params.useMultiSource = 0;
		if ( 0 == thisRank )
			printf( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"	
					"You set \"useMultiSource\" and \"useSingleSource\" at the same time. We set \"useSingleSource:\" = 0.\n"
					"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" );	
		
	}
	if ( !( params.useMultiSource || params.useSingleSource ) )
	{

		if ( 0 == thisRank )
			printf( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"	
				    "You did not set any source. The program will abort!\n" 
				    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" );	
		MPI_Barrier( MPI_COMM_WORLD );
		MPI_Abort( MPI_COMM_WORLD, 180 );
	}
	if ( params.ShenModel && params.Crust_1Medel )
	{

		if ( 0 == thisRank )
			printf( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"	
				    "You set ShenModel and Crust_1Medel both to be 1. We will use homogenourse model!\n" 
				    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" );	

		params.ShenModel = 0;
		params.Crust_1Medel = 0;
		MPI_Barrier( MPI_COMM_WORLD );
	}

	if ( params.useTerrain && params.gauss_hill )
	{
		printf( 
		"!!!!!!!!!!!!!!!!!!!!!!!WORNING!!!!!!!!!!!!!!!!!!!!!!!"
		"You set \"useTerrain\" and \"gauss_hill\" can not be 1 at the same time. You should configue your \"params.json\" file. Or, we will use \"gauss_hill\" model to run the simulation\n"
		"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		params.gauss_hill = 1;
		params.useTerrain = 0;
		params.Depth = params.DH * params.NZ * 1e-3;
	}

	if ( params.useTerrain == 0 )
	{

		if ( 0 == thisRank )
			printf( "No SRTM90 Terrain model is used!\n" );
	}
	if ( params.useTerrain == 1 )
	{

		if ( 0 == thisRank )
			printf( "SRTM90 Terrain model is used!\n" );
	}

	if ( params.useMultiSource == 0 )
	{

		if ( 0 == thisRank )
			printf( "No Multi-source model is used!\n" );
	}
	if ( params.useMultiSource == 1 )
	{
		if ( 0 == thisRank )
			printf( "Multi-source model is used!\n" );
	}

	if ( params.ShenModel == 0 )
	{

		if ( 0 == thisRank )
			printf( "No ShenModel is used!\n" );
	}
	if ( params.ShenModel == 1 )
	{

		if ( 0 == thisRank )
			printf( "ShenModel is used!\n" );
	}
	if ( params.Crust_1Medel == 0 )
	{

		if ( 0 == thisRank )
			printf( "No Crust_1Medel is used!\n" );
	}
	if ( params.Crust_1Medel == 1 )
	{

		if ( 0 == thisRank )
			printf( "Crust_1Medel is used!\n" );
	}
	MPI_Barrier( MPI_COMM_WORLD );

}



void run( MPI_Comm comm_cart, MPI_COORD thisMPICoord, MPI_NEIGHBOR mpiNeighbor, GRID grid, PARAMS params )
{
	int thisRank;
	MPI_Barrier( MPI_COMM_WORLD );
	MPI_Comm_rank( MPI_COMM_WORLD, &thisRank );
	if ( thisRank == 0 )
		printInfo( grid );
	MPI_Barrier( MPI_COMM_WORLD );

	modelChecking( params );


	SEND_RECV_DATA sr;
	allocSendRecv( grid, mpiNeighbor, &sr );


	SLICE_DATA sliceData, sliceDataCpu;
	SLICE slice = { 0 };
	locateSlice( params, grid, &slice );
	allocSliceData( grid, slice, &sliceData );
#ifdef GPU_CUDA
	allocSliceData_cpu( grid, slice, &sliceDataCpu );
#endif

	
	//construct coordinate
	if ( thisRank == 0 )
		printf( "Construct coordinate including precessing terrian data...\n"  );
	MPI_Barrier( MPI_COMM_WORLD );
	COORD coord, cpu_coord;
	allocCoord( grid, &coord );
#ifdef GPU_CUDA
	allocCoord_cpu( grid, &cpu_coord );
	constructCoord(comm_cart, thisMPICoord, grid, params, coord, cpu_coord );
#else                                                           
	constructCoord(comm_cart, thisMPICoord, grid, params, coord );
#endif 

	//construct medium
	if ( thisRank == 0 )
		printf( "Construct medium including precessing Vs Vp Rho...\n"  );
	MPI_Barrier( MPI_COMM_WORLD );
	MEDIUM medium;
	allocMedium( grid, &medium );
#ifdef GPU_CUDA
	MEDIUM cpu_medium;
	allocMedium_cpu( grid, &cpu_medium );
	constructMedium( thisMPICoord, params, grid, cpu_coord, medium, cpu_medium );
	freeMedium_cpu( cpu_medium );
#else
	constructMedium( thisMPICoord, params, grid, coord, medium );
#endif



	//read multisource model
	long long * srcIndex;
	MOMENT_RATE momentRate, momentRateSlice;
	SOURCE_FILE_INPUT src_in;

	if ( params.useMultiSource )
	{
#ifdef GPU_CUDA
		init_MultiSource( params, grid, thisMPICoord, cpu_coord, &srcIndex, &momentRate, &momentRateSlice, &src_in );
#else
		init_MultiSource( params, grid, thisMPICoord, coord, &srcIndex, &momentRate, &momentRateSlice, &src_in );
#endif
	}

	//calculate CFL condition
	calc_CFL( grid, coord, medium, params );

	if ( thisRank == 0 )
		printf( "Slice Position Coordinate(x, y, z) and Medium(Vp, Vs, Rho) data output...\n"  );
	MPI_Barrier( MPI_COMM_WORLD );
	data2D_Model_out( thisMPICoord, params, grid, coord, medium, slice, sliceData, sliceDataCpu );



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
#ifdef GPU_CUDA
	freeCoord_cpu( cpu_coord );
#endif




		
	MPI_Barrier( MPI_COMM_WORLD );

	WAVE h_W, W, t_W, m_W;
	allocWave( grid, &h_W, &W, &t_W, &m_W );

	if ( thisRank == 0 )
		printf( "Start calculating Wave Field:\n"  );
	MPI_Barrier( MPI_COMM_WORLD );
	propagate( comm_cart, thisMPICoord, mpiNeighbor,
			   grid, params, sr,
			   h_W, W, t_W, m_W,
			   _rDZ_DX, _rDZ_DY,
			   contravariant, Jac, medium,
			   src_in, srcIndex, momentRate, momentRateSlice,
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
