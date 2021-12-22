/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:propagate.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-03
*   Discription:
*
================================================================*/

#include "header.h"


void isMPIBorder( GRID grid, MPI_COORD thisMPICoord, MPI_BORDER * border )
{

	if ( 0 == thisMPICoord.X ) border->isx1 = 1;   if ( ( grid.PX - 1 ) == thisMPICoord.X ) border->isx2 = 1;
	if ( 0 == thisMPICoord.Y ) border->isy1 = 1;   if ( ( grid.PY - 1 ) == thisMPICoord.Y ) border->isy2 = 1;
	if ( 0 == thisMPICoord.Z ) border->isz1 = 1;   if ( ( grid.PZ - 1 ) == thisMPICoord.Z ) border->isz2 = 1;

}


void propagate( 
MPI_Comm comm_cart, MPI_COORD thisMPICoord, MPI_NEIGHBOR mpiNeighbor,
GRID grid, PARAMS params, SEND_RECV_DATA sr,
WAVE h_W, WAVE W, WAVE t_W, WAVE m_W,
Mat3x3 _rDZ_DX, Mat3x3 _rDZ_DY,
CONTRAVARIANT contravariant, float * Jac, MEDIUM medium,
SLICE slice, SLICE_DATA sliceData, SLICE_DATA sliceDataCpu )
{
	float DT = params.DT;
	int NT = params.TMAX / DT;
	float DH = grid.DH;


	int IT_SKIP = params.IT_SKIP;
	int sliceFreeSurf = params.sliceFreeSurf;
	SLICE freeSurfSlice;
	locateFreeSurfSlice( grid, &freeSurfSlice );
	
	
	
	SLICE_DATA freeSurfData, freeSurfDataCpu;
	PGV pgv, cpuPgv;

	int thisRank;
	MPI_Comm_rank( MPI_COMM_WORLD, &thisRank );


	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;
	

	
	int IsFreeSurface = 0;
#ifdef FREE_SURFACE
	if ( thisMPICoord.Z == grid.PZ - 1 )
		IsFreeSurface = 1;
#endif
	if ( IsFreeSurface )
	{
		allocatePGV( grid, &pgv );
		if ( sliceFreeSurf )
			allocSliceData( grid, freeSurfSlice, &freeSurfData );
#ifdef GPU_CUDA
		allocatePGV_cpu( grid, &cpuPgv );
		if ( sliceFreeSurf )
			allocSliceData_cpu( grid, freeSurfSlice, &freeSurfDataCpu );
#endif
	}
	

	MPI_BORDER border = { 0 };
	isMPIBorder( grid, thisMPICoord, &border );
	
	AUX4 Aux4_1, Aux4_2;


#ifdef PML
	allocPML( grid, &Aux4_1, &Aux4_2, border );

	PML_ALPHA pml_alpha;
	PML_BETA pml_beta;
	PML_D pml_d;

	allocPMLParameter( grid, &pml_alpha, &pml_beta, &pml_d );
	init_pml_parameter( grid, border, pml_alpha, pml_beta, pml_d );

#endif


	SOURCE S = { 0 } ;//{ _nx_ / 2, _ny_ / 2, _nz_ / 2 };
	locateSource( params, grid, &S );
	


	char VxFileName[256], VyFileName[256], VzFileName[256];
	char TxxFileName[256], TyyFileName[256], TzzFileName[256];


	int it = 0, irk = 0;
	int FB1 = 0;	int FB2 = 0;	int FB3 = 0;

	int FB[8][3] =
	{
		{ -1, -1, -1 },
		{  1,  1, -1 },
		{  1,  1,  1 },
		{ -1, -1,  1 },
		{ -1,  1, -1 },
		{  1, -1, -1 },
		{  1, -1,  1 },
		{ -1,  1,  1 },
	};// F = 1, B = -1

	/*
	int FB[8][3] =
	{
		{ -1, -1, -1 },
		{  1,  1,  1 },
		{ -1, -1, -1 },
		{  1,  1,  1 },
		{ -1, -1, -1 },
		{  1,  1,  1 },
		{ -1, -1, -1 },
		{  1,  1,  1 },
	};// F = 1, B = -1
	*/

	GaussField( grid, W );

	for ( it = 0; it < NT; it ++ )
	{
		if ( ( 0 == thisRank ) && ( it % 10 == 0 ) )
			printf( "it = %d\n", it );

		//loadPointSource( S, W, _nx_, _ny_, _nz_, Jac, it, 0, DT, DH );
	
		FB1 = FB[it % 8][0]; FB2 = FB[it % 8][1]; FB3 = FB[it % 8][2];
		for ( irk = 0; irk < 4; irk ++ )
		{
			MPI_Barrier( comm_cart );
			mpiSendRecv( grid, comm_cart, mpiNeighbor, W, sr );
#ifdef PML
			waveDeriv( grid, h_W, W, 
					   contravariant, medium, 
					   pml_beta, 
					   FB1, FB2, FB3);	

			//loadPointSource( S, h_W, _nx_, _ny_, _nz_, Jac, it, irk, DT, DH );
			
			if ( IsFreeSurface ) freeSurfaceDeriv ( grid, h_W, W, 
													contravariant, medium, Jac,
													_rDZ_DX, _rDZ_DY, 
													pml_beta, 
													FB1, FB2, FB3 );
			pmlDeriv( grid, h_W, W, 
					  contravariant, medium, 
					  Aux4_1, Aux4_2, 
					  pml_alpha,  pml_beta, pml_d,	
					  border, 
					  FB1, FB2, FB3 );
				
			if ( IsFreeSurface ) pmlFreeSurfaceDeriv( grid, h_W, W, 
													  contravariant, medium, 
													  Aux4_1, Aux4_2, 
													  _rDZ_DX, _rDZ_DY, 
													  pml_d, 
													  border, 
													  FB1, FB2 );

			waveRk( grid, irk, h_W.Vx, W.Vx, t_W.Vx, m_W.Vx, DT );
			pmlRk( grid, border, irk, Aux4_1, Aux4_2, DT );
#else
			waveDeriv( grid, h_W, W, 
					   contravariant, medium, 
					   FB1, FB2, FB3 );	
			
			//loadPointSource( S, h_W, _nx_, _ny_, _nz_, Jac, it, irk, DT, DH );
		if ( it % IT_SKIP == 0 && irk == 0  )
		{
			sprintf( TxxFileName, "%s/Txx_%d", params.OUT, it );
			sprintf( TyyFileName, "%s/Tyy_%d", params.OUT, it );
			sprintf( TzzFileName, "%s/Tzz_%d", params.OUT, it );

			data2D_output_bin( grid, slice, thisMPICoord, W.Txx, sliceData, sliceDataCpu, TxxFileName );
			data2D_output_bin( grid, slice, thisMPICoord, W.Tyy, sliceData, sliceDataCpu, TyyFileName );
			data2D_output_bin( grid, slice, thisMPICoord, W.Tzz, sliceData, sliceDataCpu, TzzFileName );

			memset( TxxFileName, 0, 256 );
			memset( TyyFileName, 0, 256 );
			memset( TzzFileName, 0, 256 );
		}
			
		if ( it % IT_SKIP == 0  )
		{
			sprintf( VxFileName, "%s/Vx_%d", params.OUT, it );
			sprintf( VyFileName, "%s/Vy_%d", params.OUT, it );
			sprintf( VzFileName, "%s/Vz_%d", params.OUT, it );

			data2D_output_bin( grid, slice, thisMPICoord, W.Vx, sliceData, sliceDataCpu, VxFileName );
			data2D_output_bin( grid, slice, thisMPICoord, W.Vy, sliceData, sliceDataCpu, VyFileName );
			data2D_output_bin( grid, slice, thisMPICoord, W.Vz, sliceData, sliceDataCpu, VzFileName );

			memset( VxFileName, 0, 256 );
			memset( VyFileName, 0, 256 );
			memset( VzFileName, 0, 256 );
		}
			if ( IsFreeSurface ) freeSurfaceDeriv ( grid, h_W, W, 
													contravariant, medium, Jac,
													_rDZ_DX, _rDZ_DY, FB1, FB2, FB3 );
			waveRk( grid, irk, h_W.Vx, W.Vx, t_W.Vx, m_W.Vx, DT );
#endif
			FB1 *= - 1; FB2 *= - 1; FB3 *= - 1; //reverse 

		} // for loop of irk: Range Kutta Four Step

		if ( IsFreeSurface ) 	comparePGV( grid, thisMPICoord, W, pgv );





		if ( it % IT_SKIP == 100000000  )
		{
			sprintf( VxFileName, "%s/Vx_%d", params.OUT, it );
			sprintf( VyFileName, "%s/Vy_%d", params.OUT, it );
			sprintf( VzFileName, "%s/Vz_%d", params.OUT, it );

			data2D_output_bin( grid, slice, thisMPICoord, W.Vx, sliceData, sliceDataCpu, VxFileName );
			data2D_output_bin( grid, slice, thisMPICoord, W.Vy, sliceData, sliceDataCpu, VyFileName );
			data2D_output_bin( grid, slice, thisMPICoord, W.Vz, sliceData, sliceDataCpu, VzFileName );

			memset( VxFileName, 0, 256 );
			memset( VyFileName, 0, 256 );
			memset( VzFileName, 0, 256 );
		}
		if ( it % IT_SKIP == 0 && sliceFreeSurf && IsFreeSurface ) 
		{
			sprintf( VxFileName, "%s/FreeSurfVx_%d", params.OUT, it );
			sprintf( VyFileName, "%s/FreeSurfVy_%d", params.OUT, it );
			sprintf( VzFileName, "%s/FreeSurfVz_%d", params.OUT, it );

			data2D_output_bin( grid, freeSurfSlice, thisMPICoord, W.Vx, freeSurfData, freeSurfDataCpu, VxFileName );
			data2D_output_bin( grid, freeSurfSlice, thisMPICoord, W.Vy, freeSurfData, freeSurfDataCpu, VyFileName );
			data2D_output_bin( grid, freeSurfSlice, thisMPICoord, W.Vz, freeSurfData, freeSurfDataCpu, VzFileName );

			memset( VxFileName, 0, 256 );
			memset( VyFileName, 0, 256 );
			memset( VzFileName, 0, 256 );

		}
		
	
	}// for loop of it: The time iterator of NT steps
	

#ifdef PML

	freePML( border, Aux4_1, Aux4_2 );

	freePMLParamter( pml_alpha, pml_beta, pml_d );

#endif


	if ( IsFreeSurface )
	{

		outputPGV( params, grid, thisMPICoord, pgv, cpuPgv );
		freePGV( pgv );
		if ( sliceFreeSurf ) 
			freeSliceData( grid, freeSurfSlice, freeSurfData );
#ifdef GPU_CUDA
		freePGV_cpu( cpuPgv );
		if ( sliceFreeSurf ) 
			freeSliceData_cpu( grid, freeSurfSlice, freeSurfDataCpu );
#endif
	}
}




