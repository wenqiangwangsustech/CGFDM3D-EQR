/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:readSource.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-17
*   Discription:
*
================================================================*/

#include "headers.h"


int readSourceInfo( MPI_COORD thisMPICoord, SOURCE_FILE_INPUT * src_in )
{
	char fileName[256];
	sprintf( fileName, "output/source_mpi_%d_%d_%d.bin", thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
	
	FILE * file = fopen( fileName, "rb" );

	if ( file == NULL )
	{
		return 0;
	}
	
	fread( &( src_in->npts ), sizeof( int ), 1, file );
	fread( &( src_in->nt ), sizeof( int ), 1, file );
	fread( &( src_in->dt ), sizeof( float ), 1, file );
	
	fclose( file );

	//printf("npts = %d\n", src_in->npts );
	//printf("nt = %d\n", src_in->nt );
	//printf("dt = %f\n", src_in->dt );
	return 1;


}

void allocReadMomentRate( SOURCE_FILE_INPUT src_in, MOMENT_RATE * momentRate )
{

	int size = src_in.npts * src_in.nt;
	momentRate->Mxx = ( float * ) malloc( sizeof( float ) * size );
	momentRate->Myy = ( float * ) malloc( sizeof( float ) * size );
	momentRate->Mzz = ( float * ) malloc( sizeof( float ) * size );
	momentRate->Mxy = ( float * ) malloc( sizeof( float ) * size );
	momentRate->Mxz = ( float * ) malloc( sizeof( float ) * size );
	momentRate->Myz = ( float * ) malloc( sizeof( float ) * size );

}


void freeReadMomentRate( MOMENT_RATE momentRate )
{

	free( momentRate.Mxx );
	free( momentRate.Myy );
	free( momentRate.Mzz );
	free( momentRate.Mxy );
	free( momentRate.Mxz );
	free( momentRate.Myz );

}

void allocReadIndex( SOURCE_FILE_INPUT src_in, long long ** srcIndex )
{
	int npts = src_in.npts;

	*srcIndex = ( long long * ) malloc( sizeof( long long ) * npts ); 

}

void freeReadIndex( long long * srcIndex )
{
	free( srcIndex );
}


void readIndexMomentRate( MPI_COORD thisMPICoord,  SOURCE_FILE_INPUT src_in, long long * srcIndex, MOMENT_RATE momentRate )
{
	char fileName[256];
	sprintf( fileName, "output/source_mpi_%d_%d_%d.bin", thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
	
	float * Mxx = momentRate.Mxx;
	float * Myy = momentRate.Myy;
	float * Mzz = momentRate.Mzz;
	float * Mxy = momentRate.Mxy;
	float * Mxz = momentRate.Mxz;
	float * Myz = momentRate.Myz;
	
	FILE * file = fopen( fileName, "rb" );

	fseek( file, sizeof( int ) + sizeof( int ) + sizeof( float ), SEEK_SET );

	int npts = src_in.npts;
	int nt = src_in.nt;

	int p = 0;
	for ( p = 0; p < npts; p ++ )
	{
		fread( srcIndex + p, sizeof( long long ), 1, file );
		fread( Mxx + p, sizeof( float ), nt, file );
		fread( Myy + p, sizeof( float ), nt, file );
		fread( Mzz + p, sizeof( float ), nt, file );
		fread( Mxy + p, sizeof( float ), nt, file );
		fread( Mxz + p, sizeof( float ), nt, file );
		fread( Myz + p, sizeof( float ), nt, file );
	}

	fclose( file );
}



void verifySrcIndex( MPI_COORD thisMPICoord, long long * srcIndex, int npts  )
{

	char fileName[256];
	sprintf( fileName, "output/srcIndex_%d_%d_%d.bin",  thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
	FILE * file;

	file = fopen( fileName, "wb" ); 
		
	fwrite( srcIndex, sizeof( long long ), npts, file );


	fclose( file );

}


void allocMomentRateSlice( SOURCE_FILE_INPUT src_in, MOMENT_RATE * momentRateSlice )
{
	
	int size = src_in.npts;
	momentRateSlice->Mxx = ( float * ) malloc( sizeof( float ) * size );
	momentRateSlice->Myy = ( float * ) malloc( sizeof( float ) * size );
	momentRateSlice->Mzz = ( float * ) malloc( sizeof( float ) * size );
	momentRateSlice->Mxy = ( float * ) malloc( sizeof( float ) * size );
	momentRateSlice->Mxz = ( float * ) malloc( sizeof( float ) * size );
	momentRateSlice->Myz = ( float * ) malloc( sizeof( float ) * size );

}

void freeMomentRateSlice( MOMENT_RATE momentRateSlice  )
{
	free( momentRateSlice.Mxx );
	free( momentRateSlice.Myy );
	free( momentRateSlice.Mzz );
	free( momentRateSlice.Mxy );
	free( momentRateSlice.Mxz );
	free( momentRateSlice.Myz );
	
}


	
void allocWave( GRID grid, WAVE * hW, float ** Jac )
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;
	
	long long num = _nx_ * _ny_ * _nz_; 
	
	hW->Vx = ( float * )malloc( sizeof( float ) * num );
	hW->Vy = ( float * )malloc( sizeof( float ) * num );
	hW->Vz = ( float * )malloc( sizeof( float ) * num );

	hW->Txx = ( float * )malloc( sizeof( float ) * num );
	hW->Tyy = ( float * )malloc( sizeof( float ) * num );
	hW->Tzz = ( float * )malloc( sizeof( float ) * num );
	hW->Txy = ( float * )malloc( sizeof( float ) * num );
	hW->Txz = ( float * )malloc( sizeof( float ) * num );
	hW->Tyz = ( float * )malloc( sizeof( float ) * num );
	
	*Jac = ( float * )malloc( sizeof( float ) * num );

	memset( hW->Vx , 0, sizeof( float ) * num );
	memset( hW->Vy , 0, sizeof( float ) * num );
	memset( hW->Vz , 0, sizeof( float ) * num );
                   
	memset( hW->Txx, 0, sizeof( float ) * num );
	memset( hW->Tyy, 0, sizeof( float ) * num );
	memset( hW->Tzz, 0, sizeof( float ) * num );
	memset( hW->Txy, 0, sizeof( float ) * num );
	memset( hW->Txz, 0, sizeof( float ) * num );
	memset( hW->Tyz, 0, sizeof( float ) * num );
	
	memset( *Jac, 0, sizeof( float ) * num );

}

void freeWave( WAVE hW, float * Jac  )
{
	
	
	free( hW.Vx );
	free( hW.Vy );
	free( hW.Vz );

	free( hW.Txx );
	free( hW.Tyy );
	free( hW.Tzz );
	free( hW.Txy );
	free( hW.Txz );
	free( hW.Tyz );
	
	free( Jac );


}




void readSource( PARAMS params, GRID grid,  MPI_COORD thisMPICoord )
{

	SOURCE_FILE_INPUT src_in;
	
	int ret = readSourceInfo( thisMPICoord, &src_in );


	MOMENT_RATE momentRate;
	long long * srcIndex;
	MOMENT_RATE momentRateSlice;

	int NT = int( params.TMAX / params.DT );
	float DT = params.DT;
	int nt = src_in.nt;
	int npts = src_in.npts;
	float dt = src_in.dt;

	float t, t0, t_weight;

	int it = 0, srcIt = 0;




	WAVE hW;
	float * Jac;

	allocWave( grid, &hW, &Jac );

	float DH = params.DH;



	if ( 1 == ret  )
	{
		allocReadMomentRate( src_in, &momentRate );
		allocReadIndex( src_in, &srcIndex );

		readIndexMomentRate( thisMPICoord, src_in, srcIndex, momentRate );

		verifySrcIndex(thisMPICoord, srcIndex, src_in.npts );
		allocReadMomentRate( src_in, &momentRateSlice);
	}
	for (  it = 0; it < NT; it ++ )
	{
		if ( 1 == ret )
		{
			t = it * DT;
			srcIt = int( t / src_in.dt );
		}
		if ( 1 == ret && srcIt - 1 < nt )
		{
			t0 = float( srcIt ) * dt;
			t_weight = ( t - t0 ) / dt;
			interpMomentRate( src_in, momentRate, momentRateSlice, t_weight, srcIt );
			addSource( hW, momentRateSlice, srcIndex, npts, DH, Jac );
		}
		if ( 0 == thisMPICoord.X && 0 == thisMPICoord.Y && 0 == thisMPICoord.Z )
			cout << "it = " << it << endl;

		MPI_Barrier( MPI_COMM_WORLD );

	}

	

	if ( 1 == ret )
	{
		
		freeMomentRateSlice( momentRateSlice );
		freeReadIndex( srcIndex );
		freeReadMomentRate( momentRate );
	}
		
	freeWave( hW, Jac );
	MPI_Barrier( MPI_COMM_WORLD );

}




