#include "header.h"

/*
Author Wenqiang Wang@SUSTech
Created on 2020-06-09
 */

int read_source_info( MPI_COORDINATE thisMpiCoord, SOURCE_INFO * srcInfo )
{
	FILE * fp;
	char fileName[256];

	sprintf( fileName, "source/source_mpi_%d_%d_%d", thisMpiCoord.X, thisMpiCoord.Y, thisMpiCoord.Z  );
	//printf("read Source Info\n");
	fp = fopen( fileName, "rb" );

	if ( fp == NULL )
	{
		printf("Can't open source file!\n");
		return 0;
	}

	fread( &( srcInfo->npts ), sizeof( int ), 1, fp );
	fread( &( srcInfo->nt ), sizeof( int ), 1, fp );
	fread( &( srcInfo->dt ), sizeof( float ), 1, fp );


	fclose( fp );


	printf("npts = %d\n", srcInfo->npts );
	printf("nt = %d\n", srcInfo->nt );
	printf("dt = %f\n", srcInfo->dt );
	return 1;
}


void cpu_allocate_moment_memory( SOURCE_INFO srcInfo, SOURCE_INDEX * srcIndex, MOMENT_RATE * momentRate )
{
	int len;
	
	int * pSrcIndex = NULL;
	float * pMomRate = NULL;
	int N = 0;

	len = sizeof( int ) * srcInfo.npts * 3; //X Y Z
	pSrcIndex = ( int * )malloc( len );
	memset( pSrcIndex, 0, len );
	srcIndex->X = pSrcIndex;
	srcIndex->Y = srcIndex->X + srcInfo.npts;
	srcIndex->Z = srcIndex->Y + srcInfo.npts;

	N = srcInfo.npts * srcInfo.nt;
	len = sizeof( float ) * N * 6;
	pMomRate = ( float * )malloc( len );
	memset( pMomRate, 0, len );

	momentRate->Mxx = pMomRate;
	momentRate->Myy = momentRate->Mxx + N;
	momentRate->Mzz = momentRate->Myy + N;
	momentRate->Mxy = momentRate->Mzz + N;
	momentRate->Mxz = momentRate->Mxy + N;
	momentRate->Myz = momentRate->Mxz + N;

}


void allocate_moment_memory( SOURCE_INFO srcInfo, SOURCE_INDEX * srcIndex, MOMENT_RATE * momentRate )
{
	int len;
	
	int * pSrcIndex = NULL;
	float * pMomRate = NULL;
	int N = 0;
    
	len = sizeof( int ) * srcInfo.npts * 3; //X Y Z
	cudaCheck( cudaMalloc( &pSrcIndex, len ) );

	cudaCheck( cudaMemset( pSrcIndex, 0, len ) );
	srcIndex->X = pSrcIndex;
	srcIndex->Y = srcIndex->X + srcInfo.npts;
	srcIndex->Z = srcIndex->Y + srcInfo.npts;

    printf("allloc moment==============1" );
	N = srcInfo.npts * srcInfo.nt;
	len = sizeof( float ) * N * 6;
	cudaCheck( cudaMalloc( &pMomRate, len ) );
	cudaMemset( pMomRate, 0, len );
    printf("allloc moment==============2" );

	momentRate->Mxx = pMomRate;
	momentRate->Myy = momentRate->Mxx + N;
	momentRate->Mzz = momentRate->Myy + N;
	momentRate->Mxy = momentRate->Mzz + N;
	momentRate->Mxz = momentRate->Mxy + N;
	momentRate->Myz = momentRate->Mxz + N;

}

void cpu_free_moment_rate_memory( SOURCE_INDEX srcIndex, MOMENT_RATE momentRate  )
{
	free( srcIndex.X );
	free( momentRate.Mxx );
}

void free_moment_rate_memory( SOURCE_INDEX srcIndex, MOMENT_RATE momentRate  )
{
	cudaCheck( cudaFree( srcIndex.X ) );
	cudaCheck( cudaFree( momentRate.Mxx ) );
}


void allocate_moment_slice_memory( SOURCE_INFO srcInfo, MOMENT_RATE * momentRateSlice )
{
	int len = sizeof( float ) * srcInfo.npts * 6;
	float * pMomentRateSlice = NULL;

	cudaCheck( cudaMalloc( &pMomentRateSlice, len ) );
	cudaCheck( cudaMemset( pMomentRateSlice, 0, len ) );
	momentRateSlice->Mxx = pMomentRateSlice;
	momentRateSlice->Myy = momentRateSlice->Mxx + srcInfo.npts;
	momentRateSlice->Mzz = momentRateSlice->Myy + srcInfo.npts;
	momentRateSlice->Mxy = momentRateSlice->Mzz + srcInfo.npts;
	momentRateSlice->Mxz = momentRateSlice->Mxy + srcInfo.npts;
	momentRateSlice->Myz = momentRateSlice->Mxz + srcInfo.npts;
}

void free_moment_rate_slice_memory( MOMENT_RATE momentRateSlice )
{
	cudaCheck( cudaFree( momentRateSlice.Mxx ) );
}


int read_moment_rate( MPI_COORDINATE thisMpiCoord, SOURCE_INFO srcInfo,  SOURCE_INDEX srcIndex, MOMENT_RATE momentRate )
{
	FILE * fp;
	char fileName[256];
	int i, j, pos;

	sprintf( fileName, "source/source_mpi_%d_%d_%d", thisMpiCoord.X, thisMpiCoord.Y, thisMpiCoord.Z  );

	//printf("npts = %d\n", srcInfo.npts );
	//printf("nt = %d\n", srcInfo.nt );
	//printf("dt = %f\n", srcInfo.dt );

	fp = fopen( fileName, "rb" );

	if ( fp == NULL )
	{
		return 0;
	}

	fseek( fp, sizeof( int ) + sizeof( int ) + sizeof( float ), SEEK_SET );

	for( i = 0; i < srcInfo.npts; ++ i ) 
	{
		fread( srcIndex.X + i, sizeof( int ), 1, fp );
		fread( srcIndex.Y + i, sizeof( int ), 1, fp );
		fread( srcIndex.Z + i, sizeof( int ), 1, fp );
		//printf("X = %d, Y = %d, Z = %d\n", srcIndex.X[i], srcIndex.Y[i], srcIndex.Z[i] );
		for( j = 0; j < srcInfo.nt; ++ j ) 
		{
			pos = i + j * srcInfo.npts;
			//printf("pos = %d\n", pos );
			fread( momentRate.Mxx + pos, sizeof( float ), 1, fp );
			// if ( momentRate.Mxx[pos] > 2692 || momentRate.Mxx[pos] < -1072468 )
			// {
				// printf("it = %d\n", j );
				// printf("momentRate.Mxx[%d] = %f\n", pos, momentRate.Mxx[pos] );
			// }
			//printf("momentRate.Mxx[%d] = %f\n", pos, momentRate.Mxx[pos]  );
		}
		for( j = 0; j < srcInfo.nt; ++ j ) 
		{
			pos = i + j * srcInfo.npts;
			fread( momentRate.Myy + pos, sizeof( float ), 1, fp );
		}

		for( j = 0; j < srcInfo.nt; ++ j ) 
		{
			pos = i + j * srcInfo.npts;
			fread( momentRate.Mzz + pos, sizeof( float ), 1, fp );
		}


		for( j = 0; j < srcInfo.nt; ++ j ) 
		{
			pos = i + j * srcInfo.npts;
			fread( momentRate.Mxy + pos, sizeof( float ), 1, fp );
		}


		for( j = 0; j < srcInfo.nt; ++ j ) 
		{
			pos = i + j * srcInfo.npts;
			fread( momentRate.Mxz + pos, sizeof( float ), 1, fp );
		}


		for( j = 0; j < srcInfo.nt; ++ j ) 
		{
			pos = i + j * srcInfo.npts;
			//printf("pos = %d\n", pos );
			//printf("momentRate.Myz[%d] = %f\n", pos, momentRate.Myz[pos]  );
			
			fread( momentRate.Myz + pos, sizeof( float ), 1, fp );
			// if ( momentRate.Myz[pos] > 523348 )
			// {
			// 	/* code */printf("momentRate.Myz[%d] = %f\n", pos, momentRate.Myz[pos]  );
			// }

			// if ( momentRate.Myz[pos] < -11212 )
			// {
			// 	/* code */printf("momentRate.Myz[%d] = %f\n", pos, momentRate.Myz[pos]  );
			// }

			
		}

	}
	fclose( fp );

	return 1;
}



__global__
void calculate_moment_rate( 
	SOURCE_INFO srcInfo, SOURCE_INDEX srcIndex, MOMENT_RATE momentRate,
	MEDIUM medium, int NX, int NY, int NZ )
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int index = 0, pos = 0;
	if ( i < srcInfo.npts && j < srcInfo.nt  )
	{
		index = INDEX( srcIndex.X[i], srcIndex.Y[i], srcIndex.Z[i] );
		pos = i + j * srcInfo.npts;
		momentRate.Mxx[pos] *= medium.mu[index];
		momentRate.Myy[pos] *= medium.mu[index];
		momentRate.Mzz[pos] *= medium.mu[index];
		momentRate.Mxy[pos] *= medium.mu[index];
		momentRate.Mxz[pos] *= medium.mu[index];
		momentRate.Myz[pos] *= medium.mu[index];
		// if ( medium.mu[index] > 1e7 )
		// {
		//	printf("===========\n");
		// }
		//printf("medium[%d] = %f\n", index, medium.mu[index] );
		//printf("momentRate.Mxx[%d] = %f\n", i, momentRate.Mxx[i] );
	}
}

__global__
void interpolate_moment_rate( 
	SOURCE_INFO srcInfo, MOMENT_RATE momentRate, 
	MOMENT_RATE momentRateSlice, int srcIt, float t_weight )
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int pos1 = 0, pos2 = 0;
	if( i < srcInfo.npts )
	{
		pos1 = srcIt * srcInfo.npts + i;
		pos2 = ( srcIt + 1 )* srcInfo.npts + i;
		//if ( momentRate.Mxx[pos2] > 2692 || momentRate.Mxx[pos2] < -1072468 )
		//{
		//	printf("it = %d\n", srcIt );
		//	printf("momentRate.Mxx[%d] = %f\n", pos2, momentRate.Mxx[pos2] );
		//}

		momentRateSlice.Mxx[i] = ( momentRate.Mxx[pos2] - momentRate.Mxx[pos1] ) * t_weight + momentRate.Mxx[pos1];
		momentRateSlice.Myy[i] = ( momentRate.Myy[pos2] - momentRate.Myy[pos1] ) * t_weight + momentRate.Myy[pos1];
		momentRateSlice.Mzz[i] = ( momentRate.Mzz[pos2] - momentRate.Mzz[pos1] ) * t_weight + momentRate.Mzz[pos1];
		momentRateSlice.Mxy[i] = ( momentRate.Mxy[pos2] - momentRate.Mxy[pos1] ) * t_weight + momentRate.Mxy[pos1];
		momentRateSlice.Mxz[i] = ( momentRate.Mxz[pos2] - momentRate.Mxz[pos1] ) * t_weight + momentRate.Mxz[pos1];
		momentRateSlice.Myz[i] = ( momentRate.Myz[pos2] - momentRate.Myz[pos1] ) * t_weight + momentRate.Myz[pos1];
		//printf(" momentRateSlice.Mxx[%d] = %f\n", i,  momentRateSlice.Mxx[i] );
		// if ( momentRateSlice.Mxx[i] > 0.00001 || momentRateSlice.Mxx[i] < -0.00001 )
		// {
		// 	printf("momentRateSlice.Mxx[%d] = %f\n", i, momentRateSlice.Mxx[i] );
		// }

	}
}

__global__
void add_moment_rate( 
	SOURCE_INFO srcInfo, SOURCE_INDEX srcIndex, MOMENT_RATE momentRateSlice, 
	WAVE h_W, float * Jac, int NX, int NY, int NZ, float DH3,
	int gaussI, int gaussJ, int gaussK, float factorGauss, int flagSurf )
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int index = 0;
	float V = 0.0;
	if ( i < srcInfo.npts )
	{
		if ( flagSurf == 1 && ( srcIndex.Z[i] + gaussK ) > NZ - 4 )
		{
			factorGauss = 0.0;
		}
		if ( ( flagSurf == 1 && srcIndex.Z[i] == NZ - 4 && gaussK <  0 ) ||
			 ( flagSurf == 1 && srcIndex.Z[i] == NZ - 5 && gaussK < -1 ) ||
			 ( flagSurf == 1 && srcIndex.Z[i] == NZ - 6 && gaussK < -2 ) )
		{
			factorGauss = factorGauss * 2;
		}

		index = INDEX( srcIndex.X[i] + gaussI, srcIndex.Y[i] + gaussJ, srcIndex.Z[i] + gaussK );
		V = Jac[index] * DH3;
		//if ( Jac[index] < 0.4 || Jac[index] > 2 )
		//	("Jac[%d] = %f\n", index, Jac[index] );
		//if ( DH3 + 1 < 500 * 500 * 500 || DH3 - 1 < 500 * 500 * 500 )
		//	printf("Serous error======================\n");
		V = -1.0f / V * factorGauss;
		//printf("DH3[%d] = %f,", index, DH3 );
		h_W.Txx[index] += momentRateSlice.Mxx[i] * V;
		h_W.Tyy[index] += momentRateSlice.Myy[i] * V;
		h_W.Tzz[index] += momentRateSlice.Mzz[i] * V;
		h_W.Txy[index] += momentRateSlice.Mxy[i] * V;
		h_W.Txz[index] += momentRateSlice.Mxz[i] * V;
		h_W.Tyz[index] += momentRateSlice.Myz[i] * V;
		
		//printf("i = %d, X = %d, Y = %d, Z = %d\n", i, srcIndex.X[i], srcIndex.Y[i], srcIndex.Z[i] );
		// if ( h_W.Txx[index] > 0.00001 || h_W.Txx[index] < -0.00001 )
		// {
		//	printf(" momentRateSlice.Mxx[%d] = %f\n", index,  momentRateSlice.Mxx[index] );
		// }


	}
}
