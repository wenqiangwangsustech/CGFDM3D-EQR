/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:coord.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-06
*   Discription:
*
================================================================*/
#include "header.h"

void allocCoord( GRID grid, COORD * coord )
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;
	
	long long num = _nx_ * _ny_ * _nz_; 

		
	float * pCoord = NULL;
	long long size = sizeof( float ) * num * COORDSIZE;

	CHECK( Malloc( ( void ** )&pCoord, size ) );
	CHECK( Memset(  pCoord, 0, size ) ); 
	
	
	coord->x = pCoord;
	coord->y = pCoord + num;
	coord->z = pCoord + num * 2;
}

void freeCoord( COORD coord )
{	
	Free( coord.x );
}

__GLOBAL__
void construct_flat_coord( 
				COORD coord,
				int _nx_, int _ny_, int _nz_,
				int frontNX, int frontNY, int frontNZ, 
				int originalX, int originalY,
				int NZ,
				float DH )
{
	long long index = 0;
	int I = 0, J = 0, K = 0;
#ifdef GPU_CUDA
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int k = blockIdx.z * blockDim.z + threadIdx.z;
#else
	int i = 0;
	int j = 0;
	int k = 0;
#endif

	CALCULATE3D( i, j, k, 0, _nx_, 0, _ny_, 0, _nz_ )
        index = INDEX( i, j, k );
		I = frontNX + i;
		J = frontNY + j;
		K = frontNZ + k;
		coord.x[index] = ( I - HALO ) * DH - originalX * DH;
		coord.y[index] = ( J - HALO ) * DH - originalY * DH;
		coord.z[index] = ( K - HALO ) * DH - NZ * DH ; 
    END_CALCULATE3D( )

}


__GLOBAL__
void construct_gauss_hill_surface( 
			float * DZ,
			int _nx_, int _ny_,
			int frontNX, int frontNY, 
			int originalX, int originalY,
			int NZ,
			float DH, float cal_depth )
{
#ifdef GPU_CUDA
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
#else
	int i = 0;
	int j = 0;
#endif
	
	float x, y;

	float height = 0.0;

	float h = 0.2 * cal_depth;
	float a = 0.1 * cal_depth, b = 0.1 * cal_depth;

	int I = 0, J = 0;

	long long index;

	CALCULATE2D( i, j, 0, _nx_, 0, _ny_ )         
		index = Index2D( i, j, _nx_, _ny_ );
		I = frontNX + i;
		J = frontNY + j;
		x = ( I - HALO ) * DH - originalX * DH;
		y = ( J - HALO ) * DH - originalY * DH;

		height = h * exp( - 0.5f * ( x * x / ( a * a ) + y * y / ( b * b ) ) );
		//printf( "_nx_ = %d, _ny_ = %d\n", _nx_, _ny_ );
		DZ[index] = double( height + abs( cal_depth ) ) / double( NZ - 1 );
		//if ( index == _nx_ * _ny_ - 1 )
		//{
		//	printf( "frontNX = %d\n", frontNX );

		//}
		//if ( index == _nx_ * _ny_ - 2 )
	//	if ( i == 100 )
	//	//if ( DZ[index] > -0.1 && DZ[index] < 0.1  )
	//	{
	//		printf( "1: &DZ = %p\n", DZ + index );
	//	}

	END_CALCULATE2D( )
}

__GLOBAL__
void verifyDZ( float * DZ, int _nx_, int _ny_, int _nz_ )
{
#ifdef GPU_CUDA
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int k = blockIdx.z * blockDim.z + threadIdx.z;
#else
	int i = 0;
	int j = 0;
	int k = 0;
#endif

	long long index;

	CALCULATE3D( i, j, k, 0, _nx_, 0, _ny_, 0, _nz_ )         
		index = Index2D( i, j, _nx_, _ny_ );
		if ( index == _nx_ * _ny_ - 2 )
		//if ( DZ[index] > -0.1 && DZ[index] < 0.1  )
		{
			printf( "2: &DZ = %p\n", DZ + index );
		}

	END_CALCULATE3D( )
}

__GLOBAL__
void construct_terrain_coord( 
			COORD coord, float * DZ,
			int _nx_, int _ny_, int _nz_,
			int frontNZ, 
			int NZ,
			float cal_depth )
{
	long long index = 0, pos = 0;
#ifdef GPU_CUDA
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int k = blockIdx.z * blockDim.z + threadIdx.z;
#else
	int i = 0;
	int j = 0;
	int k = 0;
#endif
	int K;
    CALCULATE3D( i, j, k, 0, _nx_, 0, _ny_, 0, _nz_ )			
		index = INDEX( i, j, k );
		pos = Index2D( i, j, _nx_, _ny_ );
		K = frontNZ + k - HALO; 
		//coord.z[index] = DZ[pos] * ( K + 3 * HALO  );
		coord.z[index] = - abs( cal_depth ) + DZ[pos] * K;
		//coord.z[index] = 2 * abs( cal_depth ) + DZ[pos] * ( K + HALO );
		//if ( DZ[pos] > 250 )
		////if ( pos > _nx_ * _ny_ )
		//{
		//	printf( "pos = %d\n", pos );

		//}
		//if ( index == 207684664 )
		//{
		//	printf( "1: coordZ = %f\n", coord.z[index] );

		//}
		//if ( coord.z[index] > 19999 )
		//{
		//	printf( "index = %d\n", index );

		//}
    END_CALCULATE3D( )

}

//void calculate_range( float * data, float range[2], long long num );


void constructCoord( GRID grid, COORD coord, PARAMS params )
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;


	int frontNX = grid.frontNX;
	int frontNY = grid.frontNY;
	int frontNZ = grid.frontNZ;

	int NZ = grid.NZ;

	int originalX = grid.originalX;
	int originalY = grid.originalY;
	


	float DH = grid.DH;

	float cal_depth = params.Depth * 1000;

#ifdef GPU_CUDA
	dim3 threads( 32, 4, 4);
	dim3 blocks;
	blocks.x = ( _nx_ + threads.x - 1 ) / threads.x;
	blocks.y = ( _ny_ + threads.y - 1 ) / threads.y;
	blocks.z = ( _nz_ + threads.z - 1 ) / threads.z;

	construct_flat_coord
	<<< blocks, threads >>>
	( coord, _nx_, _ny_, _nz_,
	frontNX, frontNY, frontNZ, originalX, originalY, NZ, DH );
#else
	construct_flat_coord
	( coord, _nx_, _ny_, _nz_,
	frontNX, frontNY, frontNZ, originalX, originalY, NZ, DH );
#endif

	int gauss_hill = params.gauss_hill;

	if ( gauss_hill )
	{
		float * DZ;
		long long size = sizeof( float ) * _nx_ * _ny_;
		CHECK( Malloc( ( void ** ) &DZ, size ) );
		Memset( DZ, 0, size );

#ifdef GPU_CUDA
		dim3 threadXY( 32, 16, 1 );
		dim3 blockXY;
		blockXY.x = ( _nx_ + threadXY.x - 1 ) / threadXY.x;
		blockXY.y = ( _ny_ + threadXY.y - 1 ) / threadXY.y;
		blockXY.z = 1;
		construct_gauss_hill_surface
		<<< blockXY, threadXY >>>
		( DZ, _nx_, _ny_, frontNX, frontNY, originalX, originalY, NZ, DH, cal_depth );
		//verifyDZ <<< blocks, threads >>> ( DZ, _nx_, _ny_, _nz_ );
		construct_terrain_coord
		<<< blocks, threads>>>
		( coord, DZ, _nx_, _ny_, _nz_, frontNZ, NZ, cal_depth );
#else
		construct_gauss_hill_surface( DZ, _nx_, _ny_, frontNX, frontNY, originalX, originalY, NZ, DH, cal_depth );
		construct_terrain_coord( coord, DZ, _nx_, _ny_, _nz_, frontNZ, NZ, cal_depth );
#endif

		Free( DZ );
	}

/*
*/


}
