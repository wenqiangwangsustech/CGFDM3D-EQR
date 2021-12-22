/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:medium.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-10-31
*   Discription:
*
================================================================*/
#include "header.h"

void allocMedium( GRID grid, MEDIUM * medium )
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;
	
	long long num = _nx_ * _ny_ * _nz_; 

	float * pMedium = NULL;
	long long size = sizeof( float ) * num * MEDIUMSIZE;

	CHECK( Malloc( ( void ** )&pMedium, size ) );
	CHECK( Memset(  pMedium, 0, size ) ); 
	
	medium->mu       = pMedium;
	medium->lambda   = pMedium + num;
	medium->buoyancy = pMedium + num * 2;
}

void freeMedium( MEDIUM medium )
{	
	Free( medium.mu );
}

//homogeneous medium
__GLOBAL__
void construct_homo_medium( STRUCTURE structure, int _nx_, int _ny_, int _nz_)
{
	long long index = 0;
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
		structure.Vs [index] = 3464.0f;
		structure.Vp [index] = 6000.0f;
		structure.rho[index] = 2670.0f;

    END_CALCULATE3D( )

}

__GLOBAL__
void calculate_medium( MEDIUM medium, int num )
{

#ifdef GPU_CUDA
	int i = threadIdx.x + blockIdx.x * blockDim.x;
#else
	int i = 0;
#endif
	float Vp, Vs, rho;

	CALCULATE1D( i, 0, num )	
		Vs  = medium.mu[i];
		Vp  = medium.lambda[i];
		rho = medium.buoyancy[i];
		//if ( i == num / 2 )
		//{
		//	printf( "Vs = %e, Vp = %e, rho = %e\n", Vs, Vp, rho );
		//}
		medium.mu[i] = rho * ( Vs * Vs );
		medium.lambda[i] = rho * ( Vp * Vp - 2.0f * Vs * Vs );
		medium.buoyancy[i] = 1.0f / rho;
		//if ( i == num / 2 )
		//{
		//	printf( "mu = %e, lambda = %e, buoyancy = %e\n", medium.mu[i], medium.lambda[i], medium.buoyancy[i] );
		//}
	END_CALCULATE1D( )
}


void constructMedium( GRID grid, MEDIUM medium, PARAMS params )
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;


	STRUCTURE structure = { 0 };

	structure.Vs  = medium.mu;
	structure.Vp  = medium.lambda;
	structure.rho = medium.buoyancy;

	//printf( "Vs: %p, Vp: %p, rho: %p\n", structure.Vs, structure.Vp, structure.rho  );
#ifdef GPU_CUDA
	dim3 threads( 32, 4, 4);
	dim3 blocks;
	blocks.x = ( _nx_ + threads.x - 1 ) / threads.x;
	blocks.y = ( _ny_ + threads.y - 1 ) / threads.y;
	blocks.z = ( _nz_ + threads.z - 1 ) / threads.z;

	construct_homo_medium <<< blocks, threads >>>( structure, _nx_, _ny_, _nz_ );
#else
	construct_homo_medium( structure, _nx_, _ny_, _nz_ );
#endif


}

void vs_vp_rho2lam_mu_bou( MEDIUM medium, GRID grid )
{
	
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	long long num = _nx_ * _ny_ * _nz_;

#ifdef GPU_CUDA

	dim3 threadX;
	dim3 blockX;
	threadX.x = 512;
	threadX.y = 1;
	threadX.z = 1;
	blockX.x = ( num + threadX.x - 1 ) / threadX.x;
	blockX.y = 1;
	blockX.z = 1;

	calculate_medium <<< blockX, threadX >>> ( medium, num );

#else
	calculate_medium( medium, num );
#endif

}

