/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:source.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-06
*   Discription:
*
================================================================*/

#include "header.h"
void locateSource( PARAMS params, GRID grid, SOURCE * source )
{
	int sourceX = params.sourceX;
	int sourceY = params.sourceY;
	int sourceZ = params.sourceZ;
	
	int frontNX = grid.frontNX;
	int frontNY = grid.frontNY;
	int frontNZ = grid.frontNZ;
	
 	source->X = sourceX - frontNX + HALO;
 	source->Y = sourceY - frontNY + HALO;
 	source->Z = sourceZ - frontNZ + HALO;
	
	//printf( "source.X = %d, source.Y = %d, source.Z = %d\n", source->X, source->Y, source->Z );
	
}
__DEVICE__
float sourceFunction( float rickerfc, int it, int irk, float DT )
{
    float t = 0.0f;
    float tdelay = 1.2f / rickerfc;
    tdelay = 1.0;

    if( 0 == irk ) 
    {
    	t = ( it + 0.0f ) * DT;
    }else if ( 1 == irk || 2 == irk )
    {
    	t = ( it + 0.5f ) * DT;
    }
    else if ( 3 == irk )
    {
    	t = ( it + 1.0f ) * DT;
    }
    float f0 = sqrt( PI ) * 0.5f;
    float r = PI * rickerfc * ( t - tdelay);
    float rr = r*r;
    float s = r*( 3.0f - 2.0f * rr ) * exp( -rr ) * f0 * PI * rickerfc;
    
    float M0 = 1e16;
    s *= M0;

    return s;
}

__GLOBAL__
void load_point_source( SOURCE S, WAVE h_W, int _nx_, int _ny_, int _nz_, float * Jac, int it, int irk, float DT, float DH )
{
	float amp = 1.0;
	float s = sourceFunction( 2.0f, it, irk, DT );
      
    float value = -1.0f * s * amp  / ( Jac[INDEX( S.X, S.Y, S.Z )] * (  DH * DH * DH ) );

    //printf("index = %d, value = %f\n",  INDEX( S.X, S.Y, S.Z ), value );
	h_W.Txx[INDEX( S.X, S.Y, S.Z )] += value;
	h_W.Tyy[INDEX( S.X, S.Y, S.Z )] += value;
	h_W.Tzz[INDEX( S.X, S.Y, S.Z )] += value;
}


__GLOBAL__
void load_smooth_source( SOURCE S, WAVE h_W, int _nx_, int _ny_, int _nz_, float * Jac, int it, int iRK, float DT, float DH, int nGauss, float rickerfc )
{
#ifdef GPU_CUDA
	int i = threadIdx.x + blockIdx.x * blockDim.x + S.X - nGauss;
	int j = threadIdx.y + blockIdx.y * blockDim.y + S.Y - nGauss;
	int k = threadIdx.z + blockIdx.z * blockDim.z + S.Z - nGauss;
#else
	int i = 0;
	int j = 0;
	int k = 0;
#endif
	long long  index = 0;
	float amp = 0.0f;
	float s = 0.0f;
	float value = 0.0f;


	CALCULATE3D( i, j, k, S.X - nGauss, S.X + nGauss + 1, S.Y - nGauss, S.Y + nGauss + 1, S.Z - nGauss, S.Z + nGauss + 1 )   
		index = INDEX( i, j, k );
		if ( iRK == 0 )
			s = 0.0f;
		s = sourceFunction( rickerfc, it, iRK, DT );
		float ra = nGauss * 0.5;
		float D1 = GAUSS_FUN(i-S.X, ra, 0.0);
		float D2 = GAUSS_FUN(j-S.Y, ra, 0.0);
		float D3 = GAUSS_FUN(k-S.Z, ra, 0.0);
		float amp = D1 * D2 * D3;

		//amp /= 0.998125703461425; // # 3
		amp /= 0.9951563131100551; // # 5


    	value = -1.0f * s * amp  / ( Jac[INDEX( S.X, S.Y, S.Z )] * (  DH * DH * DH ) );
    	//value = -1.0f * s * amp  / (  DH * DH * DH );
		
		//if ( i == S.X && j == S.Y && k == S.Z )
		//printf( "index = %d\n", index  );
		value *= Cs;

    	h_W.Txx[index] += value;
		h_W.Tyy[index] += value;
		h_W.Tzz[index] += value;

	END_CALCULATE3D( )

}



void loadPointSource( SOURCE S, WAVE h_W, int _nx_, int _ny_, int _nz_, float * Jac, int it, int irk, float DT, float DH, float rickerfc )
{
	//load_point_source<<< 1, 1 >>>( S, h_W, _nx_, _ny_, _nz_, Jac, it, irk, DT, DH );
	
	int nGauss = 3;

	int _nx = _nx_ - HALO;
	int _ny = _ny_ - HALO;
	int _nz = _nz_ - HALO;
	

	
#ifdef GPU_CUDA

	dim3 threads( 4, 4, 4);
	dim3 blocks;
	blocks.x = ( 2 * nGauss + 1 + threads.x - 1 ) / threads.x;
	blocks.y = ( 2 * nGauss + 1 + threads.y - 1 ) / threads.y;
	blocks.z = ( 2 * nGauss + 1 + threads.z - 1 ) / threads.z;
	
	if ( S.X >= HALO && S.X < _nx &&
		 S.Y >= HALO && S.Y < _ny &&
		 S.Z >= HALO && S.Z < _nz )
	{
		load_smooth_source
		<<< blocks, threads >>>
		( S, h_W, _nx_, _ny_, _nz_, Jac, it, irk, DT, DH, nGauss, rickerfc );
	}
	CHECK( cudaDeviceSynchronize( ) );
#else

	if ( S.X >= HALO && S.X < _nx &&
		 S.Y >= HALO && S.Y < _ny &&
		 S.Z >= HALO && S.Z < _nz )
	{
		load_smooth_source
		( S, h_W, _nx_, _ny_, _nz_, Jac, it, irk, DT, DH, nGauss, rickerfc );
	}

#endif
}



__GLOBAL__
void Gaussian_pluse( WAVE W, 
int _nx_,    int _ny_,    int _nz_, 
int frontNX, int frontNY, int frontNZ,
int _NX_,    int _NY_,	  int _NZ_ )
{
	float sigma = 10.0f;
#ifdef GPU_CUDA
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;
#else
	int i = 0;
	int j = 0;
	int k = 0;
#endif

	long long  index = 0;

	int _nx = _nx_ - HALO;
	int _ny = _ny_ - HALO;
	int _nz = _nz_ - HALO;

	int I = frontNX + i;
	int J = frontNY + j;
	int K = frontNZ + k;

	CALCULATE3D( i, j, k, HALO, _nx, HALO, _ny, HALO, _nz )   
		index = INDEX( i, j, k );
		float r2 =  - ( ( I - _NX_ / 2 ) * ( I - _NX_ / 2 ) + ( J - _NY_ / 2 ) * ( J - _NY_ / 2 ) + ( K - _NZ_ / 2 ) * ( K - _NZ_ / 2 ) );//powf( i - _nx_ * 0.5, 2 ) + powf( j - _ny_ * 0.5, 2 ) + powf( k - _nz_ * 0.5, 2 );
		float a2 = 2.0f * sigma * sigma;
		float M0 = 1e2f * expf( r2 / a2 );

		W.Txx[ index ] = 1e19f * ( 1.0f / ( PI * a2 ) ) * exp( r2 / a2 );
		W.Tyy[ index ] = 1e19f * ( 1.0f / ( PI * a2 ) ) * exp( r2 / a2 ); 
		W.Tzz[ index ] = 1e19f * ( 1.0f / ( PI * a2 ) ) * exp( r2 / a2 ); 
										 
	END_CALCULATE3D( )
}

void GaussField( GRID grid, WAVE W )
{

	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;


	int frontNX = grid.frontNX;
	int frontNY = grid.frontNY;
	int frontNZ = grid.frontNZ;

	int _NX_ = grid._NX_;
	int _NY_ = grid._NY_;
	int _NZ_ = grid._NZ_;

	int originalX = grid.originalX;
	int originalY = grid.originalY;

#ifdef GPU_CUDA
	
	//printf( "W = %p\n", W.Txx  );

	dim3 threads( 32, 4, 4);
	dim3 blocks;
	blocks.x = ( _nx_ + threads.x - 1 ) / threads.x;
	blocks.y = ( _ny_ + threads.y - 1 ) / threads.y;
	blocks.z = ( _nz_ + threads.z - 1 ) / threads.z;

	Gaussian_pluse<<< blocks, threads >>> ( W, _nx_, _ny_, _nz_, frontNX, frontNY, frontNZ, _NX_,    _NY_, _NZ_ );
#else
	Gaussian_pluse( W, _nx_, _ny_, _nz_, frontNX, frontNY, frontNZ, _NX_,    _NY_, _NZ_ );


#endif


}

