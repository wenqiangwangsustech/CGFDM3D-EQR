/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:propagate.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-05
*   Discription:
*
================================================================*/

#include "header.h"

#define P_alpha 1.0f
#define P_beta 2.0f
#define P_d 2.0f

#define cal_pml_alpha( l ) ( ( ( l ) < 0.0f ) ? 0.0f : ( pml_alpha0 * ( 1.0f - pow( ( l ) / L, P_alpha ) ) ) )
#define cal_pml_beta( l )  ( ( ( l ) < 0.0f ) ? 1.0f : ( 1.0f + ( pml_beta0 - 1.0f ) * pow( ( l ) / L, P_beta ) ) )
#define cal_pml_d( l )     ( ( ( l ) < 0.0f ) ? 0.0f : ( d0 * pow( ( l ) / L, P_d ) ) )

void allocPMLParameter( GRID grid, PML_ALPHA * pml_alpha, PML_BETA *pml_beta, PML_D * pml_d )
{

	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	float * pParams = NULL;
	long long size = sizeof( float ) * ( _nx_ + _ny_ + _nz_ ) * 3;

	CHECK( Malloc( ( void ** )&pParams, size ) );
	CHECK( Memset(  pParams, 0, size ) ); 
	
	
	pml_alpha->x = pParams + 0 * _nx_;
	pml_beta->x  = pParams + 1 * _nx_;
	pml_d->x     = pParams + 2 * _nx_; 

	pParams = pParams + 3 * _nx_;

	pml_alpha->y = pParams + 0 * _ny_;
	pml_beta->y  = pParams + 1 * _ny_;
	pml_d->y     = pParams + 2 * _ny_;

	pParams = pParams + 3 * _ny_;

	pml_alpha->z = pParams + 0 * _nz_;
	pml_beta->z  = pParams + 1 * _nz_;
	pml_d->z     = pParams + 2 * _nz_;

}

void freePMLParamter(
	PML_ALPHA pml_alpha, PML_BETA pml_beta, PML_D pml_d )
{
	Free( pml_alpha.x );
}

void cpu_allocate_pml_para_memory( GRID grid, PML_ALPHA * pml_alpha, PML_BETA * pml_beta, PML_D * pml_d )
{

	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;
		
	float * pParams = NULL;
	long long size = sizeof( float ) * ( _nx_ + _ny_ + _nz_ ) * 3;

	pParams = ( float * ) malloc( size );	
	memset(  pParams, 0, size ) ; 
	
	pml_alpha->x = pParams + 0 * _nx_;
	pml_beta->x  = pParams + 1 * _nx_;
	pml_d->x     = pParams + 2 * _nx_; 

	
	pParams = pParams + 3 * _nx_;

	pml_alpha->y = pParams + 0 * _ny_;
	pml_beta->y  = pParams + 1 * _ny_;
	pml_d->y     = pParams + 2 * _ny_;

	pParams = pParams + 3 * _ny_;

	pml_alpha->z = pParams + 0 * _nz_;
	pml_beta->z  = pParams + 1 * _nz_;
	pml_d->z     = pParams + 2 * _nz_;

}

void cpu_free_pml_para_memory( PML_ALPHA pml_alpha )
{

	free( pml_alpha.x );

}


int cpu_init_pml_parameter( GRID grid, MPI_BORDER border, PML_ALPHA pml_alpha, PML_BETA pml_beta, PML_D pml_d ) 
{
	int i = 0, j = 0, k = 0;


	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	int _nx = _nx_ - HALO;
	int _ny = _ny_ - HALO;
	int _nz = _nz_ - HALO;


	float DH = grid.DH;

	int nPML = grid.nPML;
	
	float Cp = 6000.0;//5716.0f;//Compressional Wave
	float pml_fc = 0.5;
	float L = DH * ( nPML - 1 );
	float R = pow( 10.0, - ( log10( ( double )nPML ) - 1.0 ) / log10( 2.0 ) - 3.0 );
	float d0 = - ( P_d + 1.0f ) * Cp / ( 2.0f * L ) * log( R );
	float pml_alpha0 = PI * pml_fc;
	float pml_beta0 = 3.0;
	//printf("R = %f, d0 = %f\n", R, d0 );

	float * Lx = ( float * )malloc( sizeof( float ) * _nx_ );
	float * Ly = ( float * )malloc( sizeof( float ) * _ny_ );
	float * Lz = ( float * )malloc( sizeof( float ) * _nz_ );

	memset( Lx, 0, sizeof( float ) * _nx_ );
	memset( Ly, 0, sizeof( float ) * _ny_ );
	memset( Lz, 0, sizeof( float ) * _nz_ );

	for ( i = 0; i < _nx_; ++i )	Lx[i] = - 1.0f;
	for ( j = 0; j < _ny_; ++j )	Ly[j] = - 1.0f;
	for ( k = 0; k < _nz_; ++k )	Lz[k] = - 1.0f;

	for ( i = HALO; i < _nx; ++ i )
	{
		
		if ( i < HALO + nPML && border.isx1 )
		{
			Lx[i] = ( HALO + nPML - 1 - i ) * DH;
		}

		if ( i >= _nx - nPML && border.isx2 )
		{
			Lx[i] = ( i - ( _nx - nPML ) ) * DH;
		}
	}

	for ( j = HALO; j < _ny; ++ j )
	{
		if ( j < HALO + nPML && border.isy1 )
		{
			Ly[j] = ( HALO + nPML - 1 - j ) * DH;
		}

		if ( j >= _ny - nPML && border.isy2 )
		{
			Ly[j] = ( j - ( _ny - nPML ) ) * DH;
		}
	}

	for ( k = HALO; k < _nz; ++ k )
	{
		if ( k < HALO + nPML && border.isz1 )
		{
			Lz[k] = ( HALO + nPML - 1 - k ) * DH;
		}
#ifndef FREE_SURFACE
		if ( k >= _nz - nPML && border.isz2 )
		{
			Lz[k] = ( k - ( _nz - nPML ) ) * DH;
		}
#endif
	}

	for ( i = 0; i < _nx_; ++i )	pml_alpha.x[i] = 0.0f;
	for ( j = 0; j < _ny_; ++j )	pml_alpha.y[j] = 0.0f;
	for ( k = 0; k < _nz_; ++k )	pml_alpha.z[k] = 0.0f;

	for ( i = 0; i < _nx_; ++i )	pml_beta.x[i] = 1.0f;
	for ( j = 0; j < _ny_; ++j )	pml_beta.y[j] = 1.0f;
	for ( k = 0; k < _nz_; ++k )	pml_beta.z[k] = 1.0f;

	for ( i = 0; i < _nx_; ++i )	pml_beta.x[i] = 1.0f;
	for ( j = 0; j < _ny_; ++j )	pml_beta.y[j] = 1.0f;
	for ( k = 0; k < _nz_; ++k )	pml_beta.z[k] = 1.0f;


	for ( i = HALO; i < _nx; ++ i )	
	{
		pml_alpha.x[i] = cal_pml_alpha( Lx[i] );
		pml_beta. x[i] = cal_pml_beta( Lx[i] );
		pml_d.	  x[i] = cal_pml_d( Lx[i] );
	}
	for ( j = HALO; j < _ny; ++ j )
	{	
		pml_alpha.y[j] = cal_pml_alpha( Ly[j] );
		pml_beta. y[j] = cal_pml_beta( Ly[j] );
		pml_d.	  y[j] = cal_pml_d( Ly[j] );
	}
	for ( k = HALO; k < _nz; ++ k )
	{
		pml_alpha.z[k] = cal_pml_alpha( Lz[k] );
		pml_beta. z[k] = cal_pml_beta( Lz[k] );
		pml_d.	  z[k] = cal_pml_d( Lz[k] );
	}

	for ( i = 0; i < _nx_; ++i )	
	{
		pml_d.x[i] /= pml_beta.x[i];
		pml_beta.x[i] = 1.0f / pml_beta.x[i];
	
	}
	for ( j = 0; j < _ny_; ++j )
	{
		pml_d.y[j] /= pml_beta.y[j];
		pml_beta.y[j] = 1.0f / pml_beta.y[j];
	}	
	for ( k = 0; k < _nz_; ++k )
	{
		pml_d.z[k] /= pml_beta.z[k];
		pml_beta.z[k] = 1.0f / pml_beta.z[k];
	}	


	free( Lx );
	free( Ly );
	free( Lz );

	return 0;
}

void init_pml_parameter( PARAMS params, GRID grid, MPI_BORDER border, PML_ALPHA pml_alpha, PML_BETA pml_beta, PML_D pml_d )
{
	
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

#ifdef GPU_CUDA
	PML_ALPHA cpu_pml_alpha;
	PML_BETA cpu_pml_beta;
	PML_D cpu_pml_d;


	cpu_allocate_pml_para_memory( grid, &cpu_pml_alpha, &cpu_pml_beta, &cpu_pml_d );


	cpu_init_pml_parameter( grid, border, cpu_pml_alpha, cpu_pml_beta, cpu_pml_d );

	//for ( int i = 0; i < 15; i ++ )
	//	printf( "pml_alpha_x = %f\n", cpu_pml_alpha.x[i]  );
	long long size = sizeof( float ) * ( _nx_ + _ny_ + _nz_ ) * 3;
	cudaMemcpy( pml_alpha.x, cpu_pml_alpha.x, size, cudaMemcpyHostToDevice );
	CHECK( cudaDeviceSynchronize( ) );

		
	//FILE * fp;
	//fp = fopen( "./output2/pmlParams", "wb" );
	//	
	//fwrite( cpu_pml_alpha.x, 4, ( _nx_ + _ny_ + _nz_ ) * 3, fp );

	//fclose( fp );
	cpu_free_pml_para_memory( cpu_pml_alpha );

#else
	
	cpu_init_pml_parameter( grid, border,pml_alpha, pml_beta, pml_d );
	FILE * fp;

	char fileName[256];
	sprintf( fileName, "%s/pmlParams", params.OUT );
	fp = fopen( fileName, "wb" );
		
	fwrite( pml_alpha.x, 4, ( _nx_ + _ny_ + _nz_ ) * 3, fp );

	fclose( fp );

#endif


}

