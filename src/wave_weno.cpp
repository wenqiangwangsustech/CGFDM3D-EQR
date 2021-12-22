/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:wave.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-03
*   Discription:
*
================================================================*/

#include "header.h"

#ifdef PML
#define TIMES_PML_BETA_X * pml_beta_x
#define TIMES_PML_BETA_Y * pml_beta_y
#define TIMES_PML_BETA_Z * pml_beta_z
#else
#define TIMES_PML_BETA_X 
#define TIMES_PML_BETA_Y 
#define TIMES_PML_BETA_Z 
#endif


void allocWave( GRID grid, WAVE * h_W, WAVE * W, WAVE * t_W, WAVE * m_W )
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;
	
	//printf( "_nx_ = %d, _ny_ = %d, _nz_ = %d\n", _nx_, _ny_, _nz_  );
	long long num = _nx_ * _ny_ * _nz_; 
		
	float * pWave = NULL;
	long long size = sizeof( float ) * num * WAVESIZE * 4;

	CHECK( Malloc( ( void ** )&pWave, size ) );
	CHECK( Memset(  pWave, 0, size ) ); 

	h_W->Vx  = pWave + 0 * num;
	h_W->Vy  = pWave + 1 * num;
	h_W->Vz  = pWave + 2 * num;
	h_W->Txx = pWave + 3 * num;
	h_W->Tyy = pWave + 4 * num;
	h_W->Tzz = pWave + 5 * num;
	h_W->Txy = pWave + 6 * num;
	h_W->Txz = pWave + 7 * num;
	h_W->Tyz = pWave + 8 * num;

	pWave 	 = pWave + 9 * num;
	

	W->Vx  = pWave + 0 * num;
	W->Vy  = pWave + 1 * num;
	W->Vz  = pWave + 2 * num;
	W->Txx = pWave + 3 * num;
	W->Tyy = pWave + 4 * num;
	W->Tzz = pWave + 5 * num;
	W->Txy = pWave + 6 * num;
	W->Txz = pWave + 7 * num;
	W->Tyz = pWave + 8 * num;

	pWave    = pWave + 9 * num;

	t_W->Vx  = pWave + 0 * num;
	t_W->Vy  = pWave + 1 * num;
	t_W->Vz  = pWave + 2 * num;
	t_W->Txx = pWave + 3 * num;
	t_W->Tyy = pWave + 4 * num;
	t_W->Tzz = pWave + 5 * num;
	t_W->Txy = pWave + 6 * num;
	t_W->Txz = pWave + 7 * num;
	t_W->Tyz = pWave + 8 * num;

	pWave 	 = pWave + 9 * num;

	m_W->Vx  = pWave + 0 * num;
	m_W->Vy  = pWave + 1 * num;
	m_W->Vz  = pWave + 2 * num;
	m_W->Txx = pWave + 3 * num;
	m_W->Tyy = pWave + 4 * num;
	m_W->Tzz = pWave + 5 * num;
	m_W->Txy = pWave + 6 * num;
	m_W->Txz = pWave + 7 * num;
	m_W->Tyz = pWave + 8 * num;



//	printf( "Alloc: h_W = %p,", h_W->Vx );
//	printf( " W = %p,",			  W->Vx );
//	printf( " t_W = %p,",		t_W->Vx );
//	printf( " m_W = %p\n",		m_W->Vx );

}


void freeWave( WAVE h_W, WAVE W, WAVE t_W, WAVE m_W )
{

	Free( h_W.Vx );
}


#define gama1  0.1
#define gama2  0.6
#define gama3  0.3

#define GET_VAR_I( var, VAR )  for ( n = -3; n < 3; n ++ ) { var[n + 3] = VAR[INDEX(i + n, j, k)]; } 
#define GET_VAR_J( var, VAR )  for ( n = -3; n < 3; n ++ ) { var[n + 3] = VAR[INDEX(i, j + n, k)]; } 
#define GET_VAR_K( var, VAR )  for ( n = -3; n < 3; n ++ ) { var[n + 3] = VAR[INDEX(i, j, k + n)]; } 

/*
 0  1  2  3  4  5 
				
-3 -2 -1  0  1  2
*/
#define WENOL( f ) \
u1 = U[J - 2] - 2.0 * U[J - 1] + U[J]; \
u2 = U[J - 2] - 4.0 * U[J - 1] + 3 * U[J]; \
Beta1 = 13.0 / 12.0 * u1 * u1 + 1.0 / 4.0 * u2 * u2; \
u1 = U[J - 1] - 2.0 * U[J] + U[J + 1]; \
u2 = U[J - 1] - U[J + 1]; \
Beta2 = 13.0 / 12.0 * u1 * u1 + 1.0 / 4.0 * u2 * u2; \
u1 = U[J] - 2.0 * U[J + 1] + U[J + 2]; \
u2 = 3 * U[J] - 4.0 * U[J + 1] + U[J + 2]; \
Beta3 = 13.0 / 12.0 * u1 * u1 + 1.0 / 4.0 * u2 * u2; \
_w1 = gama1 / ( ( 1e-6 + Beta1 ) * ( 1e-6 + Beta1 ) ); \
_w2 = gama2 / ( ( 1e-6 + Beta2 ) * ( 1e-6 + Beta2 ) ); \
_w3 = gama3 / ( ( 1e-6 + Beta3 ) * ( 1e-6 + Beta3 ) ); \
w1 = _w1 / ( _w1 + _w2 + _w3 ); \
w2 = _w2 / ( _w1 + _w2 + _w3 ); \
w3 = _w3 / ( _w1 + _w2 + _w3 ); \
f1 =   2.0  * U[J - 2] - 7.0 * U[J - 1] + 11.0 * U[J]; \
f2 = - 1.0  * U[J - 1] + 5.0 * U[J] + 2.0 * U[J + 1];  \
f3 =   2.0  * U[J] + 5.0 * U[J + 1] - 1.0 * U[J + 2];  \
f = w1 * f1 + w2 * f2 + w3 * f3; \
f /= 6.0;

#define WENOR( f ) \
u1 = U[J - 2] - 2.0 * U[J - 1] + U[J]; \
u2 = U[J - 2] - 4.0 * U[J - 1] + 3 * U[J]; \
Beta1 = 13.0 / 12.0 * u1 * u1 + 1.0 / 4.0 * u2 * u2; \
u1 = U[J - 1] - 2.0 * U[J] + U[J + 1]; \
u2 = U[J - 1] - U[J + 1]; \
Beta2 = 13.0 / 12.0 * u1 * u1 + 1.0 / 4.0 * u2 * u2; \
u1 = U[J] - 2.0 * U[J + 1] + U[J + 2]; \
u2 = 3 * U[J] - 4.0 * U[J + 1] + U[J + 2]; \
Beta3 = 13.0 / 12.0 * u1 * u1 + 1.0 / 4.0 * u2 * u2; \
_w1 = gama1 / ( ( 1e-6 + Beta1 ) * ( 1e-6 + Beta1 ) ); \
_w2 = gama2 / ( ( 1e-6 + Beta2 ) * ( 1e-6 + Beta2 ) ); \
_w3 = gama3 / ( ( 1e-6 + Beta3 ) * ( 1e-6 + Beta3 ) ); \
w1 = _w1 / ( _w1 + _w2 + _w3 ); \
w2 = _w2 / ( _w1 + _w2 + _w3 ); \
w3 = _w3 / ( _w1 + _w2 + _w3 ); \
f1 = - 2.0  * U[J - 2] + 5.0 * U[J - 1] + 2.0 * U[J]; \
f2 =   2.0  * U[J - 1] + 5.0 * U[J] - 2.0 * U[J + 1];  \
f3 =  11.0  * U[J] - 7.0 * U[J + 1] + 2.0 * U[J + 2];  \
f = w1 * f1 + w2 * f2 + w3 * f3; \
f /= 6.0;


__GLOBAL__															
void wave_deriv( WAVE h_W, WAVE W, CONTRAVARIANT contravariant, MEDIUM medium, 
#ifdef PML
	PML_BETA pml_beta,
#endif
	int _nx_, int _ny_, int _nz_, float rDH, int FB1, int FB2, int FB3 )	
{																	
	int _nx = _nx_ - HALO;
	int _ny = _ny_ - HALO;
	int _nz = _nz_ - HALO;

#ifdef GPU_CUDA
	int i = threadIdx.x + blockIdx.x * blockDim.x + HALO; 
	int j = threadIdx.y + blockIdx.y * blockDim.y + HALO;
	int k = threadIdx.z + blockIdx.z * blockDim.z + HALO;
#else
	int i = HALO;
	int j = HALO;
	int k = HALO;
#endif
	long long index;

#ifdef PML
	float pml_beta_x;
	float pml_beta_y;
	float pml_beta_z;
#endif
																	
	float mu;												
	float lambda;											
	float buoyancy;											
																	
	float xi_x;   float xi_y; 	float xi_z;		
	float et_x;   float et_y; 	float et_z;		
	float zt_x;   float zt_y; 	float zt_z;		
																	
	float Txx_xi; float Tyy_xi; float Txy_xi; 	
	float Txx_et; float Tyy_et; float Txy_et; 	
	float Txx_zt; float Tyy_zt; float Txy_zt; 	
	float Txz_xi; float Tyz_xi; float Tzz_xi;	
	float Txz_et; float Tyz_et; float Tzz_et;	
	float Txz_zt; float Tyz_zt; float Tzz_zt;	
	float Vx_xi ; float Vx_et ; float Vx_zt ;	
	float Vy_xi ; float Vy_et ; float Vy_zt ;	
	float Vz_xi ; float Vz_et ; float Vz_zt ;
	float Vx1;	  float Vx2;	float Vx3;
	float Vy1;	  float Vy2;	float Vy3;
	float Vz1;	  float Vz2;	float Vz3;
	float Txx1;	  float Txx2;   float Txx3;
	float Tyy1;	  float Tyy2;	float Tyy3;
	float Tzz1;   float Tzz2;	float Tzz3;
	float Txy1;   float Txy2;	float Txy3;
	float Txz1;   float Txz2;	float Txz3;
	float Tyz1;   float Tyz2;	float Tyz3;


	float u1;	float u2;
	float Beta1; float Beta2; float Beta3;
	float _w1; float _w2; float _w3;
	float w1; float w2; float w3;
	float f1; float f2; float f3;
	float FL1, FL2, FR1, FR2;

	float U[5]; int n, J;

	int i0, j0, k0;

	CALCULATE3D( i, j, k, HALO, _nx, HALO, _ny, HALO, _nz )  
		index = INDEX( i, j, k ); 									
		mu = medium.mu[index]; 											
		lambda = medium.lambda[index];										
		buoyancy = medium.buoyancy[index];

#ifdef PML
		pml_beta_x = pml_beta.x[i];
		pml_beta_y = pml_beta.y[j];
		pml_beta_z = pml_beta.z[k];
#endif

		xi_x = contravariant.xi_x[index]; xi_y = contravariant.xi_y[index]; xi_z = contravariant.xi_z[index];	
		et_x = contravariant.et_x[index]; et_y = contravariant.et_y[index]; et_z = contravariant.et_z[index];	
		zt_x = contravariant.zt_x[index]; zt_y = contravariant.zt_y[index]; zt_z = contravariant.zt_z[index];	

		if ( i == _nx_ / 2 && j == _ny_ / 2 && k == _nz_ / 2 )
		{
			for ( n = -3; n < 3; n ++ ) { 
				U[n + 3] = W.Txx[INDEX(i + n, j, k)]; 
				printf( "U[%d] = %f\n", n + 3, U[n+3]);
			}
		}
		

 
		GET_VAR_I( U, W.Vx  ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 );  Vx_xi = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_I( U, W.Vy  ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 );  Vy_xi = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_I( U, W.Vz  ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 );  Vz_xi = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_I( U, W.Txx ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 ); Txx_xi = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_I( U, W.Tyy ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 ); Tyy_xi = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_I( U, W.Tzz ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 ); Tzz_xi = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_I( U, W.Txy ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 ); Txy_xi = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_I( U, W.Txz ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 ); Txz_xi = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_I( U, W.Tyz ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 ); Tyz_xi = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
                                                                                                                                        
		GET_VAR_J( U, W.Vx  ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 );  Vx_et = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_J( U, W.Vy  ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 );  Vy_et = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_J( U, W.Vz  ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 );  Vz_et = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_J( U, W.Txx ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 ); Txx_et = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_J( U, W.Tyy ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 ); Tyy_et = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_J( U, W.Tzz ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 ); Tzz_et = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_J( U, W.Txy ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 ); Txy_et = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_J( U, W.Txz ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 ); Txz_et = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_J( U, W.Tyz ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 ); Tyz_et = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
                                                                                                                                        
		GET_VAR_K( U, W.Vx  ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 );  Vx_zt = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_K( U, W.Vy  ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 );  Vy_zt = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_K( U, W.Vz  ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 );  Vz_zt = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_K( U, W.Txx ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 ); Txx_zt = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_K( U, W.Tyy ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 ); Tyy_zt = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_K( U, W.Tzz ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 ); Tzz_zt = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_K( U, W.Txy ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 ); Txy_zt = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_K( U, W.Txz ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 ); Txz_zt = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;
		GET_VAR_K( U, W.Tyz ); J = 3; WENOL( FL1 ); J = 2; WENOR( FL2 ); J = 3; WENOL( FR2 ); J = 2; WENOR( FR1 ); Tyz_zt = ( FL1 - FL2 + FR1 - FR2 ) * rDH TIMES_PML_BETA_X;

		Vx1  = DOT_PRODUCT3D( xi_x, xi_y, xi_z, Txx_xi, Txy_xi, Txz_xi ) * buoyancy;
		Vx2  = DOT_PRODUCT3D( et_x, et_y, et_z, Txx_et, Txy_et, Txz_et ) * buoyancy;
		Vx3  = DOT_PRODUCT3D( zt_x, zt_y, zt_z, Txx_zt, Txy_zt, Txz_zt ) * buoyancy;
		Vy1  = DOT_PRODUCT3D( xi_x, xi_y, xi_z, Txy_xi, Tyy_xi, Tyz_xi ) * buoyancy;
		Vy2  = DOT_PRODUCT3D( et_x, et_y, et_z, Txy_et, Tyy_et, Tyz_et ) * buoyancy;
		Vy3  = DOT_PRODUCT3D( zt_x, zt_y, zt_z, Txy_zt, Tyy_zt, Tyz_zt ) * buoyancy;
		Vz1  = DOT_PRODUCT3D( xi_x, xi_y, xi_z, Txz_xi, Tyz_xi, Tzz_xi ) * buoyancy;
		Vz2  = DOT_PRODUCT3D( et_x, et_y, et_z, Txz_et, Tyz_et, Tzz_et ) * buoyancy;
		Vz3  = DOT_PRODUCT3D( zt_x, zt_y, zt_z, Txz_zt, Tyz_zt, Tzz_zt ) * buoyancy;

		Txx1 = DOT_PRODUCT3D( xi_x, xi_y, xi_z, Vx_xi, Vy_xi, Vz_xi ) * lambda + 2.0f * mu * ( xi_x * Vx_xi );
		Txx2 = DOT_PRODUCT3D( et_x, et_y, et_z, Vx_et, Vy_et, Vz_et ) * lambda + 2.0f * mu * ( et_x * Vx_et );
		Txx3 = DOT_PRODUCT3D( zt_x, zt_y, zt_z, Vx_zt, Vy_zt, Vz_zt ) * lambda + 2.0f * mu * ( zt_x * Vx_zt );
		Tyy1 = DOT_PRODUCT3D( xi_x, xi_y, xi_z, Vx_xi, Vy_xi, Vz_xi ) * lambda + 2.0f * mu * ( xi_y * Vy_xi );
		Tyy2 = DOT_PRODUCT3D( et_x, et_y, et_z, Vx_et, Vy_et, Vz_et ) * lambda + 2.0f * mu * ( et_y * Vy_et );
		Tyy3 = DOT_PRODUCT3D( zt_x, zt_y, zt_z, Vx_zt, Vy_zt, Vz_zt ) * lambda + 2.0f * mu * ( zt_y * Vy_zt );
		Tzz1 = DOT_PRODUCT3D( xi_x, xi_y, xi_z, Vx_xi, Vy_xi, Vz_xi ) * lambda + 2.0f * mu * ( xi_z * Vz_xi );
		Tzz2 = DOT_PRODUCT3D( et_x, et_y, et_z, Vx_et, Vy_et, Vz_et ) * lambda + 2.0f * mu * ( et_z * Vz_et );
		Tzz3 = DOT_PRODUCT3D( zt_x, zt_y, zt_z, Vx_zt, Vy_zt, Vz_zt ) * lambda + 2.0f * mu * ( zt_z * Vz_zt );

		Txy1 = DOT_PRODUCT2D( xi_y, xi_x, Vx_xi, Vy_xi ) * mu;
		Txy2 = DOT_PRODUCT2D( et_y, et_x, Vx_et, Vy_et ) * mu;
		Txy3 = DOT_PRODUCT2D( zt_y, zt_x, Vx_zt, Vy_zt ) * mu;
		Txz1 = DOT_PRODUCT2D( xi_z, xi_x, Vx_xi, Vz_xi ) * mu;
		Txz2 = DOT_PRODUCT2D( et_z, et_x, Vx_et, Vz_et ) * mu;
		Txz3 = DOT_PRODUCT2D( zt_z, zt_x, Vx_zt, Vz_zt ) * mu;
		Tyz1 = DOT_PRODUCT2D( xi_z, xi_y, Vy_xi, Vz_xi ) * mu;
		Tyz2 = DOT_PRODUCT2D( et_z, et_y, Vy_et, Vz_et ) * mu;
		Tyz3 = DOT_PRODUCT2D( zt_z, zt_y, Vy_zt, Vz_zt ) * mu;

		h_W.Vx [index] 	= Vx1  + Vx2  + Vx3 ;							
		h_W.Vy [index] 	= Vy1  + Vy2  + Vy3 ;							
		h_W.Vz [index] 	= Vz1  + Vz2  + Vz3 ;							
		h_W.Txx[index] 	= Txx1 + Txx2 + Txx3;						
		h_W.Tyy[index] 	= Tyy1 + Tyy2 + Tyy3;						
		h_W.Tzz[index] 	= Tzz1 + Tzz2 + Tzz3;						
		h_W.Txy[index] 	= Txy1 + Txy2 + Txy3;						
		h_W.Txz[index] 	= Txz1 + Txz2 + Txz3;						
		h_W.Tyz[index] 	= Tyz1 + Tyz2 + Tyz3;						
																				
																	
	END_CALCULATE3D( )												
}			


void waveDeriv( GRID grid, WAVE h_W, WAVE W, CONTRAVARIANT contravariant, MEDIUM medium, 
#ifdef PML
	PML_BETA pml_beta,
#endif
	int FB1, int FB2, int FB3 )	
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	float rDH = grid.rDH;

	//printf( "Deriv: h_W = %p,", h_W.Vx );
	//printf( " W = %p\n",	  W.Vx );
	//printf( " t_W = %p,",	t_W );
	//printf( " m_W = %p\n",	m_W );

#ifdef GPU_CUDA

	int nx = _nx_ - 2 * HALO;
	int ny = _ny_ - 2 * HALO;
	int nz = _nz_ - 2 * HALO;

	dim3 threads( 32, 4, 4);
	dim3 blocks;
	blocks.x = ( nx + threads.x - 1 ) / threads.x;
	blocks.y = ( ny + threads.y - 1 ) / threads.y;
	blocks.z = ( nz + threads.z - 1 ) / threads.z;

	wave_deriv <<< blocks, threads >>>
	( h_W, W, contravariant, medium, 
#ifdef PML
	pml_beta,
#endif  //PML
	_nx_, _ny_, _nz_, rDH, FB1, FB2, FB3 );	

	CHECK( cudaDeviceSynchronize( ) );

#else   //GPU_CUDA

	wave_deriv 
	( h_W, W, contravariant, medium, 
#ifdef PML
	pml_beta,
#endif //PML
	_nx_, _ny_, _nz_, rDH, FB1, FB2, FB3);	


#endif //GPU_CUDA
}




