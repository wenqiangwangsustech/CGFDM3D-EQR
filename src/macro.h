/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:macro.h
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-04
*   Discription:
*
================================================================*/
#ifndef __MACRO__
#define __MACRO__

#define HALO 3
#define MIN( x, y ) ( ( x ) < ( y ) ? ( x ) : ( y )  )
#define MAX( x, y ) ( ( x ) > ( y ) ? ( x ) : ( y )  )

#define WAVESIZE 9 //9 Wave components: Vx Vy Vz Txx Tyy Tzz Txy Txz Tyz
#define COORDSIZE 3 //3 coordinate components: coordX coordY coordZ
#define CONTRASIZE 9 //contravariant components
#define MEDIUMSIZE 3 //3 medium components: Vs Vp rho ( lam mu bouyancy )
#define MOMENTSIZE 6


#define PI 3.141592657f
//#define PI 3.1415926535898
#define RADIAN2DEGREE ( 180.0 / PI )
#define DEGREE2RADIAN ( PI / 180.0 )

#define Cv 1.0
#define Cs 1.0


//#define X_FAST

#define FOR_LOOP2D( i, j, startI, endI, startJ, endJ ) \
for ( j = startJ; j < endJ; j ++ ) {    \
for ( i = startI; i < endI; i ++ ) {    \


#define END_LOOP2D( ) }}

#define FOR_LOOP3D( i, j, k, startI, endI, startJ, endJ, startK, endK ) \
for ( k = startK; k < endK; k ++ ) {     \
for ( j = startJ; j < endJ; j ++ ) {     \
for ( i = startI; i < endI; i ++ ) {     \

#define END_LOOP3D( ) }}}
#define LON_FAST



#ifdef GPU_CUDA

#define Malloc cudaMalloc
#define Memset cudaMemset
#define Free cudaFree
#define Memcpy cudaMemcpy


#define __DEVICE__ __device__
#define __GLOBAL__ __global__



#define CALCULATE1D( i, startI, endI )         \
if ( i >= startI && i < endI ) {
#define END_CALCULATE1D( ) } 



#define CALCULATE2D( i, j, startI, endI, startJ, endJ )         \
if ( i >= ( startI ) && i < ( endI ) && j >= ( startJ ) && j < ( endJ ) ) { 
#define END_CALCULATE2D( ) }      


#define CALCULATE3D( i, j, k, startI, endI, startJ, endJ, startK, endK )  \
if ( i >= ( startI ) && i < ( endI ) && j >= ( startJ ) && j < ( endJ ) && k >= ( startK ) && k < ( endK ) ) {
#define END_CALCULATE3D( ) }





#define POW2(x) ( (x) * (x) )
#define GAUSS_FUN(t,a,t0) (exp(-POW2( ( (t) - (t0) ) / (a) )) / (a*1.772453850905516))



#define CHECK(call) {                                  \
  const cudaError_t error = call;                          \
  if (error != cudaSuccess) {                              \
    fprintf(stderr, "Error: %s:%d, ", __FILE__, __LINE__); \
    fprintf(stderr, "code: %d, reason: %s\n",              \
        error, cudaGetErrorString(error));                 \
  }                                                        \
}

#else
int Malloc( void ** mem, long long size  );

#define Memset memset
#define Free free
#define Memcpy memcpy
#define __DEVICE__     
#define __GLOBAL__     

#define CHECK( call ) call

#define CALCULATE2D( i, j, startI, endI, startJ, endJ ) \
for ( j = startJ; j < endJ; j ++ ) {    \
for ( i = startI; i < endI; i ++ ) {    

#define END_CALCULATE2D( ) }}

#define CALCULATE3D( i, j, k, startI, endI, startJ, endJ, startK, endK ) \
for ( k = startK; k < endK; k ++ ) {     \
for ( j = startJ; j < endJ; j ++ ) {     \
for ( i = startI; i < endI; i ++ ) {     

#define END_CALCULATE3D( ) }}}

#define CALCULATE1D( i, startI, endI )         \
    for ( i = ( startI ); i < endI; i ++ ) {
#define END_CALCULATE1D( ) } 


#endif




#define THREE 3

//forward difference coefficient
#define af_1 (-0.30874f)
#define af0  (-0.6326f )
#define af1  ( 1.2330f )
#define af2  (-0.3334f )
#define af3  ( 0.04168f)
//backward difference coefficient
#define ab_1 ( 0.30874f)
#define ab0  ( 0.6326f )
#define ab1  (-1.2330f )
#define ab2  ( 0.3334f )
#define ab3  (-0.04168f)

#define alpha1 0.0f
#define alpha2 0.5f
#define alpha3 0.5f
#define alpha4 1.0f

//#define beta1 1.0f//0.16666667f
#define beta1 0.16666667f
#define beta2 0.33333333f
#define beta3 0.33333333f
#define beta4 0.16666667f

#define Cf1 ( - 1.16666667f )
#define Cf2 (   1.33333333f )
#define Cf3 ( - 0.16666667f )

#define Cb1 (   1.16666667f )
#define Cb2 ( - 1.33333333f )
#define Cb3 (   0.16666667f )


// #define Cf1 ( - 7.0 / 6.0 )
// #define Cf2 (   4.0 / 3.0 )
// #define Cf3 ( - 1.0 / 6.0 )

// #define Cb1 (   7.0 / 6.0 )
// #define Cb2 ( - 4.0 / 3.0 )
// #define Cb3 (   1.0 / 6.0 )



//grid: x, y, z
#define Index3D( i, j, k, nx, ny, nz ) ( ( i ) + ( j ) * ( nx ) + ( k ) * ( ( nx ) * ( ny ) ) )
#define Index2D( i, j, nx, ny ) ( ( i ) + ( j ) * ( nx ) )


#define INDEX( i, j, k ) ( ( i ) + ( j ) * _nx_ + ( k ) * ( _nx_ * _ny_ ) )

//generate index of adjacent point: Up/Down/Left/Right

#define INDEX_xi( i, j, k, offset ) ( ( i  offset ) + ( j ) * _nx_ + ( k ) * ( _nx_ * _ny_ ) )
#define INDEX_et( i, j, k, offset ) ( ( i ) + ( j  offset ) * _nx_ + ( k ) * ( _nx_ * _ny_ ) )
#define INDEX_zt( i, j, k, offset ) ( ( i ) + ( j ) * _nx_ + ( k  offset ) * ( _nx_ * _ny_ ) )


#define L( W, FB, SUB ) (  FB * ( af_1 * W[INDEX_##SUB( i, j, k, - FB * 1 )] + 				\
						 		   af0 * W[index] + 										\
						 		   af1 * W[INDEX_##SUB( i, j, k, + FB * 1 )] + 				\
						 		   af2 * W[INDEX_##SUB( i, j, k, + FB * 2 )] + 				\
						 		   af3 * W[INDEX_##SUB( i, j, k, + FB * 3 )] ) * rDH )

#define L_J_T( J_T, FB ) ( FB * ( af_1 * J_T[3 - FB * 1] + \
					 			  af0  * J_T[3] + 		   \
					 			  af1  * J_T[3 + FB * 1] + \
					 			  af2  * J_T[3 + FB * 2] + \
					 			  af3  * J_T[3 + FB * 3] ) * rDH )

#define L2( W, FB, SUB ) ( FB * ( W[INDEX_##SUB( i, j, k, + FB * 1 )] - W[index] ) * rDH )
#define L3( W, FB, SUB ) ( FB * ( Cf1 * W[index] + Cf2 * W[INDEX_##SUB( i, j, k, + FB * 1 )] + Cf3 * W[INDEX_##SUB( i, j, k, + FB * 2 )] ) * rDH ) 


//#define DOT_PRODUCT( A1, A2, A3, B1, B2, B3 ) ( ( A1 ) * ( B1 ) + ( A2 ) * ( B2 ) + ( A3 ) * ( B3 ) )
#define DOT_PRODUCT3D( A1, A2, A3, B1, B2, B3 ) ( ( A1 ) * ( B1 ) + ( A2 ) * ( B2 ) + ( A3 ) * ( B3 ) )
#define DOT_PRODUCT2D( A1, A2, B1, B2 ) ( ( A1 ) * ( B1 ) + ( A2 ) * ( B2 ) )

#endif //__MACRO__
