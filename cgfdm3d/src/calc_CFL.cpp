/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:calc_CFL.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-01
*   Discription:
*
================================================================*/


#include "header.h"

__DEVICE__
double distance_point2plane( double P[3], double A[3], double B[3], double C[3] ) 
{
	double AB[3] = {
		B[0] - A[0],
		B[1] - A[1],
		B[2] - A[2]
	};

	double AC[3] = {
		C[0] - A[0],
		C[1] - A[1],
		C[2] - A[2]
	};

	//double BC[3] = {
	//	C[0] - B[0],
	//	C[1] - B[1],
	//	C[2] - B[2]
	//};

	
	double n [3] = {
		AB[1] * AC[2] - AB[2] * AC[1],
		AB[2] * AC[0] - AB[0] * AC[2],
		AB[0] * AC[1] - AB[1] * AC[0]
	};


	double PA[3] = { 
		P[0] - A[0],
		P[1] - A[1],
		P[2] - A[2]
	};

	double n_dis = sqrt( n[0] * n[0] + n[1] * n[1] + n[2] * n[2] );

	n[0] = n[0] / n_dis;
	n[1] = n[1] / n_dis;
	n[2] = n[2] / n_dis;

	 
	double d = abs( PA[0] * n[0] + PA[1] * n[1] + PA[2] * n[2] );
	return d;
}

__GLOBAL__
void calculate_DH_Range( 
	COORD coord, 
#ifdef GPU_CUDA
	float * DHRange,
#else
	float h_min_max[2], 
#endif
	int _nx_, int _ny_, int _nz_ )
{
#ifdef GPU_CUDA
	int i = blockIdx.x * blockDim.x + threadIdx.x + HALO;
	int j = blockIdx.y * blockDim.y + threadIdx.y + HALO;
	int k = blockIdx.z * blockDim.z + threadIdx.z + HALO;
#else
	int i = HALO;
	int j = HALO;
	int k = HALO;
#endif
	long long index = 0;

	int _nx = _nx_ - HALO;
	int _ny = _ny_ - HALO;
	int _nz = _nz_ - HALO;

#ifdef GPU_CUDA
	int nx = _nx - HALO;
	int ny = _ny - HALO;
	int nz = _nz - HALO;
#endif

	int ii, jj, kk;

	double P[3], A[3], B[3], C[3];
	long long posA = 0;
	long long posB = 0;
	long long posC = 0;

	double d = 0.0;
	double hmin = 1.0e20;
	double hmax = 0.0;



	CALCULATE3D( i, j, k, HALO, _nx, HALO, _ny, HALO, _nz )
		index = INDEX( i, j, k );

		P[0] = coord.x[index];
		P[1] = coord.y[index];
		P[2] = coord.z[index];

		for ( ii = - 1; ii <= 1; ii += 2 ) 
			for ( jj = - 1; jj <= 1; jj += 2 ) 
				for( kk = - 1; kk <= 1; kk += 2 ) 
				{
					posA = INDEX( i - ii, j, k );//前后
					posB = INDEX( i, j - jj, k );//左右
					posC = INDEX( i, j, k - kk );//上下

					A[0] = coord.x[posA];
					A[1] = coord.y[posA];
					A[2] = coord.z[posA];

					B[0] = coord.x[posB];
					B[1] = coord.y[posB];
					B[2] = coord.z[posB];

					C[0] = coord.x[posC];
					C[1] = coord.y[posC];
					C[2] = coord.z[posC];


					d = distance_point2plane( P, A, B, C );
					//d = 100.0f / sqrt( 3.0 );
				 	hmin = MIN( hmin, d );
				 	hmax = MAX( hmax, d );
				}
#ifdef GPU_CUDA
		index = Index3D( i - HALO, j - HALO, k - HALO, nx, ny, nz );	
		//printf( "hmax = %f\n", hmax );
		DHRange[index] = hmin;
		DHRange[index + nx * ny * nz] = hmax;

#endif

	END_CALCULATE3D( )

#ifndef GPU_CUDA
	h_min_max[0] = hmin;
	h_min_max[1] = hmax;
#endif

}


#ifdef GPU_CUDA
void allocDHRange( float ** DHRange, int nx, int ny, int nz )
{
	long long size = sizeof( float ) * nx * ny * nz * 2;

	CHECK( Malloc( DHRange, size ) );
	
}


void freeDHRange( float * DHRange )
{

	Free( DHRange );
}

//The integer number can not represent array size. 
//Especially, the type of last parameter of cublasIsamin/cublasIsamax is int(interger). 
//Then, we must segment the array(data) into numBlocks + 1 parts;
//cublasIsamin/cublasIsamax solve only absolute max/min value.
//
void calculate_range( float * data, float range[2], long long num )
{
    cublasHandle_t handle;
    cublasStatus_t handleStat, stat;
    handleStat = cublasCreate( &handle );

		//system( "sleep 100s"  );
	if ( handleStat != CUBLAS_STATUS_SUCCESS )
	{
		printf(  "CUBLAS initialization failed\n"  );
	}
	

    long long indexMin, indexMax;
	int t_index;
	float maxValue = -1.0e20, minValue = 1.0e20;
	float tmpMax = -1.0e20, tmpMin = 1.0e20;

	int blockSize = 1024 * 1024 * 1024;

	int	numBlocks = num / blockSize + 1;
	int retSize   = num % blockSize;

	//printf( "numBlocks = %d, retSize = %d\n", numBlocks, retSize );
	
	int cnt;

	int i = 0;
	for ( i = 0; i < numBlocks; i ++ )
	{

		long long stride = i * blockSize;
		cnt = blockSize;
		if ( i == numBlocks - 1 )
			cnt = retSize;
	
		//printf( "cnt = %d\n", cnt  );

		stat = cublasIsamin( handle, cnt, data + stride, 1, &t_index );
		//printf( "t_index = %d\n", t_index  );
		indexMin = stride + t_index - 1;
		cudaMemcpy( &tmpMin, &data[indexMin], sizeof( float ), cudaMemcpyDeviceToHost  );
		//printf( "tmpMin = %f\n",  tmpMin );
		minValue = MIN( minValue, tmpMin );

		if ( stat != CUBLAS_STATUS_SUCCESS )
		    printf("CUDA find Max failed!\n");
		
		stat = cublasIsamax( handle, cnt, data + stride, 1, &t_index );
		indexMax = stride + t_index - 1;
		cudaMemcpy( &tmpMax, &data[indexMax], sizeof( float ), cudaMemcpyDeviceToHost  );
		//printf( "tmpMax = %f\n", tmpMax  );
		maxValue = MAX( maxValue, tmpMax );

		if ( stat != CUBLAS_STATUS_SUCCESS )
		    printf("CUDA find Max failed!\n");
	}
//	printf( "num1 = %ld, num2= %ld\n", num, index + retSize  );

	range[0] = minValue;
	range[1] = maxValue;
    
	if ( handleStat == CUBLAS_STATUS_SUCCESS ) 
		cublasDestroy( handle );
	//printf( "indexMin = %ld, indexMax = %ld\n", indexMin, indexMax );
	//printf( "min = %7.f, max = %7.f\n", range[0], range[1] );
	
	
}
#else
void calculate_range( float * data, float range[2], long long num )
{
	long long i = 0;
	float maxValue = -1.0e20, minValue = 1.0e20;
	for ( i = 0; i < num; i ++ )
	{
			
		if ( data[i] < minValue )
			minValue = data[i];

		if ( data[i] > maxValue )
			maxValue = data[i];

	}
	range[0] = minValue;
	range[1] = maxValue;

}
	

#endif

void calc_CFL( GRID grid, COORD coord, MEDIUM medium, PARAMS params )
{

	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;

	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	STRUCTURE structure = { 0 };

	structure.Vs  = medium.mu;
	structure.Vp  = medium.lambda;
	structure.rho = medium.buoyancy;

	//printf( "Vs: %p, Vp: %p, rho: %p\n", structure.Vs, structure.Vp, structure.rho  );
	

	float h_min_max[2];

	float Vs_min_max[2];
	float Vp_min_max[2];
	float rho_min_max[2];



	long long num = nx * ny * nz;

#ifdef GPU_CUDA
	float * DHRange;
	allocDHRange( &DHRange, nx, ny, nz );
#endif
	
#ifdef GPU_CUDA
	dim3 threads( 32, 4, 4);
	dim3 blocks;
	
	blocks.x = ( nx + threads.x - 1 ) / threads.x;
	blocks.y = ( ny + threads.y - 1 ) / threads.y;
	blocks.z = ( nz + threads.z - 1 ) / threads.z;

	calculate_DH_Range
	<<< blocks, threads >>>
	( coord, DHRange, _nx_, _ny_, _nz_ );

	calculate_range( DHRange, h_min_max, num * 2 );



	//MPI_Barrier( MPI_COMM_WORLD );
	//calculate_range( coord.z, h_min_max, _nx_ * _ny_ * _nz_ );

	//float coordZ;
	//cudaMemcpy( &coordZ, coord.z + 207684664, sizeof( float ), cudaMemcpyDeviceToHost  );
	//printf( "2: coordZ = %f\n", coordZ );
	//calculate_range( coord.z + 207684664, h_min_max, 1 );


	freeDHRange( DHRange );
	

#else
	calculate_DH_Range( coord, h_min_max, _nx_, _ny_, _nz_ );
#endif

	//printf( "min = %f, max = %f\n", h_min_max[0], h_min_max[1] );

	calculate_range( structure.Vs, Vs_min_max, num );
	calculate_range( structure.Vp, Vp_min_max, num );
	calculate_range( structure.rho, rho_min_max, num );



	MPI_Barrier( MPI_COMM_WORLD );

	int thisRank;
	MPI_Comm_rank( MPI_COMM_WORLD, &thisRank );
	//printf( "thisRank = %d, H   Range: %5.f ~ %5.f\n", thisRank, h_min_max[0], h_min_max[1] );
	//printf( "thisRank = %d, Vs  Range: %5.f ~ %5.f\n", thisRank, Vs_min_max[0], Vs_min_max[1] );
	//printf( "thisRank = %d, Vp  Range: %5.f ~ %5.f\n", thisRank, Vp_min_max[0], Vp_min_max[1] );
	//printf( "thisRank = %d, rho Range: %5.f ~ %5.f\n", thisRank, rho_min_max[0], rho_min_max[1] );

	float H_min,   H_max;
	float Vs_min,  Vs_max;
	float Vp_min,  Vp_max;
	float rho_min, rho_max;
	

  	MPI_Allreduce( &h_min_max[0],   &H_min,   1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD );
  	MPI_Allreduce( &h_min_max[1],   &H_max,   1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD );
                                              
  	MPI_Allreduce( &Vs_min_max[0],  &Vs_min,  1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD );
  	MPI_Allreduce( &Vs_min_max[1],  &Vs_max,  1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD );
                                              
  	MPI_Allreduce( &Vp_min_max[0],  &Vp_min,  1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD );
  	MPI_Allreduce( &Vp_min_max[1],  &Vp_max,  1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD );
                                              
  	MPI_Allreduce( &rho_min_max[0], &rho_min, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD );
  	MPI_Allreduce( &rho_min_max[1], &rho_max, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD );

	float dtmax = 1.32 * H_min / Vp_max;

	

	if ( 0 == thisRank )
	{
		printf( "H   Range: %5.2e ~ %5.2e\n", H_min,   H_max	 );
		printf( "Vs  Range: %5.2e ~ %5.2e\n", Vs_min,  Vs_max  );
		printf( "Vp  Range: %5.2e ~ %5.2e\n", Vp_min,  Vp_max  );
		printf( "rho Range: %5.2e ~ %5.2e\n", rho_min, rho_max );
		
		printf( "dtmax = %5.2e, DT = %5.2e\n", dtmax, params.DT );
		if ( params.DT > dtmax )
		{
			printf( "The parameters can afford the CFL condition!\n" );
			MPI_Abort( MPI_COMM_WORLD, 110 );
		}
	}

	

	MPI_Barrier( MPI_COMM_WORLD );

}




