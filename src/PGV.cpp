/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:propagate.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-20
*   Discription:
*
================================================================*/
#include "header.h"


void allocatePGV( GRID grid, PGV * pgv )
{
	int nx = grid.nx;
	int ny = grid.ny;


	int len = sizeof( float ) * nx * ny * 2;
	float * pPgv = NULL;

	CHECK( Malloc( ( void ** )&pPgv, len) );
	CHECK( Memset( pPgv, 0, len ));

	int num = nx * ny;

	pgv->pgvh = pPgv;
	pgv->pgv  = pgv->pgvh + num;
	//printf( "pgvh = %p, pgv = %p\n",  pgv->pgvh, pgv->pgv );
}

void freePGV( PGV pgv )
{
	Free( pgv.pgvh );
}


void allocatePGV_cpu( GRID grid, PGV * cpuPgv )
{
	int nx = grid.nx;
	int ny = grid.ny;

	int len = sizeof( float ) * nx * ny * 2;
	float * pCpuPgv = NULL;

	pCpuPgv = ( float * )malloc( len );
	memset( pCpuPgv, 0, len );
	int num = nx * ny;

	cpuPgv->pgvh = pCpuPgv;
	cpuPgv->pgv  = cpuPgv->pgvh + num;
}

void freePGV_cpu( PGV cpuPgv )
{
	free( cpuPgv.pgvh );
}


__GLOBAL__	
void compare_pgv( PGV pgv, WAVE W, int nx, int ny, int nz )							
{
	int _nx_ = nx + 2 * HALO;
	int _ny_ = ny + 2 * HALO;
	int _nz_ = nz + 2 * HALO;

#ifdef GPU_CUDA
	int i0 = threadIdx.x + blockIdx.x * blockDim.x;
	int j0 = threadIdx.y + blockIdx.y * blockDim.y;
#else
	int i0 = 0;
	int j0 = 0;
#endif


	//printf( "==========compare pgv===============\n"  );
	
	int i = i0 + HALO;
	int j = j0 + HALO;
	int k = _nz_ - HALO - 1;

	long long index, pos;

	float Vx = 0.0f, Vy = 0.0f, Vz = 0.0f, Vh = 0.0f, V = 0.0f;

	double c = 1.0 / Cv;
	CALCULATE2D( i0, j0, 0, nx, 0, ny )
		i = i0 + HALO;
		j = j0 + HALO;
		index = INDEX( i, j, k );	
		pos = Index2D( i0, j0, nx, ny );
			
		//if ( i0 == nx - 1 && j0 == ny - 1 )
		//	printf( "nx = %d, ny = %d, pos = %d, pgvh = %p, pgv = %p\n", nx, ny, pos, pgv.pgvh, pgv.pgv );

		Vx = W.Vx[index] * c;
		Vy = W.Vy[index] * c;
		Vz = W.Vz[index] * c;
		
		Vh = sqrtf( Vx * Vx + Vy * Vy );
		V  = sqrtf( Vx * Vx + Vy * Vy + Vz * Vz);
		
		if ( pgv.pgvh[pos] < Vh )
		{
			pgv.pgvh[pos] = Vh;
		}
		if ( pgv.pgv[pos] < V )
		{
			pgv.pgv[pos] = V;
		}
		/*
		*/
	END_CALCULATE2D( )
}

void outputPgvData( PARAMS params, MPI_COORD thisMPICoord, PGV cpuPgv, int nx, int ny )
{
	char fileName1[1024] = { 0 };
	char fileName2[1024] = { 0 };
	FILE * filePgvh = NULL;
	FILE * filePgv  = NULL;

	int i, j;

	sprintf( fileName1, "%s/PGVh_Z_mpi_%d_%d_%d.bin", params.OUT, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
	sprintf( fileName2, "%s/PGV_Z_mpi_%d_%d_%d.bin",  params.OUT, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
	filePgvh = fopen( fileName1, "wb" );
	filePgv  = fopen( fileName2, "wb" );

	//printf( fileName1 );
	//printf( fileName2 );
	float vm = 0.0;

	fwrite( cpuPgv.pgvh, sizeof( float ), nx * ny, filePgvh );
	fwrite( cpuPgv.pgv , sizeof( float ), nx * ny, filePgv  );
	

	fclose( filePgvh );
	fclose( filePgv  );
}



void comparePGV( GRID grid, MPI_COORDINATE thisMPICoord, WAVE W, PGV pgv )
{
	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;
#ifdef GPU_CUDA
	dim3 threads( 32, 16, 1 );
	dim3 blocks;
	blocks.x = ( nx + threads.x - 1 ) / threads.x;
	blocks.y = ( ny + threads.y - 1 ) / threads.y;
	blocks.z = 1;

	compare_pgv <<< blocks, threads >>> ( pgv, W, nx, ny, nz );
	CHECK( cudaDeviceSynchronize( ) );
#else
	compare_pgv ( pgv, W, nx, ny, nz );
#endif 



}

void outputPGV( PARAMS params, GRID grid, MPI_COORDINATE thisMPICoord, PGV pgv, PGV cpuPgv )
{

	int nx = grid.nx;
	int ny = grid.ny;

	int size;
	size = sizeof( float ) * nx * ny * 2;
#ifdef GPU_CUDA
	CHECK( cudaMemcpy( cpuPgv.pgvh, pgv.pgvh, size, cudaMemcpyDeviceToHost ) );
	outputPgvData( params, thisMPICoord, cpuPgv, nx, ny );
#else
	outputPgvData( params, thisMPICoord, pgv, nx, ny );
#endif
}	

