/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:data_io.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-06
*   Discription:
*
================================================================*/
#include "header.h"

void locateSlice( PARAMS params, GRID grid, SLICE * slice )
{
	int sliceX = params.sliceX;
	int sliceY = params.sliceY;
	int sliceZ = params.sliceZ;
	
	int frontNX = grid.frontNX;
	int frontNY = grid.frontNY;
	int frontNZ = grid.frontNZ;
	
 	slice->X = sliceX - frontNX + HALO;
 	slice->Y = sliceY - frontNY + HALO;
 	slice->Z = sliceZ - frontNZ + HALO;
	
	//printf( "slice.X = %d, slice.Y = %d, slice.Z = %d\n", slice->X, slice->Y, slice->Z );
	
}

void locateFreeSurfSlice( GRID grid, SLICE * slice )
{
	int sliceX = -1;
	int sliceY = -1;
	int sliceZ = grid.NZ - 1;
	
	int frontNX = grid.frontNX;
	int frontNY = grid.frontNY;
	int frontNZ = grid.frontNZ;
	
 	slice->X = sliceX - frontNX + HALO;
 	slice->Y = sliceY - frontNY + HALO;
 	slice->Z = sliceZ - frontNZ + HALO;
	
	//printf( "slice.X = %d, slice.Y = %d, slice.Z = %d\n", slice->X, slice->Y, slice->Z );
	
}


void allocDataout( GRID grid, char XYZ, float ** dataout )
{
	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;

	long long num = 0;
	
	switch( XYZ )
	{
		case 'X':
			num = ny * nz;
			break;
		case 'Y':
			num = nx * nz;
			break;
		case 'Z':
			num = nx * ny;
			break;
	}
		
	float * pData = NULL;
	long long size = sizeof( float ) * num;

	CHECK( Malloc( ( void ** )&pData, size ) );
	CHECK( Memset(  pData, 0, size ) ); 
	
	*dataout = pData;
}

void freeDataout( float * dataout  )
{

	Free( dataout  );
}


__GLOBAL__
void pack_iodata_x( float * datain, float * dataout, int nx, int ny, int nz, int I  )
{
	int _nx_ = nx + 2 * HALO;
	int _ny_ = ny + 2 * HALO;
	//int _nz_ = nz + 2 * HALO;

	//printf( "datain = %p\n", datain  );

#ifdef GPU_CUDA
	int j0 = threadIdx.x + blockIdx.x * blockDim.x;
	int k0 = threadIdx.y + blockIdx.y * blockDim.y;
#else
	int j0 = 0;
	int k0 = 0;
#endif
	
	int i = I;
	int j = j0 + HALO;
	int k = k0 + HALO;

	long long index, pos;
	
	CALCULATE2D( j0, k0, 0, ny, 0, nz )
		j = j0 + HALO;
		k = k0 + HALO;
		index = INDEX( i, j, k );	
		pos = Index2D( j0, k0, ny, nz );
		dataout[pos] = datain[index];
		//printf( "1:datain = %e\n", datain[pos]  );
	END_CALCULATE2D( )

}


__GLOBAL__
void pack_iodata_y( float * datain, float * dataout, int nx, int ny, int nz, int J  )
{
	int _nx_ = nx + 2 * HALO;
	int _ny_ = ny + 2 * HALO;
	//int _nz_ = nz + 2 * HALO;

#ifdef GPU_CUDA
	int i0 = threadIdx.x + blockIdx.x * blockDim.x;
	int k0 = threadIdx.y + blockIdx.y * blockDim.y;
#else
	int i0 = 0;
	int k0 = 0;
#endif
	
	int i = i0 + HALO;
	int j = J;
	int k = k0 + HALO;

	long long index, pos;
	
	CALCULATE2D( i0, k0, 0, nx, 0, nz )
		i = i0 + HALO;
		k = k0 + HALO;
		index = INDEX( i, j, k );	
		pos = Index2D( i0, k0, nx, nz );
		dataout[pos] = datain[index];
		//printf( "2:datain = %e\n", datain[pos]  );
	END_CALCULATE2D( )

}


__GLOBAL__
void pack_iodata_z( float * datain, float * dataout, int nx, int ny, int nz, int K  )
{
	int _nx_ = nx + 2 * HALO;
	int _ny_ = ny + 2 * HALO;
	//int _nz_ = nz + 2 * HALO;

#ifdef GPU_CUDA
	int i0 = threadIdx.x + blockIdx.x * blockDim.x;
	int j0 = threadIdx.y + blockIdx.y * blockDim.y;
#else
	int i0 = 0;
	int j0 = 0;
#endif
	
	int i = i0 + HALO;
	int j = j0 + HALO;
	int k = K;

	long long index, pos;
	
	CALCULATE2D( i0, j0, 0, nx, 0, ny )
		i = i0 + HALO;
		j = j0 + HALO;
		index = INDEX( i, j, k );	
		pos = Index2D( i0, j0, nx, ny );
		dataout[pos] = datain[index];
		//printf( "3:datain = %e\n", datain[index]  );
	END_CALCULATE2D( )


}

void allocDataout_cpu(  GRID grid, char XYZ, float ** dataout )
{
	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;

	long long num = 0;
	
	switch( XYZ )
	{
		case 'X':
			num = ny * nz;
			break;
		case 'Y':
			num = nx * nz;
			break;
		case 'Z':
			num = nx * ny;
			break;
	}
		
	long long size = sizeof( float ) * num;

	*dataout = (float * )malloc( size );
	Memset(  *dataout, 0, size );
	
}

void freeDataout_cpu( float * dataout  )
{
	Free( dataout  );
}



void allocSliceData( GRID grid, SLICE slice, SLICE_DATA * sliceData )
{
	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;

	if ( slice.X >= HALO && slice.X < _nx )
		allocDataout( grid, 'X', &( sliceData->x ) );
	if ( slice.Y >= HALO && slice.Y < _ny )
		allocDataout( grid, 'Y', &( sliceData->y ) );
	if ( slice.Z >= HALO && slice.Z < _nz )
		allocDataout( grid, 'Z', &( sliceData->z ) );
}

void allocSliceData_cpu( GRID grid, SLICE slice, SLICE_DATA * sliceDataCpu )
{
	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;

	if ( slice.X >= HALO && slice.X < _nx )
		allocDataout_cpu( grid, 'X', &( sliceDataCpu->x ) );
	if ( slice.Y >= HALO && slice.Y < _ny )
		allocDataout_cpu( grid, 'Y', &( sliceDataCpu->y ) );
	if ( slice.Z >= HALO && slice.Z < _nz )
		allocDataout_cpu( grid, 'Z', &( sliceDataCpu->z ) );
}


void freeSliceData( GRID grid, SLICE slice, SLICE_DATA sliceData )
{
	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;

	if ( slice.X >= HALO && slice.X < _nx )
		freeDataout( sliceData.x );
	if ( slice.Y >= HALO && slice.Y < _ny )
		freeDataout( sliceData.y );
	if ( slice.Z >= HALO && slice.Z < _nz )
		freeDataout( sliceData.z );
}

void freeSliceData_cpu( GRID grid, SLICE slice, SLICE_DATA sliceDataCpu )
{
	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;

	if ( slice.X >= HALO && slice.X < _nx )
		freeDataout_cpu( sliceDataCpu.x );
	if ( slice.Y >= HALO && slice.Y < _ny )
		freeDataout_cpu( sliceDataCpu.y );
	if ( slice.Z >= HALO && slice.Z < _nz )
		freeDataout_cpu( sliceDataCpu.z );
}

void data2D_output_bin( GRID grid, SLICE slice, 
						MPI_COORD thisMPICoord, 
						float * datain, SLICE_DATA sliceData, SLICE_DATA sliceDataCpu,	
						const char * name )
{

	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;

	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;


	if ( slice.X >= HALO && slice.X < _nx )
	{

#ifdef GPU_CUDA
		dim3 threads( 32, 16, 1 );
		dim3 blocks;
		blocks.x = ( ny + threads.x - 1 ) / threads.x;
		blocks.y = ( nz + threads.y - 1 ) / threads.y;
		blocks.z = 1;
		pack_iodata_x
		<<< blocks, threads >>>
		( datain, sliceData.x, nx, ny, nz, slice.X  );
		long long size = sizeof( float ) * ny * nz;
		CHECK( cudaMemcpy( sliceDataCpu.x, sliceData.x, size, cudaMemcpyDeviceToHost ) );

#else
		pack_iodata_x( datain, sliceData.x, nx, ny, nz, slice.X  );

		sliceDataCpu.x = sliceData.x;

#endif

		FILE * fp;
		char fileName[256];
		sprintf( fileName, "%s_X_mpi_%d_%d_%d.bin", name, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );

		fp = fopen( fileName, "wb" ); 
		
		fwrite( sliceDataCpu.x, sizeof( float ), ny * nz, fp );
		
		fclose( fp );
	}

	if ( slice.Y >= HALO && slice.Y < _ny )
	{

#ifdef GPU_CUDA
		dim3 threads( 32, 16, 1 );
		dim3 blocks;
		blocks.x = ( nx + threads.x - 1 ) / threads.x;
		blocks.y = ( nz + threads.y - 1 ) / threads.y;
		blocks.z = 1;
		pack_iodata_y
		<<< blocks, threads >>>
		( datain, sliceData.y, nx, ny, nz, slice.Y  );
		long long size = sizeof( float ) * nx * nz;
		CHECK( cudaMemcpy( sliceDataCpu.y, sliceData.y, size, cudaMemcpyDeviceToHost ) );
#else
		pack_iodata_y( datain, sliceData.y, nx, ny, nz, slice.Y  );
		sliceDataCpu.y = sliceData.y;
#endif

		FILE * fp;
		char fileName[256];
		sprintf( fileName, "%s_Y_mpi_%d_%d_%d.bin", name, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );

		fp = fopen( fileName, "wb" ); 

		fwrite( sliceDataCpu.y, sizeof( float ), nx * nz, fp );
		
		fclose( fp );
	}

	if ( slice.Z >= HALO && slice.Z < _nz )
	{

#ifdef GPU_CUDA
		dim3 threads( 32, 16, 1 );
		dim3 blocks;
		blocks.x = ( nx + threads.x - 1 ) / threads.x;
		blocks.y = ( ny + threads.y - 1 ) / threads.y;
		blocks.z = 1;
		pack_iodata_z
		<<< blocks, threads >>>
		( datain, sliceData.z, nx, ny, nz, slice.Z  );
		long long size = sizeof( float ) * nx * ny;
		CHECK( cudaMemcpy( sliceDataCpu.z, sliceData.z, size, cudaMemcpyDeviceToHost ) );
#else
		pack_iodata_z( datain, sliceData.z, nx, ny, nz, slice.Z  );
		sliceDataCpu.z = sliceData.z;
#endif

		FILE * fp;
		char fileName[256];
		sprintf( fileName, "%s_Z_mpi_%d_%d_%d.bin", name, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );

		fp = fopen( fileName, "wb" ); 
		
		fwrite( sliceDataCpu.z, sizeof( float ), nx * ny, fp );

		fclose( fp );
	}

}



