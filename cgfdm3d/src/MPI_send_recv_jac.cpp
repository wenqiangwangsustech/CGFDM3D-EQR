/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:MPI_send_recv.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-17
*   Discription:
*
================================================================*/
#include "header.h"

typedef void (*PACK_UNPACK_FUNC)( float * jac, float * thisSend, int xStartHalo, int _nx_, int _ny_, int _nz_ );


__GLOBAL__
void packJacX( float * jac, float * thisSend, 
	int xStartHalo, int _nx_, int _ny_, int _nz_ )
{
	//printf("packJac_MPI_x\n");
#ifdef GPU_CUDA
	int i0 = threadIdx.x + blockIdx.x * blockDim.x;
	int j0 = threadIdx.y + blockIdx.y * blockDim.y;
	int k0 = threadIdx.z + blockIdx.z * blockDim.z;
#else
	int i0 = 0;
	int j0 = 0;
	int k0 = 0;
#endif

	long long index;
	long long pos;

	int i = 0;
	int j = 0;
	int k = 0;
	
	CALCULATE3D( i0, j0, k0, 0, HALO, 0, _ny_, 0, _nz_ )
		i = i0 + xStartHalo;
		j = j0;
		k = k0;
		//k = k0 / _nz_ * _nz_ + k0 % _nz_;
		index = INDEX( i, j, k );
		pos = Index3D( i0, j0, k0, HALO, _ny_, _nz_ );//i - xStartHalo + j * HALO + k * HALO * _ny_;
		thisSend[pos] = jac[index];
	END_CALCULATE3D( )
}

__GLOBAL__
void unpackJacX( float * jac, float * thisRecv,  
	int xStartHalo, int _nx_, int _ny_, int _nz_ )
{
#ifdef GPU_CUDA
	int i0 = threadIdx.x + blockIdx.x * blockDim.x;
	int j0 = threadIdx.y + blockIdx.y * blockDim.y;
	int k0 = threadIdx.z + blockIdx.z * blockDim.z;
#else
	int i0 = 0;
	int j0 = 0;
	int k0 = 0;
#endif

	long long index;
	long long pos;

	int i = 0;
	int j = 0;
	int k = 0;
	
	CALCULATE3D( i0, j0, k0, 0, HALO, 0, _ny_, 0, _nz_ )
		i = i0 + xStartHalo;
		j = j0;
		k = k0;
		//k = k0 / _nz_ * _nz_ + k0 % _nz_;
		index = INDEX( i, j, k );
		pos = Index3D( i0, j0, k0, HALO, _ny_, _nz_ );//i - xStartHalo + j * HALO + k * HALO * _ny_;
		jac[index] = thisRecv[pos];
	END_CALCULATE3D( )
}

void PackUnpackJacX( float * jac, float * thisSendRecv, 
	int xStartHalo, int _nx_, int _ny_, int _nz_, PACK_UNPACK_FUNC packJac_unpackJac_func )
{
#ifdef GPU_CUDA
	dim3 threads( 4, 8, 16);
	dim3 blocks;
	blocks.x = ( HALO + threads.x - 1 ) / threads.x;
	blocks.y = ( _ny_ + threads.y - 1 ) / threads.y;
	blocks.z = ( _nz_ + threads.z - 1 ) / threads.z;
	packJac_unpackJac_func<<< blocks, threads >>>
	( jac, thisSendRecv, xStartHalo, _nx_, _ny_, _nz_ );
	CHECK( cudaDeviceSynchronize( ) );
#else
	packJac_unpackJac_func
	( jac, thisSendRecv, xStartHalo, _nx_, _ny_, _nz_ );
#endif
}


__GLOBAL__
void packJacY( float * jac, float * thisSend, 
	int yStartHalo, int _nx_, int _ny_, int _nz_ )
{
	//printf("packJac_MPI_y\n");
#ifdef GPU_CUDA
	int i0 = threadIdx.x + blockIdx.x * blockDim.x;
	int j0 = threadIdx.y + blockIdx.y * blockDim.y;
	int k0 = threadIdx.z + blockIdx.z * blockDim.z;
#else
	int i0 = 0;
	int j0 = 0;
	int k0 = 0;
#endif

	long long index;
	long long pos;

	int i = 0;
	int j = 0;
	int k = 0;

	CALCULATE3D( i0, j0, k0, 0, _nx_, 0, HALO, 0, _nz_ )
		i = i0;
		j = j0 + yStartHalo;
		k = k0;
		index = INDEX( i, j, k );
		pos = Index3D( i0, j0, k0, _nx_, HALO, _nz_ );//i + ( j - yStartHalo ) * _nx_ + k * HALO * _nx_;
		thisSend[pos] = jac[index];
	END_CALCULATE3D( )

}

__GLOBAL__
void unpackJacY( float * jac, float * thisRecv, 
	int yStartHalo, int _nx_, int _ny_, int _nz_ )
{
	//printf("unpackJac_MPI_y\n");
#ifdef GPU_CUDA
	int i0 = threadIdx.x + blockIdx.x * blockDim.x;
	int j0 = threadIdx.y + blockIdx.y * blockDim.y;
	int k0 = threadIdx.z + blockIdx.z * blockDim.z;
#else
	int i0 = 0;
	int j0 = 0;
	int k0 = 0;
#endif

	long long index;
	long long pos;

	int i = 0;
	int j = 0;
	int k = 0;

	CALCULATE3D( i0, j0, k0, 0, _nx_, 0, HALO, 0, _nz_ )
	    i = i0;
	    j = j0 + yStartHalo;
	    k = k0;
		index = INDEX( i, j, k );
		pos = Index3D( i0, j0, k0, _nx_, HALO, _nz_ );//i + ( j - yStartHalo ) * _nx_ + k * HALO * _nx_;
		jac[index] = thisRecv[pos];
	END_CALCULATE3D( )

}
void PackUnpackJacY( float * jac, float * thisSendRecv, 
	int yStartHalo, int _nx_, int _ny_, int _nz_, PACK_UNPACK_FUNC packJac_unpackJac_func )
{
#ifdef GPU_CUDA
	dim3 threads( 8, 4, 16);
	dim3 blocks;
	blocks.x = ( _nx_ + threads.x - 1 ) / threads.x;
	blocks.y = ( HALO + threads.y - 1 ) / threads.y;
	blocks.z = ( _nz_ + threads.z - 1 ) / threads.z;
	packJac_unpackJac_func<<< blocks, threads >>>
	( jac, thisSendRecv, yStartHalo, _nx_, _ny_, _nz_ );
	CHECK( cudaDeviceSynchronize( ) );
#else
	packJac_unpackJac_func
	( jac, thisSendRecv, yStartHalo, _nx_, _ny_, _nz_ );
#endif
}

void packJacZ( float * jac, float * thisSend, int zStartHalo, int _nx_, int _ny_, int _nz_ )
{
	long long blockLen = _nx_ * _ny_ * HALO;
	long long size = sizeof( float ) * blockLen; 
	long long zStartStride = zStartHalo * _nx_ * _ny_;
	long long num = _nx_ * _ny_ * _nz_; 

#ifdef GPU_CUDA
	CHECK( Memcpy( thisSend + 0 * blockLen, jac + zStartStride, size, cudaMemcpyDeviceToDevice ) );
#else
	Memcpy( thisSend + 0 * blockLen, jac + zStartStride, size );
#endif
}

void unpackJacZ( float * jac, float * thisRecv, int zStartHalo, int _nx_, int _ny_, int _nz_ )
{
	long long blockLen = _nx_ * _ny_ * HALO;
	long long size = sizeof( float ) * blockLen; 
	long long zStartStride = zStartHalo * _nx_ * _ny_;
	long long num = _nx_ * _ny_ * _nz_; 

#ifdef GPU_CUDA
	CHECK( Memcpy( jac + zStartStride, thisRecv + 0 * blockLen, size, cudaMemcpyDeviceToDevice ) );
#else
	Memcpy( jac  + zStartStride, thisRecv + 0 * blockLen, size );
#endif
}
void PackUnpackJacZ( float * jac, float * thisSendRecv, int zStartHalo, int _nx_, int _ny_, int _nz_,  PACK_UNPACK_FUNC packJac_unpackJac_func )
{

	packJac_unpackJac_func( jac, thisSendRecv, zStartHalo, _nx_, _ny_, _nz_ );

}


void mpiSendRecvJac( GRID grid, MPI_Comm comm_cart, MPI_NEIGHBOR mpiNeighbor, float * jac, SEND_RECV_DATA sr )
{

	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;

	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;

	long long num = 0;

	float * thisXSend1 = sr.thisXSend1;
	float * thisXRecv1 = sr.thisXRecv1;
	float * thisYSend1 = sr.thisYSend1;
	float * thisYRecv1 = sr.thisYRecv1;
	float * thisZSend1 = sr.thisZSend1;
	float * thisZRecv1 = sr.thisZRecv1;
                                 
	float * thisXSend2 = sr.thisXSend2;
	float * thisXRecv2 = sr.thisXRecv2;
	float * thisYSend2 = sr.thisYSend2;
	float * thisYRecv2 = sr.thisYRecv2;
	float * thisZSend2 = sr.thisZSend2;
	float * thisZRecv2 = sr.thisZRecv2;

	int xStartHalo, yStartHalo, zStartHalo;


	MPI_Status stat;

//x direction data exchange

	xStartHalo = nx;
	if ( mpiNeighbor.X2 >= 0 ) PackUnpackJacX( jac, thisXSend2, xStartHalo, _nx_, _ny_, _nz_,   packJacX );

	//printf( "X1 = %d, X2 = %d\n", mpiNeighbor.X1, mpiNeighbor.X2  );
	//printf( "======================================================\n"  );

	num = HALO * _ny_ * _nz_;
	MPI_Sendrecv( sr.thisXSend2, num, MPI_FLOAT, mpiNeighbor.X2, 101,
				  sr.thisXRecv1, num, MPI_FLOAT, mpiNeighbor.X1, 101,
				  comm_cart, &stat );

	//printf( "X1 = %d, X2 = %d\n", mpiNeighbor.X1, mpiNeighbor.X2  );
	xStartHalo = 0;
	if ( mpiNeighbor.X1 >= 0 ) PackUnpackJacX( jac, thisXRecv1, xStartHalo, _nx_, _ny_, _nz_, unpackJacX );

	
	xStartHalo = HALO;
	if ( mpiNeighbor.X1 >= 0 ) PackUnpackJacX( jac, thisXSend1, xStartHalo, _nx_, _ny_, _nz_,   packJacX );

	num = HALO * _ny_ * _nz_;
	MPI_Sendrecv( sr.thisXSend1, num, MPI_FLOAT, mpiNeighbor.X1, 102,
				  sr.thisXRecv2, num, MPI_FLOAT, mpiNeighbor.X2, 102,
				  comm_cart, &stat );

	xStartHalo = _nx;
	if ( mpiNeighbor.X2 >= 0 ) PackUnpackJacX( jac, thisXRecv2, xStartHalo, _nx_, _ny_, _nz_, unpackJacX );

//y direction data exchange
	yStartHalo = ny;
	if ( mpiNeighbor.Y2 >= 0 ) PackUnpackJacY( jac, thisYSend2, yStartHalo, _nx_, _ny_, _nz_,   packJacY );

	num = HALO * _nx_ * _nz_;
	MPI_Sendrecv( sr.thisYSend2, num, MPI_FLOAT, mpiNeighbor.Y2, 103,
				  sr.thisYRecv1, num, MPI_FLOAT, mpiNeighbor.Y1, 103,
				  comm_cart, &stat );

	yStartHalo = 0;
	if ( mpiNeighbor.Y1 >= 0 ) PackUnpackJacY( jac, thisYRecv1, yStartHalo, _nx_, _ny_, _nz_, unpackJacY );

	
	yStartHalo = HALO;
	if ( mpiNeighbor.Y1 >= 0 ) PackUnpackJacY( jac, thisYSend1, yStartHalo, _nx_, _ny_, _nz_,   packJacY );

	num = HALO * _nx_ * _nz_;
	MPI_Sendrecv( sr.thisYSend1, num, MPI_FLOAT, mpiNeighbor.Y1, 104,
				  sr.thisYRecv2, num, MPI_FLOAT, mpiNeighbor.Y2, 104,
				  comm_cart, &stat );

	yStartHalo = _ny;
	if ( mpiNeighbor.Y2 >= 0 ) PackUnpackJacY( jac, thisYRecv2, yStartHalo, _nx_, _ny_, _nz_, unpackJacY );

//z direction data exchange
	zStartHalo = nz;
	if ( mpiNeighbor.Z2 >= 0 ) PackUnpackJacZ( jac, thisZSend2, zStartHalo, _nx_, _ny_, _nz_,   packJacZ );

	//char fileName[256];

	//if ( mpiNeighbor.Z2 > 0 )
	//{
	//	sprintf( fileName, "./output2/thisZSend2_%d_%d_%d", thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
	//	FILE * fp = fopen( fileName, "wb" );	
	//	fwrite( thisZSend2, sizeof( float ), ny * nz, fp );
	//
	//	fclose( fp );
	//
	//}	
	
	num = HALO * _nx_ * _ny_;
	MPI_Sendrecv( sr.thisZSend2, num, MPI_FLOAT, mpiNeighbor.Z2, 105,
				  sr.thisZRecv1, num, MPI_FLOAT, mpiNeighbor.Z1, 105,
				  comm_cart, &stat );

	zStartHalo = 0;
	if ( mpiNeighbor.Z1 >= 0 ) PackUnpackJacZ( jac, thisZRecv1, zStartHalo, _nx_, _ny_, _nz_, unpackJacZ );

	
	zStartHalo = HALO;
	if ( mpiNeighbor.Z1 >= 0 ) PackUnpackJacZ( jac, thisZSend1, zStartHalo, _nx_, _ny_, _nz_,   packJacZ );

	num = HALO * _nx_ * _ny_;
	MPI_Sendrecv( sr.thisZSend1, num, MPI_FLOAT, mpiNeighbor.Z1, 106,
				  sr.thisZRecv2, num, MPI_FLOAT, mpiNeighbor.Z2, 106,
				  comm_cart, &stat );

	zStartHalo = _nz;
	if ( mpiNeighbor.Z2 >= 0 ) PackUnpackJacZ( jac, thisZRecv2, zStartHalo, _nx_, _ny_, _nz_, unpackJacZ );
}




