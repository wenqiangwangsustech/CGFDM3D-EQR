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

typedef void (*PACK_UNPACK_FUNC)( float * waveData, float * thisSend, int xStartHalo, int _nx_, int _ny_, int _nz_ );


__GLOBAL__
void packX( float * waveData, float * thisSend, 
	int xStartHalo, int _nx_, int _ny_, int _nz_ )
{
	//printf("pack_MPI_x\n");
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
	
	CALCULATE3D( i0, j0, k0, 0, HALO, 0, _ny_, 0, _nz_ * WAVESIZE )
		i = i0 + xStartHalo;
		j = j0;
		k = k0;
		//k = k0 / _nz_ * _nz_ + k0 % _nz_;
		index = INDEX( i, j, k );
		pos = Index3D( i0, j0, k0, HALO, _ny_, _nz_ * WAVESIZE );//i - xStartHalo + j * HALO + k * HALO * _ny_;
		thisSend[pos] = waveData[index];
	END_CALCULATE3D( )
}

__GLOBAL__
void unpackX( float * waveData, float * thisRecv,  
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
	
	CALCULATE3D( i0, j0, k0, 0, HALO, 0, _ny_, 0, _nz_ * WAVESIZE )
		i = i0 + xStartHalo;
		j = j0;
		k = k0;
		//k = k0 / _nz_ * _nz_ + k0 % _nz_;
		index = INDEX( i, j, k );
		pos = Index3D( i0, j0, k0, HALO, _ny_, _nz_ * WAVESIZE );//i - xStartHalo + j * HALO + k * HALO * _ny_;
		waveData[index] = thisRecv[pos];
	END_CALCULATE3D( )
}

void PackUnpackX( float * waveData, float * thisSendRecv, 
	int xStartHalo, int _nx_, int _ny_, int _nz_, PACK_UNPACK_FUNC pack_unpack_func )
{
#ifdef GPU_CUDA
	dim3 threads( 4, 8, 16);
	dim3 blocks;
	blocks.x = ( HALO + threads.x - 1 ) / threads.x;
	blocks.y = ( _ny_ + threads.y - 1 ) / threads.y;
	blocks.z = ( _nz_ * WAVESIZE + threads.z - 1 ) / threads.z;
	pack_unpack_func<<< blocks, threads >>>
	( waveData, thisSendRecv, xStartHalo, _nx_, _ny_, _nz_ );
	CHECK( cudaDeviceSynchronize( ) );
#else
	pack_unpack_func
	( waveData, thisSendRecv, xStartHalo, _nx_, _ny_, _nz_ );
#endif
}


__GLOBAL__
void packY( float * waveData, float * thisSend, 
	int yStartHalo, int _nx_, int _ny_, int _nz_ )
{
	//printf("pack_MPI_y\n");
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

	CALCULATE3D( i0, j0, k0, 0, _nx_, 0, HALO, 0, _nz_ * WAVESIZE )
		i = i0;
		j = j0 + yStartHalo;
		k = k0;
		index = INDEX( i, j, k );
		pos = Index3D( i0, j0, k0, _nx_, HALO, _nz_ * WAVESIZE );//i + ( j - yStartHalo ) * _nx_ + k * HALO * _nx_;
		thisSend[pos] = waveData[index];
	END_CALCULATE3D( )

}

__GLOBAL__
void unpackY( float * waveData, float * thisRecv, 
	int yStartHalo, int _nx_, int _ny_, int _nz_ )
{
	//printf("unpack_MPI_y\n");
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

	CALCULATE3D( i0, j0, k0, 0, _nx_, 0, HALO, 0, _nz_ * WAVESIZE )
	    i = i0;
	    j = j0 + yStartHalo;
	    k = k0;
		index = INDEX( i, j, k );
		pos = Index3D( i0, j0, k0, _nx_, HALO, _nz_ * WAVESIZE );//i + ( j - yStartHalo ) * _nx_ + k * HALO * _nx_;
		waveData[index] = thisRecv[pos];
	END_CALCULATE3D( )

}
void PackUnpackY( float * waveData, float * thisSendRecv, 
	int yStartHalo, int _nx_, int _ny_, int _nz_, PACK_UNPACK_FUNC pack_unpack_func )
{
#ifdef GPU_CUDA
	dim3 threads( 8, 4, 16);
	dim3 blocks;
	blocks.x = ( _nx_ + threads.x - 1 ) / threads.x;
	blocks.y = ( HALO + threads.y - 1 ) / threads.y;
	blocks.z = ( _nz_ * WAVESIZE + threads.z - 1 ) / threads.z;
	pack_unpack_func<<< blocks, threads >>>
	( waveData, thisSendRecv, yStartHalo, _nx_, _ny_, _nz_ );
	CHECK( cudaDeviceSynchronize( ) );
#else
	pack_unpack_func
	( waveData, thisSendRecv, yStartHalo, _nx_, _ny_, _nz_ );
#endif
}

void packZ( float * waveData, float * thisSend, int zStartHalo, int _nx_, int _ny_, int _nz_ )
{
	long long blockLen = _nx_ * _ny_ * HALO;
	long long size = sizeof( float ) * blockLen; 
	long long zStartStride = zStartHalo * _nx_ * _ny_;
	long long num = _nx_ * _ny_ * _nz_; 

	WAVE W = { 0 };

	W.Vx  = waveData + 0 * num;
	W.Vy  = waveData + 1 * num;
	W.Vz  = waveData + 2 * num;
	W.Txx = waveData + 3 * num;
	W.Tyy = waveData + 4 * num;
	W.Tzz = waveData + 5 * num;
	W.Txy = waveData + 6 * num;
	W.Txz = waveData + 7 * num;
	W.Tyz = waveData + 8 * num;

#ifdef GPU_CUDA
	CHECK( Memcpy( thisSend + 0 * blockLen, W.Vx  + zStartStride, size, cudaMemcpyDeviceToDevice ) );
	CHECK( Memcpy( thisSend + 1 * blockLen, W.Vy  + zStartStride, size, cudaMemcpyDeviceToDevice ) );
	CHECK( Memcpy( thisSend + 2 * blockLen, W.Vz  + zStartStride, size, cudaMemcpyDeviceToDevice ) );
	CHECK( Memcpy( thisSend + 3 * blockLen, W.Txx + zStartStride, size, cudaMemcpyDeviceToDevice ) );
	CHECK( Memcpy( thisSend + 4 * blockLen, W.Tyy + zStartStride, size, cudaMemcpyDeviceToDevice ) );
	CHECK( Memcpy( thisSend + 5 * blockLen, W.Tzz + zStartStride, size, cudaMemcpyDeviceToDevice ) );
	CHECK( Memcpy( thisSend + 6 * blockLen, W.Txy + zStartStride, size, cudaMemcpyDeviceToDevice ) );
	CHECK( Memcpy( thisSend + 7 * blockLen, W.Txz + zStartStride, size, cudaMemcpyDeviceToDevice ) );
	CHECK( Memcpy( thisSend + 8 * blockLen, W.Tyz + zStartStride, size, cudaMemcpyDeviceToDevice ) );
#else
	Memcpy( thisSend + 0 * blockLen, W.Vx  + zStartStride, size );
	Memcpy( thisSend + 1 * blockLen, W.Vy  + zStartStride, size );
	Memcpy( thisSend + 2 * blockLen, W.Vz  + zStartStride, size );
	Memcpy( thisSend + 3 * blockLen, W.Txx + zStartStride, size );
	Memcpy( thisSend + 4 * blockLen, W.Tyy + zStartStride, size );
	Memcpy( thisSend + 5 * blockLen, W.Tzz + zStartStride, size );
	Memcpy( thisSend + 6 * blockLen, W.Txy + zStartStride, size );
	Memcpy( thisSend + 7 * blockLen, W.Txz + zStartStride, size );
	Memcpy( thisSend + 8 * blockLen, W.Tyz + zStartStride, size );
#endif
}

void unpackZ( float * waveData, float * thisRecv, int zStartHalo, int _nx_, int _ny_, int _nz_ )
{
	long long blockLen = _nx_ * _ny_ * HALO;
	long long size = sizeof( float ) * blockLen; 
	long long zStartStride = zStartHalo * _nx_ * _ny_;
	long long num = _nx_ * _ny_ * _nz_; 

	WAVE W = { 0 };

	W.Vx  = waveData + 0 * num;
	W.Vy  = waveData + 1 * num;
	W.Vz  = waveData + 2 * num;
	W.Txx = waveData + 3 * num;
	W.Tyy = waveData + 4 * num;
	W.Tzz = waveData + 5 * num;
	W.Txy = waveData + 6 * num;
	W.Txz = waveData + 7 * num;
	W.Tyz = waveData + 8 * num;

#ifdef GPU_CUDA
	CHECK( Memcpy( W.Vx  + zStartStride, thisRecv + 0 * blockLen, size, cudaMemcpyDeviceToDevice ) );
	CHECK( Memcpy( W.Vy  + zStartStride, thisRecv + 1 * blockLen, size, cudaMemcpyDeviceToDevice ) );
	CHECK( Memcpy( W.Vz  + zStartStride, thisRecv + 2 * blockLen, size, cudaMemcpyDeviceToDevice ) );
	CHECK( Memcpy( W.Txx + zStartStride, thisRecv + 3 * blockLen, size, cudaMemcpyDeviceToDevice ) );
	CHECK( Memcpy( W.Tyy + zStartStride, thisRecv + 4 * blockLen, size, cudaMemcpyDeviceToDevice ) );
	CHECK( Memcpy( W.Tzz + zStartStride, thisRecv + 5 * blockLen, size, cudaMemcpyDeviceToDevice ) );
	CHECK( Memcpy( W.Txy + zStartStride, thisRecv + 6 * blockLen, size, cudaMemcpyDeviceToDevice ) );
	CHECK( Memcpy( W.Txz + zStartStride, thisRecv + 7 * blockLen, size, cudaMemcpyDeviceToDevice ) );
	CHECK( Memcpy( W.Tyz + zStartStride, thisRecv + 8 * blockLen, size, cudaMemcpyDeviceToDevice ) );
#else
	Memcpy( W.Vx  + zStartStride, thisRecv + 0 * blockLen, size );
	Memcpy( W.Vy  + zStartStride, thisRecv + 1 * blockLen, size );
	Memcpy( W.Vz  + zStartStride, thisRecv + 2 * blockLen, size );
	Memcpy( W.Txx + zStartStride, thisRecv + 3 * blockLen, size );
	Memcpy( W.Tyy + zStartStride, thisRecv + 4 * blockLen, size );
	Memcpy( W.Tzz + zStartStride, thisRecv + 5 * blockLen, size );
	Memcpy( W.Txy + zStartStride, thisRecv + 6 * blockLen, size );
	Memcpy( W.Txz + zStartStride, thisRecv + 7 * blockLen, size );
	Memcpy( W.Tyz + zStartStride, thisRecv + 8 * blockLen, size );
#endif
}
void PackUnpackZ( float * waveData, float * thisSendRecv, int zStartHalo, int _nx_, int _ny_, int _nz_,  PACK_UNPACK_FUNC pack_unpack_func )
{

	pack_unpack_func( waveData, thisSendRecv, zStartHalo, _nx_, _ny_, _nz_ );

}

void alloc_send_recv( GRID grid, float ** send, float ** recv, char XYZ )
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	long long num = 0;
	
	switch( XYZ )
	{
		case 'X':
			num = _ny_ * _nz_* HALO * WAVESIZE;
			break;       
		case 'Y':         
			num = _nx_ * _nz_* HALO * WAVESIZE;
			break;       
		case 'Z':         
			num = _nx_ * _ny_* HALO * WAVESIZE;
			break;
	}
		
	float * pSend = NULL;
	float * pRecv = NULL;
	long long size = sizeof( float ) * num;

	CHECK( Malloc( ( void ** )&pSend, size ) );
	CHECK( Memset(  pSend, 0, size ) ); 
	
	CHECK( Malloc( ( void ** )&pRecv, size ) );
	CHECK( Memset(  pRecv, 0, size ) ); 

	*send = pSend;
	*recv = pRecv;

}


void allocSendRecv( GRID grid, MPI_NEIGHBOR mpiNeighbor, SEND_RECV_DATA * sr )
{
	memset( sr, 0, sizeof( SEND_RECV_DATA ) );
	
	/*if ( mpiNeighbor.X1 > 0 )*/	alloc_send_recv( grid, &( sr->thisXSend1 ), &( sr->thisXRecv1 ), 'X' );
	/*if ( mpiNeighbor.Y1 > 0 )*/	alloc_send_recv( grid, &( sr->thisYSend1 ), &( sr->thisYRecv1 ), 'Y' );
	/*if ( mpiNeighbor.Z1 > 0 )*/	alloc_send_recv( grid, &( sr->thisZSend1 ), &( sr->thisZRecv1 ), 'Z' ); 

	/*if ( mpiNeighbor.X2 > 0 )*/	alloc_send_recv( grid, &( sr->thisXSend2 ), &( sr->thisXRecv2 ), 'X' );
	/*if ( mpiNeighbor.Y2 > 0 )*/	alloc_send_recv( grid, &( sr->thisYSend2 ), &( sr->thisYRecv2 ), 'Y' );
	/*if ( mpiNeighbor.Z2 > 0 )*/	alloc_send_recv( grid, &( sr->thisZSend2 ), &( sr->thisZRecv2 ), 'Z' );
	
}

void freeSendRecv( MPI_NEIGHBOR mpiNeighbor, SEND_RECV_DATA sr )
{
	/*if ( mpiNeighbor.X1 > 0 )*/	{ Free( sr.thisXSend1);	Free( sr.thisXRecv1 ); };
	/*if ( mpiNeighbor.Y1 > 0 )*/	{ Free( sr.thisYSend1);	Free( sr.thisYRecv1 ); };
	/*if ( mpiNeighbor.Z1 > 0 )*/	{ Free( sr.thisZSend1);	Free( sr.thisZRecv1 ); }; 
                                 
	/*if ( mpiNeighbor.X2 > 0 )*/	{ Free( sr.thisXSend2);	Free( sr.thisXRecv2 ); };
	/*if ( mpiNeighbor.Y2 > 0 )*/	{ Free( sr.thisYSend2);	Free( sr.thisYRecv2 ); };
	/*if ( mpiNeighbor.Z2 > 0 )*/	{ Free( sr.thisZSend2);	Free( sr.thisZRecv2 ); };
	
}


void mpiSendRecv( GRID grid, MPI_Comm comm_cart, MPI_NEIGHBOR mpiNeighbor, WAVE W, SEND_RECV_DATA sr )
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

	float * waveData = W.Vx;

	MPI_Status stat;

//x direction data exchange

	xStartHalo = nx;
	if ( mpiNeighbor.X2 >= 0 ) PackUnpackX( waveData, thisXSend2, xStartHalo, _nx_, _ny_, _nz_,   packX );

	//printf( "X1 = %d, X2 = %d\n", mpiNeighbor.X1, mpiNeighbor.X2  );
	//printf( "======================================================\n"  );

	num = HALO * _ny_ * _nz_ * WAVESIZE;
	MPI_Sendrecv( sr.thisXSend2, num, MPI_FLOAT, mpiNeighbor.X2, 101,
				  sr.thisXRecv1, num, MPI_FLOAT, mpiNeighbor.X1, 101,
				  comm_cart, &stat );

	//printf( "X1 = %d, X2 = %d\n", mpiNeighbor.X1, mpiNeighbor.X2  );
	xStartHalo = 0;
	if ( mpiNeighbor.X1 >= 0 ) PackUnpackX( waveData, thisXRecv1, xStartHalo, _nx_, _ny_, _nz_, unpackX );

	
	xStartHalo = HALO;
	if ( mpiNeighbor.X1 >= 0 ) PackUnpackX( waveData, thisXSend1, xStartHalo, _nx_, _ny_, _nz_,   packX );

	num = HALO * _ny_ * _nz_ * WAVESIZE;
	MPI_Sendrecv( sr.thisXSend1, num, MPI_FLOAT, mpiNeighbor.X1, 102,
				  sr.thisXRecv2, num, MPI_FLOAT, mpiNeighbor.X2, 102,
				  comm_cart, &stat );

	xStartHalo = _nx;
	if ( mpiNeighbor.X2 >= 0 ) PackUnpackX( waveData, thisXRecv2, xStartHalo, _nx_, _ny_, _nz_, unpackX );

//y direction data exchange
	yStartHalo = ny;
	if ( mpiNeighbor.Y2 >= 0 ) PackUnpackY( waveData, thisYSend2, yStartHalo, _nx_, _ny_, _nz_,   packY );

	num = HALO * _nx_ * _nz_ * WAVESIZE;
	MPI_Sendrecv( sr.thisYSend2, num, MPI_FLOAT, mpiNeighbor.Y2, 103,
				  sr.thisYRecv1, num, MPI_FLOAT, mpiNeighbor.Y1, 103,
				  comm_cart, &stat );

	yStartHalo = 0;
	if ( mpiNeighbor.Y1 >= 0 ) PackUnpackY( waveData, thisYRecv1, yStartHalo, _nx_, _ny_, _nz_, unpackY );

	
	yStartHalo = HALO;
	if ( mpiNeighbor.Y1 >= 0 ) PackUnpackY( waveData, thisYSend1, yStartHalo, _nx_, _ny_, _nz_,   packY );

	num = HALO * _nx_ * _nz_ * WAVESIZE;
	MPI_Sendrecv( sr.thisYSend1, num, MPI_FLOAT, mpiNeighbor.Y1, 104,
				  sr.thisYRecv2, num, MPI_FLOAT, mpiNeighbor.Y2, 104,
				  comm_cart, &stat );

	yStartHalo = _ny;
	if ( mpiNeighbor.Y2 >= 0 ) PackUnpackY( waveData, thisYRecv2, yStartHalo, _nx_, _ny_, _nz_, unpackY );

//z direction data exchange
	zStartHalo = nz;
	if ( mpiNeighbor.Z2 >= 0 ) PackUnpackZ( waveData, thisZSend2, zStartHalo, _nx_, _ny_, _nz_,   packZ );

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
	
	num = HALO * _nx_ * _ny_ * WAVESIZE;
	MPI_Sendrecv( sr.thisZSend2, num, MPI_FLOAT, mpiNeighbor.Z2, 105,
				  sr.thisZRecv1, num, MPI_FLOAT, mpiNeighbor.Z1, 105,
				  comm_cart, &stat );

	zStartHalo = 0;
	if ( mpiNeighbor.Z1 >= 0 ) PackUnpackZ( waveData, thisZRecv1, zStartHalo, _nx_, _ny_, _nz_, unpackZ );

	
	zStartHalo = HALO;
	if ( mpiNeighbor.Z1 >= 0 ) PackUnpackZ( waveData, thisZSend1, zStartHalo, _nx_, _ny_, _nz_,   packZ );

	num = HALO * _nx_ * _ny_ * WAVESIZE;
	MPI_Sendrecv( sr.thisZSend1, num, MPI_FLOAT, mpiNeighbor.Z1, 106,
				  sr.thisZRecv2, num, MPI_FLOAT, mpiNeighbor.Z2, 106,
				  comm_cart, &stat );

	zStartHalo = _nz;
	if ( mpiNeighbor.Z2 >= 0 ) PackUnpackZ( waveData, thisZRecv2, zStartHalo, _nx_, _ny_, _nz_, unpackZ );
}




