/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:source.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-14
*   Discription:
*
================================================================*/
#include "header.h"



void readSourceInfo  ( SOURCE_FILE_INPUT * src_in,  char * sourceFileName )
{
	FILE * sourceFile = fopen( sourceFileName, "rb" );
	
#ifndef SOURCE_NPTS_LONG_LONG
	int nnnnnn;
	fread( &( nnnnnn ), sizeof( int ), 1, sourceFile ); 
	src_in->npts = nnnnnn;
#else
	fread( &( src_in->npts ), sizeof( long long  ), 1, sourceFile ); 
#endif
	fread( &( src_in->nt ),   sizeof( int ), 1, sourceFile ); 
	fread( &( src_in->dt ),   sizeof( float ), 1, sourceFile ); 
	
	fclose( sourceFile );

}
void allocSourceLocation( SOURCE_FILE_INPUT * src_in  )
{
	long long npts = src_in->npts;
	
	long long size = sizeof( float ) * npts * COORDSIZE;

	float * pCoord = NULL;
	pCoord = ( float * )malloc( size );
	memset( pCoord,   0, size );

	src_in->lon =		pCoord + 0 * npts;
	src_in->lat =		pCoord + 1 * npts;
	src_in->coordZ =    pCoord + 2 * npts;
	
}


void freeSourceLocation( SOURCE_FILE_INPUT src_in  )
{
	free( src_in.lon );
}

void readSourceLocation( SOURCE_FILE_INPUT src_in,  char * sourceFileName )
{
	FILE * sourceFile = fopen( sourceFileName, "rb" );

	if ( NULL == sourceFile )
		cout << sourceFileName << ", is source file" << endl;
	fseek( sourceFile, sizeof( int ) + sizeof( int ) + sizeof( float ), SEEK_CUR );

	long long npts = src_in.npts;
	int nt = src_in.nt;
	
	long long i = 0;
	for ( i = 0; i < npts; i ++ )
	{
		fread( &( src_in.lon[i]    ), 1, sizeof( float ), sourceFile ); 
		fread( &( src_in.lat[i]    ), 1, sizeof( float ), sourceFile ); 
		fread( &( src_in.coordZ[i] ), 1, sizeof( float ), sourceFile ); 
		
		fseek( sourceFile, 3 * sizeof( float ) + 2 * nt * sizeof( float ), SEEK_CUR );
	}

	//int nnnn = 1000;

	//cout << "lon[" << nnnn  <<"] = " << src_in.lon[nnnn] << endl;
	//cout << "npts = " << npts << ", nt = " << nt << ", dt = " << src_in.dt << endl;
	fclose( sourceFile );
}


void LonLat2cart( PJ * P,  PJ_DIRECTION PD, SOURCE_FILE_INPUT src_in )
{
	long long npts = src_in.npts;

	PJ_COORD * pj_coord;


	pj_coord = ( PJ_COORD * )malloc( sizeof( PJ_COORD ) * src_in.npts );

	long long i = 0;
	for ( i = 0; i < npts; i ++ )
	{
		pj_coord[i].lp.lam = src_in.lon[i] * DEGREE2RADIAN;
		pj_coord[i].lp.phi = src_in.lat[i] * DEGREE2RADIAN;
	}
	

	proj_trans_array( P, PD, npts, pj_coord );


	for ( i = 0; i < npts; i ++ )
	{
		src_in.lon[i] = pj_coord[i].xy.x;
		src_in.lat[i] = pj_coord[i].xy.y;
	}
	

	free( pj_coord );
}

void projTrans( PARAMS params, SOURCE_FILE_INPUT src_in )
{
	float lon_0 = params.centerLongitude;
	float lat_0 = params.centerLatitude;

	long long i = 0;
	float maxLon = - 180.0, minLon = 180.0, maxLat = -90.0, minLat = 90.0;
	for ( i = 0; i < src_in.npts; i ++ )
	{	
		if ( maxLon < src_in.lon[i] )
		{
			maxLon = src_in.lon[i];
		}
		if ( minLon > src_in.lon[i] )
		{
			minLon = src_in.lon[i];
		}

		if ( maxLat < src_in.lat[i] )
		{
			maxLat = src_in.lat[i];
		}
		if ( minLat > src_in.lat[i] )
		{
			minLat = src_in.lat[i];
		}
	}

	
	//cout << "maxLon = " << maxLon << ", minLon = " << minLon << ", maxLat = " << maxLat << ", minLat = " << minLat << endl;
	
	PJ_CONTEXT * C;
	PJ * P;
	
	C = proj_context_create( );
	
	char projStr[256];//""
	sprintf( projStr, "+proj=aeqd +lon_0=%lf +lat_0=%lf +x_0=0.0 +y_0=0.0 +ellps=WGS84", lon_0, lat_0 );

	//printf( projStr  );
	//printf( "\n"  );
	P = proj_create( C, projStr );
	if ( NULL == P )
	{
		printf( "Failed to create projection\n"  );
	}


	LonLat2cart( P, PJ_FWD, src_in );
	proj_destroy( P );
	proj_context_destroy( C );
}



long long srcCoord2Index( GRID grid, float Depth, COORD coord, SOURCE_FILE_INPUT src_in, map< long long, long long > & point_index )//map<int, POINT_INDEX> pos_pointIndex )
{
	long long i = 0;
	long long npts = src_in.npts;
	int srcX, srcY, srcZ;

	float DH = grid.DH;
	double DZ = 0.0;
	double terrain;
	

	float * x = src_in.lon;
	float * y = src_in.lat;
	float * z = src_in.coordZ;

	int originalX = grid.originalX;
	int originalY = grid.originalY;

	int frontNX = grid.frontNX;
	int frontNY = grid.frontNY;
	int frontNZ = grid.frontNZ;

	int halo = grid.halo;
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;

	int NZ = grid.NZ;

	long long  pointNum = 0;

	//cout << "x = " << x[npts-1] << ", y = " << y[npts-1] << ", z = " << z[npts-1] << endl;
	float maxCoordZ = -Depth;
	float minCoordZ = Depth;
	for ( int jj = 0; jj < npts; jj ++  )
	{
		if( maxCoordZ < z[jj] )
			maxCoordZ = z[jj];

		if( minCoordZ > z[jj] )
			minCoordZ = z[jj];
	}

	for ( int jj = 0; jj < npts; jj ++  )
	{
		z[jj] = z[jj] - maxCoordZ;
	}
	//cout << "maxCoordZ = " << maxCoordZ << ", minCoordZ = " << minCoordZ  << endl;
	
	int extendNZ = 0;
	
	//cout << npts << endl;
	int rank;

	int k = 0;

	long long index1 = 0, index2 = 0, index = 0;
	double dz1 = 0.0, dz2 = 0.0;
	for ( i = 0; i < npts; i ++ )
	{
		srcX = int( x[i] / DH + 0.5 ) + originalX - frontNX + halo; 
		srcY = int( y[i] / DH + 0.5 ) + originalY - frontNY + halo; 

#define SourceTerrain
#ifdef SourceTerrain
		DZ = coord.z[INDEX( srcX, srcY, 4)] - coord.z[INDEX( srcX, srcY, 3)];
		terrain = DZ * double( NZ - 1 ) - abs( Depth );
		//if ( z[i] > terrain )
		//{
		//	z[i] = terrain;
		//}
		z[i] = terrain + z[i];

#endif //SourceTerrain




#define NO_SOURCE_SMOOTH

#ifdef NO_SOURCE_SMOOTH
		if ( srcX >= halo && srcX < _nx && srcY >= halo && srcY < _ny )
		{
			for ( k = halo; k < _nz; k ++ )
			{
				index1 = INDEX( srcX, srcY, k - 1 );
				index2 = INDEX( srcX, srcY, k + 1 );
				index  = INDEX( srcX, srcY, k );
				
				dz1 = ( coord.z[index] - coord.z[index1] ) * 0.5;
				dz2 = ( coord.z[index2] - coord.z[index] ) * 0.5;
				
				if ( coord.z[index] - dz1 <= z[i] && z[i] <  dz2 + coord.z[index] )
				{	
					srcZ = k;
					point_index[i] = INDEX( srcX, srcY, srcZ );
					pointNum ++;					
				}
			}
		}
#else
		
		if ( srcX >= 0 && srcX < _nx_ && srcY >= 0 && srcY < _ny_ )
		{
			for ( k = 0; k < _nz_; k ++ )
			{
				if ( k - 1 == -1 )
				{	
					index  = INDEX( srcX, srcY, 0 );
					index2 = INDEX( srcX, srcY, 1 );
					dz2 = ( coord.z[index2] - coord.z[index] ) * 0.5;
					dz1 = dz2;
				}
				if ( k + 1 == _nz_ )
				{
					index1 = INDEX( srcX, srcY, _nz_ - 2 );
					index  = INDEX( srcX, srcY, _nz_ - 1 );
					dz1 = ( coord.z[index] - coord.z[index1] ) * 0.5;
					dz2 = dz1;
					
				}
				if ( k - 1 != -1 && k + 1 != _nz_ )
				{
					index1 = INDEX( srcX, srcY, k - 1 );
					index  = INDEX( srcX, srcY, k );
					index2 = INDEX( srcX, srcY, k + 1 );
					dz1 = ( coord.z[index] - coord.z[index1] ) * 0.5;
					dz2 = ( coord.z[index2] - coord.z[index] ) * 0.5;
					
				}
				
				
				if ( coord.z[index] - dz1 <= z[i] && z[i] <  dz2 + coord.z[index] )
				{	
					srcZ = k;
					point_index[i] = INDEX( srcX, srcY, srcZ );
					pointNum ++;					

				}
			}
		}
#endif // NO_SOURCE_SMOOTH



	}

	return pointNum;
}






void allocSourceParams( SOURCE_FILE_INPUT * src_in, long long pointNum  )
{
	if ( 0 == pointNum )
		return;
	int nt = src_in->nt;

	float * pTmp = NULL;
	long long num = pointNum * ( 3 + 2 * nt );
	long long size = num * sizeof( float );

	pTmp = ( float * ) malloc( size );
	memset( pTmp, 0, size );
	
	src_in->area   = pTmp + 0 * pointNum;
	src_in->strike = pTmp + 1 * pointNum;
	src_in->dip    = pTmp + 2 * pointNum;

	pTmp = pTmp + 3 * pointNum;

	src_in->rake = pTmp + 0  * pointNum;
	src_in->rate = pTmp + nt * pointNum;


}
void freeSourceParams( SOURCE_FILE_INPUT src_in, long long pointNum  )
{
	if ( 0 == pointNum )
		return;
	free( src_in.area   );
}


void readSourceParams( SOURCE_FILE_INPUT src_in,  char * sourceFileName, map<long long, long long> & point_index  )
{

	int size = point_index.size( );
	if ( 0 == size )
	{
		return;

	}
	FILE * sourceFile = fopen( sourceFileName, "rb" );

	if ( NULL == sourceFile )
		cout << sourceFileName << ", is source file" << endl;
	
	int npts = src_in.npts;
	int nt = src_in.nt;
	long long i = 0;
	long long nptsIndex; 




	int headSize = sizeof( int ) + sizeof( int ) + sizeof( int ); // npts nt dt
	long long byteSize = 6 * sizeof( float ) + //lon lat coordZ area strike dip
						 2 * sizeof( float ) * nt;

	for ( map< long long, long long>::iterator it = point_index.begin( ); it != point_index.end( ); it ++ )
	{
		nptsIndex = it->first;
		//cout << nptsIndex << endl;
		fseek( sourceFile, headSize + byteSize * nptsIndex + 3 * sizeof( float ), SEEK_SET )  ;  
		//	
		fread( src_in.area + i,   sizeof( float ), 1, sourceFile );
		fread( src_in.strike + i, sizeof( float ), 1, sourceFile );
		fread( src_in.dip + i,    sizeof( float ), 1, sourceFile );

		fread( src_in.rake + i * nt, sizeof( float ), nt, sourceFile );
		fread( src_in.rate + i * nt, sizeof( float ), nt, sourceFile );
		i ++;
	}

//	for ( i = 0; i < size; i ++ )
//	{
//		cout << "strike = " << src_in.strike[i] << endl;
//	}
	

	//int nnnn = 1000;
	//cout <<  size   << endl;

	//cout << "lon[" << nnnn  <<"] = " << src_in.lon[nnnn] << endl;
	//cout << "npts = " << npts << ", nt = " << nt << ", dt = " << src_in.dt << endl;
	fclose( sourceFile );
}

void allocMomentRate( MOMENT_RATE * momentRate, long long pointNum, int nt )
{
	if ( 0 == pointNum )
	{
		return;
	}
	long long num = pointNum * nt;
	long long size = sizeof( float ) * num * MOMENTSIZE;

	float * pMomentRate = NULL;
	CHECK( Malloc( ( void ** )&pMomentRate, size ) );
	CHECK( Memset( pMomentRate, 0, size ) ); 
	
	momentRate->Mxx = pMomentRate + 0 * num;
	momentRate->Myy = pMomentRate + 1 * num;
	momentRate->Mzz = pMomentRate + 2 * num;
	momentRate->Mxy = pMomentRate + 3 * num;
	momentRate->Mxz = pMomentRate + 4 * num;
	momentRate->Myz = pMomentRate + 5 * num;
}

void allocMomentRate_cpu( MOMENT_RATE * momentRate, long long pointNum, int nt )
{
	if ( 0 == pointNum )
		return;

	long long num = pointNum * nt;
	long long size = sizeof( float ) * num * MOMENTSIZE;

	float * pMomentRate = NULL;
	pMomentRate = ( float * )malloc( size );
	memset( pMomentRate, 0, size ); 

	momentRate->Mxx = pMomentRate + 0 * num;
	momentRate->Myy = pMomentRate + 1 * num;
	momentRate->Mzz = pMomentRate + 2 * num;
	momentRate->Mxy = pMomentRate + 3 * num;
	momentRate->Mxz = pMomentRate + 4 * num;
	momentRate->Myz = pMomentRate + 5 * num;
}

void freeMomentRate( MOMENT_RATE momentRate, long long pointNum )
{

	if ( 0 == pointNum )
		return;

	Free( momentRate.Mxx );

}
void freeMomentRate_cpu( MOMENT_RATE momentRate, long long pointNum )
{

	if ( 0 == pointNum )
		return;

	free( momentRate.Mxx );

}

void solveMomentRate( SOURCE_FILE_INPUT src_in, MOMENT_RATE momentRate, long long pointNum )
{
	if ( 0 == pointNum )
		return;
	
	float s, d, r, a, rt;
	float * strike = src_in.strike, * dip = src_in.dip, * area = src_in.area;
	float * rake = src_in.rake, * rate = src_in.rate;
		
	long long p = 0, index = 0;

	int it = 0, nt = src_in.nt;
	
	float M11 = 0.0f;
	float M22 = 0.0f;
	float M33 = 0.0f;
	float M12 = 0.0f;
	float M13 = 0.0f;
	float M23 = 0.0f;


	float * Mxx = momentRate.Mxx;
	float * Myy = momentRate.Myy;
	float * Mzz = momentRate.Mzz;
	float * Mxy = momentRate.Mxy;
	float * Mxz = momentRate.Mxz;
	float * Myz = momentRate.Myz;

	for ( p = 0; p < pointNum; p ++ )
	{
		s = strike[p] * DEGREE2RADIAN;
		d = dip[p] * DEGREE2RADIAN;
		a = area[p];

		for ( it = 0; it < nt; it ++ )
		{
			index = it + p * nt;
			r = rake[index] * DEGREE2RADIAN;
			rt = rate[index];

			/*
			M11 = -(sin(d) * cos(r) * sin(s * 2.0) + sin(d * 2.0) * sin(r) * sin(s) * sin(s));   
			M22 = sin(d) * cos(r) * sin(s * 2.0) - sin(d * 2.0) * sin(r) * cos(s) * cos(s);    
			M33 = -(M11+M22);                                                                                 
			M12 = sin(d) * cos(r) * cos(s * 2.0) + 0.5 * sin(d * 2.0) * sin(r) * sin(s * 2.0); 
			M13 = -(cos(d) * cos(r) * cos(s) + cos(d * 2.0) * sin(r) * sin(s));                   
			M23 = -(cos(d) * cos(r) * sin(s) - cos(d * 2.0) * sin(r) * cos(s));                   

			Mxx[index] =  M22 * a * rt;
			Myy[index] =  M11 * a * rt;
			Mzz[index] =  M33 * a * rt;
			Mxy[index] =  M12 * a * rt;
			Mxz[index] = -M23 * a * rt;
			Myz[index] = -M13 * a * rt;
			*/

			M11 = -( sin(d)* cos(r)* sin(2.0*s) +  sin(2.0*d)* sin(r)* sin(s)* sin(s));
			M22 =  sin(d)* cos(r)* sin(2.0*s) -  sin(2.0*d)* sin(r)* cos(s)* cos(s);
			M33 = -(M11+M22);
			M12 =  sin(d)* cos(r)* cos(2.0*s) + 0.5* sin(2.0*d)* sin(r)* sin(2.0*s);
			M13 = -( cos(d)* cos(r)* cos(s) +  cos(2.0*d)* sin(r)* sin(s));
			M23 = -( cos(d)* cos(r)* sin(s) -  cos(2.0*d)* sin(r)* cos(s));

			Mxx[index] =   M22 * a * rt;
			Myy[index] =   M11 * a * rt;
			Mzz[index] =   M33 * a * rt;
			Mxy[index] =   M12 * a * rt;
			Mxz[index] = - M23 * a * rt;
			Myz[index] = - M13 * a * rt;

		}
	}
	
	
/*
    for (int i = 0; i < nt; i++)
    {
        M11[i] = -(sin(d) * cos(r[i]) * sin(s * 2.0) + sin(d * 2.0) * sin(r[i]) * sin(s) * sin(s)) * area * rate[i];
        M22[i] = (sin(d) * cos(r[i]) * sin(s * 2.0) - sin(d * 2.0) * sin(r[i]) * cos(s) * cos(s)) * area * rate[i];
        M33[i] = (-(M11[i] + M22[i]));
        M12[i] = (sin(d) * cos(r[i]) * cos(s * 2.0) + 0.5 * sin(d * 2.0) * sin(r[i]) * sin(s * 2.0)) * area * rate[i];
        M13[i] = (cos(d) * cos(r[i]) * cos(s) + cos(d * 2.0) * sin(r[i]) * sin(s)) * area * rate[i];
        M23[i] = (cos(d) * cos(r[i]) * sin(s) - cos(d * 2.0) * sin(r[i]) * cos(s)) * area * rate[i];
    }
*/
}




void dealDuplicateIndex( MOMENT_RATE momentRate, 
						 map< long long, long long > & point_index, 
						 map< long long, long long > & index_point, 
						 SOURCE_FILE_INPUT src_in )
{

	long long point_num = point_index.size( );
	if ( 0 == point_num  )
	{
		return;

	}
	long long pnt = 0, idx = 0;

	//map < long long, long long > index_point;

	int i = 0;
	for ( map < long long, long long >::iterator it = point_index.begin( ); it != point_index.end( ); it ++ )
	{	
		idx = it->second;
		index_point[idx] = i;//It means that the hashmap index_point stores the max pnt since pnt is in a high order.
		i ++;	
	}


	int t = 0;
	int nt = src_in.nt;
	int ii, indexII, indexI;

	float * Mxx = momentRate.Mxx;
	float * Myy = momentRate.Myy;
	float * Mzz = momentRate.Mzz;
	float * Mxy = momentRate.Mxy;
	float * Mxz = momentRate.Mxz;
	float * Myz = momentRate.Myz;


	i = 0;
	for ( map < long long, long long >::iterator it = point_index.begin( ); it != point_index.end( ); it ++ )
	{	
		idx = it->second;
		ii = index_point[idx];
		if ( ii > i )
		{
			for ( t = 0; t < nt; t ++ )
			{
				indexI =  t + i  * nt;
				indexII = t + ii * nt;

				Mxx[indexII] += Mxx[indexI]; 
				Myy[indexII] += Myy[indexI]; 
				Mzz[indexII] += Mzz[indexI]; 
				Mxy[indexII] += Mxy[indexI]; 
				Mxz[indexII] += Mxz[indexI]; 
				Myz[indexII] += Myz[indexI]; 
			}
		}
		i ++;

	}

	
	//cout << "point_index size = " << point_index.size( )  << ", index_point size = " << index_point.size()  << endl;

}

void outputSourceData( PARAMS params,  SOURCE_FILE_INPUT src_in,
					   map< long long, long long> & index_point, 
					   MOMENT_RATE momentRate, MPI_COORD thisMPICoord )
{

	
	int size = index_point.size( );
	if ( 0 == size )
	{
		return;

	}
	char fileName[256];
	sprintf( fileName, "%s/source_mpi_%d_%d_%d.bin", params.OUT,  thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
	
	FILE * file = fopen( fileName, "wb" );
	
	if ( NULL == file )
		printf( "The file %s can not be opened\n", fileName );
	
	int npts = src_in.npts;
	int nt = src_in.nt;
	float dt = src_in.dt;

	fwrite( &size, sizeof( int ), 1, file );
	fwrite( &nt, sizeof( int ), 1, file );
	fwrite( &dt, sizeof( float ), 1, file );

	long long index = 0;
	long long ii = 0;
	long long pos = 0;

	float * Mxx = momentRate.Mxx;
	float * Myy = momentRate.Myy;
	float * Mzz = momentRate.Mzz;
	float * Mxy = momentRate.Mxy;
	float * Mxz = momentRate.Mxz;
	float * Myz = momentRate.Myz;

	for( map<long long, long long>::iterator it = index_point.begin( ); it != index_point.end( ); it ++ )
	{
		index = it->first;
		ii = it->second;
		pos = ii * nt;
		fwrite( &index, sizeof( long long ), 1, file);
		fwrite( Mxx + pos, sizeof( float ), nt, file);
		fwrite( Myy + pos, sizeof( float ), nt, file);
		fwrite( Mzz + pos, sizeof( float ), nt, file);
		fwrite( Mxy + pos, sizeof( float ), nt, file);
		fwrite( Mxz + pos, sizeof( float ), nt, file);
		fwrite( Myz + pos, sizeof( float ), nt, file);
			
	}
	

	fclose( file );



}

void verifyLocation( PARAMS params, SOURCE_FILE_INPUT src_in )
{
	FILE * file[3];

	char fileNameX[256], fileNameY[256], fileNameZ[256];
	
	sprintf( fileNameX, "%s/source_coord_X.bin", params.OUT  );
	sprintf( fileNameY, "%s/source_coord_Y.bin", params.OUT  );
	sprintf( fileNameZ, "%s/source_coord_Z.bin", params.OUT  );

	file[0] = fopen( fileNameX, "wb" );
	file[1] = fopen( fileNameY, "wb" );
	file[2] = fopen( fileNameZ, "wb" );

	long long npts = src_in.npts;

	float * x = src_in.lon;
	float * y = src_in.lat;
	float * z = src_in.coordZ;



	fwrite( x, sizeof( float ), npts, file[0] );
	fwrite( y, sizeof( float ), npts, file[1] );
	fwrite( z, sizeof( float ), npts, file[2] );

	fclose( file[0] );
	fclose( file[1] );
	fclose( file[2] );

}

void allocSrcIndex( long long ** srcIndex, long long npts )
{
	if ( npts == 0 )
		return;
		
	long long size = sizeof( long long ) * npts;
	
	CHECK( Malloc( ( void ** )srcIndex, size ) );
	CHECK( Memset( *srcIndex, 0, size ) ); 
}


void freeSrcIndex( long long * srcIndex, long long npts )
{
	if ( npts == 0 )
		return;
	Free( srcIndex );
}

void allocSrcIndex_cpu( long long ** srcIndex, long long npts )
{
	if ( npts == 0 )
		return;
		
	long long size = sizeof( long long ) * npts;
	
	*srcIndex = ( long long * )malloc( size );
	memset( *srcIndex, 0, size ); 
}


void freeSrcIndex_cpu( long long * srcIndex, long long npts )
{
	if ( npts == 0 )
		return;
	free( srcIndex );
}

void changeStorageOrder( SOURCE_FILE_INPUT src_in, 
					     map< long long, long long> & index_point, 
					     long long * srcIndex,
					     MOMENT_RATE momentRateOld, 
					     MOMENT_RATE momentRateNew )
{
	int nt = src_in.nt;
	int npts = index_point.size( );
	if ( 0 == npts )
	{
		return;
	}


	float * Mxx0 = momentRateOld.Mxx;
	float * Myy0 = momentRateOld.Myy;
	float * Mzz0 = momentRateOld.Mzz;
	float * Mxy0 = momentRateOld.Mxy;
	float * Mxz0 = momentRateOld.Mxz;
	float * Myz0 = momentRateOld.Myz;

	float * Mxx1 = momentRateNew.Mxx;
	float * Myy1 = momentRateNew.Myy;
	float * Mzz1 = momentRateNew.Mzz;
	float * Mxy1 = momentRateNew.Mxy;
	float * Mxz1 = momentRateNew.Mxz;
	float * Myz1 = momentRateNew.Myz;

	int t = 0;
	long long ii = 0;
	long long pos0 = 0, pos1 = 0;

	long long J = 0;
	for( map<long long, long long>::iterator it = index_point.begin( ); it != index_point.end( ); it ++ )
	{
		srcIndex[J] = it->first;
		ii = it->second;
		for ( t = 0; t < nt; t ++ )
		{
			pos0 = ii * nt + t;
			pos1 = J + t * npts;
			Mxx1[pos1] = Mxx0[pos0];
			Myy1[pos1] = Myy0[pos0];
			Mzz1[pos1] = Mzz0[pos0];
			Mxy1[pos1] = Mxy0[pos0];
			Mxz1[pos1] = Mxz0[pos0];
			Myz1[pos1] = Myz0[pos0];
		}
		J ++;
	}
}

void allocMomentRateSlice( MOMENT_RATE * momentRateSlice, long long npts )
{
	if ( 0 == npts )
		return;

	long long num = npts;
	long long size = sizeof( float ) * num * MOMENTSIZE;

	float * pMomentRate = NULL;
	CHECK( Malloc( ( void ** )&pMomentRate, size ) );
	CHECK( Memset( pMomentRate, 0, size ) ); 

	momentRateSlice->Mxx = pMomentRate + 0 * num;
	momentRateSlice->Myy = pMomentRate + 1 * num;
	momentRateSlice->Mzz = pMomentRate + 2 * num;
	momentRateSlice->Mxy = pMomentRate + 3 * num;
	momentRateSlice->Mxz = pMomentRate + 4 * num;
	momentRateSlice->Myz = pMomentRate + 5 * num;

}

void freeMomentRateSlice( MOMENT_RATE momentRateSlice, long long npts  )
{
	if ( 0 == npts )
		return;
	Free( momentRateSlice.Mxx );
}

//This function includes allocating srcIndex and momentRate memory.
void init_MultiSource( 
		PARAMS params, GRID grid, MPI_COORD thisMPICoord, COORD coord, 
		long long ** srcIndex, MOMENT_RATE * momentRate, MOMENT_RATE * momentRateSlice, SOURCE_FILE_INPUT * ret_src_in )
{
	int thisRank;
	MPI_Barrier( MPI_COMM_WORLD );
	MPI_Comm_rank( MPI_COMM_WORLD, &thisRank );
	if ( thisRank == 0 )
		printf( "Precessing fualt(Source data precessing) data...\n"  );
	MPI_Barrier( MPI_COMM_WORLD );



	SOURCE_FILE_INPUT src_in;
	char sourceFileName[256];
	sprintf( sourceFileName, "%s/%s", params.sourceDir, params.sourceFile );
	
	//cout << sourceFileName << ", is source file" << endl;

	readSourceInfo( &src_in, sourceFileName );
	allocSourceLocation( &src_in ); //A1
	readSourceLocation( src_in, sourceFileName ); 	

	
	projTrans( params, src_in );	
	


	//allocSrcIndex( srcIndex, src_in.npts );
	
	float Depth = params.Depth * 1000;

	map< long long, long long >  point_index;
	
	map<long long, long long >  index_point;

	long long pointNum = srcCoord2Index( grid, Depth, coord, src_in, point_index );
	
	if( 0 == thisMPICoord.X && 0 == thisMPICoord.Y && 0 == thisMPICoord.Z )
	{
		verifyLocation( params, src_in );
	}
	

	freeSourceLocation( src_in ); //F1
	MPI_Barrier( MPI_COMM_WORLD );

	allocSourceParams( &src_in, pointNum ); //A2
	readSourceParams( src_in, sourceFileName, point_index );
/**/
	MOMENT_RATE momentRateOld;
	allocMomentRate_cpu( &momentRateOld, pointNum, src_in.nt ); //A3

	solveMomentRate( src_in, momentRateOld, pointNum );
	freeSourceParams( src_in, pointNum ); //F2

	dealDuplicateIndex( momentRateOld, point_index, index_point, src_in );
	point_index.clear( );
	outputSourceData( params, src_in, index_point, momentRateOld, thisMPICoord );



	long long npts = index_point.size( );
	int nt = src_in.nt;
	MOMENT_RATE momentRateNew;
	allocMomentRate_cpu( &momentRateNew, npts, nt ); //A4

	long long * srcIndexNew;
	allocSrcIndex_cpu( &srcIndexNew, npts ); //A5
	changeStorageOrder( src_in, index_point, srcIndexNew, momentRateOld, momentRateNew );
	freeMomentRate_cpu( momentRateOld, pointNum ); //F3

	index_point.clear( );
	allocMomentRateSlice( momentRateSlice, npts );
#ifdef GPU_CUDA
	allocMomentRate( momentRate, npts, nt );
	allocSrcIndex( srcIndex, npts );
	long long mr_size = npts * nt * sizeof( float ) * MOMENTSIZE;
	CHECK( cudaMemcpy( momentRate->Mxx, momentRateNew.Mxx, mr_size, cudaMemcpyHostToDevice ));

	long long si_size = npts * sizeof( long long );
	CHECK( cudaMemcpy( *srcIndex, srcIndexNew, si_size, cudaMemcpyHostToDevice ));

	freeMomentRate_cpu( momentRateNew, npts ); //F4
	freeSrcIndex_cpu( srcIndexNew, npts ); //F5

#else
	//memcpy( *momentRate, momentRateNew, sizeof( MOMENT_RATE ));
	momentRate->Mxx = momentRateNew.Mxx;
	momentRate->Myy = momentRateNew.Myy;
	momentRate->Mzz = momentRateNew.Mzz;
	momentRate->Mxy = momentRateNew.Mxy;
	momentRate->Mxz = momentRateNew.Mxz;
	momentRate->Myz = momentRateNew.Myz;

	*srcIndex = srcIndexNew;
#endif
    ret_src_in->nt = src_in.nt;
    ret_src_in->dt = src_in.dt;
    ret_src_in->npts = npts;
	//MPI_Barrier( MPI_COMM_WORLD );
	//cout << "npts = "<< npts << endl;
	//cout << "1111111111111111111111111111111111111111111"<< endl;

}



void finish_MultiSource( long long * srcIndex, MOMENT_RATE momentRate, MOMENT_RATE momentRateSlice, long long npts )
{
	freeSrcIndex( srcIndex, npts );
	freeMomentRate( momentRate, npts );	
	freeMomentRateSlice( momentRateSlice, npts  );
}


__GLOBAL__
void calculate_MomentRate( long long npts, long long nt, MEDIUM medium, float * Jac, MOMENT_RATE momentRate, long long * srcIndex, float DH )
{

#ifdef GPU_CUDA
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
	long long j = threadIdx.y + blockIdx.y * blockDim.y;
#else
	long long i = 0;
	long long j = 0;
#endif
	long long idx = 0;
	long long pos = 0;
	
	float V = 1.0;

	CALCULATE2D( i, j, 0, npts, 0, nt )
		idx = srcIndex[i];
		pos = i + j * npts;
		V = Jac[idx] * DH * DH * DH;
		V = 1.0 / V;
		V = V * Cv;
		momentRate.Mxx[pos] *= medium.mu[idx] * V;
		momentRate.Myy[pos] *= medium.mu[idx] * V;
		momentRate.Mzz[pos] *= medium.mu[idx] * V;
		momentRate.Mxy[pos] *= medium.mu[idx] * V;
		momentRate.Mxz[pos] *= medium.mu[idx] * V;
		momentRate.Myz[pos] *= medium.mu[idx] * V;
	END_CALCULATE2D( )


}

void calculateMomentRate( SOURCE_FILE_INPUT src_in, MEDIUM medium, float * Jac, MOMENT_RATE momentRate, long long * srcIndex, float DH )
{

	long long npts = src_in.npts;
	long long nt = src_in.nt;

	if ( 0 == npts )
		return;

	
#ifdef GPU_CUDA
	dim3 threads( 16, 16, 1 );
	dim3 blocks;
	blocks.x = ( npts + threads.x - 1 ) / threads.x;
	blocks.y = ( nt + threads.y - 1 ) / threads.y;
	blocks.z = 1;
	calculate_MomentRate
	<<< blocks, threads >>>
	( npts, nt, medium, Jac, momentRate, srcIndex, DH );

#else
	calculate_MomentRate
	( npts, nt, medium, Jac, momentRate, srcIndex, DH );
#endif



}



__GLOBAL__
void interp_momentRate( long long npts, 
					   MOMENT_RATE momentRate, 
					   MOMENT_RATE momentRateSlice, 
					   float t_weight, int srcIt )
{
#ifdef GPU_CUDA
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
#else
	long long i = 0;
#endif
	long long pos0 = 0, pos1 = 0;

	CALCULATE1D( i, 0, npts )
	//for ( i = 0; i < npts; i ++ )
		pos0 = srcIt * npts + i;
		pos1 = ( srcIt + 1 ) * npts + i;
		momentRateSlice.Mxx[i] = ( momentRate.Mxx[pos1] - momentRate.Mxx[pos0] ) * t_weight + momentRate.Mxx[pos0];
		momentRateSlice.Myy[i] = ( momentRate.Myy[pos1] - momentRate.Myy[pos0] ) * t_weight + momentRate.Myy[pos0];
		momentRateSlice.Mzz[i] = ( momentRate.Mzz[pos1] - momentRate.Mzz[pos0] ) * t_weight + momentRate.Mzz[pos0];
		momentRateSlice.Mxy[i] = ( momentRate.Mxy[pos1] - momentRate.Mxy[pos0] ) * t_weight + momentRate.Mxy[pos0];
		momentRateSlice.Mxz[i] = ( momentRate.Mxz[pos1] - momentRate.Mxz[pos0] ) * t_weight + momentRate.Mxz[pos0];
		momentRateSlice.Myz[i] = ( momentRate.Myz[pos1] - momentRate.Myz[pos0] ) * t_weight + momentRate.Myz[pos0];
	END_CALCULATE1D( )
}


__GLOBAL__
void addSource( WAVE hW,  
				MOMENT_RATE momentRateSlice, 
				long long * srcIndex, int npts,
				int gaussI, int gaussJ, int gaussK, float factorGauss, int _nx_, int _ny_, int _nz_, int flagSurf )
{
	
#ifdef GPU_CUDA
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
#else
	long long i = 0;
#endif
	long long idx = 0;
	//float V = 1.0;

	int Z = 0;
	CALCULATE1D( i, 0, npts )
		idx = srcIndex[i];

		if ( flagSurf == 1 )
		{
			Z = idx / _nx_ / _ny_;
			if ( Z + gaussK > _nz_ - 4 )
				factorGauss = 0.0;

			if ( ( Z == _nz_ - 4 && gaussK <  0 ) ||
				 ( Z == _nz_ - 5 && gaussK < -1 ) ||
				 ( Z == _nz_ - 6 && gaussK < -2 ) )
			{
				factorGauss = factorGauss * 2;
			}
			
		}
		//V = Jac[idx] * DH * DH * DH;
		//V = -1.0 / V;
		hW.Txx[idx] -= momentRateSlice.Mxx[i] * factorGauss;
		hW.Tyy[idx] -= momentRateSlice.Myy[i] * factorGauss;
		hW.Tzz[idx] -= momentRateSlice.Mzz[i] * factorGauss;
		hW.Txy[idx] -= momentRateSlice.Mxy[i] * factorGauss;
		hW.Txz[idx] -= momentRateSlice.Mxz[i] * factorGauss;
		hW.Tyz[idx] -= momentRateSlice.Myz[i] * factorGauss;
	END_CALCULATE1D( ) 
	
}


__GLOBAL__
void addSource1( WAVE hW,  
				MOMENT_RATE momentRateSlice, 
				long long * srcIndex, int npts )
{
	
#ifdef GPU_CUDA
	long long i = threadIdx.x + blockIdx.x * blockDim.x;
#else
	long long i = 0;
#endif
	long long idx = 0;
	//float V = 1.0;

	CALCULATE1D( i, 0, npts )
		idx = srcIndex[i];
		//V = Jac[idx] * DH * DH * DH;
		//V = -1.0 / V;
		hW.Txx[idx] -= momentRateSlice.Mxx[i];
		hW.Tyy[idx] -= momentRateSlice.Myy[i];
		hW.Tzz[idx] -= momentRateSlice.Mzz[i];
		hW.Txy[idx] -= momentRateSlice.Mxy[i];
		hW.Txz[idx] -= momentRateSlice.Mxz[i];
		hW.Tyz[idx] -= momentRateSlice.Myz[i];
	END_CALCULATE1D( ) 
	
}


void addMomenteRate(   
		GRID grid, 
		SOURCE_FILE_INPUT src_in, 
		WAVE hW, float * Jac, 
		long long * srcIndex, MOMENT_RATE momentRate, MOMENT_RATE momentRateSlice, 
		int it, int irk, float DT, float DH, float * gaussFactor, int nGauss, int flagSurf )

{

	MPI_Barrier( MPI_COMM_WORLD );

		
	long long npts = src_in.npts;
	int nt = src_in.nt;
	float dt = src_in.dt;

	float tmpT, t1, t_weight;
	
	if( 0 == irk ) 
	{
		tmpT = ( it + 0.0f ) * DT;
	}else if ( 1 == irk || 2 == irk )
	{
		tmpT = ( it + 0.5f ) * DT;
	}
	else if ( 3 == irk )
	{
		tmpT = ( it + 1.0f ) * DT;
	}
	
	
	int srcIt = int( tmpT / dt );

	if ( ( srcIt + 1 ) >= nt )
	{

		//printf( "srcIt = %d\n", srcIt );
		finish_MultiSource( srcIndex, momentRate, momentRateSlice, src_in.npts );
		return;
	}

	t1 = float( srcIt ) * dt;
	t_weight = ( tmpT - t1 ) / dt;

	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;
	int gaussI, gaussJ, gaussK, gPos;
	int lenGauss = nGauss * 2 + 1;
	float factorGauss = 0.0f;
#ifdef GPU_CUDA
	long long num = npts;
	dim3 threads( 256, 1, 1 );
	dim3 blocks;
	blocks.x = ( num + threads.x - 1 ) / threads.x;
	blocks.y = 1;
	blocks.z = 1;
	CHECK( cudaDeviceSynchronize( ) );
	interp_momentRate <<< blocks, threads >>>
	( npts, momentRate, momentRateSlice, t_weight, srcIt );
	CHECK( cudaDeviceSynchronize( ) );
#ifdef NO_SOURCE_SMOOTH
	addSource1  <<< blocks, threads >>>
	( hW, momentRateSlice, srcIndex, npts );
	CHECK( cudaDeviceSynchronize( ) );
	
#else

	for( gaussK = - nGauss ; gaussK < nGauss + 1; gaussK ++ )
	{
		for( gaussJ = - nGauss; gaussJ < nGauss + 1; gaussJ ++ )
		{
			for( gaussI = - nGauss; gaussI < nGauss + 1; gaussI ++ )
			{
				gPos = ( gaussI + nGauss ) + ( gaussJ + nGauss ) * lenGauss + ( gaussK + nGauss ) * lenGauss * lenGauss;
				factorGauss = gaussFactor[gPos];
				addSource  <<< blocks, threads >>>
				( hW, momentRateSlice, srcIndex, npts, gaussI, gaussJ, gaussK, 
				_nx_, _ny_, _nz_, factorGauss, flagSurf );
				CHECK( cudaDeviceSynchronize( ) );
			}
		}
	}
#endif


#else
	interp_momentRate 
	( npts, momentRate, momentRateSlice, t_weight, srcIt );
#ifdef NO_SOURCE_SMOOTH
	addSource1 
	( hW, momentRateSlice, srcIndex, npts );
	CHECK( cudaDeviceSynchronize( ) );
	
#else

	for( gaussK = - nGauss ; gaussK < nGauss + 1; gaussK ++ )
	{
		for( gaussJ = - nGauss; gaussJ < nGauss + 1; gaussJ ++ )
		{
			for( gaussI = - nGauss; gaussI < nGauss + 1; gaussI ++ )
			{
				gPos = ( gaussI + nGauss ) + ( gaussJ + nGauss ) * lenGauss + ( gaussK + nGauss ) * lenGauss * lenGauss;
				factorGauss = gaussFactor[gPos];
				addSource 
				( hW, momentRateSlice, srcIndex, npts, gaussI, gaussJ, gaussK, 
				_nx_, _ny_, _nz_, factorGauss, flagSurf );
				CHECK( cudaDeviceSynchronize( ) );
			}
		}
	}
#endif


#endif
	MPI_Barrier( MPI_COMM_WORLD );

}

