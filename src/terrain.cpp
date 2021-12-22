/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:terrain.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-05
*   Discription:
*
================================================================*/
#include "header.h"

#define degreePerBlcokX 5
#define degreePerBlcokY 5

#define pointsPerBlockX 6000
#define pointsPerBlockY 6000
#define extendPoint 50


#define errorNegativeElevation -32768
#define errorPositiveElevation 32768
#define ERROR_ELEVATION

void calculateCutLonLat( GRID grid, PJ * P,  PJ_DIRECTION PD, PJ_COORD pj_LonLat[4] )
{		
	
	int _originalX = grid._originalX;
	int _originalY = grid._originalY;

	int _NX_ = grid._NX_;
	int _NY_ = grid._NY_;

	


	double leftDis  = - ( _originalX + extendPoint ) * grid.DH;
	double rightDis =   ( ( _NX_ - _originalX ) + extendPoint ) * grid.DH;
	
	double downDis  = - ( _originalY + extendPoint ) * grid.DH;
	double upDis    =   ( ( _NY_ - _originalY ) + extendPoint ) * grid.DH; 
	
	//PJ_XY xy = { leftDis, downDis };
	//printf( "leftDis: %f, downDis: %f\n", leftDis, downDis );

	PJ_XY leftDownXY, leftUpXY, rightDown, rightUpXY;

	leftDownXY.x = leftDis;	leftDownXY.y = downDis;
	leftUpXY.x   = leftDis;	leftUpXY.y   = upDis;

	rightDown.x = rightDis;		rightDown.y  = downDis;
	rightUpXY.x   = rightDis;	rightUpXY.y	 = upDis;

	PJ_XY xy[4] = { leftDownXY, leftUpXY, rightDown, rightUpXY};
	
	int i = 0;
	for ( i = 0; i < 4; i ++ )
		pj_LonLat[i].xy = xy[i];

	//PJ_COORD leftDownLonLat, leftUpLonLat, rightDownLonLat, rightUpLonLat;
	
	proj_trans_array( P, PD, 4, pj_LonLat );
	
	for ( i = 0; i < 4; i ++ )
		printf( "longitude: %lf, latitude: %lf\n", pj_LonLat[i].lp.lam * RADIAN2DEGREE, pj_LonLat[i].lp.phi * RADIAN2DEGREE );
	
		
}

//Eastern Hemisphere
void readSRTM90( PARAMS params, float * totalTerrain )
{
	int lon, lat;

	int blockX = params.blockX;	
	int blockY = params.blockY;

	int lonStart = params.lonStart;
	int latStart = params.latStart;
	
	int lonEnd = lonStart + ( blockX - 1 ) * degreePerBlcokX;
	int latEnd = latStart + ( blockY - 1 ) * degreePerBlcokY;


	char srtm90FileName[512];

	

	int totalPointX = ( pointsPerBlockX - 1 ) * blockX + 1;
	int totalPointY = ( pointsPerBlockY - 1 ) * blockY + 1;


	float * oneBlockTerrain = ( float * )malloc( pointsPerBlockX * pointsPerBlockY * sizeof( float ) );
	memset( oneBlockTerrain, 0, pointsPerBlockX * pointsPerBlockY * sizeof( float ) );
	
	int i = 0, j = 0;

//#ifdef LON_FAST

	for ( lat = latStart; lat <= latEnd; lat += 5 )
	{
		
		i = 0;
		for ( lon = lonStart; lon <= lonEnd; lon += 5 )
		{

			memset( srtm90FileName, 0, sizeof( char ) * 512 );
			sprintf( srtm90FileName, "%s/srtm_%dN%dE.bin", params.TerrainDir, lat, lon );
			FILE * fp = fopen( srtm90FileName, "rb"  );
			if ( NULL == fp )
			{	
				printf( "There is no such file %s\n", srtm90FileName );
				exit( 1 );
			}
			int startI = i * ( pointsPerBlockX - 1);	
			int startJ = j * ( pointsPerBlockY - 1);	
			
			int endX = ( i + 1 ) * ( pointsPerBlockX - 1 ) + 1;
			int endY = ( j + 1 ) * ( pointsPerBlockY - 1 ) + 1;
			int ii, jj;
			long long pos, index;
			//int pos = ( i + j * blockX ) * pointsPerBlockX * pointsPerBlockY;
			//fread( &totalTerrain[pos], sizeof( float ), pointsPerBlockX * pointsPerBlockY,fp);			
			fread( oneBlockTerrain, sizeof( float ), pointsPerBlockX * pointsPerBlockY, fp );
			//for ( jj = 0; jj < ( pointsPerBlockY - 1 ); jj ++  )
			//{
			//	for ( ii = 0; ii < ( pointsPerBlockX - 1 ); ii ++ )
			//	{
			//		index = ii + jj * pointsPerBlockX;
			//		pos = ( startI + ii ) + ( startJ + jj ) * totalPointX;
			//		totalTerrain[pos] = oneBlockTerrain[index];
			//	}
			//}

			FOR_LOOP2D( ii, jj, 0, pointsPerBlockX - 1, 0, pointsPerBlockY - 1 ) 

					index = Index2D( ii, jj, pointsPerBlockX, pointsPerBlockY );
					//pos = ( startI + ii ) + ( startJ + jj ) * totalPointX;
					pos = Index2D( startI + ii, startJ + jj,  totalPointX, totalPointY );
					totalTerrain[pos] = oneBlockTerrain[index];
			END_LOOP2D( )

			//for ( int j = 0; j < pointsPerBlockY; j ++  )
			//{
			//	printf( "\n"  );
			//	for ( int i = 0; i < pointsPerBlockX; i ++ )
			//	{
			//		printf( "%e ", totalTerrain[i + j * pointsPerBlockX]  );
			//	}
			//}
			fclose( fp );
			i ++;
		}
		j ++;
	}


//#else
    
/*
	for ( lon = lonStart; lon <= lonEnd; lon += 5 )
	{
		
		j = 0;
	    for ( lat = latStart; lat <= latEnd; lat += 5 )
		{

			memset( srtm90FileName, 0, sizeof( char ) * 512 );
			sprintf( srtm90FileName, "%s/srtm_%dN%dE.bin", params.TerrainDir, lat, lon );
			FILE * fp = fopen( srtm90FileName, "rb"  );
			if ( NULL == fp )
			{	
				printf( "There is no such file %s\n", srtm90FileName );
				exit( 1 );
			}
			int startI = i * ( pointsPerBlockX - 1);	
			int startJ = j * ( pointsPerBlockY - 1);	
			
			int endX = ( i + 1 ) * ( pointsPerBlockX - 1 ) + 1;
			int endY = ( j + 1 ) * ( pointsPerBlockY - 1 ) + 1;
			int ii, jj;
			long long pos, index;
			//int pos = ( i + j * blockX ) * pointsPerBlockX * pointsPerBlockY;
			//fread( &totalTerrain[pos], sizeof( float ), pointsPerBlockX * pointsPerBlockY,fp);			
			fread( oneBlockTerrain, sizeof( float ), pointsPerBlockX * pointsPerBlockY, fp );

			FOR_LOOP2D( ii, jj, 0, pointsPerBlockX - 1, 0, pointsPerBlockY - 1 ) 

					index = Index2D( ii, jj, pointsPerBlockX, pointsPerBlockY );
					//pos = ( startI + ii ) + ( startJ + jj ) * totalPointX;
					pos = Index2D( startI + ii, startJ + jj,  totalPointX, totalPointY );
					totalTerrain[pos] = oneBlockTerrain[index];
			END_LOOP2D( )


			//for ( int j = 0; j < pointsPerBlockY; j ++  )
			//{
			//	printf( "\n"  );
			//	for ( int i = 0; i < pointsPerBlockX; i ++ )
			//	{
			//		printf( "%e ", totalTerrain[i + j * pointsPerBlockX]  );
			//	}
			//}
			fclose( fp );
			j ++;
		}
		i ++;
	}
*/

	
//#endif //LON_FAST

	free( oneBlockTerrain );
	//printf( "i = %d, j = %d\n", i, j );

}

void resetTotalTerrain( PARAMS params, float * totalTerrain  )
{
	int totalPointX = ( pointsPerBlockX - 1 ) * params.blockX + 1;
	int totalPointY = ( pointsPerBlockY - 1 ) * params.blockY + 1;
	int i = 0, j = 0;
	long long index;

	FOR_LOOP2D( i, j, 0, totalPointX, 0, totalPointY )
        index = Index2D( i, j, totalPointX, totalPointY );
		if ( totalTerrain[index] < errorNegativeElevation + 1.0 || totalTerrain[index] > errorPositiveElevation - 1.0 ) 
		{
			totalTerrain[index] = 0.0;
		}
	END_LOOP2D( )
    

}


void callibrateTotalTerrain( PARAMS params, float * totalTerrain )
{

	int totalPointX = ( pointsPerBlockX - 1 ) * params.blockX + 1;
	int totalPointY = ( pointsPerBlockY - 1 ) * params.blockY + 1;
	int i = 0, j = 0;
	long long index, index1, index2;

	float * gradTotalTerrain = ( float * )malloc( sizeof( float ) * totalPointX * totalPointY );
	float * newTotalTerrain = ( float * )malloc( sizeof( float ) * totalPointX * totalPointY );
	float gradX = 0.0, gradY = 0.0;
	
	long long cnt = 0;

	float a, b, c, d;

	memcpy( newTotalTerrain, totalTerrain, sizeof( float ) * totalPointX * totalPointY );

	FOR_LOOP2D( i, j, 1, totalPointX - 1, 1, totalPointY - 1 )
        index = Index2D( i, j, totalPointX, totalPointY );

        index1 = Index2D( i - 1, j, totalPointX, totalPointY );
        index2 = Index2D( i + 1, j, totalPointX, totalPointY );
		a = totalTerrain[index2];
		b = totalTerrain[index1];
		gradX = a - b;

        index1 = Index2D( i, j - 1, totalPointX, totalPointY );
        index2 = Index2D( i, j + 1, totalPointX, totalPointY );
		c = totalTerrain[index2];
		d = totalTerrain[index1];
		gradY = c - d;

		gradTotalTerrain[index] = sqrt( gradX * gradX + gradY * gradY );
		
		if ( gradTotalTerrain[index] > 100.0 )
		{
			newTotalTerrain[index] = ( a + b + c + d ) * 0.25;
		}
		

	END_LOOP2D( )

	int thisRank;
	MPI_Comm_rank( MPI_COMM_WORLD, &thisRank );
//	if ( 0 == thisRank  )
//	{
//		char gradFileName[256];
//		sprintf( gradFileName, "%s/gradTotalTerrain.bin", params.TerrainDir );
//		FILE * gradTerrainFile = fopen( gradFileName, "wb"  );
//		fwrite( gradTotalTerrain, sizeof( float ), totalPointX * totalPointY, gradTerrainFile );
//		fclose( gradTerrainFile );
//	}
	
	memcpy( totalTerrain, newTotalTerrain, sizeof( float ) * totalPointX * totalPointY );

	free( gradTotalTerrain );
	free( newTotalTerrain );

}


void cart2LonLat( GRID grid, PJ * P,  PJ_DIRECTION PD, COORD coord, LONLAT LonLat )
{
	PJ_COORD * pj_coord;

	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	int _nz  = grid._nz;


	long long num = _nx_ * _ny_; 
	

	int i = 0, j = 0;
	long long index = 0, pos = 0;
	int k = _nz - 1;  

	pj_coord = ( PJ_COORD * )malloc( sizeof( PJ_COORD ) * num );

	FOR_LOOP2D( i, j, 0, _nx_, 0, _ny_ )
			index = INDEX( i, j, k );//i + j * _nx_ + k *_nx_ * _ny_;
			pos = Index2D( i, j, _nx_, _ny_ );
			pj_coord[pos].xy.x = coord.x[index];
			pj_coord[pos].xy.y = coord.y[index];

	END_LOOP2D( )
	//cout << "========================" << endl;
	

	proj_trans_array( P, PD, num, pj_coord );

	FOR_LOOP2D( i, j, 0, _nx_, 0, _ny_ )
		pos = Index2D( i, j, _nx_, _ny_ );
		LonLat.lon[pos] = pj_coord[pos].lp.lam * RADIAN2DEGREE;
		LonLat.lat[pos] = pj_coord[pos].lp.phi * RADIAN2DEGREE;
	END_LOOP2D( )

	

	free( pj_coord );
}


void projTrans( double lon_0, double lat_0, GRID grid, COORD coord, LONLAT LonLat  )
{
	
	//#include "/public/software/proj-8.1.0/include/proj.h"
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

	PJ_COORD pj_LonLat[4];

	//calculateCutLonLat( grid, P, PJ_INV, pj_LonLat );
	
	cart2LonLat(grid, P, PJ_INV, coord, LonLat );

	proj_destroy( P );
	proj_context_destroy( C );

}

/*   
 *	C----------D
 *	|		   |		
 *	|		   |
 *	|		   |
 *	|		   |
 *	A----------B
 */



double bilinear( double x, double y, double x1, double x2, double y1, double y2, double f11, double f12, double f21, double f22 ) {
    //return ( f11 + f21 + f12 + f22 ) * 0.25;//(f11*(x2 - x)*(y2 - y) + f21*(x - x1)*(y2 - y) + f12*(x2 - x)*(y - y1) + f22*(x - x1)*(y - y1));///((x2 - x1)*(y2 - y1));
    return ((f11*(x2 - x)*(y2 - y) + f21*(x - x1)*(y2 - y) + f12*(x2 - x)*(y - y1) + f22*(x - x1)*(y - y1))/((x2 - x1)*(y2 - y1)) );
}
double bilinearInterp( double x, double y, double x1, double y1, double x2, double y2, double A, double B, double C, double D )
{
	double AB = ( x2 - x ) / ( x2 - x1 ) * A + ( x - x1 ) / ( x2 - x1 ) * B;
	double CD = ( x2 - x ) / ( x2 - x1 ) * C + ( x - x1 ) / ( x2 - x1 ) * D; 
	double R  = ( y2 - y ) / ( y2 - y1 ) * AB+ ( y - y1 ) / ( y2 - y1 ) * CD;
	//R = ( A + B + C + D ) * 0.25;
	
	return R;

}

double interp2d(double x[2], double y[2], double z[4], double x_, double y_ )
{
	int i, j;
	double Li = 0.0;
	double Lx[2], Ly[2];

	for( i = 0; i < 2; i ++ ) 
	{
		Lx[i] = 1;
		for( j = 0; j < 2; j ++ ) 
		{
			if( i == j ) continue;
			Lx[i] = Lx[i] * (x_-x[j]) / (x[i]-x[j]);
		}
	}

	for( i = 0; i < 2; i ++ ) 
	{
		Ly[i] = 1;
		for( j = 0; j < 2; j ++ ) 
		{
			if( i == j ) continue;
			Ly[i] = Ly[i] * (y_-y[j]) / (y[i]-y[j]);
		}
	}

	for( j = 0; j < 2; j ++ ) 
	{
		for( i = 0; i < 2; i ++ ) 
		{
			Li = Li + Lx[i] * Ly[j] * z[i*2+j];
		}
	}

  return Li;
}


void terrainInterp( PARAMS params, GRID grid, float * terrain, float * totalTerrain, LONLAT LonLat )
{			
	int totalPointX = ( pointsPerBlockX - 1 ) * params.blockX + 1;
	int totalPointY = ( pointsPerBlockY - 1 ) * params.blockY + 1;

	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;

	int i = 0, j = 0;
	long long index = 0;

	double lonStart = params.lonStart;
	double latStart = params.latStart;

	double deltaLon = double( degreePerBlcokX *  params.blockX ) / double( totalPointX );
	double deltaLat = double( degreePerBlcokY *  params.blockY ) / double( totalPointY );

	//cout << "blockX = " << params.blockX << " ";
	//cout << "blockY = " << params.blockY << " ";
	//cout << "deltaLon = " << deltaLon << " ";
	//cout << "deltaLat = " << deltaLat << endl;

	double x_; double y_; double x1; double y1; double x2; double y2; 
	//double f11; double f12; double f21; double f22;

	double x[2] = { 0.0 };
	double y[2] = { 0.0 };
	double z[4] = { 0.0 };
	
	int I = 0, J = 0, pos = 0;

	FOR_LOOP2D( i, j, 0, _nx_, 0, _ny_ )
		index = Index2D( i, j, _nx_, _ny_ );

		x_ = LonLat.lon[index];		
		y_ = LonLat.lat[index];		

		I = int( ( x_ - lonStart ) / deltaLon ); 
		J = int( ( y_ - latStart ) / deltaLat ); 
			
		x1 = I * deltaLon + lonStart;
		y1 = J * deltaLat + latStart;

		x2 = x1 + deltaLon;
		y2 = y1 + deltaLat;

		x[0] = x1;
		x[1] = x2;
		y[0] = y1;
		y[1] = y2;

		
		pos = Index2D( I, J, totalPointX, totalPointY );
		//f11 = totalTerrain[pos]; //A
		z[0] = totalTerrain[pos]; //A

		pos = Index2D( I + 1, J, totalPointX, totalPointY );
		//f12 = totalTerrain[pos]; //B
		z[2] = totalTerrain[pos]; //B

		pos = Index2D( I + 1, J + 1, totalPointX, totalPointY );
		//f21 = totalTerrain[pos]; //B
		z[1] = totalTerrain[pos]; //B

		pos = Index2D( I + 1, J + 1, totalPointX, totalPointY );
		//f22 = totalTerrain[pos]; //D
		z[3] = totalTerrain[pos]; //D
						
		terrain[index] = interp2d( x, y, z, x_, y_ );
	END_LOOP2D( )

}

void preprocessTerrain( PARAMS params, MPI_Comm comm_cart, MPI_COORD thisMPICoord, GRID grid, COORD coord )
{
	double lon_0 = params.centerLongitude;
	double lat_0 = params.centerLatitude;
	
	int halo = grid.halo;

	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;

	float * terrain;
	terrain = ( float * ) malloc( sizeof( float ) * _nx_ * _ny_ );
	
	int i = 0, j = 0, k = 0;
	long long index = 0;
	
	//if ( 1 == params.SRTM90 && thisMPICoord.Z == grid.PZ - 1 )
	if ( 1 == params.SRTM90 )
	{
		float *  totalTerrain;
		int totalPointX = ( pointsPerBlockX - 1 ) * params.blockX + 1;
		int totalPointY = ( pointsPerBlockY - 1 ) * params.blockY + 1;

		totalTerrain = ( float * ) malloc( totalPointX * totalPointY * sizeof( float )  );
		memset( totalTerrain, 0, totalPointX * totalPointY * sizeof( float ) );
		
		readSRTM90( params, totalTerrain );
#ifdef ERROR_ELEVATION
		resetTotalTerrain( params, totalTerrain );
#endif
		callibrateTotalTerrain( params, totalTerrain );

		LONLAT LonLat;
		LonLat.lon = ( double * ) malloc( sizeof( double ) * _nx_ * _ny_  );
		LonLat.lat = ( double * ) malloc( sizeof( double ) * _nx_ * _ny_  );
		memset( LonLat.lon, 0, sizeof( double ) * _nx_ * _ny_ );
		memset( LonLat.lat, 0, sizeof( double ) * _nx_ * _ny_ );


		projTrans( lon_0, lat_0, grid, coord, LonLat );

		//cout << "==================" << endl;
		//MPI_Barrier( comm_cart );	
	//	cout << "========================" << endl;
		terrainInterp( params, grid, terrain, totalTerrain, LonLat );

		/*
		int maxValue = -100000;
	
		
		for ( j = 0; j < _ny_; j ++ )
			for ( i = 0; i < _nx_; i ++ )
			{
				index = i + j * _nx_;
				if ( terrain[index] > maxValue )
					maxValue = terrain[index];
				//if ( DZ[index] > maxValue )
				//	maxValue = DZ[index];
				//if ( DZ[index] < minValue )
				//	minValue = DZ[index];
			}
		
		cout << "terrain max = " << maxValue << endl;
		*/

		if ( thisMPICoord.X == 0 && thisMPICoord.Y == 0 && thisMPICoord.Z == grid.PZ - 1 )
		{
			char totalTerrainFileName[256];
			sprintf( totalTerrainFileName, "%s/totalTerrain.bin", params.TerrainDir );
			FILE * totalTerrainFile = fopen( totalTerrainFileName, "wb"  );

		
			fwrite( totalTerrain, sizeof( float ), totalPointX * totalPointY, totalTerrainFile );

			fclose( totalTerrainFile );
		}
		free( totalTerrain );
		//MPI_Barrier( comm_cart );	
		if ( thisMPICoord.Z == grid.PZ - 1 )
		{

			if ( 0 == thisMPICoord.X && 0 == thisMPICoord.Y )
				printf( "ouput \"projection data\" including longitude, latitude and terrain on the gound of the calculation area.\n"  );

			FILE * lonFile, * latFile, * terrainFile, * coordxFile, * coordyFile;

			char lonFileName[256], latFileName[256], terrainFileName[256];

			sprintf( lonFileName,     "%s/lon_mpi_%d_%d_%d.bin", params.OUT, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
			sprintf( latFileName,	  "%s/lat_mpi_%d_%d_%d.bin", params.OUT, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
			sprintf( terrainFileName, "%s/terrain_mpi_%d_%d_%d.bin", params.OUT, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );

			lonFile = fopen( lonFileName, "wb" );
			latFile = fopen( latFileName, "wb" );
			terrainFile = fopen( terrainFileName, "wb" );
			
			float lon, lat;
	    
            FOR_LOOP2D( i, j, halo, _nx, halo, _ny )
	            index = Index2D( i, j, _nx_, _ny_ );
		    	lon = LonLat.lon[index];
		    	lat = LonLat.lat[index];
		    	fwrite( &(lon), sizeof( float ), 1, lonFile );
		    	fwrite( &(lat), sizeof( float ), 1, latFile );
		    	fwrite( terrain + index, sizeof( float ),  1, terrainFile );
	        END_LOOP2D( )
				
			fclose( lonFile );
			fclose( latFile );
			fclose( terrainFile );
		}
		/*
		*/
		
		free( LonLat.lon );
		free( LonLat.lat );

	}

	MPI_Barrier( MPI_COMM_WORLD );
	//float * DZ = ( float * ) malloc( sizeof( float ) * _nx_ * _ny_ );
	
	
	double Depth = params.Depth * 1e3;
	int NZ = grid.NZ;
	int pos = 0;

	int K = 0;
	int frontNZ = grid.frontNZ;
	//float minValue = 1000000;
	//float maxValue = -1000000;
	//for ( j = 0; j < _ny_; j ++ )
	//	for ( i = 0; i < _nx_; i ++ )
	//	{
	//		index = i + j * _nx_;
	//		DZ[index] = ( terrain[index] + abs( Depth ) ) / ( NZ - 1 );
	////		if ( terrain[index] > maxValue )
	////			maxValue = terrain[index];
	////		if ( DZ[index] > maxValue )
	////			maxValue = DZ[index];
	////		if ( DZ[index] < minValue )
	////			minValue = DZ[index];
	//	}
		
	//cout << "DZ max = " << maxValue << ", DZ min = " << minValue << endl;
	double DZ = 0.0;
	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;
    
    FOR_LOOP3D( i, j, k, 0, _nx_, 0, _ny_, 0, _nz_ )			
		index = INDEX( i, j, k );
		pos = Index2D( i, j, _nx_, _ny_ );
		K = frontNZ + k - halo; 
		DZ = double( terrain[pos] + abs( Depth ) ) / double( NZ - 1 );
		coord.z[index] = - abs( Depth ) + DZ  * K;
    END_LOOP3D( )

	free( terrain );	
	

}

