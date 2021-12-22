/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:structure.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-23
*   Discription:
*
================================================================*/
#include "header.h"


void verifyInterpVs( PARAMS params, GRID grid, MPI_COORD thisMPICoord, float * Vs, int k )
{

	char fileName[256];

	sprintf( fileName, "%s/InterpVs%d_mpi_%d_%d_%d.bin", params.OUT, k, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );

	int halo = grid.halo;

	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;

	FILE * file = fopen( fileName, "wb" );
	
	int i, j, index = 0;

	//cout << "k = " << k << endl;
    FOR_LOOP2D( i, j, halo, _nx, halo, _ny )
		//index = i + j * _nx_ + k * _nx_ * _ny_;
		index = k * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_);
		fwrite( Vs + index, sizeof( float ), 1, file );
    END_LOOP2D( )

	fclose( file );

}

void structureInterp(	PARAMS params, GRID grid, 
					float * VModel, int nLon, int nLat, int NLayer, 
					STRUCTURE structure, LONLAT LonLat, COORD coord, MPI_COORD thisMPICoord )
{			

	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	int i = 0, j = 0;
	long long index = 0;

	double LonStart = params.MLonStart;
	double LatStart = params.MLatStart;
	double LonEnd = params.MLonEnd;
	double LatEnd = params.MLatEnd;

	float LonStep  = params.MLonStep ;
	float LatStep  = params.MLatStep ;

	double x_; double y_; double x1; double y1; double x2; double y2; 

	double x[2] = { 0.0 };
	double y[2] = { 0.0 };
	double z[4] = { 0.0 };

	float * Vs = ( float * ) malloc( sizeof( float ) * _nx_ * _ny_ * NLayer );
	
	float * coordX = coord.x;	
	float * coordY = coord.y;
	float * coordZ = coord.z;	

	int I = 0, J = 0, pos = 0;

	int k = 0, K = 0;

	int IMax = 0, JMax = 0;
	int IMin = 0, JMin = 0;


    for ( K = 0; K < NLayer; K ++ )
    {
        FOR_LOOP2D( i, j, 0, _nx_, 0, _ny_ )
			//index = i + j * _nx_;
			index = Index2D( i, j, _nx_, _ny_);

			x_ = LonLat.lon[index];		
			y_ = LonLat.lat[index];		

			I = int( ( x_ - LonStart ) / LonStep ); 
			J = int( ( y_ - LatStart ) / LatStep ); 
				
			x1 = I * LonStep + LonStart;
			y1 = J * LatStep + LatStart;
			

			x2 = x1 + LonStep;
			y2 = y1 + LatStep;

			x[0] = x1;	x[1] = x2;	y[0] = y1;	y[1] = y2;

			//pos = K + I * NLayer + J * NLayer * nLon;
			pos = K + Index2D( I, J, nLon, nLat ) * NLayer;
			z[0] = VModel[pos]; //A

			//pos = K + ( I + 1 ) * NLayer + J * NLayer * nLon;
			pos = K + Index2D( I + 1, J, nLon, nLat ) * NLayer;
			z[2] = VModel[pos]; //B

			//pos = K + I * NLayer + ( J + 1 ) * NLayer * nLon;
			pos = K + Index2D( I, J + 1, nLon, nLat ) * NLayer;
			z[1] = VModel[pos]; //C

			//pos = K + ( I + 1 ) * NLayer + ( J + 1 ) * NLayer * nLon;
			pos = K + Index2D( I + 1, J + 1, nLon, nLat ) * NLayer;
			z[3] = VModel[pos]; //D
			
			//pos = i + j * _nx_ + K * _nx_ * _ny_;
			pos = K * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_);
			Vs[pos] = interp2d( x, y, z, x_, y_ );
        END_LOOP2D( )
    }


	//cout << "IMax = " << IMax << ", JMax = " << JMax;
	//cout << ", IMin = " << IMin << ", JMin = " << JMin << endl;
	if ( grid.PZ - 1 == thisMPICoord.Z )
		verifyInterpVs( params, grid, thisMPICoord, Vs, 40 );


#define DealWithFirstLayer
#ifdef DealWithFirstLayer
	int pos0 = 0;
	int pos1 = 0;
	
    FOR_LOOP2D( i, j, 0, _nx_, 0, _ny_ ) 
		pos0 = 0 * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_);
		pos1 = 1 * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_);
		Vs[pos0] = Vs[pos1];
    END_LOOP2D( )

#endif

	
	double DZ = 0.0;
	int NZ = grid.NZ;
	double Depth = params.Depth * 1e3;
	double terrain = 0.0;
	float MVeticalStep = params.MVeticalStep;
	int index1 = 0, index2 = 0;


	float VS = 0.0, VP = 0.0, RHO = 0.0;
    float LAM = 0.0, MU = 0.0;

	FOR_LOOP3D( i, j, k, 0, _nx_, 0, _ny_, 0, _nz_ )			
        index = INDEX( i, j, k );

#define StructureTerrain
#ifdef StructureTerrain
		index1 = INDEX( i, j, 3);
		index2 = INDEX( i, j, 4);
		DZ = coordZ[index2] - coordZ[index1];
		terrain = DZ * double( NZ - 1 ) - abs( Depth );
		
		K = int(-( coordZ[index] - terrain)) / MVeticalStep;
#else
	    K = int( - coordZ[index] / MVeticalStep );				
#endif
		if ( K <= 0 )
			K = 0;
		
		pos = K * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_);




		VS = Vs[pos] * 1e-3;
		VP = ( 0.9409 + 2.0947 * VS - 0.8206 * ( pow(VS, 2) ) + 0.2683 * ( pow(VS, 3) ) - 0.0251 * ( pow(VS, 4) ) );
		RHO = ( 1.6612 * VP - 0.4721 * ( pow(VP, 2) ) + 0.0671 * ( pow(VP, 3) ) - 0.0043 * ( pow(VP, 4) ) + 0.000106 * ( pow(VP, 5)) );




		structure.Vs[index] = VS * 1e3;//LAM;
		structure.Vp[index] = VP * 1e3;//MU;
		structure.rho[index] = RHO * 1e3;//RHO;
 
 





    END_LOOP3D( )
	
	free( Vs );
	
}



#define Deeper 1e4 
#define MUnit 1e3// 1km/s = 1e3m
void readWeisenShenModel( PARAMS params, GRID grid, MPI_COORD thisMPICoord, COORD coord, STRUCTURE structure )
{
	
	float LonStart = params.MLonStart;
	float LatStart = params.MLatStart;
	float LonEnd = params.MLonEnd;
	float LatEnd = params.MLatEnd;



	float LonStep  = params.MLonStep ;
	float LatStep  = params.MLatStep ;

	
	//cout << "LonStart = " << LonStart << ", LatStart = " << LatStart << endl;
	//cout << "LonEnd = " << LonEnd << ", LatEnd = " << LatEnd << endl;
	//cout << "LonStep = " << LonStep << ", LatStep = " << LatStep << endl;
	
	
	char fileName[256];

	int i, j, nLon, nLat;

	nLon = int( ( LonEnd - LonStart ) / LonStep + 0.5 + 1); 
	nLat = int( ( LatEnd - LatStart ) / LatStep + 0.5 + 1); 

	//cout << "========================" << endl;
	char lonName[64];
	char latName[64];

	//float * MLon = ( float * )malloc( nLat );
	//float * VModel = ( float * )malloc( sizeof( float ) * nLat * nLon );
	

	float Depth = params.Depth * 1e3 + Deeper;
	int NLayer = int( Depth / params.MVeticalStep );

	//vector< list<float> >  VModel( nLon * nLat );
	float * VModel = ( float * ) malloc( nLon * nLat * NLayer * sizeof( float ) );
	int pos = 0;
	double dep, vs, err;
	

    FOR_LOOP2D( i, j, 0, nLon, 0, nLat )
		//cout << params.StructureDir << endl;

		sprintf( lonName, "%.1f", LonStart + i * LonStep );
		sprintf( latName, "%.1f", LatStart + j * LatStep );
		
		if( lonName[strlen( lonName ) - 1] == '0' )
		{
			lonName[strlen( lonName ) - 1] = 0;
			lonName[strlen( lonName ) - 1] = 0;
		}
		if( latName[strlen( latName ) - 1] == '0' )
		{
			latName[strlen( latName ) - 1] = 0;
			latName[strlen( latName ) - 1] = 0;
		}

		sprintf( fileName, "%s/%s_%s.mod", params.MediumDir, lonName, latName );
		
		
		//if ( thisMPICoord.X == 0 && thisMPICoord.Y == 0 && thisMPICoord.Z == 0 )
		//	cout << fileName << endl;

		FILE * file = fopen( fileName, "r" );
		

		if ( NULL == file )
		{
			cout << "Can't open file: " << fileName << endl;
		}

		

		for ( int K = 0; K < NLayer; K ++ )
		{
			fscanf( file, "%lf %lf %lf", &dep, &vs, &err );

			pos = K + Index2D( i, j, nLon, nLat ) * NLayer;
			VModel[pos] = vs * MUnit;
		}


		fclose( file );

		memset( fileName, 0, 256 );
		
    END_LOOP2D( )

	if ( thisMPICoord.X == 0 && thisMPICoord.Y == 0 && thisMPICoord.Z == 0 )
	{
		char VModelFileName[256];
		sprintf( VModelFileName, "%s/VModel.txt", params.MediumDir);
		FILE * fileTxt = fopen( VModelFileName, "w" );
		

		for ( int K = 0; K < NLayer; K ++ )
		{
			for ( int j = 0; j < nLat; j ++ )
			{
				for ( int i = 0; i < nLon; i ++ )
				{
								
				    //pos = k + i * NLayer + j * NLayer * nLon;
			        pos = K + Index2D( i, j, nLon, nLat ) * NLayer;
					fprintf( fileTxt, "%.0f ", VModel[pos] );
				}

				fprintf( fileTxt, "\n" );
			}
			fprintf( fileTxt, "K = %d\n", K );
		}


		fclose( fileTxt );
	}
	

	
	double lon_0 = params.centerLongitude;
	double lat_0 = params.centerLatitude;
	
	int halo = grid.halo;

	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;
	

	LONLAT LonLat;
	LonLat.lon = ( double * ) malloc( sizeof( double ) * _nx_ * _ny_  );
	LonLat.lat = ( double * ) malloc( sizeof( double ) * _nx_ * _ny_  );
	memset( LonLat.lon, 0, sizeof( double ) * _nx_ * _ny_ );
	memset( LonLat.lat, 0, sizeof( double ) * _nx_ * _ny_ );

	projTrans( lon_0, lat_0, grid, coord, LonLat );

	structureInterp( params, grid, VModel, nLon, nLat, NLayer, structure, LonLat, coord, thisMPICoord );



		if ( thisMPICoord.Z == grid.PZ - 1 && !params.useTerrain )
		{
			FILE * lonFile, * latFile, * terrainFile, * coordxFile, * coordyFile;

			char lonFileName[256], latFileName[256], terrainFileName[256];

			sprintf( lonFileName,     "%s/lon_mpi_%d_%d_%d.bin", params.OUT, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
			sprintf( latFileName,	  "%s/lat_mpi_%d_%d_%d.bin", params.OUT, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );
			sprintf( terrainFileName, "%s/terrain_mpi_%d_%d_%d.bin", params.OUT, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );

			lonFile = fopen( lonFileName, "wb" );
			latFile = fopen( latFileName, "wb" );
			terrainFile = fopen( terrainFileName, "wb" );

				
			float lon, lat;
			long long index;
	    
            FOR_LOOP2D( i, j, halo, _nx, halo, _ny )
	            index = Index2D( i, j, _nx_, _ny_ );
		    	lon = LonLat.lon[index];
		    	lat = LonLat.lat[index];
		    	fwrite( &(lon), sizeof( float ), 1, lonFile );
		    	fwrite( &(lat), sizeof( float ), 1, latFile );
	        END_LOOP2D( )
				
			fclose( lonFile );
			fclose( latFile );
			fclose( terrainFile );
		}

	
	free( LonLat.lon );
	free( LonLat.lat );


	free( VModel );	
	

}

void dealWeisenShenStructure( PARAMS params, GRID grid, COORD coord, MPI_COORD thisMPICoord )
{


	STRUCTURE structure;
	MEDIUM medium;
	allocMedium_cpu( grid, &medium );
	structure.Vs  = medium.mu;
	structure.Vp  = medium.lambda;
	structure.rho = medium.buoyancy;
	
	readWeisenShenModel( params, grid, thisMPICoord, coord, structure );
	//data2D_output_bin( params, grid, thisMPICoord, structure.Vs, "Vs" );

	freeMedium_cpu(  medium );


}


