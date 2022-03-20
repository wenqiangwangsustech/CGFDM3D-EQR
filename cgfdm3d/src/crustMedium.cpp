/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:crustMedium.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-12-20
*   Discription:
*
================================================================*/
#include "header.h"



#define CRUSTSIZE 4
#define CUnit 1e3// 1km/s = 1e3m
#define BoundaryDifference 1000.0


void verifyInterpVs( PARAMS params, GRID grid, MPI_COORD thisMPICoord, float * Vs, int k );


typedef struct CRUST
{
	float * vs		;
	float * vp		;
	float * rho		;
	float * bnds	;
}CRUST;


void allocCrustForFile( int nLon, int nLat, int nLayer, CRUST * crust)
{
	int num = nLayer * nLon * nLat;
	int size = num * sizeof( float ) * CRUSTSIZE;

	float * pCrust = ( float * ) malloc( size );
	
	crust->vs	= pCrust + 0 * num;
	crust->vp	= pCrust + 1 * num;
	crust->rho	= pCrust + 2 * num;
	crust->bnds = pCrust + 3 * num;
	
}

void allocCrustInterp( GRID grid, int nLayer, CRUST * crust )
{
	
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;

	int num = nLayer * grid._nx_ * grid._ny_;
	int size = num * sizeof( float ) * CRUSTSIZE;

	float * pCrust = ( float * ) malloc( size );
	
	crust->vs	= pCrust + 0 * num;
	crust->vp	= pCrust + 1 * num;
	crust->rho	= pCrust + 2 * num;
	crust->bnds = pCrust + 3 * num;
	
}

void freeCrust( CRUST crust )
{
	free( crust.vs );
}

void crustInterpFunc(	PARAMS params, GRID grid, int nLon, int nLat, int nLayer, 
					CRUST crustFile, CRUST crustInterp, 
					LONLAT LonLat, COORD coord, MPI_COORD thisMPICoord )
{			

	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	int i = 0, j = 0;
	long long index = 0;

	float LonStart = -179.5;
	float LatStart =   89.5;

	float LonStep  = params.CrustLonStep;
	float LatStep  = params.CrustLatStep;

	double x_; double y_; double x1; double y1; double x2; double y2; 

	double x[2] = { 0.0 };
	double y[2] = { 0.0 };
	double z[4] = { 0.0 };
	double vs  [4] = { 0.0 };
	double vp  [4] = { 0.0 };
	double rho [4] = { 0.0 };
	double bnds[4] = { 0.0 };

	int I = 0, J = 0, pos = 0;
	int k = 0, K = 0;

	int IMax = 0, JMax = 0;
	int IMin = 0, JMin = 0;

	int nLonLat = nLon * nLat;



    for ( K = 0; K < nLayer; K ++ )
    {
		//printf( "K = %d\n", K  );
        FOR_LOOP2D( i, j, 0, _nx_, 0, _ny_ )
			//index = i + j * _nx_;
			index = Index2D( i, j, _nx_, _ny_);

			x_ = LonLat.lon[index];		
			y_ = LonLat.lat[index];		

			I = int( ( x_ - LonStart ) / LonStep ); 
			J = int( ( LatStart - y_ ) / LatStep ) + 1; 
				
			x1 = I * LonStep + LonStart;
			y1 = LatStart - J * LatStep ;
			

			x2 = x1 + LonStep;
			y2 = y1 + LatStep;

			x[0] = x1;	x[1] = x2;	y[0] = y1;	y[1] = y2;

			//pos = K + I * NLayer + J * NLayer * nLon;
			pos = Index2D( I, J, nLon, nLat ) + K * nLonLat;
			//A
			vs  [0] = crustFile.vs  [pos];
        	vp  [0] = crustFile.vp  [pos];
        	rho [0] = crustFile.rho [pos];
			bnds[0] = crustFile.bnds[pos];

			//pos = K + ( I + 1 ) * NLayer + J * NLayer * nLon;
			pos = Index2D( I + 1, J, nLon, nLat ) + K * nLonLat;
			//B
			vs  [2] = crustFile.vs  [pos];
        	vp  [2] = crustFile.vp  [pos];
        	rho [2] = crustFile.rho [pos];
			bnds[2] = crustFile.bnds[pos];

			//pos = K + I * NLayer + ( J + 1 ) * NLayer * nLon;
			pos = Index2D( I, J - 1, nLon, nLat ) + K * nLonLat;
			//C
			vs  [1] = crustFile.vs  [pos];
        	vp  [1] = crustFile.vp  [pos];
        	rho [1] = crustFile.rho [pos];
			bnds[1] = crustFile.bnds[pos];

			//pos = K + ( I + 1 ) * NLayer + ( J + 1 ) * NLayer * nLon;
			pos = Index2D( I + 1, J - 1, nLon, nLat ) + K * nLonLat;
			//D
			vs  [3] = crustFile.vs  [pos];
        	vp  [3] = crustFile.vp  [pos];
        	rho [3] = crustFile.rho [pos];
			bnds[3] = crustFile.bnds[pos];
			
			//pos = i + j * _nx_ + K * _nx_ * _ny_;
			pos = K * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_);
			crustInterp.vs  [pos] = interp2d( x, y, vs  , x_, y_ );
			crustInterp.vp  [pos] = interp2d( x, y, vp  , x_, y_ );
			crustInterp.rho [pos] = interp2d( x, y, rho , x_, y_ );
			crustInterp.bnds[pos] = interp2d( x, y, bnds, x_, y_ );
        END_LOOP2D( )
	if ( grid.PZ - 1 == thisMPICoord.Z )
		verifyInterpVs( params, grid, thisMPICoord, crustInterp.bnds, K );
    }
	
}


void solveStructure( GRID grid, STRUCTURE structure, COORD coord, CRUST crustInterp, int nLayer  )
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	long long index = 0, 
	posl = 0, posl_1 = 0,
	posnLayer_1 = 0,
	pos0, posm, posm1;

	float vs, vp, rho, bnds, bnds1;
	int i, j, k, m, l;	
	
	float coordZ = 0.0f;

	int startLayer = 0;


	float maxBnd = -10000000;
	float minBnd = 10000000;
	int K = 0;
	long long pos = 0;
    for ( K = 0; K < nLayer; K ++ )
    {
        FOR_LOOP2D( i, j, 0, _nx_, 0, _ny_ )
			index = Index2D( i, j, _nx_, _ny_);
			pos = K * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_);
			crustInterp.bnds[pos] *= CUnit;
			if (  crustInterp.bnds[pos] >= maxBnd )
			{
				maxBnd = crustInterp.bnds[pos];
			}
			if (  crustInterp.bnds[pos] < minBnd )
			{
				minBnd = crustInterp.bnds[pos];
			}
        END_LOOP2D( )

	}


	//printf( "maxBnd = %f, minBnd = %f\n", maxBnd, minBnd );

	FOR_LOOP3D( i, j, k, 0, _nx_, 0, _ny_, 0, _nz_ )			
        index = INDEX( i, j, k );

		coordZ = coord.z[index];
		
		pos0 = 0 * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_ ); 
		if ( coordZ > crustInterp.bnds[pos0] )
		{
			l = 0;
			posl = l * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_ ); 
			while ( abs( crustInterp.bnds[pos0] - crustInterp.bnds[posl] ) < BoundaryDifference )
			{
				l += 1;
				posl = l * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_ ); 
			}

			posl_1 = ( l - 1 ) * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_ ); 
			vs   = crustInterp.vs  [posl_1];
			vp   = crustInterp.vp  [posl_1];
			rho  = crustInterp.rho [posl_1];
			
		}

		posnLayer_1 = ( nLayer - 1 ) * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_ ); 
		if ( coordZ <= crustInterp.bnds[posnLayer_1] )
		{
			vs   = crustInterp.vs  [posnLayer_1];
			vp   = crustInterp.vp  [posnLayer_1];
			rho  = crustInterp.rho [posnLayer_1];
		}
		
		pos0 = 0 * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_ ); 
		posnLayer_1 = ( nLayer - 1 ) * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_ ); 
		if ( coordZ <= crustInterp.bnds[pos0] && coordZ > crustInterp.bnds[posnLayer_1] )
		{
			for ( m = 0; m < nLayer - 1; m ++ )
			{
				posm = m * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_ ); 
				posm1 = ( m + 1 ) * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_ ); 
				if ( coordZ <= crustInterp.bnds[posm] && coordZ > crustInterp.bnds[posm1]  )
				{
					l = m;
					posl = l * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_ );
					while ( abs( crustInterp.bnds[posm] - crustInterp.bnds[posl] ) < BoundaryDifference )
					{
						l += 1;
						posl = l * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_ );
					}
					posl_1 = ( l - 1 ) * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_ );
					vs   = crustInterp.vs  [posl_1];
					vp   = crustInterp.vp  [posl_1];
					rho  = crustInterp.rho [posl_1];
				}

			}

		}


		/*

		if ( vp > 8.04 ) vp = 8.04;
		if ( vs > 4.47 ) vs = 4.47;
		if ( rho > 3.31 ) rho = 3.31;

		if ( vp < 1.80 )
		{
			vp  = 1.50;
			vs  = 0.00;
			rho = 1.02;
		}
		*/
		structure.Vs[index]  = vs  * CUnit;//LAM;
		structure.Vp[index]  = vp  * CUnit;//MU;
		structure.rho[index] = rho * CUnit;//RHO;
		/*
		for ( K = startLayer; K < nLayer - 1; K ++ )
		{
			pos = K * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_ );
			pos1 = ( K + 1 ) * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_ );
			vs   = crustInterp.vs  [pos];
			vp   = crustInterp.vp  [pos];
			rho  = crustInterp.rho [pos];
			bnds = crustInterp.bnds[pos] * CUnit;
			bnds1 = crustInterp.bnds[pos1] * CUnit;

			coordZ = coord.z[index];

			if ( coordZ >= bnds1 && coordZ < bnds )
			{
				structure.Vs[index]  = vs;//LAM;
				structure.Vp[index]  = vp;//MU;
				structure.rho[index] = rho;//RHO;
			}


			if( K == startLayer )
			{	
				if ( coordZ >= bnds )
				{
					pos = 2 * _nx_ * _ny_ + Index2D( i, j, _nx_, _ny_ );
					vs   = crustInterp.vs  [pos];
					vp   = crustInterp.vp  [pos];
					rho  = crustInterp.rho [pos];
					structure.Vs[index]  = vs;//LAM;
					structure.Vp[index]  = vp;//MU;
					structure.rho[index] = rho;//RHO;
				}
			}
			if( K + 1 == nLayer - 1 )
			{	
				if ( coordZ < bnds1 )
				{
					vs   = crustInterp.vs  [pos1];
					vp   = crustInterp.vp  [pos1];
					rho  = crustInterp.rho [pos1];
					structure.Vs[index]  = vs;//LAM;
					structure.Vp[index]  = vp;//MU;
					structure.rho[index] = rho;//RHO;
				}
			}


		}
		structure.Vs[index]  = vs  * CUnit;//LAM;
		structure.Vp[index]  = vp  * CUnit;//MU;
		structure.rho[index] = rho * CUnit;//RHO;
		*/
		

    END_LOOP3D( )



}


 

void readCrustal_1( PARAMS params, GRID grid, MPI_COORD thisMPICoord, COORD coord, STRUCTURE structure )
{
	FILE * vsFile, * vpFile, * rhoFile, * bndsFile;

	char vsName[256], vpName[256], rhoName[256], bndsName[256];

	sprintf( vsName  , "%s/crust1.vs",   params.crustDir );   
	sprintf( vpName  , "%s/crust1.vp",   params.crustDir ); 
	sprintf( rhoName , "%s/crust1.rho",  params.crustDir );   
	sprintf( bndsName, "%s/crust1.bnds", params.crustDir ); 
	

	vsFile   = fopen( vsName  , "rb" );
	vpFile   = fopen( vpName  , "rb" );
	rhoFile  = fopen( rhoName , "rb" );
	bndsFile = fopen( bndsName, "rb" );
	
	int nLon = 360, nLat = 180;

	int nLonLat = nLon * nLat;
	int nLayer = 9;
	float vs, vp, rho, bnds;

	CRUST crustFile, crustInterp;

	allocCrustForFile( nLon, nLat, nLayer, &crustFile);
	allocCrustInterp( grid, nLayer, &crustInterp );



	int i, n = 0;
	FILE * txtFile = fopen( "crustTxt.txt", "w"  );
	for ( i = 0; i < nLonLat; i ++ )
	{
		for ( n = 0; n < nLayer; n ++ )
		{		
			fscanf( vsFile  , "%f", &vs );
			fscanf( vpFile  , "%f", &vp );
			fscanf( rhoFile , "%f", &rho );
			fscanf( bndsFile, "%f", &bnds );
			crustFile.vs  [i + n * nLonLat ] = vs  ;
        	crustFile.vp  [i + n * nLonLat ] = vp  ;
        	crustFile.rho [i + n * nLonLat ] = rho ;
			crustFile.bnds[i + n * nLonLat ] = bnds;
		}
		
		fprintf( txtFile,
			"vs = %e, vp = %e, rho = %e, bnds = %e\n", 
			crustFile.vs  [i + 3 * nLonLat ],	
			crustFile.vp  [i + 3 * nLonLat ],	
			crustFile.rho [i + 3 * nLonLat ],	
			crustFile.bnds[i + 3 * nLonLat ]	
		);
	}
	fclose( txtFile );
	

	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	double lon_0 = params.centerLongitude;
	double lat_0 = params.centerLatitude;
	LONLAT LonLat;
	LonLat.lon = ( double * ) malloc( sizeof( double ) * _nx_ * _ny_  );
	LonLat.lat = ( double * ) malloc( sizeof( double ) * _nx_ * _ny_  );
	memset( LonLat.lon, 0, sizeof( double ) * _nx_ * _ny_ );
	memset( LonLat.lat, 0, sizeof( double ) * _nx_ * _ny_ );

	projTrans( lon_0, lat_0, grid, coord, LonLat );

	crustInterpFunc(	params, grid, nLon, nLat, nLayer, 
					crustFile, crustInterp, 
					LonLat, coord, thisMPICoord );

	
	solveStructure( grid, structure, coord, crustInterp, nLayer  );




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
	    
			int i, j;
			int halo = grid.halo;
			int _nx = grid._nx;
			int _ny = grid._ny;
			

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



	freeCrust( crustFile );
	freeCrust( crustInterp );

	free( LonLat.lon );
	free( LonLat.lat );

	fclose( vsFile   );
	fclose( vpFile   );
	fclose( rhoFile  );
	fclose( bndsFile );

}





