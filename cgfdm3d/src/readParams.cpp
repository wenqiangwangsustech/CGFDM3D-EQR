/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:readJson.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-04
*   Discription:
*
================================================================*/
#include "header.h"

void getParams( PARAMS * params)
{
	//printf( "********************\n" );
	char jsonFile[1024] = { 0 };
	strcpy( jsonFile, "params.json" );
	FILE * fp;
	fp = fopen( jsonFile, "r" );

	if ( NULL == fp )
	{
		printf( "There is not %s file!\n", jsonFile );
		MPI_Abort( MPI_COMM_WORLD, 100  );//exit( 1 );
	}
	
	fseek( fp, 0, SEEK_END );
	int len = ftell( fp );
	
	fseek( fp, 0, SEEK_SET );

	char * jsonStr = ( char * ) malloc( len * sizeof( char ) );

	if ( NULL == jsonStr )
	{
		printf( "Can't allocate json string memory\n" );
	}

	fread( jsonStr, sizeof( char ), len, fp );
	
	//printf( "%s\n", jsonStr );
	cJSON * object, *item;
	object = cJSON_Parse( jsonStr );
	if ( NULL == object )
	{
		printf( "Can't parse json file!\n");
		//exit( 1 );	
		return;
	}

	fclose( fp );

	
	
	if (item = cJSON_GetObjectItem(object, "TMAX"))
		params->TMAX = item->valuedouble;
    
	if (item = cJSON_GetObjectItem(object, "DT"))
		params->DT = item->valuedouble;
	
	if (item = cJSON_GetObjectItem(object, "DH"))
		params->DH = item->valuedouble;
    
	if (item = cJSON_GetObjectItem(object, "NX"))
		params->NX = item->valueint;
	if (item = cJSON_GetObjectItem(object, "NY"))
		params->NY = item->valueint;
	if (item = cJSON_GetObjectItem(object, "NZ"))
		params->NZ = item->valueint;

	if (item = cJSON_GetObjectItem(object, "PX"))
		params->PX = item->valueint;
	if (item = cJSON_GetObjectItem(object, "PY"))
		params->PY = item->valueint;
	if (item = cJSON_GetObjectItem(object, "PZ"))
		params->PZ = item->valueint;

	if (item = cJSON_GetObjectItem(object, "centerX"))
		params->centerX = item->valueint;
	if (item = cJSON_GetObjectItem(object, "centerY"))
		params->centerY = item->valueint;

	if (item = cJSON_GetObjectItem(object, "centerLongitude"))
		params->centerLongitude = item->valuedouble;
	if (item = cJSON_GetObjectItem(object, "centerLatitude"))
		params->centerLatitude = item->valuedouble;

	if (item = cJSON_GetObjectItem(object, "IT_SKIP"))
		params->IT_SKIP = item->valueint;
    
	if (item = cJSON_GetObjectItem(object, "sliceX"))
		params->sliceX = item->valueint;
	if (item = cJSON_GetObjectItem(object, "sliceY"))
		params->sliceY = item->valueint;
	if (item = cJSON_GetObjectItem(object, "sliceZ"))
		params->sliceZ = item->valueint;
	
	if (item = cJSON_GetObjectItem(object, "sliceFreeSurf"))
		params->sliceFreeSurf = item->valueint;

	if (item = cJSON_GetObjectItem(object, "nPML"))
		params->nPML = item->valueint;

	if (item = cJSON_GetObjectItem(object, "gauss_hill"))
		params->gauss_hill = item->valueint;

	if (item = cJSON_GetObjectItem(object, "itSlice"))
		params->itSlice = item->valueint;
	if (item = cJSON_GetObjectItem(object, "itStep"))
		params->itStep = item->valueint;

	if (item = cJSON_GetObjectItem(object, "waveOutput"))
		strcpy( params->waveOutput, item->valuestring);
	if (item = cJSON_GetObjectItem(object, "sliceName"))
		strcpy( params->sliceName , item->valuestring);


	if (item = cJSON_GetObjectItem(object, "itStart"))
		params->itStart = item->valueint;
	if (item = cJSON_GetObjectItem(object, "itEnd"))
		params->itEnd = item->valueint;

	if (item = cJSON_GetObjectItem(object, "igpu"))
		params->igpu = item->valueint;


	if (item = cJSON_GetObjectItem(object, "OUT"))
		strcpy( params->OUT, item->valuestring);
	
	if (item = cJSON_GetObjectItem(object, "TerrainDir"))
		strcpy( params->TerrainDir, item->valuestring);

	if (item = cJSON_GetObjectItem(object, "SRTM90"))
		params->SRTM90 = item->valueint;
	
	if (item = cJSON_GetObjectItem(object, "lonStart"))
		params->lonStart = item->valueint;
	if (item = cJSON_GetObjectItem(object, "latStart"))
		params->latStart = item->valueint;
	
	
	if (item = cJSON_GetObjectItem(object, "blockX"))
		params->blockX = item->valueint;
	if (item = cJSON_GetObjectItem(object, "blockY"))
		params->blockY = item->valueint;
	

	
	
	
	if (item = cJSON_GetObjectItem(object, "Depth(km)"))
		params->Depth = item->valuedouble;


	if (item = cJSON_GetObjectItem(object, "MLonStart"))
		params->MLonStart = item->valuedouble;
	if (item = cJSON_GetObjectItem(object, "MLatStart"))
		params->MLatStart = item->valuedouble;
	if (item = cJSON_GetObjectItem(object, "MLonEnd"))
		params->MLonEnd = item->valuedouble;
	if (item = cJSON_GetObjectItem(object, "MLatEnd"))
		params->MLatEnd = item->valuedouble;
	if (item = cJSON_GetObjectItem(object, "MLonStep"))
		params->MLonStep  = item->valuedouble;
	if (item = cJSON_GetObjectItem(object, "MLatStep"))
	 	params->MLatStep  = item->valuedouble;

	if (item = cJSON_GetObjectItem(object, "MVeticalStep"))
		params->MVeticalStep = item->valuedouble;

	if (item = cJSON_GetObjectItem(object, "MediumDir"))
		strcpy( params->MediumDir, item->valuestring);

	//cout << "LonStart = " << params->MLonStart << ", LatStart = " << params->MLatStart << endl;
	//cout << "LonEnd = "   << params->MLonEnd   << ", LatEnd = "   << params->MLatEnd << endl;
	//cout << "LonStep = "  << params->MLonStep  << ", LatStep = "  << params->MLatStep << endl;
	
	if (item = cJSON_GetObjectItem(object, "sourceX"))
		params->sourceX = item->valueint;
	if (item = cJSON_GetObjectItem(object, "sourceY"))
		params->sourceY = item->valueint;
	if (item = cJSON_GetObjectItem(object, "sourceZ"))
		params->sourceZ = item->valueint;

	
	if (item = cJSON_GetObjectItem(object, "sourceFile"))
		strcpy( params->sourceFile, item->valuestring);
	
	
	free( jsonStr );

	return;
}


/*
int main( int argc, char ** argv )
{
	
	PARAMS params;
	getParam( &params );
	
    printf( "%lf\n", params. TMAX );
    
	printf( "%lf\n", params. DT );
    printf( "%lf\n", params. DH );
    
	printf( "%d\n", params.NX );
    printf( "%d\n", params.NY );
    printf( "%d\n", params.NZ );

    printf( "%d\n", params.PX );
    printf( "%d\n", params.PY );
    printf( "%d\n", params.PZ );

	printf( "%d\n", params.centerX );
	printf( "%d\n", params.centerY );

	printf( "%lf\n", params. centerLatitude  ); 
	printf( "%lf\n", params. centerLongitude );
    
	printf( "%d\n", params.SourceX ); 
	printf( "%d\n", params.SourceY );
	printf( "%d\n", params.SourceZ );

	
	printf( "%d\n", params.nPML );
	
	printf( "%d\n", params.SliceX ); 
	printf( "%d\n", params.SliceY );
	printf( "%d\n", params.SliceZ );

	printf( "%d\n", params.itSlice			 );
	printf( "%d\n", params.itStep			 );
	printf( "%s\n", params.waveOutput	 );
	printf( "%s\n", params.sliceName	 );

	printf( "%d\n", params.itStart			 );
	printf( "%d\n", params.itEnd			 );
	printf( "%d\n", params.igpu			 );
	printf( "%s\n", params.OUT	 );
}
*/
