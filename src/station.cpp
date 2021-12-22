/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:station.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-22
*   Discription:
*
================================================================*/

void allocStation( STATION_IDX stationIndex, int pointNum )
{
	
	long long size = sizeof( float ) * pointNum * 3;
	
	float * pIndex = NULL;

	CHECH( Malloc( ( void ** )&pIndex, size ) );
	CHECH( Memset( pIndex, 0, size );

	stationIndex.X = pIndex;
	stationIndex.Y = pIndex + pointNum;
	stationIndex.Z = pIndex + pointNum * 2;

}

void freeStationIndex( STATION_IDX stationIndex )
{
	Free( stationIndex.X );
}


void readStation(    )
{




}

