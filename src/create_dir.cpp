/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:create_dir.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-30
*   Discription:
*
================================================================*/

#include "header.h"


void createDir( PARAMS params )
{
	int thisRank;
	MPI_Comm_rank( MPI_COMM_WORLD, &thisRank );
	
	if ( 0 == thisRank )
	{
#if __GNUC__
 		mkdir( params.OUT, 0777 );
#elif _MSC_VER
		_mkdir( params.OUT );
#endif 

	}

	MPI_Barrier( MPI_COMM_WORLD );
}

