/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:printInfo.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-03
*   Discription:
*
================================================================*/

#include "header.h"
void printInfo( GRID grid )
{

	
	printf( 
	"=============================================\n"
	"MPI:  PX = %5d, PY = %5d, PZ = %5d\n"
	"GRID: NX = %5d, NY = %5d, NZ = %5d\n"
	"DH = %5.2e\n"
	"=============================================\n", 
	grid.PX, grid.PY, grid.PZ,
	grid.NX, grid.NY, grid.NZ,
	grid.DH 
	);



}

