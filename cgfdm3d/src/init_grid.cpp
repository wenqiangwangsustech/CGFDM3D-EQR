/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:init_grid.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-06
*   Discription:
*
================================================================*/
#include "header.h"

void init_grid( PARAMS params, GRID * grid, MPI_COORD thisMPICoord )
{
	int resX = 0;	
	int resY = 0;
	int resZ = 0;	
		
	grid->PX = params.PX;
	grid->PY = params.PY;
	grid->PZ = params.PZ;

	grid->_NX_ = params.NX + 2 * HALO;
	grid->_NY_ = params.NY + 2 * HALO;
	grid->_NZ_ = params.NZ + 2 * HALO;
	
	grid->_NX = params.NX + HALO;
	grid->_NY = params.NY + HALO;
	grid->_NZ = params.NZ + HALO;

	grid->NX = params.NX;
	grid->NY = params.NY;
	grid->NZ = params.NZ;

	grid->nx = params.NX / params.PX;
	grid->ny = params.NY / params.PY;
	grid->nz = params.NZ / params.PZ;
	
	resX = params.NX % params.PX;
	resY = params.NY % params.PY;
	resZ = params.NZ % params.PZ;

	if ( thisMPICoord.X < resX )
	{
		grid->nx ++;
		grid->frontNX = thisMPICoord.X * grid->nx;
	}	
	else
	{	
		grid->frontNX = resX * ( grid->nx + 1 ) + ( thisMPICoord.X - resX ) * grid->nx;
	}

	if ( thisMPICoord.Y < resY )
	{
		grid->ny ++;
		grid->frontNY = thisMPICoord.Y * grid->ny;
	}	
	else
	{	
		grid->frontNY = resY * ( grid->ny + 1 ) + ( thisMPICoord.Y - resY ) * grid->ny;
	}

	if ( thisMPICoord.Z < resZ )
	{
		grid->nz ++;
		grid->frontNZ = thisMPICoord.Z * grid->nz;
	}	
	else
	{	
		grid->frontNZ = resZ * ( grid->nz + 1 ) + ( thisMPICoord.Z - resZ ) * grid->nz;
	}

	MPI_Barrier( MPI_COMM_WORLD );
	//cout << "X: " << thisMPICoord.X << ", Y: " << thisMPICoord.Y << ", Z: " << thisMPICoord.Z << " ====> frontNX = " << grid->frontNX << ", frontNY = " << grid->frontNY << ", frontNZ = " << grid->frontNZ << endl;
	//cout << "X: " << thisMPICoord.X << ", Y: " << thisMPICoord.Y << ", Z: " << thisMPICoord.Z << " ====> nx = " << grid->nx << ", ny = " << grid->ny << ", nz = " << grid->nz << endl;

	grid->_frontNX = grid->frontNX + HALO;
	grid->_frontNY = grid->frontNY + HALO;
	grid->_frontNZ = grid->frontNZ + HALO;


	grid->_nx = grid->nx + HALO;
	grid->_ny = grid->ny + HALO;
	grid->_nz = grid->nz + HALO;

	grid->_nx_ = grid->nx + 2 * HALO;
	grid->_ny_ = grid->ny + 2 * HALO;
	grid->_nz_ = grid->nz + 2 * HALO;
	
	grid->originalX = params.centerX;
	grid->originalY = params.centerY;

	grid->_originalX = grid->originalX + HALO;
	grid->_originalY = grid->originalY + HALO;

	grid->halo = HALO;

	grid->DH = params.DH;
	grid->rDH = 1.0 / grid->DH;
	//printf( "dh = %e\n", grid->DH );
	

	grid->nPML = params.nPML;
}
