/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:cpu_Malloc.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-02
*   Discription:
*
================================================================*/

#ifndef GPU_CUDA

#include <stdio.h>
#include <stdlib.h>
int Malloc( void ** mem, long long  size  )
{
	*mem = malloc( size );
	if ( *mem == NULL )
	{
		printf( "can not malloc, Error: %s:%d\n", __FILE__, __LINE__ );
	}
	return 0;

}


#endif

