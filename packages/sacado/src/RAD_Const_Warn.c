#include <stdio.h>

int RAD_Const_Warn_verbose = 0;

 int
RAD_Const_Warn(void *v)
{
	if (RAD_Const_Warn_verbose)
		printf("RAD_Const_Warn(%lx)\n", (unsigned long)v);
	}
