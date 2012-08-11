

/* Wrapper for MPI routines that take string arguments, to avoid the
   difference between how Fortran compilers pass strings */
#include <mpi.h>
#include <stdlib.h>

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

void my_get_processor_name_(int *int_name, int *namelen, int *ierr)
{
   int i, fnamelen;
   char *name;
   fnamelen = *namelen;
   name = (char *)malloc(MPI_MAX_PROCESSOR_NAME*sizeof(char));
   *ierr = MPI_Get_processor_name(name, namelen);
   if (*namelen > fnamelen) *namelen=fnamelen;
   for (i=0; i<*namelen; i++) int_name[i] = (int)name[i];
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
