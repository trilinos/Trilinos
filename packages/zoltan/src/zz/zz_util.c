/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>

#include "zz_util_const.h"
#include "zoltan_mem.h"
#include "zz_const.h"

extern int fsync(int fd);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* Remove leading & trailing white space and convert to upper case. */

int Zoltan_Clean_String(
const char *string1,			/* original string */
char **pstring2) 		/* cleaned string to return */
{

    char     *string2;		/* cleaned up string */
    int       start, end;	/* indices bounding true string */
    int       length1;		/* length of string 1 */
    int       i;		/* loop counter */

    length1 = strlen(string1);
    start = 0;
    end = length1;
    while (start < length1 && isspace((int)(string1[start])))
	start++;
    while (end > start && isspace((int)(string1[end])))
	end--;

    string2 = (char *) ZOLTAN_MALLOC((end - start + 1) * sizeof(char));
    *pstring2 = string2;

    if (string2 == NULL)
	return (ZOLTAN_MEMERR);

    for (i = start; i < end; i++) {
	*string2++ = toupper(string1[i]);
    }
    *string2 = '\0';

    return (ZOLTAN_OK);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* Version of strdup that uses Zoltan_Malloc */

char *Zoltan_Strdup(const char *str)
{
  /* char *yo = "Zoltan_Strdup"; */
  char *c = NULL;
  if (!str){
    return c;
  }

  c = (char *)ZOLTAN_MALLOC(strlen(str) + 1);

  if (c){
    strcpy(c, str);
  }

  return c;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
 /*
  * These transformations are used by the box assign functions when 
  * coordinates have been transformed due to degenerate geometries.
  */
void Zoltan_Transform_Point(
  double *p,         /* point to transform */
  double (*m)[3],    /* linear transformation */
  int *permute,      /* simplified transformation, permute coordinates */
  int d,             /* dimension of input (2 or 3) */
  int ndims,         /* dimension of output (1, 2 or 3) */
  double *v)         /* output */
{
  double tmp[3];

  tmp[0] = p[0];
  tmp[1] = p[1];
  tmp[2] = p[2];

  if (permute[0] >= 0){
    v[0] = p[permute[0]];
    v[1] = (ndims < 2) ? 0.0 : p[permute[1]];

    if (d == 3){
      v[2] = (ndims < 3) ? 0.0 : p[permute[2]];
    }
  }
  else{
    if (d == 2){  
      v[0] = m[0][0] * tmp[0]  +  m[0][1] * tmp[1];
   
      v[1] = (ndims < 2) ? 
              0.0 :
              m[1][0] * tmp[0]  +  m[1][1] * tmp[1];
    }
    else if (d == 3) {
      v[0] = m[0][0]*tmp[0] + m[0][1]*tmp[1] + m[0][2]*tmp[2];
   
      v[1] = (ndims < 2) ? 
              0.0 :
              m[1][0]*tmp[0] + m[1][1]*tmp[1] + m[1][2]*tmp[2];
   
      v[2] = (ndims < 3) ?
              0.0 :
              m[2][0]*tmp[0] + m[2][1]*tmp[1] + m[2][2]*tmp[2];
     
    }
  }
}

/* 
 * Input: the bounds of an axis-aligned box, a linear transformation, the
 * dimension (2 or 3) of the box, and the dimension of the transformed box
 * (1, 2 or 3).
 *
 * Output: the 4 or 8 vertices of the box obtained by applying the transformation
 * followed by:
 *  ndims is 1: projecting to the X-axis
 *  ndims is 2: projecting to the XY plane
 *  ndims is 3: no projection
 */

void Zoltan_Transform_Box_Points(
  double *lo, double *hi,    /* input: bounds of 2D or 3D axis-aligned box */
  double (*m)[3],            /* input: linear transformation     */
  int *perm,                 /* input: coordinate permutation */
  int d,                     /* dimension of box (2 or 3) */
  int ndims,                 /* dimension of transformed box (1, 2 or 3) */
  double (*v)[3])            /* output: 4 or 8 vertices of resulting box    */
{
     int i;

     if (d == 2){  
       v[0][0] = lo[0]; v[0][1] = lo[1];
       v[1][0] = lo[0]; v[1][1] = hi[1];
       v[2][0] = hi[0]; v[2][1] = hi[1]; 
       v[3][0] = hi[0]; v[3][1] = lo[1]; 

       for (i=0; i<4; i++){
         Zoltan_Transform_Point(v[i], m, perm, 2, ndims, v[i]);
       }
     }
     else if (d == 3) {
       v[0][0] = lo[0]; v[0][1] = lo[1]; v[0][2] = lo[2];
       v[1][0] = lo[0]; v[1][1] = hi[1]; v[1][2] = lo[2];
       v[2][0] = lo[0]; v[2][1] = hi[1]; v[2][2] = hi[2];
       v[3][0] = lo[0]; v[3][1] = lo[1]; v[3][2] = hi[2];
       v[4][0] = hi[0]; v[4][1] = lo[1]; v[4][2] = lo[2];
       v[5][0] = hi[0]; v[5][1] = hi[1]; v[5][2] = lo[2];
       v[6][0] = hi[0]; v[6][1] = hi[1]; v[6][2] = hi[2];
       v[7][0] = hi[0]; v[7][1] = lo[1]; v[7][2] = hi[2];

       for (i=0; i<8; i++){
         Zoltan_Transform_Point(v[i], m, perm, 3, ndims, v[i]);
       }
     }
}

/* 
 * Given the bounds of a 2D or 3D axis-aligned box and a linear
 * transformation, return the bounds of the transformed box,  Set
 * ndims as described above if the box is to be projected to a
 * lower dimension.
 *
 * This is used by the box assign functions when coordinates have been
 * transformed due to degenerate geometries.
 */

void Zoltan_Transform_Box( 
  double *lo, double *hi,  /* input: box bounds, output: bounds of transformed box */
  double (*m)[3],          /* 3x3 transformation */
  int *perm,               /* if transformation is simple coordinate permutation */
  int d,                   /* dimension of box */
  int ndims)               /* dimension of transformed box */
{          
     double v[8][3];
     int i, npoints;

     npoints = ((d == 2) ? 4 : 8);

     Zoltan_Transform_Box_Points(lo, hi, m, perm, d, ndims, v);

     lo[0] = hi[0] = v[0][0];
     lo[1] = hi[1] = 0.0;
     lo[2] = hi[2] = 0.0;

     if (ndims > 1){
       lo[1] = hi[1] = v[0][1];
       if (ndims > 2){
         lo[2] = hi[2] = v[0][2];
       }
     }

     for (i=1; i<npoints; i++){
       if (v[i][0] < lo[0] )     lo[0] = v[i][0];
       else if (v[i][0] > hi[0]) hi[0] = v[i][0];
       if (ndims > 1){
         if (v[i][1] < lo[1])      lo[1] = v[i][1];
         else if (v[i][1] > hi[1]) hi[1] = v[i][1];
         if (ndims > 2){
           if (v[i][2] < lo[2])      lo[2] = v[i][2];
           else if (v[i][2] > hi[2]) hi[2] = v[i][2];
         }
       }
     }
}


int
Zoltan_AllReduceInPlace(void *sndrcvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  int ierr;

#ifndef MPI_IN_PLACE
  void * dummy;
  int size;

  MPI_Type_size(datatype, &size);

  dummy = ZOLTAN_MALLOC(size*count);
  if (dummy == NULL)
    return ZOLTAN_MEMERR;
  memcpy (dummy, sndrcvbuf, size*count);
  ierr = MPI_Allreduce(dummy, sndrcvbuf, count, datatype, op, comm);
  ZOLTAN_FREE(&dummy);
#else /* MPI_IN_PLACE */
  ierr = MPI_Allreduce(MPI_IN_PLACE, sndrcvbuf, count, datatype, op, comm);
#endif /* MPI_IN_PLACE */
  return (ierr);
}

/* The MPI_Datatype for a ZOLTAN_GNO_TYPE can only be determined at run time */

static MPI_Datatype zz_mpi_gno_type;
static char* zz_mpi_gno_name=NULL;
static char* zz_mpi_datatype_names[5] = 
{
"MPI_SHORT",
"MPI_INT",
"MPI_LONG",
"MPI_LONG_LONG_INT",     /* not standard, but shows up */
"MPI_LONG_LONG"             /* standard */
};

MPI_Datatype Zoltan_mpi_gno_type()
{
  int size_short=0, size_int=0, size_long=0, size_long_long=0;

  if (zz_mpi_gno_name){
    return zz_mpi_gno_type;
  }

  MPI_Type_size(MPI_SHORT, &size_short);
  MPI_Type_size(MPI_INT, &size_int);
  MPI_Type_size(MPI_LONG, &size_long);
#ifdef MPI_LONG_LONG
  MPI_Type_size(MPI_LONG_LONG, &size_long_long);
#else
  #ifdef MPI_LONG_LONG_INT
    MPI_Type_size(MPI_LONG_LONG_INT, &size_long_long);
  #endif
#endif

  if (sizeof(ssize_t) == size_short){
    zz_mpi_gno_type = MPI_SHORT;
    zz_mpi_gno_name=zz_mpi_datatype_names[0];
  }
  else if (sizeof(ssize_t) == size_int){
    zz_mpi_gno_type = MPI_INT;
    zz_mpi_gno_name=zz_mpi_datatype_names[1];
  }
  else if (sizeof(ssize_t) == size_long){
    zz_mpi_gno_type = MPI_LONG;
    zz_mpi_gno_name=zz_mpi_datatype_names[2];
  }
  else if (sizeof(ssize_t) == size_long_long){

#ifdef MPI_LONG_LONG
    zz_mpi_gno_type = MPI_LONG_LONG;
    zz_mpi_gno_name=zz_mpi_datatype_names[4];
#else
  #ifdef MPI_LONG_LONG_INT
    zz_mpi_gno_type = MPI_LONG_LONG_INT;
    zz_mpi_gno_name=zz_mpi_datatype_names[3];
  #endif
#endif

  }

  if (!zz_mpi_gno_name){
    /* should never happen */
    fprintf(stderr,"Zoltan_mpi_gno_type: It happened\n");
    MPI_Abort(MPI_COMM_WORLD,99);
  }

  return zz_mpi_gno_type;
}

char *Zoltan_mpi_gno_name()
{
  Zoltan_mpi_gno_type();
  return zz_mpi_gno_name;
}

/* On a linux node, try to write the contents of /proc/meminfo to a file.
 * If committedOnly, then only write the Committed_AS line.  This is the
 * amount of memory that has been granted for memory allocation requests.
 * It may exceed the amount of physical memory, which will cause a fault
 * on a system that doesn't swap if that memory is written to.
 */
void Zoltan_write_linux_meminfo(int append, char *msg, int committedOnly)
{
int rank;
int f, n;
size_t fsize, rc;
char *c=NULL, *next=NULL, *c_end;
char fbuf[64],buf[2048],label[64],value[64],units[64];

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  f = open("/proc/meminfo", O_RDONLY);
  if (f == -1) return;

  c = buf;
  rc = read(f, (void *)c++, 1);

  while ((rc == 1) && (c - buf < 2047)){
    rc = read(f, (void *)c++, 1);
  }

  fsize = c-buf-1;

  close(f);

  sprintf(fbuf,"meminfo_%d.txt",rank);

  if (append){
    f = open(fbuf,O_WRONLY | O_APPEND | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  }
  else{
    f = open(fbuf,O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  }

  if (f == -1) return;

  if (committedOnly){

    c = buf;
    c_end = buf + fsize;

    while( c < c_end){
      next = strchr(c, '\n');
      *next = 0;
      n = sscanf(c, "%s %s %s", label, value, units);
      if (n == 3){
        if (strcmp(label, "Committed_AS:") == 0){
          if (msg != NULL) sprintf(buf,"%s: \t%s \t%s %s\n",msg,label,value,units);
          else             sprintf(buf,"%s %s %s\n",label,value,units);

          fsize = strlen(buf);
          break;
        }
      }
      c = next + 1;
    }
  }
  else{
    if (msg != NULL){
      write(f, msg, strlen(msg));
    }
  }

  write(f,buf,fsize);

  fsync(f);
  close(f);
}

int Zoltan_get_global_id_type(char **name)
{
  if (name){
    *name = zoltan_id_datatype_name;
  }
  return sizeof(ZOLTAN_ID_TYPE);
}

int Zoltan_overflow_test(size_t val)
{
ssize_t mask;

  /* is value too large to store an int */

  if (sizeof(size_t) <= sizeof(int))
    return 0;

  mask = 0xffffffff00000000;

  if ((val & mask) != 0x0000000000000000){
    return 1;
  }

  return 0;
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
