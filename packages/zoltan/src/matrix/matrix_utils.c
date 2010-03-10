/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2009 Sandia National Laboratories.                          *
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

#include <math.h>
#include "zz_const.h"
#include "zz_util_const.h"
#include "zoltan_dd.h"
#include "phg.h"
#include "matrix.h"


static int zoltan_matrix_msg_tag = 19570222;

void
Zoltan_Matrix_Free(Zoltan_matrix *m)
{
  if (m->yend != m->ystart + 1)
    ZOLTAN_FREE(&m->yend);
  ZOLTAN_FREE(&m->ystart);
  ZOLTAN_FREE(&m->yGNO);
  ZOLTAN_FREE(&m->pinGNO);
  ZOLTAN_FREE(&m->pinwgt);
  ZOLTAN_FREE(&m->yGID);
  ZOLTAN_FREE(&m->ywgt);
  if (m->ypid != m->xpid)
    ZOLTAN_FREE(&m->ypid);
  ZOLTAN_FREE(&m->xpid);
  if (m->yoffset != m->xoffset)
    ZOLTAN_FREE(&m->yoffset);
  ZOLTAN_FREE(&m->xoffset);
  ZOLTAN_FREE(&m->ybipart);

  memset (m, 0, sizeof(Zoltan_matrix));
}

void
Zoltan_Matrix2d_Free(Zoltan_matrix_2d *m)
{
  Zoltan_Matrix_Free (&m->mtx);

  ZOLTAN_FREE(&m->dist_x);
  ZOLTAN_FREE(&m->dist_y);

  Zoltan_PHGComm_Destroy(m->comm);
}

void
Zoltan_Matrix_Reset(Zoltan_matrix* m)
{
  m->yGNO = NULL;
  m->ystart = NULL;
  m->yend = NULL;
  m->pinGNO = NULL;
  m->pinwgt = NULL;
  m->ywgt = NULL;
  m->yGID = NULL;
  m->ypid = NULL;
  m->yoffset = NULL;
  m->ybipart = NULL;
  m->xpid = NULL;
  m->xoffset = NULL;
}


void
Zoltan_Matrix2d_Init(Zoltan_matrix_2d *m)
{
  memset(m, 0, sizeof(Zoltan_matrix_2d));

  Zoltan_Distribute_Set(m, &Zoltan_Distribute_Linear, (void*)m);
}

int
Zoltan_Distribute_Set(Zoltan_matrix_2d* mat,
		      distFnct *hashDistFct, void * hashDistData)
{
  mat->hashDistData = hashDistData;
  mat->hashDistFct = hashDistFct;
  return (ZOLTAN_OK);
}


int
Zoltan_Matrix_Complete(ZZ* zz,Zoltan_matrix* m)
{
  static char *yo = "Zoltan_Matrix_Complete";
  int ierr = ZOLTAN_OK;

  if(m->completed)
    return (ZOLTAN_OK);

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (m->yend != m->ystart +1 ) {/* Not in compact mode yet */
    int y;
    /* I have to also rewrites all the pinarrays */

    for (y = 1 ; y <= m->nY ; ++y) {
      int length;
      if (m->ystart[y] == m->yend[y-1]) /* No hole */
	continue;
      length = m->yend[y]-m->ystart[y];
      memcpy(m->pinGNO+m->yend[y-1], m->pinGNO+m->ystart[y], length*sizeof(int));
      memcpy(m->pinwgt+m->yend[y-1]*m->pinwgtdim,
	     m->pinGNO+m->ystart[y]*m->pinwgtdim, length*sizeof(float)*m->pinwgtdim);
      m->ystart[y] = m->yend[y-1];
      m->yend[y] = m->ystart[y] + length;
    }

    ZOLTAN_FREE(&m->yend);
    m->yend = m->ystart + 1;
  }

  m->completed = 1;
 End:
  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}


/* Create necessary data to do projection from user distributed data to zoltan ones
 */
static int
Zoltan_Matrix_Project_Create_Kernel(ZOLTAN_COMM_OBJ **plan, MPI_Comm communicator,
				    int n, int* inpid, int *inoffset,
				    int *outsize, int **output)
{
  int ierr = ZOLTAN_OK;

  zoltan_matrix_msg_tag --;

  ierr = Zoltan_Comm_Create(plan, n, inpid, communicator, zoltan_matrix_msg_tag, outsize);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
    return ierr;

  *output = (int*)ZOLTAN_MALLOC((*outsize)*sizeof(int));
  if ((*outsize) > 0 && (*output) == NULL) return ZOLTAN_MEMERR;

  ierr = Zoltan_Comm_Do(*plan, zoltan_matrix_msg_tag, (char *) inoffset, sizeof(int),
		 (char *) (*output));
  return (ierr);
}


int
Zoltan_Matrix_Project_Create(ZZ* zz, Zoltan_matrix* m)
{
  int ierr = ZOLTAN_OK;
/* Build comm_plan to be able to import/export objs informations */
  ierr = Zoltan_Matrix_Project_Create_Kernel(&m->planX, zz->Communicator,
				      m->nX, m->xpid, m->xoffset,
				      &m->xNSend, &m->xsend);
  if (m->ypid != m->xpid) {
    /* Do specific plan */
    ierr = Zoltan_Matrix_Project_Create_Kernel(&m->planY, zz->Communicator,
					m->nY, m->ypid, m->yoffset,
					&m->yNSend, &m->ysend);
  }
  else {
    m->yNSend = m->xNSend;
    m->planY = m->planX;
    m->ysend = m->xsend;
  }

  return (ierr);
}

/* Project obj informations from the original (queries) distribution to the actual matrix */
int
Zoltan_Matrix_Project_Forward(ZZ *zz, Zoltan_matrix_2d *m, char* inputX, char* inputY, int elemsize,
			      char *outputX, char *outputY)
{
  int i,j;
  char *tmp;
  int ierr = ZOLTAN_OK;

  if (m->mtx.planX == NULL)
    Zoltan_Matrix_Project_Create(zz, &m->mtx);

  for (j = 0 ; j < 2 ; ++j) {
    char *in, *out;
    ZOLTAN_COMM_OBJ *plan;
    int nSend;
    int *send_ptr;
    if (j == 0){
      in = inputX;
      out = outputX;
      plan = m->mtx.planX;
      nSend = m->mtx.xNSend;
      send_ptr = m->mtx.xsend;
    }
    else {
      if (m->mtx.planY == m->mtx.planX) /* same thing than X ! */
	break;
      in = inputY;
      out = outputY;
      plan = m->mtx.planY;
      nSend = m->mtx.yNSend;
      send_ptr = m->mtx.ysend;
    }
    if (in == NULL || out == NULL || plan == NULL)
      continue;

    tmp = (char*) ZOLTAN_MALLOC(nSend*elemsize);
    if (nSend >0 && tmp == NULL) return (ZOLTAN_MEMERR);
    for (i = 0 ; i < nSend ; ++i) {
      memcpy(tmp+i*elemsize, in+send_ptr[i]*elemsize, elemsize);
    }
    zoltan_matrix_msg_tag --;
    ierr = Zoltan_Comm_Do_Reverse(plan, zoltan_matrix_msg_tag, (char *) tmp, sizeof(int), NULL,
				  (char *) out);
    ZOLTAN_FREE(&tmp);
  }
  return (ierr);
}

int
Zoltan_Matrix_Project_Backward(ZZ *zz, Zoltan_matrix_2d *m, char* infoX, char* infoY, int elemsize,
			       char *outputX, char *outputY)
{
  int i,j;
  char *tmp;
  int ierr = ZOLTAN_OK;

  if (m->mtx.planX == NULL)
    Zoltan_Matrix_Project_Create(zz, &m->mtx);

  for (j = 0 ; j < 2 ; ++j) {
    char *in, *out;
    ZOLTAN_COMM_OBJ *plan;
    int nSend;
    int *send_ptr;
    if (j == 0){
      in = infoX;
      out = outputX;
      plan = m->mtx.planX;
      nSend = m->mtx.xNSend;
      send_ptr = m->mtx.xsend;
    }
    else {
      in = infoY;
      out = outputY;
      plan = m->mtx.planY;
      nSend = m->mtx.yNSend;
      send_ptr = m->mtx.ysend;
    }
    if (in == NULL || out == NULL || plan == NULL)
      continue;

    tmp = (char*) ZOLTAN_MALLOC(nSend*elemsize);
    if (nSend >0 && tmp == NULL) return (ZOLTAN_MEMERR);
    zoltan_matrix_msg_tag --;
    ierr = Zoltan_Comm_Do(plan, zoltan_matrix_msg_tag, (char *) in, sizeof(int),
			  (char *) tmp);
    for (i = 0 ; i < nSend ; ++i) {
      memcpy(out+send_ptr[i]*elemsize, tmp+i*elemsize, elemsize);
    }

    ZOLTAN_FREE(&tmp);
  }
  return (ierr);
}


int
Zoltan_Matrix2d_adjproc (ZZ* zz, const Zoltan_matrix_2d * const mat, int **adjproc)
{
  static char *yo = "Zoltan_Matrix2d_adjproc";
  int ierr = ZOLTAN_OK;
  int i;
  ZOLTAN_TRACE_ENTER(zz, yo);

  *adjproc = (int*) ZOLTAN_MALLOC(mat->mtx.nPins*sizeof(int));
  if (mat->mtx.nPins && (*adjproc == NULL))
    MEMORY_ERROR;

  for (i = 0 ; i < mat->mtx.nPins ; ++i ) {
    (*adjproc)[i] = EDGE_TO_PROC_Y(mat, mat->mtx.pinGNO[i]);
  }

 End:
  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}



#ifdef __cplusplus
}
#endif
