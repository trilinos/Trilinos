/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
#include "mpi.h"

_MPI_COMM_IMPL* _MPI_COMM_LIST = 0;
_MPI_DATA_ENTRY* _MPI_DATA_BUFF = 0;
_MPI_TYPE_DES* _MPI_TYPE_LIST = 0;
_MPI_OP_TYPE* _MPI_OP_LIST = 0;
_MPI_REQUEST_OBJECT** _MPI_REQ_LIST_OF_LISTS = 0;

int _MPI_INIT_STATUS = 0;
int _MPI_FINALIZED_FLAG = 0;
int _MPI_INITIALIZED_FLAG = 0;
int _MPI_TYPE_COUNT = 0;
int _MPI_DATA_BUFF_COUNT = 0;
int _MPI_COMM_COUNT = 0;
int _MPI_OP_COUNT = 0;
int _MPI_REQ_COUNT = 0;

int _MPI_COMM_ARRAY_SIZE = 0;
int _MPI_DATA_ARRAY_SIZE = 0;
int _MPI_TYPE_ARRAY_SIZE = 0;
int _MPI_OP_ARRAY_SIZE = 0;
int _MPI_REQ_ARRAY_SIZE = 0;

MPI_Request* _MPI_REQNULL;
