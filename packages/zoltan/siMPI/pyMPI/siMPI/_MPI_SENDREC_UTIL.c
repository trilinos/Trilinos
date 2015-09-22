/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:54 $
 *    Revision: 1.3 $
 ****************************************************************************/
/***********************************************************************************************/
/* FILE  **************************     _MPI_SENDREC_UTIL.c    *********************************/
/***********************************************************************************************/
/* Author : Lisa Alano June 26 2002                                                            */
/* Copyright (c) 2002 University of California Regents                                         */
/***********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"


/*=============================================================================================*/
int _MPI_checks (void *message, int count, MPI_Datatype data, int dest, int tag, MPI_Comm comm)
{
   int buf,cnt,dat,des,tad,com;
  _MPI_COVERAGE();
   buf=_MPI_checkBuffer(message);
   cnt=_MPI_checkCount(count);
   dat=_MPI_checkDatatype(data);
   des=_MPI_checkDestination(dest);
   tad=_MPI_checkTag(tag);
   com=_MPI_checkCommunicator(comm);

   if (buf == MPI_SUCCESS) {
     if (cnt == MPI_SUCCESS) {
       if (dat == MPI_SUCCESS) {
         if (des == MPI_SUCCESS) {
           if (tad == MPI_SUCCESS) {
             if (com == MPI_SUCCESS) {
               return (int)MPI_SUCCESS;
             }
             return MPI_ERR_COMM;
           }
           return MPI_ERR_TAG;
         }
         return MPI_ERR_RANK;
       }
       return MPI_ERR_TYPE;
     }
     return MPI_ERR_COUNT;
   }
   return MPI_ERR_BUFFER;
}
/*=============================================================================================*/
int _MPI_checkBuffer(void* message)
{
  if(message == (void *)NULL)
     return _MPI_NOT_OK;
  return MPI_SUCCESS;
}
/*=============================================================================================*/
int _MPI_checkCount(int count)
{
  if(count >= 0)
     return MPI_SUCCESS;
  return MPI_ERR_COUNT;
}
/*=============================================================================================*/

int _MPI_checkDatatype(MPI_Datatype type)
{
  if ( (type>=_MPI_TYPE_START) && (type<=_MPI_TYPE_END+_MPI_TYPE_COUNT) )
    return MPI_SUCCESS;
  return MPI_ERR_TYPE;
}

/*=============================================================================================*/
int _MPI_checkDestination(int dest)
{
  if ( (dest==_MPI_RANK)||dest==(MPI_ANY_SOURCE) )
    return MPI_SUCCESS;
  return MPI_ERR_RANK;
}
/*=============================================================================================*/
int _MPI_checkTag(int tag)
{
  if ( (tag == MPI_ANY_TAG) || (tag>=0) )
    return MPI_SUCCESS;
  return MPI_ERR_TAG;
}
/*=============================================================================================*/
int _MPI_checkCommunicator(MPI_Comm comm)
{
int index;
  return _MPI_Comm_check_legal(comm, &index);
}
/*=============================================================================================*/
int _MPI_calculateSize(int count, MPI_Datatype type)
{
  int size;
  if (count == 0)
    return 0;
  size = _MPI_getSize(type);
  return size*count;
}
/*=============================================================================================*/
/* EDIT for other types ========================= */
/* Not included: MPI_PACKED and FORTRAN DATATYPES */
/* ============================================== */
int _MPI_getSize(MPI_Datatype type)
{
  switch(type)
  {
    case MPI_CHAR:
      return sizeof(char);
    case MPI_BYTE:
      return sizeof(char);
    case MPI_SHORT:
      return sizeof(short);
    case MPI_INT:
      return sizeof(int);
    case MPI_LONG:
      return sizeof(long);
    case MPI_FLOAT:
      return sizeof(float);
    case MPI_DOUBLE:
      return sizeof(double);
    case MPI_UNSIGNED_CHAR:
      return sizeof(unsigned char);
    case MPI_UNSIGNED_SHORT:
      return sizeof(unsigned short);
    case MPI_UNSIGNED:
      return sizeof(unsigned);
    case MPI_UNSIGNED_LONG:
      return sizeof(unsigned long);
    case MPI_LONG_DOUBLE:
      return sizeof(long double);
    case MPI_UB:
      return sizeof(MPI_Aint);
    case MPI_LB:
      return sizeof(MPI_Aint);
    case MPI_FLOAT_INT:
      return sizeof(_MPI_FLOAT_INT);
    case MPI_LONG_INT:
      return sizeof(_MPI_LONG_INT);
    case MPI_DOUBLE_INT:
      return sizeof(_MPI_DOUBLE_INT);
    case MPI_2INT:
      return sizeof(_MPI_2INT);
    case MPI_LONG_DOUBLE_INT:
      return sizeof(_MPI_LONG_DOUBLE_INT);
    case MPI_LONG_LONG_INT:
      return sizeof(_MPI_LONG_LONG_INT);
    default:
      return _MPI_calculateStructureSize(type);
  }
}
/*=============================================================================================*/
int _MPI_calculateStructureSize (MPI_Datatype type)
{
  int index, size;
  _MPI_TYPE_DES* datatype;

  index = _MPI_FindType(type);
  if (index == _MPI_NOT_OK)
  {
    _MPI_ERR_ROUTINE (MPI_ERR_TYPE, "MPI: Datatype error. No such datatype exists.");
    MPI_Abort (MPI_COMM_NULL, MPI_ERR_TYPE);
  }
  datatype = &_MPI_TYPE_LIST[index];
  size = datatype->size;
  /*EDIT should be extent */

  return size;
}
/*=============================================================================================*/
int _MPI_checkRequest (MPI_Request request)
{
  if (request == MPI_REQUEST_NULL) return MPI_ERR_REQUEST;
  return MPI_SUCCESS;
}
/*=============================================================================================*/
/* Same as _MPI_Buff_Insert below. */
/* This is for Buffered Send or Bsend - where you use the user's buffer*/
int _MPI_Bbuff_Insert(void *message,int count,MPI_Datatype datatype,int tag,MPI_Comm comm)
{
  int index;
  for (index=0; index < _MPI_DATA_ARRAY_SIZE ; index++) {
    if (_MPI_DATA_BUFF[index].valid == _MPI_NOT_VALID) {
      _MPI_DATA_BUFF[index].valid = _MPI_VALID;
      _MPI_DATA_BUFF[index].buffer = message;
      _MPI_DATA_BUFF[index].count = count;
      _MPI_DATA_BUFF[index].type = datatype;
      _MPI_DATA_BUFF[index].tag = tag;
      _MPI_DATA_BUFF[index].comm = comm;       /* ----------------------------------- */
      _MPI_DATA_BUFF[index].user = _MPI_TRUE;  /* this is where the DIFFERENCE occurs */
      _MPI_DATA_BUFF_COUNT++;
      return MPI_SUCCESS;
    }
  }

  /* ----------------------------------------------- */
  /* This realloc is OK because internal references  */
  /* to the list never escape.  We just search the   */
  /* whole list for matches.                         */
  /* ----------------------------------------------- */
  _MPI_DATA_BUFF = (_MPI_DATA_ENTRY*) _MPI_safeRealloc
    (_MPI_DATA_BUFF,
     (_MPI_DATA_ARRAY_SIZE+_MPI_PREALLOCATION_SIZE)*sizeof(_MPI_DATA_ENTRY),
     "Error in _MPI_Buff_Insert for reallocation");
  _MPI_DATA_ARRAY_SIZE+=_MPI_PREALLOCATION_SIZE;
  _MPI_DATA_BUFF[index].valid = _MPI_VALID;
  _MPI_DATA_BUFF[index].buffer = message;
  _MPI_DATA_BUFF[index].count = count;
  _MPI_DATA_BUFF[index].type = datatype;
  _MPI_DATA_BUFF[index].tag = tag;
  _MPI_DATA_BUFF[index].comm = comm;
  _MPI_DATA_BUFF[index].user = _MPI_TRUE;
  _MPI_DATA_BUFF_COUNT++;
  return MPI_SUCCESS;
}
/*=============================================================================================*/
/* Same as _MPI_Bbuff_Insert above. */
int _MPI_Buff_Insert(void *message,int count,MPI_Datatype datatype,int tag,MPI_Comm comm) {
  int index;
  for (index=0; index < _MPI_DATA_ARRAY_SIZE ; index++) {
    if (_MPI_DATA_BUFF[index].valid == _MPI_NOT_VALID) {
      _MPI_DATA_BUFF[index].valid = _MPI_VALID;
      _MPI_DATA_BUFF[index].buffer = message;
      _MPI_DATA_BUFF[index].count = count;
      _MPI_DATA_BUFF[index].type = datatype;
      _MPI_DATA_BUFF[index].tag = tag;
      _MPI_DATA_BUFF[index].comm = comm;
      _MPI_DATA_BUFF[index].user = _MPI_FALSE;
      _MPI_DATA_BUFF_COUNT++;
      return MPI_SUCCESS;
    }
  }

  _MPI_DATA_BUFF = (_MPI_DATA_ENTRY*) _MPI_safeRealloc
    (_MPI_DATA_BUFF,
     (_MPI_DATA_ARRAY_SIZE+_MPI_PREALLOCATION_SIZE)*sizeof(_MPI_DATA_ENTRY), 
     "Error in _MPI_Buff_Insert for reallocation");
  _MPI_DATA_ARRAY_SIZE+=_MPI_PREALLOCATION_SIZE;
  _MPI_DATA_BUFF[index].valid = _MPI_VALID;
  _MPI_DATA_BUFF[index].buffer = message;
  _MPI_DATA_BUFF[index].count = count;
  _MPI_DATA_BUFF[index].type = datatype;
  _MPI_DATA_BUFF[index].tag = tag;
  _MPI_DATA_BUFF[index].comm = comm;
  _MPI_DATA_BUFF[index].user = _MPI_FALSE;
  _MPI_DATA_BUFF_COUNT++;
  return MPI_SUCCESS;
}

/*=============================================================================================*/
int _MPI_Data_Invalid(int index) {
  if ( (index>=0)&&(index<_MPI_DATA_ARRAY_SIZE) ) {
    if(_MPI_DATA_BUFF[index].user==_MPI_TRUE) {
      ;
    } else {
      _MPI_safeFree(_MPI_DATA_BUFF[index].buffer,"BUFF buffer"); 
      _MPI_DATA_BUFF[index].buffer = 0;
    }

    _MPI_DATA_BUFF[index].valid = _MPI_NOT_VALID;
    _MPI_DATA_BUFF[index].buffer = (void *) NULL;
    _MPI_DATA_BUFF[index].count = _MPI_NOT_VALID;
    _MPI_DATA_BUFF[index].type = _MPI_NOT_VALID;
    _MPI_DATA_BUFF[index].tag = _MPI_NOT_VALID;
    _MPI_DATA_BUFF[index].comm = MPI_COMM_NULL;
    _MPI_DATA_BUFF[index].user = _MPI_NOT_VALID;

    return MPI_SUCCESS;
  }

  return _MPI_NOT_OK;
}
/*=============================================================================================*/
int _MPI_Type_Invalid(int index)
{
  if ( (index>=0)&&(index<_MPI_TYPE_ARRAY_SIZE) )
  {
    _MPI_TYPE_LIST[index].id = _MPI_NOT_VALID;
    _MPI_TYPE_LIST[index].size = 0;
    _MPI_TYPE_LIST[index].extent = 0;
    _MPI_TYPE_LIST[index].next = (_MPI_TYPE_DES *) NULL;
    return MPI_SUCCESS;
  }
  return _MPI_NOT_OK;
}
/*=============================================================================================*/
int _MPI_Buff_Ifind(int tag, MPI_Comm comm) {
  int index;
  for (index=0; index<_MPI_DATA_BUFF_COUNT; index++) {
    if(_MPI_DATA_BUFF[index].valid == _MPI_VALID) {
      if (_MPI_DATA_BUFF[index].comm == comm) {
        if ( (_MPI_DATA_BUFF[index].tag == tag)||(tag == MPI_ANY_TAG) ) {
          return index;
        }
      }
    }
  }
  return _MPI_NOT_OK;
}

/**************************************************************************/
/* GLOBAL **************       _MPI_Buff_Find      ************************/
/**************************************************************************/
/* Match the tag to a current message                                     */
/**************************************************************************/
int _MPI_Buff_Find(int tag, MPI_Comm comm) {
  int index;
  for (index=0; index<_MPI_DATA_BUFF_COUNT; index++) {
    if(_MPI_DATA_BUFF[index].valid == _MPI_VALID) {
      if (_MPI_DATA_BUFF[index].comm == comm) {
        if ( (_MPI_DATA_BUFF[index].tag == tag)||(tag == MPI_ANY_TAG) ) {
#ifdef _MPI_DEBUG
          ;/*printf("Found at index=%d\n",index);*/
#endif
          return index;
        }
      }
    }
  }

  return _MPI_NOT_OK;
}
/*==========================================================================*/

