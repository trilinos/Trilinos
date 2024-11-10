// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <MFLOCA.H>

extern "C" {

#include <MFNVector.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//#include <MFTime.h>

void MFLOCANVDiff(void*,void*,void*,MFErrorHandler);
void MFLOCANVAdd(void*,void*,void*,MFErrorHandler);
void MFLOCANVPrint(FILE*, void*,MFErrorHandler);
MFNVector MFCloneLOCANVector(void*,MFErrorHandler);
void MFFreeLOCANVectorData(void*,MFErrorHandler);
MFNVector MFCreateLOCANVectorWithData(const Teuchos::RCP<LMCEV>&,MFErrorHandler);

static char MFNVectorErrorMsg[256]="";

void MFLOCANVDiff(void *adata, void *bdata, void *cdata, MFErrorHandler err)
{
  static char RoutineName[]={"MFLOCANVDiff"};
  LOCANVectorData *a_vec_data;
  LOCANVectorData *b_vec_data;
  LOCANVectorData *c_vec_data;
  Teuchos::RCP<LMCEV> a;
  Teuchos::RCP<LMCEV> b;
  Teuchos::RCP<LMCEV> c;

  a_vec_data = (LOCANVectorData*) adata;
  if(a_vec_data==(LOCANVectorData*)NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for a (argument 1) is NULL");
    printf("%s -- %s\n",RoutineName,MFNVectorErrorMsg);fflush(stdout);
    MFSetError(err, 12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  b_vec_data = (LOCANVectorData*) bdata;
  if(b_vec_data==(LOCANVectorData*)NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for b (argument 2) is NULL");
    printf("%s -- %s\n",RoutineName,MFNVectorErrorMsg);fflush(stdout);
    MFSetError(err, 12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  c_vec_data = (LOCANVectorData*) cdata;
  if(c_vec_data==(LOCANVectorData*)NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for c (argument 3) is NULL");
    printf("%s -- %s\n",RoutineName,MFNVectorErrorMsg);fflush(stdout);
    MFSetError(err, 12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  a = a_vec_data->u_ptr;
  b = b_vec_data->u_ptr;
  c = c_vec_data->u_ptr;

  c->update(1.0, *a, -1.0, *b, 0.0);

  return;
}

void MFLOCANVAdd(void *adata,void *bdata,void *cdata, MFErrorHandler err)
{
  static char RoutineName[]={"MFLOCANVAdd"};

  LOCANVectorData *a_vec_data;
  LOCANVectorData *b_vec_data;
  LOCANVectorData *c_vec_data;
  Teuchos::RCP<LMCEV> a;
  Teuchos::RCP<LMCEV> b;
  Teuchos::RCP<LMCEV> c;

  a_vec_data = (LOCANVectorData*) adata;
  if(a_vec_data==(LOCANVectorData*)NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for a (argument 1) is NULL");
    printf("%s -- %s\n",RoutineName,MFNVectorErrorMsg);fflush(stdout);
    MFSetError(err, 12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  b_vec_data = (LOCANVectorData*) bdata;
  if(b_vec_data==(LOCANVectorData*)NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for b (argument 2) is NULL");
    printf("%s -- %s\n",RoutineName,MFNVectorErrorMsg);fflush(stdout);
    MFSetError(err, 12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  c_vec_data = (LOCANVectorData*) cdata;
  if(c_vec_data==(LOCANVectorData*)NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for c (argument 3) is NULL");
    printf("%s -- %s\n",RoutineName,MFNVectorErrorMsg);fflush(stdout);
    MFSetError(err, 12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  a = a_vec_data->u_ptr;
  b = b_vec_data->u_ptr;
  c = c_vec_data->u_ptr;

  c->update(1.0, *a, 1.0, *b, 0.0);

  return;
}

void MFLOCANVPrint(FILE *ifp, void *u, MFErrorHandler err)
{
}

MFNVector MFCloneLOCANVector(void *data, MFErrorHandler err)
{
  LOCANVectorData *vec_data = (LOCANVectorData*) data;
  Teuchos::RCP<LMCEV> u2 =
    Teuchos::rcp_dynamic_cast<LMCEV>(vec_data->u_ptr->clone(),true);

  return MFCreateLOCANVectorWithData(u2, err);
}

void MFFreeLOCANVectorData(void *data, MFErrorHandler err)
{
  LOCANVectorData *vec_data = (LOCANVectorData*) data;
  delete vec_data;
}

MFNVector MFCreateLOCANVectorWithData(const Teuchos::RCP<LMCEV>& u,
                      MFErrorHandler err)
{
  MFNVector cthis;

#ifdef MFTIMINGS
  clock_t starttime;

  MFCalledMFCreateNVector++;
  starttime=clock();
#endif

  cthis=MFCreateNVectorBaseClass("LOCA",err);
  LOCANVectorData *vec_data = new LOCANVectorData(u);

  MFNVectorSetData(cthis,vec_data,err);

  MFNVectorSetDiff(cthis,MFLOCANVDiff,err);
  MFNVectorSetAdd(cthis,MFLOCANVAdd,err);
  MFNVectorSetClone(cthis,MFCloneLOCANVector,err);
  MFNVectorSetPrint(cthis,MFLOCANVPrint,err);
  MFNVectorSetFreeData(cthis,MFFreeLOCANVectorData,err);

#ifdef MFTIMINGS
    MFTimeMFCreateNVector+=clock()-starttime;
#endif
  return cthis;
}

}  /*extern C*/
