// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

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
