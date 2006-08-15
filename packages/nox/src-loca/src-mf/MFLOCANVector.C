// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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

#include <MFError.h>
#include <MFNVector.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <MFTime.h>

void MFLOCANVDiff(void*,void*,void*);
void MFLOCANVAdd(void*,void*,void*);
void MFLOCANVPrint(FILE*, void*);
MFNVector MFCloneLOCANVector(void*);
void MFFreeLOCANVectorData(void*);
MFNVector MFCreateLOCANVectorWithData(const Teuchos::RefCountPtr<LMCEV>&);
void MFSetError(int,char*,char*,int,char*);

static char MFNVectorErrorMsg[256]="";

void MFLOCANVDiff(void *adata, void *bdata, void *cdata)
{
  static char RoutineName[]={"MFLOCANVDiff"};
  LOCANVectorData *a_vec_data;
  LOCANVectorData *b_vec_data;
  LOCANVectorData *c_vec_data;
  Teuchos::RefCountPtr<LMCEV> a;
  Teuchos::RefCountPtr<LMCEV> b;
  Teuchos::RefCountPtr<LMCEV> c;

  a_vec_data = (LOCANVectorData*) adata;
  if(a_vec_data==(LOCANVectorData*)NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for a (argument 1) is NULL");
    printf("%s -- %s\n",RoutineName,MFNVectorErrorMsg);fflush(stdout);
    MFSetError(12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  b_vec_data = (LOCANVectorData*) bdata;
  if(b_vec_data==(LOCANVectorData*)NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for b (argument 2) is NULL");
    printf("%s -- %s\n",RoutineName,MFNVectorErrorMsg);fflush(stdout);
    MFSetError(12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  c_vec_data = (LOCANVectorData*) cdata;
  if(c_vec_data==(LOCANVectorData*)NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for c (argument 3) is NULL");
    printf("%s -- %s\n",RoutineName,MFNVectorErrorMsg);fflush(stdout);
    MFSetError(12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  a = a_vec_data->u_ptr;
  b = b_vec_data->u_ptr;
  c = c_vec_data->u_ptr;

  c->update(1.0, *a, -1.0, *b, 0.0);

  return;
}

void MFLOCANVAdd(void *adata,void *bdata,void *cdata)
{
  static char RoutineName[]={"MFLOCANVAdd"};

  LOCANVectorData *a_vec_data;
  LOCANVectorData *b_vec_data;
  LOCANVectorData *c_vec_data;
  Teuchos::RefCountPtr<LMCEV> a;
  Teuchos::RefCountPtr<LMCEV> b;
  Teuchos::RefCountPtr<LMCEV> c;

  a_vec_data = (LOCANVectorData*) adata;
  if(a_vec_data==(LOCANVectorData*)NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for a (argument 1) is NULL");
    printf("%s -- %s\n",RoutineName,MFNVectorErrorMsg);fflush(stdout);
    MFSetError(12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  b_vec_data = (LOCANVectorData*) bdata;
  if(b_vec_data==(LOCANVectorData*)NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for b (argument 2) is NULL");
    printf("%s -- %s\n",RoutineName,MFNVectorErrorMsg);fflush(stdout);
    MFSetError(12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  c_vec_data = (LOCANVectorData*) cdata;
  if(c_vec_data==(LOCANVectorData*)NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for c (argument 3) is NULL");
    printf("%s -- %s\n",RoutineName,MFNVectorErrorMsg);fflush(stdout);
    MFSetError(12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  a = a_vec_data->u_ptr;
  b = b_vec_data->u_ptr;
  c = c_vec_data->u_ptr;

  c->update(1.0, *a, 1.0, *b, 0.0);

  return;
}

void MFLOCANVPrint(FILE *ifp, void *u)
{
}

MFNVector MFCloneLOCANVector(void *data)
{
  LOCANVectorData *vec_data = (LOCANVectorData*) data;
  Teuchos::RefCountPtr<LMCEV> u2 = 
    Teuchos::rcp_dynamic_cast<LMCEV>(vec_data->u_ptr->clone(),true);
  
  return  MFCreateLOCANVectorWithData(u2);
}

void MFFreeLOCANVectorData(void *data)
{
  LOCANVectorData *vec_data = (LOCANVectorData*) data;
  delete vec_data;
}

MFNVector MFCreateLOCANVectorWithData(const Teuchos::RefCountPtr<LMCEV>& u)
{
  MFNVector cthis;

#ifdef MFTIMINGS
  clock_t starttime;

  MFCalledMFCreateNVector++;
  starttime=clock();
#endif

  cthis=MFCreateNVectorBaseClass("LOCA");
  LOCANVectorData *vec_data = new LOCANVectorData(u);

  MFNVectorSetData(cthis,vec_data);

  MFNVectorSetDiff(cthis,MFLOCANVDiff);
  MFNVectorSetAdd(cthis,MFLOCANVAdd);
  MFNVectorSetClone(cthis,MFCloneLOCANVector);
  MFNVectorSetPrint(cthis,MFLOCANVPrint);
  MFNVectorSetFreeData(cthis,MFFreeLOCANVectorData);

#ifdef MFTIMINGS
    MFTimeMFCreateNVector+=clock()-starttime;
#endif
  return cthis;
 }

}  /*extern C*/
