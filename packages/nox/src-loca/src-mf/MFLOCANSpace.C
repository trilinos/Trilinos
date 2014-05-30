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

#include <MFNSpace.h>
#include <MFNVector.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//#include <MFTime.h>

double MFLOCANSpaceDistance(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler);
void MFLOCANSpaceDirection(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler);
void MFLOCANSpaceAdd(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler);
void MFLOCANSpaceScale(MFNSpace,double,MFNVector,MFNVector,void*,MFErrorHandler);
double MFLOCANSpaceInner(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler);

MFNSpace MFCreateLOCANSpace(LOCAData* data)
{
  MFNSpace cthis;

#ifdef MFTIMINGS
  clock_t starttime;

  MFCalledMFCreateNSpace++;
  starttime=clock();
#endif

  cthis=MFCreateNSpaceBaseClass("LOCANSpace", data->mfErrorHandler);
  MFNSpaceSetData(cthis, (void*)data, data->mfErrorHandler);
  MFNSpaceSetDistance(cthis,MFLOCANSpaceDistance, data->mfErrorHandler);
  MFNSpaceSetInnerProduct(cthis,MFLOCANSpaceInner, data->mfErrorHandler);
  MFNSpaceSetDirection(cthis,MFLOCANSpaceDirection, data->mfErrorHandler);
  MFNSpaceSetAdd(cthis,MFLOCANSpaceAdd, data->mfErrorHandler);
  MFNSpaceSetScale(cthis,MFLOCANSpaceScale, data->mfErrorHandler);

#ifdef MFTIMINGS
    MFTimeMFCreateNSpace+=clock()-starttime;
#endif
  return cthis;
}

double MFLOCANSpaceDistance(MFNSpace cthis,MFNVector v0,MFNVector v1,void *d,
                MFErrorHandler err)
{
  double result;
  MFNVector dv;

#ifdef MFTIMINGS
  clock_t starttime;

  MFCalledMFNSpaceDistance++;
  starttime=clock();
#endif
  dv=MFCloneNVector(v0,err);
  MFLOCANSpaceDirection(cthis,v0,v1,dv,d,err);
  result=MFLOCANSpaceInner(cthis,dv,dv,d,err);
  MFFreeNVector(dv,err);

#ifdef MFTIMINGS
    MFTimeMFNSpaceDistance+=clock()-starttime;
#endif
  return sqrt(result);
}

void MFLOCANSpaceDirection(MFNSpace cthis,MFNVector v0,MFNVector v1,MFNVector diff,void *d, MFErrorHandler err)
{
#ifdef MFTIMINGS
  clock_t starttime;

  MFCalledMFNSpaceDirection++;
  starttime=clock();
#endif

  MFNVDiff(v1,v0,diff,err);

#ifdef MFTIMINGS
    MFTimeMFNSpaceDirection+=clock()-starttime;
#endif
  return;
}

double MFLOCANSpaceInner(MFNSpace cthis,MFNVector v0,MFNVector v1,void *d, MFErrorHandler err)
{

#ifdef MFTIMINGS
  clock_t starttime;

  MFCalledMFNSpaceInner++;
  starttime=clock();
#endif

  LOCAData* data = (LOCAData*) d;

  LOCANVectorData *v0_data = (LOCANVectorData *) MFNVectorGetData(v0,err);
  LOCANVectorData *v1_data = (LOCANVectorData *) MFNVectorGetData(v1,err);

  double dotp = data->grp->computeScaledDotProduct(*(v0_data->u_ptr),
                           *(v1_data->u_ptr));

#ifdef MFTIMINGS
    MFTimeMFNSpaceInner+=clock()-starttime;
#endif
  return dotp;
}

void MFLOCANSpaceAdd(MFNSpace cthis,MFNVector v0,MFNVector v1,MFNVector sum,
             void *d, MFErrorHandler err)
{

#ifdef MFTIMINGS
  clock_t starttime;

  MFCalledMFNSpaceAdd++;
  starttime=clock();
#endif

  MFNVAdd(v0,v1,sum,err);

#ifdef MFTIMINGS
    MFTimeMFNSpaceAdd+=clock()-starttime;
#endif
  return;
}

void MFLOCANSpaceScale(MFNSpace cthis,double s, MFNVector v,MFNVector w,
               void *d, MFErrorHandler err)
{
#ifdef MFTIMINGS
  clock_t starttime;

  MFCalledMFNSpaceScale++;
  starttime=clock();
#endif

  LOCANVectorData* u_data = (LOCANVectorData *) MFNVectorGetData(v,err);
  LOCANVectorData* ur_data = (LOCANVectorData *) MFNVectorGetData(w,err);
  ur_data->u_ptr->update(s, *(u_data->u_ptr), 0.0);

#ifdef MFTIMINGS
    MFTimeMFNSpaceScale+=clock()-starttime;
#endif
  return;
}

}
