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
