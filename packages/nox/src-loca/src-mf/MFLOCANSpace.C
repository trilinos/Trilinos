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

#include <MFNSpace.h>
#include <MFNVector.h>
#include <MFError.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <MFTime.h>

void MFSetError(int,char*,char*,int,char*);

double MFLOCANSpaceDistance(MFNSpace,MFNVector,MFNVector,void*);
void MFLOCANSpaceDirection(MFNSpace,MFNVector,MFNVector,MFNVector,void*);
void MFLOCANSpaceAdd(MFNSpace,MFNVector,MFNVector,MFNVector,void*);
void MFLOCANSpaceScale(MFNSpace,double,MFNVector,MFNVector,void*);
double MFLOCANSpaceInner(MFNSpace,MFNVector,MFNVector,void*);

MFNSpace MFCreateLOCANSpace(LOCAData* data)
{
  MFNSpace cthis;

#ifdef MFTIMINGS
  clock_t starttime;

  MFCalledMFCreateNSpace++;
  starttime=clock();
#endif

  cthis=MFCreateNSpaceBaseClass("LOCANSpace");
  MFNSpaceSetData(cthis, (void*)data);
  MFNSpaceSetDistance(cthis,MFLOCANSpaceDistance);
  MFNSpaceSetInnerProduct(cthis,MFLOCANSpaceInner);
  MFNSpaceSetDirection(cthis,MFLOCANSpaceDirection);
  MFNSpaceSetAdd(cthis,MFLOCANSpaceAdd);
  MFNSpaceSetScale(cthis,MFLOCANSpaceScale);

#ifdef MFTIMINGS
    MFTimeMFCreateNSpace+=clock()-starttime;
#endif
  return cthis;
}

double MFLOCANSpaceDistance(MFNSpace cthis,MFNVector v0,MFNVector v1,void *d)
{
  double result;
  MFNVector dv;

#ifdef MFTIMINGS
  clock_t starttime;

  MFCalledMFNSpaceDistance++;
  starttime=clock();
#endif
  dv=MFCloneNVector(v0);
  MFLOCANSpaceDirection(cthis,v0,v1,dv,d);
  result=MFLOCANSpaceInner(cthis,dv,dv,d);
  MFFreeNVector(dv);

#ifdef MFTIMINGS
    MFTimeMFNSpaceDistance+=clock()-starttime;
#endif
  return sqrt(result);
}

void MFLOCANSpaceDirection(MFNSpace cthis,MFNVector v0,MFNVector v1,MFNVector diff,void *d)
{
#ifdef MFTIMINGS
  clock_t starttime;

  MFCalledMFNSpaceDirection++;
  starttime=clock();
#endif

  MFNVDiff(v1,v0,diff);

#ifdef MFTIMINGS
    MFTimeMFNSpaceDirection+=clock()-starttime;
#endif
  return;
}

double MFLOCANSpaceInner(MFNSpace cthis,MFNVector v0,MFNVector v1,void *d)
{

#ifdef MFTIMINGS
  clock_t starttime;

  MFCalledMFNSpaceInner++;
  starttime=clock();
#endif

  LOCAData* data = (LOCAData*) d;

  LOCANVectorData *v0_data = (LOCANVectorData *) MFNVectorGetData(v0);
  LOCANVectorData *v1_data = (LOCANVectorData *) MFNVectorGetData(v1);

  double dotp = data->grp->computeScaledDotProduct(*(v0_data->u_ptr),
						   *(v1_data->u_ptr));

#ifdef MFTIMINGS
    MFTimeMFNSpaceInner+=clock()-starttime;
#endif
  return dotp;
}

void MFLOCANSpaceAdd(MFNSpace cthis,MFNVector v0,MFNVector v1,MFNVector sum,
		     void *d)
{

#ifdef MFTIMINGS
  clock_t starttime;

  MFCalledMFNSpaceAdd++;
  starttime=clock();
#endif

  MFNVAdd(v0,v1,sum);

#ifdef MFTIMINGS
    MFTimeMFNSpaceAdd+=clock()-starttime;
#endif
  return;
}

void MFLOCANSpaceScale(MFNSpace cthis,double s, MFNVector v,MFNVector w,
		       void *d)
{
#ifdef MFTIMINGS
  clock_t starttime;

  MFCalledMFNSpaceScale++;
  starttime=clock();
#endif

  LOCANVectorData* u_data = (LOCANVectorData *) MFNVectorGetData(v);
  LOCANVectorData* ur_data = (LOCANVectorData *) MFNVectorGetData(w);
  ur_data->u_ptr->update(s, *(u_data->u_ptr), 0.0);

#ifdef MFTIMINGS
    MFTimeMFNSpaceScale+=clock()-starttime;
#endif
  return;
}

}
