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
#include <stdio.h>
#include <MFNRegion.h>
#include <stdlib.h>

int LOCATest(MFNVector, void *, MFErrorHandler);
void MFFreeLOCAData(void*, MFErrorHandler);

void MFFreeLOCAData(void* data, MFErrorHandler err)
{
  LOCAData* locaData = (LOCAData*) data;
  delete locaData;
}

MFNRegion MFNRegionCreateLOCA(LOCAData* data)
{
  MFNRegion loca;

  loca=MFNRegionCreateBaseClass("LOCA", data->mfErrorHandler);
  MFNRegionSetTest(loca,LOCATest, data->mfErrorHandler);
  MFNRegionSetData(loca,(void *)data, data->mfErrorHandler);
  MFNRegionSetFreeData(loca,MFFreeLOCAData, data->mfErrorHandler);

  return(loca);
}

int LOCATest(MFNVector u, void *d, MFErrorHandler err)
{
   
  LOCANVectorData* v_data = (LOCANVectorData *)MFNVectorGetData(u,err);
  LOCAData* data = (LOCAData*) d;
  
  list<ParamData>::iterator it = data->paramData->begin();
  for (unsigned int i=0; i<data->paramData->size(); i++) {
    
    if (v_data->u_ptr->getScalar(i) < it->minValue)
      return 0;
    
    if (v_data->u_ptr->getScalar(i) > it->maxValue)
      return 0;
    
    ++it;
    
  }

  if (v_data->u_ptr->getXVec()->norm(NOX::Abstract::Vector::MaxNorm) > 
      data->solutionMax)
    return 0;
  
  return 1;
}


}
