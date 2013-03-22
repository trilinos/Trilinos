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
  
  std::list<ParamData>::iterator it = data->paramData->begin();
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
