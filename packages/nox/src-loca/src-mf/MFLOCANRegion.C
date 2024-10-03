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
