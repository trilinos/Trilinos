/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <snl_fei_Broker_FEData.hpp>

//----------------------------------------------------------------------------
snl_fei::Broker_FEData::Broker_FEData(fei::SharedPtr<FiniteElementData> feData,
			      fei::SharedPtr<fei::MatrixGraph> matrixGraph,
				      int nodeIDType)
  : feData_(feData),
    matrixGraph_(matrixGraph),
    nodeIDType_(nodeIDType),
    setStructure_(false),
    setMatrixMatrixGraph_(false),
    lookup_(NULL)
{
}

//----------------------------------------------------------------------------
snl_fei::Broker_FEData::~Broker_FEData()
{
  delete lookup_;
}

