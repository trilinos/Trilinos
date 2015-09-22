/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <snl_fei_Broker_LinSysCore.hpp>

#undef fei_file
#define fei_file "snl_fei_Broker_LinSysCore.cpp"
#include <fei_ErrMacros.hpp>

//----------------------------------------------------------------------------
snl_fei::Broker_LinSysCore::Broker_LinSysCore(fei::SharedPtr<LinearSystemCore> lsc,
			      fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                              fei::SharedPtr<fei::Reducer> reducer,
                              bool blockMatrix)
  : linsyscore_(lsc),
    matrixGraph_(matrixGraph),
    reducer_(reducer),
    lookup_(NULL),
    setGlobalOffsets_(false),
    numLocalEqns_(0),
    setMatrixStructure_(false),
    blockMatrix_(blockMatrix)
{
  int dummyID = -1;
  lsc->setNumRHSVectors(1, &dummyID);
}

//----------------------------------------------------------------------------
snl_fei::Broker_LinSysCore::~Broker_LinSysCore()
{
  delete lookup_;
}

