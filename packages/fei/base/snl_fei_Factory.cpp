/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef _snl_fei_Factory_cpp_
#define _snl_fei_Factory_cpp_

#include <fei_macros.hpp>

#include <snl_fei_Factory.hpp>

//----------------------------------------------------------------------------
snl_fei::Factory::Factory(MPI_Comm comm,
                          fei::SharedPtr<LibraryWrapper> wrapper)
  : fei::Factory(comm),
    comm_(comm),
    broker_(),
    matrixGraph_(),
    nodeIDType_(0),
    lsc_(),
    feData_(),
    wrapper_(wrapper),
    outputLevel_(0)
{
  if (wrapper_.get() != NULL) {
    lsc_ = wrapper->getLinearSystemCore();
    feData_ = wrapper->getFiniteElementData();
  }
}

//----------------------------------------------------------------------------
snl_fei::Factory::Factory(MPI_Comm comm,
                          fei::SharedPtr<LinearSystemCore> lsc)
  : fei::Factory(comm),
    comm_(comm),
    broker_(),
    matrixGraph_(),
    nodeIDType_(0),
    lsc_(lsc),
    feData_(),
    wrapper_(),
    outputLevel_(0)
{
}

//----------------------------------------------------------------------------
snl_fei::Factory::Factory(MPI_Comm comm,
                          fei::SharedPtr<FiniteElementData> feData, int nodeIDType)
  : fei::Factory(comm),
    comm_(comm),
    broker_(),
    matrixGraph_(),
    nodeIDType_(nodeIDType),
    lsc_(),
    feData_(feData),
    wrapper_(NULL),
    outputLevel_(0)
{
}

//----------------------------------------------------------------------------
snl_fei::Factory::~Factory()
{
}

#endif

