/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_macros.hpp"
#include "fei_LibraryWrapper.hpp"
#include "fei_LinearSystemCore.hpp"
#include "fei_FiniteElementData.hpp"
#include <cstdlib>

LibraryWrapper::LibraryWrapper(fei::SharedPtr<LinearSystemCore> lsc)
 : haveLinearSystemCore_(true),
   haveFiniteElementData_(false),
   lsc_(lsc),
   feData_()
{
}

LibraryWrapper::LibraryWrapper(fei::SharedPtr<FiniteElementData> feData)
 : haveLinearSystemCore_(false),
   haveFiniteElementData_(true),
   lsc_(NULL),
   feData_(feData)
{
}

LibraryWrapper::~LibraryWrapper()
{
}

