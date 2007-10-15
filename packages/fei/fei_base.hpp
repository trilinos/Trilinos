/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_base_hpp_
#define _fei_base_hpp_

#include <base/fei_macros.hpp>

#include <base/fei_ParameterSet.hpp>

// snl_fei::Factory header in turn includes fei::Factory header
#include <base/snl_fei_Factory.hpp>

#include <base/fei_VectorSpace.hpp>
#include <base/fei_MatrixGraph.hpp>

// fei::Vector_Impl header in turn includes fei::Vector header
#include <base/fei_Vector_Impl.hpp>

// fei::Matrix_Impl header in turn includes fei::Matrix header
#include <base/fei_Matrix_Impl.hpp>

#include <base/fei_LinearSystem.hpp>

#include <base/fei_utils.hpp>

//fei::FEI_Impl header in turn includes FEI.hpp header
#include <base/fei_FEI_Impl.hpp>

#endif // _fei_base_hpp_

