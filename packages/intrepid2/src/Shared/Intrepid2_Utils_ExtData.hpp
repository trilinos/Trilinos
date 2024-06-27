// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_Utils_ExtData.hpp
    \brief  Header file for external data utility functions.
    \author Created by P. Bochev and D. Ridzal and Kyungjoo Kim.
*/

#ifndef __INTREPID2_UTILS_EXTDATA_HPP__
#define __INTREPID2_UTILS_EXTDATA_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Types.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

namespace Intrepid2 {

  /**
     \brief Declarations of templated utility functions

     Note: 
     - NOT compiled on device (host only)
   */
  enum TypeOfExactData{
    INTREPID2_UTILS_FRACTION = 0,
    INTREPID2_UTILS_SCALAR
  };


  /** \brief  Compares the values in the test matrix <var><b>testMat</b></var> to precomputed
      analytic values stored in a file, where the input matrix is an array of arrays.

      \param  inputFile        [in]     -  input file
      \param  testMat          [in]     -  test matrix
      \param  reltol           [in]     -  relative tolerance for equality comparisons
      \param  iprint           [in]     -  if 0, no output; if 1, details are printed
      \param  analyticDataType [in]     -  type of analytic data for comparison:
      \li if INTREPID2_UTILS_FRACTION, analytic fractions are parsed and computed
      \li if INTREPID2_UTILS_SCALAR, high-precision scalar data is read
      \return 0 if pass; error code otherwise
  */
  template<typename ValueType,
           class ...testMatProperties>
  ordinal_type compareToAnalytic( std::ifstream &inputFile,
                                  const Kokkos::DynRankView<ValueType,testMatProperties...> testMat,
                                  const ValueType reltol,
                                  const ordinal_type iprint,
                                  const TypeOfExactData analyticDataType = INTREPID2_UTILS_FRACTION );

  /** \brief  Loads analytic values stored in a file into the matrix <var><b>testMat</b></var>,
      where the output matrix is an array of arrays.

      \param  testMat          [in]     -  test matrix
      \param  inputFile        [in]     -  input file
      \param  analyticDataType [in]     -  type of analytic data for comparison:
      \li if INTREPID2_UTILS_FRACTION, analytic fractions are parsed and computed
      \li if INTREPID2_UTILS_SCALAR, high-precision scalar data is read
  */
  template<typename ValueType,
           class ...testMatProperties>
  void getAnalytic( Kokkos::DynRankView<ValueType,testMatProperties...> testMat,
                    std::ifstream &inputFile,
                    const TypeOfExactData analyticDataType = INTREPID2_UTILS_FRACTION );

} // end namespace Intrepid2

#include "Intrepid2_Utils_ExtDataDef.hpp"

#endif
