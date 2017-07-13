// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_Utils_ExtData.hpp
    \brief  Intrepid utilities for external data.
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
     Note: 
     - NOT compiled on device (host only)
   */

  /***************************************************************************************************
   **                                                                                               **
   **                      Declarations of templated utility functions                              **
   **                                                                                               **
   ***************************************************************************************************/

  enum TypeOfExactData{
    INTREPID2_UTILS_FRACTION = 0,
    INTREPID2_UTILS_SCALAR
  };

  /***************************************************************************************************
   *                                                                                                 *
   *               Utility functions for handling external data in tests                             *
   *                                                                                                 *
   ***************************************************************************************************/

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
