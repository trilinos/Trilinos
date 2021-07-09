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

/** \file   Intrepid2_IntegrationTools.hpp
    \brief  Header file for the Intrepid2::IntegrationTools class; provides support for structure-aware integration.
    \author Created by Nathan V. Roberts.
*/

#ifndef __INTREPID2_INTEGRATIONTOOLS_HPP__
#define __INTREPID2_INTEGRATIONTOOLS_HPP__

#include "Intrepid2_ConfigDefs.hpp"

#include "Shards_CellTopology.hpp"
#include "Shards_BasicTopologies.hpp"

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_Kernels.hpp"

#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_RealSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"

#include "Intrepid2_Data.hpp"
#include "Intrepid2_TransformedVectorData.hpp"
#include "Intrepid2_VectorData.hpp"

#include "Kokkos_Core.hpp"


namespace Intrepid2 {

  /** \class Intrepid2::IntegrationTools
      \brief Provides support for structure-aware integration.
  */
  template<typename DeviceType = void>
  class IntegrationTools {
  public:
    /** \brief   Allocates storage for the contraction of \a <b>vectorDataLeft</b> and \a <b>vectorDataRight</b> containers on
                 point and space dimensions, weighting each point according to <b>cellMeasures</b>.

        \param  vectorDataRight      [in] - Left input container, with logical shape (C,F,P,D)
        \param  cellMeasures             [in] - Point weight container, with logical shape (C,P)
        \param  vectorDataLeft        [in] - Right input container with logical shape (C,F,P,D)
        \param  sumInto                        [in] - If TRUE, sum into given output array, otherwise overwrite it. Default: FALSE.

        \return <b>integrals</b>, a container with logical shape (C,F,F), suitable for passing as the first argument to the integrate() variant that takes an Intrepid2::Data object as its first, <b>integrals</b>, argument.
    */
    template<class Scalar>
    static Data<Scalar,DeviceType> allocateIntegralData(const TransformedVectorData<Scalar,DeviceType> vectorDataLeft,
                                                        const TensorData<Scalar,DeviceType> cellMeasures,
                                                        const TransformedVectorData<Scalar,DeviceType> vectorDataRight);
    
    /** \brief   Contracts \a <b>vectorDataLeft</b> and \a <b>vectorDataRight</b> containers on
        point and space dimensions, weighting each point according to <b>cellMeasures</b>,
        and stores the result in \a <b>outputValues</b>.  The <b>integrals</b> container can be constructed using allocateIntegralData().

        \param  outputValues          [out] - Output array, with logical shape (C,F,F)
        \param  vectorDataRight      [in] - Left input container, with logical shape (C,F,P,D)
        \param  cellMeasures             [in] - Point weight container, with logical shape (C,P)
        \param  vectorDataLeft        [in] - Right input container with logical shape (C,F,P,D)
        \param  sumInto                        [in] - If TRUE, sum into given output array, otherwise overwrite it. Default: FALSE.
        \param  approxFlops               [in] - if not NULL, the double pointed to will be set with an estimated number of floating point operations.  Intended for performance assessment purposes.
    */
    template<class Scalar>
    static void integrate(Data<Scalar,DeviceType> integrals, const TransformedVectorData<Scalar,DeviceType> &vectorDataLeft,
                          const TensorData<Scalar,DeviceType> &cellMeasures,
                          const TransformedVectorData<Scalar,DeviceType> &vectorDataRight, const bool sumInto = false,
                          double* approximateFlops = NULL);
  }; // end IntegrationTools class

} // end namespace Intrepid2

// include templated definitions
#include <Intrepid2_IntegrationToolsDef.hpp>

#endif
