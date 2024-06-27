// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "Intrepid2_TransformedBasisValues.hpp"
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
    static Data<Scalar,DeviceType> allocateIntegralData(const TransformedBasisValues<Scalar,DeviceType> vectorDataLeft,
                                                        const TensorData<Scalar,DeviceType> cellMeasures,
                                                        const TransformedBasisValues<Scalar,DeviceType> vectorDataRight);
    
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
    static void integrate(Data<Scalar,DeviceType> integrals, const TransformedBasisValues<Scalar,DeviceType> &vectorDataLeft,
                          const TensorData<Scalar,DeviceType> &cellMeasures,
                          const TransformedBasisValues<Scalar,DeviceType> &vectorDataRight, const bool sumInto = false,
                          double* approximateFlops = NULL);
  }; // end IntegrationTools class

} // end namespace Intrepid2

// include templated definitions
#include <Intrepid2_IntegrationToolsDef.hpp>

#endif
