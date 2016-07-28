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

/** \file   Intrepid_CubaturePolylib.hpp
    \brief  Header file for the Intrepid2::CubaturePolylib class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CUBATURE_POLYLIB_HPP__
#define __INTREPID2_CUBATURE_POLYLIB_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Cubature.hpp"

#include "Intrepid2_Polylib.hpp"


namespace Intrepid2 {

  namespace Impl {
    class CubaturePolylib {
    public:

      static ordinal_type
      getNumPoints( const ordinal_type degree,
                    const EPolyType polytype );

      template<typename cubPointValueType,  class ...cubPointProperties,
               typename cubWeightValueType, class ...cubWeightProperties>
      static void
      getCubature( Kokkos::DynRankView<cubPointValueType, cubPointProperties...>  cubPoints,
                   Kokkos::DynRankView<cubWeightValueType,cubWeightProperties...> cubWeights,
                   const ordinal_type degree,
                   const double alpha,
                   const double beta,
                   const EPolyType polytype );
    };
  }

  /** \class Intrepid2::CubaturePolylib
      \brief Utilizes cubature (integration) rules contained in the library Polylib
      (Spencer Sherwin, Aeronautics, Imperial College London) within Intrepid.

      They are based on zeros of Jacobi polynomials,
      e.g. Legendre (alpha=beta=0, default), Chebyshev (alpha=beta=-0.5), etc.
      They are given on the interval [-1,1] and are optimal with respect to
      the following requirements, yielding 4 subclasses:
      \li Gauss             - no restrictions (all integration points are contained in the open interval (-1,1))
      \li Gauss-Radau-Left  - left-most integration point is fixed at -1
      \li Gauss-Radau-Right - right-most integration point is fixed at +1
      \li Gauss-Lobatto     - left-most and right-most integration points are fixed at -1 and +1, respectively
  */
  template<typename ExecSpaceType = void,
           typename pointValueType = double,
           typename weightValueType = double>
  class CubaturePolylib
    : public Cubature<ExecSpaceType,pointValueType,weightValueType> {
  private:

    /** \brief The degree of polynomials that are integrated
        exactly by this cubature rule.
    */
    ordinal_type degree_;

    /** \brief Dimension of integration domain.
     */
    ordinal_type dimension_;

    /** \brief Type of integration points.
     */
    EPolyType polytype_;

    /** \brief Jacobi parameter alpha and beta.
     */
    double alpha_, beta_;

  public:

    typedef typename Cubature<ExecSpaceType,pointValueType,weightValueType>::pointViewType  pointViewType;
    typedef typename Cubature<ExecSpaceType,pointValueType,weightValueType>::weightViewType weightViewType;

    virtual
    void
    getCubature( pointViewType  cubPoints,
                 weightViewType cubWeights ) const {
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_EXCEPTION( cubPoints.dimension(0)  < this->getNumPoints() ||
                                    cubPoints.dimension(1)  < this->getDimension() ||
                                    cubWeights.dimension(0) < this->getNumPoints(), std::out_of_range,
                                    ">>> ERROR (CubaturePolylib): Insufficient space allocated for cubature points or weights.");
#endif
      Impl::CubaturePolylib::getCubature( cubPoints,
                                          cubWeights,
                                          degree_,
                                          alpha_,
                                          beta_,
                                          polytype_
                                          );
    }

    /** \brief Returns the number of cubature points.
     */
    virtual
    ordinal_type
    getNumPoints() const {
      return Impl::CubaturePolylib::getNumPoints( degree_, polytype_ );
    }

    /** \brief Returns dimension of integration domain.
     */
    virtual
    ordinal_type
    getDimension() const {
      return dimension_;
    }

    /** \brief Returns cubature name.
     */
    virtual
    const char*
    getName() const {
      return "CubaturePolylib";
    }

    /** \brief Returns max. degree of polynomials that are integrated exactly.
        The return vector has size 1.
    */
    ordinal_type
    getAccuracy() const {
      return degree_;
    }

    CubaturePolylib( const ordinal_type degree,
                     const EPolyType polytype = POLYTYPE_GAUSS,
                     const double alpha = 0.0,
                     const double beta = 0.0 );

  };

}

#include "Intrepid2_CubaturePolylibDef.hpp"

#endif
