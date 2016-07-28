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

/** \file   Intrepid_CubatureTensor.hpp
    \brief  Header file for the Intrepid2::CubatureTensor class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CUBATURE_TENSOR_HPP__
#define __INTREPID2_CUBATURE_TENSOR_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Cubature.hpp"
#include "Intrepid2_CubatureDirect.hpp"

namespace Intrepid2 {

  /** \class Intrepid2::CubatureTensor
      \brief Defines tensor-product cubature (integration) rules in Intrepid.
  */
  template<typename ExecSpaceType = void,
           typename pointValueType = double,
           typename weightValueType = double>
  class CubatureTensor
    : public Cubature<ExecSpaceType,pointValueType,weightValueType> {
  private:

    /** \brief Array of cubature rules.
     */
    ordinal_type numCubatures_;

    CubatureDirect<ExecSpaceType,pointValueType,weightValueType> cubatures_[Parameters::MaxDimension];

    /** \brief Dimension of integration domain.
     */
    ordinal_type dimension_;

  public:

    class Internal {
    private:
      CubatureTensor *obj_;

    public:
      Internal(CubatureTensor *obj)
        : obj_(obj) {}

      /** \brief Returns cubature points and weights
          (return arrays must be pre-sized/pre-allocated).

          \param cubPoints       [out]     - Vector containing the cubature points.
          \param cubWeights      [out]     - Vector of corresponding cubature weights.
      */

      template<typename cubPointValueType,  class ...cubPointProperties,
               typename cubWeightValueType, class ...cubWeightProperties>
      void
      getCubature( Kokkos::DynRankView<cubPointValueType, cubPointProperties...>  cubPoints,
                   Kokkos::DynRankView<cubWeightValueType,cubWeightProperties...> cubWeights ) const;

    };
    Internal impl_;

    typedef typename Cubature<ExecSpaceType,pointValueType,weightValueType>::pointViewType  pointViewType;
    typedef typename Cubature<ExecSpaceType,pointValueType,weightValueType>::weightViewType weightViewType;

    virtual
    void
    getCubature( pointViewType  cubPoints,
                 weightViewType cubWeights ) const {
      impl_.getCubature( cubPoints,
                         cubWeights );
    }

    /** \brief Returns the number of cubature points.
     */
    virtual
    ordinal_type
    getNumPoints() const {
      ordinal_type numCubPoints = 1;
      for (auto i=0;i<numCubatures_;++i)
        numCubPoints *= cubatures_[i].getNumPoints();
      return numCubPoints;
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
      return "CubatureTensor";
    }

    virtual
    ordinal_type 
    getAccuracy() const {
      ordinal_type r_val = 0;
      for (auto i=0;i<numCubatures_;++i)
        r_val = Util::max(r_val, cubatures_[i].getAccuracy());
      return r_val;
    }

    /** \brief Return the number of cubatures.
     */
    ordinal_type getNumCubatures() const {
      return numCubatures_;
    }

    /** \brief Returns max. degree of polynomials that are integrated exactly.
     */
    void getAccuracy( ordinal_type *accuracy ) const {
      for (auto i=0;i<numCubatures_;++i)
        accuracy[i] = cubatures_[i].getAccuracy();
    }

    CubatureTensor() 
      : numCubatures_(0),
        dimension_(0),
        impl_(this) {}

    CubatureTensor(const CubatureTensor &b)
      : numCubatures_(b.numCubatures_),
        dimension_(b.dimension_),
        impl_(this) {
      for (auto i=0;i<numCubatures_;++i) 
        cubatures_[i] = b.cubatures_[i];
    }

    /** \brief Constructor.

        \param cubature1        [in]     - First direct cubature rule.
        \param cubature2        [in]     - Second direct cubature rule.
    */
    template<typename CubatureType0,
             typename CubatureType1>
    CubatureTensor( const CubatureType0 cubature0,
                    const CubatureType1 cubature1 )
      : numCubatures_(2),
        dimension_(cubature0.getDimension()+cubature1.getDimension()),
        impl_(this) {
      cubatures_[0] = cubature0;
      cubatures_[1] = cubature1;
    }

    /** \brief Constructor.

        \param cubature1        [in]     - First direct cubature rule.
        \param cubature2        [in]     - Second direct cubature rule.
        \param cubature3        [in]     - Third direct cubature rule.
    */
    template<typename CubatureType0,
             typename CubatureType1,
             typename CubatureType2>
    CubatureTensor( const CubatureType0 cubature0,
                    const CubatureType1 cubature1,
                    const CubatureType2 cubature2 ) 
      : numCubatures_(3),
        dimension_(cubature0.getDimension()+cubature1.getDimension()+cubature2.getDimension()),
        impl_(this) {
      cubatures_[0] = cubature0;
      cubatures_[1] = cubature1;
      cubatures_[2] = cubature2;
    }

    CubatureTensor& operator=(const CubatureTensor &b) {
      if (this != &b) {
        Cubature<ExecSpaceType,pointValueType,weightValueType>::operator= (b);
        numCubatures_ = b.numCubatures_;
        for (auto i=0;i<numCubatures_;++i)
          cubatures_[i] = b.cubatures_[i];
        dimension_ = b.dimension_;
        // do not copy impl
      }
      return *this;
    }

  };


} // end namespace Intrepid2


// include templated definitions
#include <Intrepid2_CubatureTensorDef.hpp>

#endif
