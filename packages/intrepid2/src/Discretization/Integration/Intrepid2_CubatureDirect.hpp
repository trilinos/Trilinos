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

/** \file   Intrepid2_CubatureDirect.hpp
    \brief  Header file for the Intrepid2::CubatureDirect class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CUBATURE_DIRECT_HPP__
#define __INTREPID2_CUBATURE_DIRECT_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Cubature.hpp"

namespace Intrepid2 {

  /** \class Intrepid2::CubatureDirect
      \brief Defines direct cubature (integration) rules in Intrepid.

      Cubature template (rule) consists of cubature points and cubature weights.
      Intrepid provides a small collection of frequently used cubature rule templates
      for FEM reconstructions on simplices (edge, tri, tet) and pyramid cells.

      For quad, hex, and triprism rules, see tensor-product rules
      defined in the class CubatureTensor, and its derived classes.

      Cubature rules for simplices and the pyramid are stored in the
      <var>cubature_data_</var> array.

      All templates are defined on a reference cell and can be mapped to physical space
      cells by the methods available in the MultiCell class.
  */
  template<typename ExecSpaceType = void,
           typename pointValueType = double,
           typename weightValueType = double>
  class CubatureDirect
    : public Cubature<ExecSpaceType,pointValueType,weightValueType> {
  protected:

    /**
     \brief Cubature data is defined on the host space and is static
    */
    struct CubatureDataStatic {

      /** \brief  Number of cubature points stored in the template.
       */
      ordinal_type numPoints_;

      /** \brief  Array with the (X,Y,Z) coordinates of the cubature points.
       */
      pointValueType points_[Parameters::MaxIntegrationPoints][Parameters::MaxDimension];

      /** \brief  Array with the associated cubature weights.
       */
      weightValueType weights_[Parameters::MaxIntegrationPoints];
    };

    /**
     \brief Cubature data is defined on exec space and deep-copied when an object is created
    */
    struct CubatureData {

      /** \brief  Number of cubature points stored in the template.
       */
      ordinal_type numPoints_;

      /** \brief  Array with the (X,Y,Z) coordinates of the cubature points.
       */
      Kokkos::DynRankView<pointValueType,ExecSpaceType> points_;

      /** \brief  Array with the associated cubature weights.
       */
      Kokkos::DynRankView<weightValueType,ExecSpaceType> weights_;
    };

    /** \brief The degree of polynomials that are integrated
        exactly by this cubature rule.
    */
    ordinal_type degree_;

    /** \brief Dimension of integration domain.
     */
    ordinal_type dimension_;

    /** \brief Cubature data on device
     */
    CubatureData cubatureData_;

    /** \brief Returns cubature points and weights

        \param cubPoints       [out]     - Array containing the cubature points.
        \param cubWeights      [out]     - Array of corresponding cubature weights.
        \param cubData          [in]     - Cubuture data object
    */
    template<typename cubPointValueType,  class ...cubPointProperties,
             typename cubWeightValueType, class ...cubWeightProperties>
    void
    getCubatureFromData( Kokkos::DynRankView<cubPointValueType, cubPointProperties...>  cubPoints,
                         Kokkos::DynRankView<cubWeightValueType,cubWeightProperties...> cubWeights,
                         const CubatureData cubData) const {
#ifdef HAVE_INTREPID2_DEBUG
      // check size of cubPoints and cubWeights
      INTREPID2_TEST_FOR_EXCEPTION( cubPoints.rank() != 2, std::invalid_argument,
                                    ">>> ERROR (CubatureDirect): cubPoints must be rank 2." );

      INTREPID2_TEST_FOR_EXCEPTION( cubWeights.rank() != 1, std::invalid_argument,
                                    ">>> ERROR (CubatureDirect): cubPoints must be rank 1." );

      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(cubPoints.extent(0))  < this->getNumPoints() ||
                                    static_cast<ordinal_type>(cubPoints.extent(1))  < this->getDimension(), std::out_of_range,
                                    ">>> ERROR (CubatureDirect): Insufficient space allocated for cubature points.");

      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(cubWeights.extent(0)) < this->getNumPoints(), std::out_of_range,
                                    ">>> ERROR (CubatureDirect): Insufficient space allocated for cubature weights.");
#endif
      // need subview here
      typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

      range_type pointRange(0, this->getNumPoints());
      range_type dimRange  (0, this->getDimension());
      {
        const auto src = Kokkos::subdynrankview(cubData.points_, pointRange, dimRange);
              auto dst = Kokkos::subdynrankview(cubPoints,       pointRange, dimRange);

        Kokkos::deep_copy( dst, src );
      }
      {
        const auto src = Kokkos::subdynrankview(cubData.weights_, pointRange);
              auto dst = Kokkos::subdynrankview(cubWeights,       pointRange);

        Kokkos::deep_copy(dst ,src);
      }
    }

  public:

    //
    // Cubature public functions
    //
    typedef typename Cubature<ExecSpaceType,pointValueType,weightValueType>::pointViewType  pointViewType;
    typedef typename Cubature<ExecSpaceType,pointValueType,weightValueType>::weightViewType weightViewType;

    using Cubature<ExecSpaceType,pointValueType,weightValueType>::getCubature;

    virtual
    void
    getCubature( pointViewType  cubPoints,
                 weightViewType cubWeights ) const {
      this->getCubatureFromData(cubPoints, cubWeights, this->cubatureData_);
    }

    /** \brief Returns the number of cubature points.
     */
    virtual
    ordinal_type
    getNumPoints() const {
      return cubatureData_.numPoints_;
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
      return "CubatureDirect";
    }

    /** \brief Returns max. degree of polynomials that are integrated exactly.
        The return vector has size 1.
    */
    virtual
    ordinal_type 
    getAccuracy() const {
      return degree_;
    }

    CubatureDirect()
      : degree_(),
        dimension_(),
        cubatureData_() {}

    CubatureDirect(const CubatureDirect &b)
      : degree_(b.degree_),
        dimension_(b.dimension_),
        cubatureData_(b.cubatureData_) {}

    CubatureDirect(const ordinal_type degree,
                   const ordinal_type dimension)
      : degree_(degree),
        dimension_(dimension),
        cubatureData_() {}

  };

} // end namespace Intrepid2


#endif
