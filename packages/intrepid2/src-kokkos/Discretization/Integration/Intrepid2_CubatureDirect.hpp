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

/** \file   Intrepid_CubatureDirect.hpp
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
  template<typename ExecSpaceType>
  class CubatureDirect : public Cubature<ExecSpaceType> {
  protected:
    
    /** \brief The degree of polynomials that are integrated
        exactly by this cubature rule.
    */
    ordinal_type degree_;
    
    /** \brief Dimension of integration domain.
     */
    ordinal_type dimension_;


    /** \brief Returns cubature points and weights
        
        \param cubPoints       [out]     - Array containing the cubature points.
        \param cubWeights      [out]     - Array of corresponding cubature weights.
        \param cubData          [in]     - Cubuture data object
    */
    CubatureData
    getCubatureData(const ordinal_type degree) const {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::logic_error,
                                    ">>> ERROR (CubatureDirect::getCubatureData): this method should be over-riden by derived classes.");
    }

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
                         const CubatureData  cubData) const;
    
  public:

    CubatureDirect() = default;    
    ~CubatureDirect() = default; 

    //
    // Cubature public functions
    //

    /** \brief Returns cubature points and weights
        
        \param cubPoints       [out]     - Array containing the cubature points.
        \param cubWeights      [out]     - Array of corresponding cubature weights.
    */
    template<typename cubPointValueType,  class ...cubPointProperties,
             typename cubWeightValueType, class ...cubWeightProperties>
    void 
    getCubature( Kokkos::DynRankView<cubPointValueType, cubPointProperties...>  cubPoints,
                 Kokkos::DynRankView<cubWeightValueType,cubWeightProperties...> cubWeights ) const;
    
    
    /** \brief Returns cubature points and weights.
               Method for physical space cubature, throws an exception.
        
        \param cubPoints             [out]        - Array containing the cubature points.
        \param cubWeights            [out]        - Array of corresponding cubature weights.
        \param cellCoords             [in]        - Array of cell coordinates
    */
    template<typename cubPointValueType,  class ...cubPointProperties,
             typename cubWeightValueType, class ...cubWeightProperties,
             typename cellCoordValueType, class ...cellCoordProperties>
    void
    getCubature( Kokkos::DynRankView<cubPointValueType, cubPointProperties...>  cubPoints,
                 Kokkos::DynRankView<cubWeightValueType,cubWeightProperties...> cubWeights,
                 Kokkos::DynRankView<cellCoordValueType,cellCoordProperties...> cellCoords ) const;
    
    /** \brief Returns the number of cubature points.
     */
    ordinal_type getNumPoints() const;
    
    /** \brief Returns dimension of integration domain.
     */
    ordinal_type getDimension() const;

    /** \brief Returns cubature name.
     */
    const char* getName() const;

    //
    // CubatureDirect specific functions
    //

    /** \brief Returns max. degree of polynomials that are integrated exactly.
        The return vector has size 1.
    */
    ordinal_type getAccuracy() const;
    
    /** \brief Returns maximum cubature accuracy.
     */
    // never used 
    // ordinal_type getMaxAccuracy() const;
  };  
  
} // end namespace Intrepid2


// include templated definitions
#include <Intrepid2_CubatureDirectDef.hpp>

#endif
