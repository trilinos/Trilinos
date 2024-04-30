// @HEADER
// ************************************************************************
//
//                           Intrepid Package
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_CubatureControlVolumeSide.hpp
    \brief  Header file for the Intrepid::CubatureControlVolumeSide class.
    \author Created by K. Peterson, P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_CUBATURE_CONTROLVOLUMESIDE_HPP
#define INTREPID_CUBATURE_CONTROLVOLUMESIDE_HPP

#include "Intrepid_Cubature.hpp"
#include "Teuchos_Assert.hpp"
#include "Shards_CellTopology.hpp"
#include "Intrepid_CellTools.hpp"

namespace Intrepid{

  /** \class Intrepid::CubatureControlVolumeSide
      \brief Defines cubature (integration) rules over control volumes.
  */
  template<class Scalar, class ArrayPoint, class ArrayWeight>
  class CubatureControlVolumeSide : public Intrepid::Cubature<Scalar,ArrayPoint,ArrayWeight>{
  public:
    
    /** brief Constructor.
	
	\param cellTopology           [in]     - The topology of the primary cell.
    */
    CubatureControlVolumeSide(const Teuchos::RCP<const shards::CellTopology> & cellTopology);

    /** \brief Returns cubature points and weights
	       Method for reference space cubature - throws an exception.
	
	\param cubPoints             [out]        - Array containing the cubature points.
        \param cubWeights            [out]        - Array of corresponding cubature weights.
    */
    void getCubature(ArrayPoint& cubPoints,
		     ArrayWeight& cubWeights) const;
    
    /** \brief Returns cubature points and weights
	       (return arrays must be pre-sized/pre-allocated).
	
	\param cubPoints             [out]        - Array containing the cubature points.
        \param cubWeights            [out]        - Array of corresponding cubature weights.
        \param cellCoords             [in]        - Array of cell coordinates
    */
    void getCubature(ArrayPoint& cubPoints,
		     ArrayWeight& cubWeights,
                     ArrayPoint& cellCoords) const;
    
    /** \brief Returns the number of cubature points.
     */
    int getNumPoints() const;
    
    /** \brief Returns dimension of integration domain.
     */
    int getDimension() const;
    
     /** \brief Returns max. degree of polynomials that are integrated exactly on each triangle.
             The return vector has size 1.
     */
    void getAccuracy(std::vector<int> & accuracy) const;

    
    virtual ~CubatureControlVolumeSide() {}
    
  private:
    
    
    /** \brief The topology of the primary cell.
     */
    Teuchos::RCP<const shards::CellTopology> primaryCellTopo_;

    /** \brief The topology of the sub-control volume.
     */
    Teuchos::RCP<const shards::CellTopology> subCVCellTopo_;

    /** \brief The degree of the polynomials that are integrated exactly.
     */
    int degree_;
    
    /** \brief The number of cubature points.
     */
    int numPoints_;
    
    /** \brief Dimension of integration domain.
     */
    int cubDimension_;
    
  }; // end class CubatureControlVolumeSide

} // end namespace Intrepid

#include "Intrepid_CubatureControlVolumeSideDef.hpp"

#endif


#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

