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

/** \file   Intrepid_CubatureControlVolumeBoundary.hpp
    \brief  Header file for the Intrepid::CubatureControlVolumeBoundary class.
    \author Created by K. Peterson, P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_CUBATURE_CONTROLVOLUMEBOUNDARY_HPP
#define INTREPID_CUBATURE_CONTROLVOLUMEBOUNDARY_HPP

#include "Intrepid_Cubature.hpp"
#include "Teuchos_Assert.hpp"
#include "Shards_CellTopology.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"

namespace Intrepid{

  /** \class Intrepid::CubatureControlVolumeBoundary
      \brief Defines cubature (integration) rules over Neumann boundaries for control volume method.

      Integration on Neumann boundaries for the control volume method requires integration points
      defined on primary cell sides. These points are not equivalent to control volume points on lower
      dimensional topologies and therefore require a separate class to define them. 
  */
  template<class Scalar, class ArrayPoint, class ArrayWeight>
  class CubatureControlVolumeBoundary : public Intrepid::Cubature<Scalar,ArrayPoint,ArrayWeight>{
  public:
    
    /** brief Constructor.
	
	\param cellTopology           [in]     - The topology of the primary cell.
	\param cellSide               [in]     - The index of the boundary side of the primary cell 
    */
    CubatureControlVolumeBoundary(const Teuchos::RCP<const shards::CellTopology>& cellTopology, int cellSide=0);

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
    
     /** \brief Returns max. degree of polynomials that are integrated exactly.
             The return vector has size 1.
     */
    void getAccuracy(std::vector<int> & accuracy) const;

    
    virtual ~CubatureControlVolumeBoundary() {}
    
  private:
    
    
    /** \brief The topology of the primary cell side.
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

    /** \brief Index of cell side
     */
    int sideIndex_;
    
  }; // end class CubatureControlVolumeBoundary

} // end namespace Intrepid

#include "Intrepid_CubatureControlVolumeBoundaryDef.hpp"

#endif


#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

