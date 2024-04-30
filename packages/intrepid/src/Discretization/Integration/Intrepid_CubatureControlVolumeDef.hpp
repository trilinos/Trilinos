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

/** \file   Intrepid_CubatureControlVolume.hpp
    \brief  Header file for the Intrepid::CubatureControlVolume class.
    \author Created by K. Peterson, P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_CUBATURE_CONTROLVOLUMEDEF_HPP
#define INTREPID_CUBATURE_CONTROLVOLUMEDEF_HPP

namespace Intrepid{

template<class Scalar, class ArrayPoint, class ArrayWeight>
CubatureControlVolume<Scalar,ArrayPoint,ArrayWeight>::CubatureControlVolume(const Teuchos::RCP<const shards::CellTopology> & cellTopology)
{
  // topology of primary cell
  primaryCellTopo_ = cellTopology;

  // get topology of sub-control volume (will be quad or hex depending on dimension)
  const CellTopologyData &myCellData =
         (primaryCellTopo_->getDimension() > 2) ? *shards::getCellTopologyData<shards::Hexahedron<8> >() :
                                                  *shards::getCellTopologyData<shards::Quadrilateral<4> >();
  subCVCellTopo_ = Teuchos::rcp(new shards::CellTopology(&myCellData));

  degree_ = 1;

  // one control volume cubature point per primary cell node
  numPoints_ = primaryCellTopo_->getNodeCount();

  cubDimension_ = primaryCellTopo_->getDimension();

}

template<class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureControlVolume<Scalar,ArrayPoint,ArrayWeight>::getCubature(ArrayPoint& cubPoints,
		                                                       ArrayWeight& cubWeights) const
{
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (CubatureControlVolume): Cubature defined in physical space calling method for reference space cubature.");
}

template<class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureControlVolume<Scalar,ArrayPoint,ArrayWeight>::getCubature(ArrayPoint& cubPoints,
		                                                       ArrayWeight& cubWeights,
                                                                       ArrayPoint& cellCoords) const
{
  // get array dimensions
  int numCells         = cellCoords.dimension(0);
  int numNodesPerCell  = cellCoords.dimension(1);
  int spaceDim         = cellCoords.dimension(2);
  int numNodesPerSubCV = subCVCellTopo_->getNodeCount();

  // get sub-control volume coordinates (one sub-control volume per node of primary cell)
  Intrepid::FieldContainer<Scalar> subCVCoords(numCells,numNodesPerCell,numNodesPerSubCV,spaceDim);
  Intrepid::CellTools<Scalar>::getSubCVCoords(subCVCoords,cellCoords,*(primaryCellTopo_));

  // Integration points and weights for calculating sub-control volumes
  Intrepid::DefaultCubatureFactory<double>  subCVCubFactory;
  int subcvCubDegree = 2;
  Teuchos::RCP<Intrepid::Cubature<double,Intrepid::FieldContainer<double>  > > subCVCubature;
  subCVCubature = subCVCubFactory.create(*(subCVCellTopo_), subcvCubDegree);

  int subcvCubDim       = subCVCubature -> getDimension();
  int numSubcvCubPoints = subCVCubature -> getNumPoints();

   // Get numerical integration points and weights
  Intrepid::FieldContainer<double> subcvCubPoints (numSubcvCubPoints, subcvCubDim);
  Intrepid::FieldContainer<double> subcvCubWeights(numSubcvCubPoints);

  subCVCubature -> getCubature(subcvCubPoints, subcvCubWeights);

  // Loop over cells
  for (int icell = 0; icell < numCells; icell++){

    // get sub-control volume centers (integration points)
     Intrepid::FieldContainer<Scalar> subCVCenter(numNodesPerCell,1,spaceDim);
     Intrepid::FieldContainer<Scalar> cellCVCoords(numNodesPerCell,numNodesPerSubCV,spaceDim);
     for (int isubcv = 0; isubcv < numNodesPerCell; isubcv++){
       for (int idim = 0; idim < spaceDim; idim++){
          for (int inode = 0; inode < numNodesPerSubCV; inode++){
              subCVCenter(isubcv,0,idim) += subCVCoords(icell,isubcv,inode,idim)/numNodesPerSubCV;
              cellCVCoords(isubcv,inode,idim) = subCVCoords(icell,isubcv,inode,idim);
          }
          cubPoints(icell,isubcv,idim) = subCVCenter(isubcv,0,idim);
        }
     }

   // calculate Jacobian and determinant for each subCV quadrature point
     Intrepid::FieldContainer<Scalar> subCVJacobian(numNodesPerCell, numSubcvCubPoints, spaceDim, spaceDim);
     Intrepid::FieldContainer<Scalar> subCVJacobDet(numNodesPerCell, numSubcvCubPoints);
     Intrepid::CellTools<Scalar>::setJacobian(subCVJacobian, subcvCubPoints, cellCVCoords, *(subCVCellTopo_));
     Intrepid::CellTools<Scalar>::setJacobianDet(subCVJacobDet, subCVJacobian );

    // fill array with sub control volumes (the sub control volume cell measure)
     for (int inode = 0; inode < numNodesPerCell; inode++){
         Scalar vol = 0;
         for (int ipt = 0; ipt < numSubcvCubPoints; ipt++){
            vol += subcvCubWeights(ipt)*subCVJacobDet(inode,ipt);
         }
         cubWeights(icell,inode) = vol;
     }

 } // end cell loop

} // end getCubature
    

template<class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureControlVolume<Scalar,ArrayPoint,ArrayWeight>::getNumPoints()const{
  return numPoints_;
} // end getNumPoints

template<class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureControlVolume<Scalar,ArrayPoint,ArrayWeight>::getDimension()const{
  return cubDimension_;
} // end getNumPoints

template<class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureControlVolume<Scalar,ArrayPoint,ArrayWeight>::getAccuracy(std::vector<int>& accuracy)const{
  accuracy.assign(1,degree_);
} // end getAccuracy

} // end namespace Intrepid

#endif


#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

