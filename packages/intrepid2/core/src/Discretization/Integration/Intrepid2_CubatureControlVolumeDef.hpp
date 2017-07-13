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

/** \file   Intrepid_CubatureControlVolume.hpp
    \brief  Header file for the Intrepid2::CubatureControlVolume class.
    \author Created by K. Peterson, P. Bochev and D. Ridzal.
*/

#ifndef INTREPID2_CUBATURE_CONTROLVOLUMEDEF_HPP
#define INTREPID2_CUBATURE_CONTROLVOLUMEDEF_HPP

namespace Intrepid2{

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
  index_type numCells         = static_cast<index_type>(cellCoords.dimension(0));
  index_type numNodesPerCell  = static_cast<index_type>(cellCoords.dimension(1));
  index_type spaceDim         = static_cast<index_type>(cellCoords.dimension(2));
  int numNodesPerSubCV = subCVCellTopo_->getNodeCount();

  // get sub-control volume coordinates (one sub-control volume per node of primary cell)
  Intrepid2::FieldContainer<Scalar> subCVCoords(numCells,numNodesPerCell,numNodesPerSubCV,spaceDim);
  Intrepid2::CellTools<Scalar>::getSubCVCoords(subCVCoords,cellCoords,*(primaryCellTopo_));

  // Integration points and weights for calculating sub-control volumes
  Intrepid2::DefaultCubatureFactory<double>  subCVCubFactory;
  int subcvCubDegree = 2;
  Teuchos::RCP<Intrepid2::Cubature<double,Intrepid2::FieldContainer<double>  > > subCVCubature;
  subCVCubature = subCVCubFactory.create(*(subCVCellTopo_), subcvCubDegree);

  int subcvCubDim       = subCVCubature -> getDimension();
  int numSubcvCubPoints = subCVCubature -> getNumPoints();

   // Get numerical integration points and weights
  Intrepid2::FieldContainer<double> subcvCubPoints (numSubcvCubPoints, subcvCubDim);
  Intrepid2::FieldContainer<double> subcvCubWeights(numSubcvCubPoints);

  subCVCubature -> getCubature(subcvCubPoints, subcvCubWeights);

  // Loop over cells
  for (index_type icell = 0; icell < numCells; icell++){

    // get sub-control volume centers (integration points)
     Intrepid2::FieldContainer<Scalar> subCVCenter(numNodesPerCell,1,spaceDim);
     Intrepid2::FieldContainer<Scalar> cellCVCoords(numNodesPerCell,numNodesPerSubCV,spaceDim);
     for (index_type isubcv = 0; isubcv < numNodesPerCell; isubcv++){
       for (index_type idim = 0; idim < spaceDim; idim++){
          for (int inode = 0; inode < numNodesPerSubCV; inode++){
              subCVCenter(isubcv,0,idim) += subCVCoords(icell,isubcv,inode,idim)/numNodesPerSubCV;
              cellCVCoords(isubcv,inode,idim) = subCVCoords(icell,isubcv,inode,idim);
          }
          cubPoints(icell,isubcv,idim) = subCVCenter(isubcv,0,idim);
        }
     }

   // calculate Jacobian and determinant for each subCV quadrature point
     Intrepid2::FieldContainer<Scalar> subCVJacobian(numNodesPerCell, numSubcvCubPoints, spaceDim, spaceDim);
     Intrepid2::FieldContainer<Scalar> subCVJacobDet(numNodesPerCell, numSubcvCubPoints);
     Intrepid2::CellTools<Scalar>::setJacobian(subCVJacobian, subcvCubPoints, cellCVCoords, *(subCVCellTopo_));
     Intrepid2::CellTools<Scalar>::setJacobianDet(subCVJacobDet, subCVJacobian );

    // fill array with sub control volumes (the sub control volume cell measure)
     for (index_type inode = 0; inode < numNodesPerCell; inode++){
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

} // end namespace Intrepid2

#endif

