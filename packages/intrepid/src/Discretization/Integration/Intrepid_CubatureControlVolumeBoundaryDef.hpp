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

/** \file   Intrepid_CubatureControlVolumeBoundaryDef.hpp
    \brief  Definition file for the Intrepid::CubatureControlVolumeBoundary class.
    \author Created by K. Peterson, P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_CUBATURE_CONTROLVOLUMEBOUNDARYDEF_HPP
#define INTREPID_CUBATURE_CONTROLVOLUMEBOUNDARYDEF_HPP

namespace Intrepid{

template<class Scalar, class ArrayPoint, class ArrayWeight>
CubatureControlVolumeBoundary<Scalar,ArrayPoint,ArrayWeight>::CubatureControlVolumeBoundary(const Teuchos::RCP<const shards::CellTopology> & cellTopology, int cellSide)
{
  // topology of primary cell 
  primaryCellTopo_ = cellTopology;

  // get topology of sub-control volume (will be quad or hex depending on dimension)
  const CellTopologyData &myCellData =
  (primaryCellTopo_->getDimension() > 2) ? *shards::getCellTopologyData<shards::Hexahedron<8> >() :
                                           *shards::getCellTopologyData<shards::Quadrilateral<4> >();

  subCVCellTopo_ = Teuchos::rcp(new shards::CellTopology(&myCellData));

  degree_ = 1;

  cubDimension_ = primaryCellTopo_->getDimension();

  sideIndex_ = cellSide;

  // one control volume boundary cubature point per subcell node (for now)
  numPoints_ = primaryCellTopo_->getNodeCount(primaryCellTopo_->getDimension()-1,cellSide);

}

template<class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureControlVolumeBoundary<Scalar,ArrayPoint,ArrayWeight>::getCubature(ArrayPoint& cubPoints,
		                                                               ArrayWeight& cubWeights) const
{
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (CubatureControlVolumeBoundary): Cubature defined in physical space calling method for reference space cubature.");
}

template<class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureControlVolumeBoundary<Scalar,ArrayPoint,ArrayWeight>::getCubature(ArrayPoint& cubPoints,
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

  // define subcontrol volume side index corresponding to primary cell side
  int numPrimarySideNodes = primaryCellTopo_->getNodeCount(spaceDim-1,sideIndex_);  
  int numCVSideNodes = subCVCellTopo_->getNodeCount(spaceDim-1,0);  
  int numPrimarySides = primaryCellTopo_->getSubcellCount(spaceDim-1);
  Intrepid::FieldContainer<int> CVSideonBoundary(numPrimarySides,numCVSideNodes);

  switch(primaryCellTopo_->getKey() ) {

      case shards::Triangle<3>::key:
      case shards::Quadrilateral<4>::key:

         for (int iside=0; iside<numPrimarySides; iside++) {
              CVSideonBoundary(iside,0) = 0; CVSideonBoundary(iside,1) = 3;
          }

      break;

      case shards::Hexahedron<8>::key:

         // sides 0-3
         for (int iside=0; iside<4; iside++) {
              CVSideonBoundary(iside,0) = 0; CVSideonBoundary(iside,1) = 3;
              CVSideonBoundary(iside,2) = 3; CVSideonBoundary(iside,3) = 0;
         }
         // side 4
         CVSideonBoundary(4,0) = 4; CVSideonBoundary(4,1) = 4;
         CVSideonBoundary(4,2) = 4; CVSideonBoundary(4,3) = 4;

         // side 5
         CVSideonBoundary(5,0) = 5; CVSideonBoundary(5,1) = 5;
         CVSideonBoundary(5,2) = 5; CVSideonBoundary(5,3) = 5;

      break;

      case shards::Tetrahedron<4>::key:

         CVSideonBoundary(0,0) = 0; CVSideonBoundary(0,1) = 3; CVSideonBoundary(0,2) = 0;
         CVSideonBoundary(1,0) = 0; CVSideonBoundary(1,1) = 3; CVSideonBoundary(1,2) = 3;
         CVSideonBoundary(2,0) = 3; CVSideonBoundary(2,1) = 4; CVSideonBoundary(2,2) = 0;
         CVSideonBoundary(3,0) = 4; CVSideonBoundary(3,1) = 4; CVSideonBoundary(3,2) = 4;

      break;

      default:
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                            ">>> ERROR (CubatureControlVolumeBoundary: invalid cell topology.");
    } // cell key

    Intrepid::FieldContainer<double> sideCenterLocal(1,spaceDim-1);
    for (int idim = 0; idim < spaceDim-1; idim++){
        sideCenterLocal(0,idim) = 0.0;
    }

    // get side cubature points
    for (int icell = 0; icell < numCells; icell++)
    {

      for (int inode=0; inode<numPrimarySideNodes; inode++){

         int cvind = primaryCellTopo_->getNodeMap(spaceDim-1,sideIndex_,inode);
         int cvside = CVSideonBoundary(sideIndex_,inode);

         Intrepid::FieldContainer<double> cubpoint(spaceDim);
         for (int idim=0; idim<spaceDim; idim++) { 
            for (int icvnode=0; icvnode<numCVSideNodes; icvnode++) { 
               int cvnode = subCVCellTopo_->getNodeMap(spaceDim-1,cvside,icvnode);
               cubpoint(idim) = cubpoint(idim) + subCVCoords(icell,cvind,cvnode,idim);       
            }
            cubPoints(icell,inode,idim) = cubpoint(idim)/numCVSideNodes; 
         }

         // map side center to reference subcell
          Intrepid::FieldContainer<Scalar> refSidePoints(1,spaceDim);
          Intrepid::CellTools<Scalar>::mapToReferenceSubcell(refSidePoints,
                                        sideCenterLocal,
                                        spaceDim-1, cvside, *(subCVCellTopo_));

         // array of sub-control volume coordinates
          Intrepid::FieldContainer<Scalar> cellCVCoords(1, numNodesPerSubCV, spaceDim);
          for (int icvnode = 0; icvnode < numNodesPerSubCV; icvnode++){
              for (int idim = 0; idim < spaceDim; idim++){
                   cellCVCoords(0,icvnode,idim) = subCVCoords(icell,cvind,icvnode,idim);
              }
          }

         // calculate Jacobian at side centers
          Intrepid::FieldContainer<Scalar> subCVsideJacobian(1, 1, spaceDim, spaceDim);
          Intrepid::FieldContainer<Scalar> subCVsideJacobianDet(1, 1);
          Intrepid::CellTools<Scalar>::setJacobian(subCVsideJacobian, refSidePoints, cellCVCoords, *(subCVCellTopo_));
          Intrepid::CellTools<Scalar>::setJacobianDet(subCVsideJacobianDet, subCVsideJacobian);

         // calculate Jacobian at side centers
          Intrepid::FieldContainer<Scalar> measure(1, 1);
          Intrepid::FieldContainer<Scalar> weights(1, 1);
          if (spaceDim == 3){
             weights(0,0) = 4.0;
             Intrepid::FunctionSpaceTools::computeFaceMeasure<Scalar>(measure,subCVsideJacobian,weights,cvside,*(subCVCellTopo_));
          }
          else if (spaceDim == 2){
             weights(0,0) = 2.0;
             Intrepid::FunctionSpaceTools::computeEdgeMeasure<Scalar>(measure,subCVsideJacobian,weights,cvside,*(subCVCellTopo_));
          }

          cubWeights(icell,inode) = measure(0,0);

       } // end loop over primary side nodes

 } // end cell loop

} // end getCubature
    

template<class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureControlVolumeBoundary<Scalar,ArrayPoint,ArrayWeight>::getNumPoints()const{
  return numPoints_;
} // end getNumPoints

template<class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureControlVolumeBoundary<Scalar,ArrayPoint,ArrayWeight>::getDimension()const{
  return cubDimension_;
} // end getDimension

template<class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureControlVolumeBoundary<Scalar,ArrayPoint,ArrayWeight>::getAccuracy(std::vector<int>& accuracy)const{
  accuracy.assign(1,degree_);
} // end getAccuracy

} // end namespace Intrepid

#endif


#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

