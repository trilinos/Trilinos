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

#ifndef INTREPID_CUBATURE_CONTROLVOLUMESIDEDEF_HPP
#define INTREPID_CUBATURE_CONTROLVOLUMESIDEDEF_HPP

namespace Intrepid{

template<class Scalar, class ArrayPoint, class ArrayWeight>
CubatureControlVolumeSide<Scalar,ArrayPoint,ArrayWeight>::CubatureControlVolumeSide(const Teuchos::RCP<const shards::CellTopology>& cellTopology)
{
  // topology of primary cell
  primaryCellTopo_ = cellTopology;

  // get topology of sub-control volume (will be quad or hex depending on dimension)
  const CellTopologyData &myCellData =
         (primaryCellTopo_->getDimension() > 2) ? *shards::getCellTopologyData<shards::Hexahedron<8> >() :
                                                  *shards::getCellTopologyData<shards::Quadrilateral<4> >();

  subCVCellTopo_ = Teuchos::rcp(new shards::CellTopology(&myCellData));

  degree_ = 1;

  numPoints_ = primaryCellTopo_->getEdgeCount();

  cubDimension_ = primaryCellTopo_->getDimension();

}

template<class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureControlVolumeSide<Scalar,ArrayPoint,ArrayWeight>::getCubature(ArrayPoint& cubPoints,
		                                                       ArrayWeight& cubWeights) const
{
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (CubatureControlVolumeSide): Cubature defined in physical space calling method for reference space cubature.");
}

template<class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureControlVolumeSide<Scalar,ArrayPoint,ArrayWeight>::getCubature(ArrayPoint& cubPoints,
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

 // num edges per primary cell
  int numEdgesPerCell = primaryCellTopo_->getEdgeCount();

  // Loop over cells
  for (int icell = 0; icell < numCells; icell++){

     // Get subcontrol volume side midpoints and normals
      int iside = 1;
      int numNodesPerSide = subCVCellTopo_->getNodeCount(spaceDim-1,iside);
      Intrepid::FieldContainer<int> sideNodes(numNodesPerSide);
      for (int i=0; i<numNodesPerSide; i++){
          sideNodes(i) = subCVCellTopo_->getNodeMap(spaceDim-1,iside,i);
      }

      // Loop over primary cell nodes and get side midpoints
      //   In each primary cell the number of control volume side integration
      //   points is equal to the number of primary cell edges. In 2d the
      //   number of edges = number of nodes and this loop defines all side
      //   points. In 3d this loop computes the side points for all
      //   subcontrol volume sides for iside = 1. Additional code below
      //   computes the remaining points for particular 3d topologies.
       for (int inode=0; inode < numNodesPerCell; inode++){
          for(int idim=0; idim < spaceDim; idim++){
             Scalar midpt = 0.0;
             for (int i=0; i<numNodesPerSide; i++){
                  midpt += subCVCoords(icell,inode,sideNodes(i),idim);
             }
             cubPoints(icell,inode,idim) = midpt/numNodesPerSide;
          }
       }

      // Map side center to reference subcell
       //Intrepid::FieldContainer<Scalar> sideCenterLocal(1,spaceDim-1);
       Intrepid::FieldContainer<double> sideCenterLocal(1,spaceDim-1);
       for (int idim = 0; idim < spaceDim-1; idim++){
          sideCenterLocal(0,idim) = 0.0;
       }

       Intrepid::FieldContainer<Scalar> refSidePoints(1,spaceDim);
       iside = 1;
       Intrepid::CellTools<Scalar>::mapToReferenceSubcell(refSidePoints,
                                    sideCenterLocal,
                                    spaceDim-1, iside, *(subCVCellTopo_));

      // Array of cell control volume coordinates
       Intrepid::FieldContainer<Scalar> cellCVCoords(numNodesPerCell, numNodesPerSubCV, spaceDim);
       for (int isubcv = 0; isubcv < numNodesPerCell; isubcv++) {
         for (int inode = 0; inode < numNodesPerSubCV; inode++){
           for (int idim = 0; idim < spaceDim; idim++){
               cellCVCoords(isubcv,inode,idim) = subCVCoords(icell,isubcv,inode,idim);
           }
         }
       }

      // calculate Jacobian at side centers
       Intrepid::FieldContainer<Scalar> subCVsideJacobian(numNodesPerCell, 1, spaceDim, spaceDim);
       Intrepid::CellTools<Scalar>::setJacobian(subCVsideJacobian, refSidePoints, cellCVCoords, *(subCVCellTopo_));

      // Get subcontrol volume side normals
       Intrepid::FieldContainer<Scalar> normals(numNodesPerCell, 1, spaceDim);
       Intrepid::CellTools<Scalar>::getPhysicalSideNormals(normals,subCVsideJacobian,iside,*(subCVCellTopo_));

       for (int inode = 0; inode < numNodesPerCell; inode++) {
          for (int idim = 0; idim < spaceDim; idim++){
             cubWeights(icell,inode,idim) = normals(inode,0,idim)*pow(2,spaceDim-1);
          }
       }

       if (primaryCellTopo_->getKey()==shards::Hexahedron<8>::key)
         {
           // first set of side midpoints and normals (above) associated with
           // primary cell edges 0-7 are obtained from side 1 of the
           // eight control volumes

           // second set of side midpoints and normals associated with
           // primary cell edges 8-11 are obtained from side 5 of the
           // first four control volumes.
           iside = 5;
           for (int i=0; i<numNodesPerSide; i++){
              sideNodes(i) = subCVCellTopo_->getNodeMap(spaceDim-1,iside,i);
           }
           int numExtraSides = numEdgesPerCell - numNodesPerCell;
             for (int icount=0; icount < numExtraSides; icount++){
                int iedge = icount + numNodesPerCell;
                for(int idim=0; idim < spaceDim; idim++){
                    Scalar midpt = 0.0;
                    for (int i=0; i<numNodesPerSide; i++){
                        midpt += subCVCoords(icell,icount,sideNodes(i),idim)/numNodesPerSide;
                    }
                    cubPoints(icell,iedge,idim) = midpt;
                }
            }

           // Map side center to reference subcell
           iside = 5;
           Intrepid::CellTools<Scalar>::mapToReferenceSubcell(refSidePoints,
                                        sideCenterLocal,
                                        spaceDim-1, iside, *(subCVCellTopo_));

           // calculate Jacobian at side centers
           Intrepid::CellTools<Scalar>::setJacobian(subCVsideJacobian, refSidePoints, cellCVCoords, *(subCVCellTopo_));

           // Get subcontrol volume side normals
           Intrepid::CellTools<Scalar>::getPhysicalSideNormals(normals,subCVsideJacobian,iside,*(subCVCellTopo_));

           for (int icount = 0; icount < numExtraSides; icount++) {
              int iedge = icount + numNodesPerCell;
              for (int idim = 0; idim < spaceDim; idim++){
                  cubWeights(icell,iedge,idim) = normals(icount,0,idim)*pow(2,spaceDim-1);
              }
           }

         } // end if Hex

        if (primaryCellTopo_->getKey()==shards::Tetrahedron<4>::key)
          {
           // first set of side midpoints and normals associated with
           // primary cell edges 0-2 are obtained from side 1 of the
           // eight control volumes (above)

           // second set of side midpoints and normals associated with
           // primary cell edges 3-5 are obtained from side 5 of the
           // first three control volumes.
           iside = 5;
           for (int i=0; i<numNodesPerSide; i++){
              sideNodes(i) = subCVCellTopo_->getNodeMap(spaceDim-1,iside,i);
           }
           for (int icount=0; icount < 3; icount++){
                int iedge = icount + 3;
                for(int idim=0; idim < spaceDim; idim++){
                    Scalar midpt = 0.0;
                    for (int i=0; i<numNodesPerSide; i++){
                        midpt += subCVCoords(icell,icount,sideNodes(i),idim)/numNodesPerSide;
                    }
                    cubPoints(icell,iedge,idim) = midpt;
                }
           }

          // Map side center to reference subcell
           iside = 5;
           Intrepid::CellTools<Scalar>::mapToReferenceSubcell(refSidePoints,
                                        sideCenterLocal,
                                        spaceDim-1, iside, *(subCVCellTopo_));

           // calculate Jacobian at side centers
           Intrepid::CellTools<Scalar>::setJacobian(subCVsideJacobian, refSidePoints, cellCVCoords, *(subCVCellTopo_));

           // Get subcontrol volume side normals
           Intrepid::CellTools<Scalar>::getPhysicalSideNormals(normals,subCVsideJacobian,iside,*(subCVCellTopo_));

           for (int icount = 0; icount < 3; icount++) {
              int iedge = icount + 3;
              for (int idim = 0; idim < spaceDim; idim++){
                  cubWeights(icell,iedge,idim) = normals(icount,0,idim)*pow(2,spaceDim-1);
              }
           }

       }// if tetrahedron

  } // end loop over cells

} // end getCubature
    

template<class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureControlVolumeSide<Scalar,ArrayPoint,ArrayWeight>::getNumPoints()const{
  return numPoints_;
} // end getNumPoints

template<class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureControlVolumeSide<Scalar,ArrayPoint,ArrayWeight>::getDimension()const{
  return cubDimension_;
} // end getNumPoints

template<class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureControlVolumeSide<Scalar,ArrayPoint,ArrayWeight>::getAccuracy(std::vector<int>& accuracy)const{
  accuracy.assign(1,degree_);
} // end getAccuracy

} // end namespace Intrepid

#endif


#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

