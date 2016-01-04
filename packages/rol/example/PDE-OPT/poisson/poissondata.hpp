// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  poissondata.hpp
    \brief Generates and manages data for the Poisson example, including
           all mesh and discretization data, matrices, etc.
*/

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"
#include "../TOOLS/dofmanager.hpp"

template<class Real>
class PoissonData {

private:
  Teuchos::RCP<MeshManager<Real> > meshMgr_;
  Teuchos::RCP<DofManager<Real> >  dofMgr_;
  std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;

  shards::CellTopology cellType_;
  int spaceDim_;
  int numNodesPerCell_;
  int numCubPoints_;

  int numCells_;

  Teuchos::RCP<Intrepid::FieldContainer<Real> > cubPoints_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cubWeights_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cellNodes_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cellJac_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cellJacInv_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cellJacDet_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cellWeightedMeasure_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradReference_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradPhysical_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradPhysicalWeighted_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > stiffMats_;

public:

  PoissonData(Teuchos::ParameterList &parlist) {
    /*** Initialize mesh / finite element fields / degree-of-freedom manager. ***/
    // Mesh manager.
    meshMgr_ = Teuchos::rcp(new MeshManager_Rectangle<Real>(parlist));
    printMeshData();
    // Finite element fields.
    Teuchos::RCP<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> > > basisPtrQ1 =
      Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >);
    basisPtrs_.resize(1, Teuchos::null);
    basisPtrs_[0] = basisPtrQ1;
    // Degree-of-freedom manager.
    dofMgr_ = Teuchos::rcp(new DofManager<Real>(meshMgr_, basisPtrs_));

    // Retrieve some basic information.
    cellType_ = (basisPtrs_[0])->getBaseCellTopology();   // get the cell type from any basis
    spaceDim_ = cellType_.getDimension();                 // retrieve spatial dimension
    numNodesPerCell_ = cellType_.getNodeCount();          // retrieve number of nodes per cell
    numCells_ = meshMgr_->getNumCells();                  // retrieve total number of cells in the mesh

    // Cubature data.
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                                          // create cubature factory
    int cubDegree = 4;                                                                          // set cubature degree, e.g. 2
    Teuchos::RCP<Intrepid::Cubature<Real> > cellCub = cubFactory.create(cellType_, cubDegree);  // create default cubature
    numCubPoints_ = cellCub->getNumPoints();                                                    // retrieve number of cubature points

    int lfs = dofMgr_->getLocalFieldSize(0);

    // Discretization data. 
    cubPoints_            = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCubPoints_, spaceDim_));
    cubWeights_           = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCubPoints_));
    cellNodes_            = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, numNodesPerCell_, spaceDim_));
    cellJac_              = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, numCubPoints_, spaceDim_, spaceDim_));
    cellJacInv_           = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, numCubPoints_, spaceDim_, spaceDim_));
    cellJacDet_           = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, numCubPoints_));
    cellWeightedMeasure_  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, numCubPoints_));
    gradReference_        = Teuchos::rcp(new Intrepid::FieldContainer<Real>(lfs, numCubPoints_, spaceDim_));  
    gradPhysical_         = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, lfs, numCubPoints_, spaceDim_));
    gradPhysicalWeighted_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, lfs, numCubPoints_, spaceDim_));
    stiffMats_            = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, lfs, lfs));

    // Geometric definition of the cells in the mesh, based on the cell-to-node map.
    Intrepid::FieldContainer<Real> &nodes = *meshMgr_->getNodes();
    Intrepid::FieldContainer<int>  &ctn   = *meshMgr_->getCellToNodeMap();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numNodesPerCell_; ++j) {
        for (int k=0; k<spaceDim_; ++k) {
          (*cellNodes_)(i, j, k) = nodes(ctn(i,j), k);
        }
      }
    }

    cellCub->getCubature(*cubPoints_, *cubWeights_);                                         // retrieve cubature points and weights
    (*basisPtrs_[0]).getValues(*gradReference_, *cubPoints_, Intrepid::OPERATOR_GRAD);       // evaluate grad operator at cubature points

    Intrepid::CellTools<Real>::setJacobian(*cellJac_, *cubPoints_, *cellNodes_, cellType_);  // compute cell Jacobians
    Intrepid::CellTools<Real>::setJacobianInv(*cellJacInv_, *cellJac_);                      // compute inverses of cell Jacobians
    Intrepid::CellTools<Real>::setJacobianDet(*cellJacDet_, *cellJac_);                      // compute determinants of cell Jacobians

    Intrepid::FunctionSpaceTools::computeCellMeasure<Real>(*cellWeightedMeasure_,            // compute weighted cell measure
                                                           *cellJacDet_,
                                                           *cubWeights_);
    Intrepid::FunctionSpaceTools::HGRADtransformGRAD<Real>(*gradPhysical_,                   // transform reference gradients into physical space
                                                           *cellJacInv_,
                                                           *gradReference_);
    Intrepid::FunctionSpaceTools::multiplyMeasure<Real>(*gradPhysicalWeighted_,              // multiply with weighted measure
                                                        *cellWeightedMeasure_,
                                                        *gradPhysical_);
    Intrepid::FunctionSpaceTools::integrate<Real>(*stiffMats_,                               // compute local stiffness matrices
                                                  *gradPhysical_,
                                                  *gradPhysicalWeighted_,
                                                  Intrepid::COMP_CPP);

  }


  void printMeshData() const {
    Teuchos::RCP<Intrepid::FieldContainer<Real> > nodesPtr = meshMgr_->getNodes();
    Teuchos::RCP<Intrepid::FieldContainer<int> >  cellToNodeMapPtr = meshMgr_->getCellToNodeMap();
    Intrepid::FieldContainer<Real> &nodes = *nodesPtr;
    Intrepid::FieldContainer<int>   &cellToNodeMap = *cellToNodeMapPtr;
    std::cout << "Number of nodes = " << meshMgr_->getNumNodes() << std::endl;
    std::cout << "Number of cells = " << meshMgr_->getNumCells() << std::endl;
    std::cout << "Number of edges = " << meshMgr_->getNumEdges() << std::endl;
    // Print mesh to file.
    std::ofstream meshfile;
    meshfile.open("mesh.txt");
    for (int i=0; i<cellToNodeMap.dimension(0); ++i) {
      meshfile << nodes(cellToNodeMap(i,0), 0) << "  " << nodes(cellToNodeMap(i,0), 1) << std::endl;
      meshfile << nodes(cellToNodeMap(i,1), 0) << "  " << nodes(cellToNodeMap(i,1), 1) << std::endl;
      meshfile << nodes(cellToNodeMap(i,2), 0) << "  " << nodes(cellToNodeMap(i,2), 1) << std::endl;
      meshfile << nodes(cellToNodeMap(i,3), 0) << "  " << nodes(cellToNodeMap(i,3), 1) << std::endl;
      meshfile << nodes(cellToNodeMap(i,0), 0) << "  " << nodes(cellToNodeMap(i,0), 1) << std::endl;
      meshfile << nodes(cellToNodeMap(i,1), 0) << "  " << nodes(cellToNodeMap(i,1), 1) << std::endl;
      meshfile << nodes(cellToNodeMap(i,2), 0) << "  " << nodes(cellToNodeMap(i,2), 1) << std::endl;
    }
    meshfile.close();
  }

}; // class Poisson_Data
