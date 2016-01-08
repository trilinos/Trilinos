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

/*! \file  data.hpp
    \brief Generates and manages data for the Poisson example, including
           all mesh and discretization data, matrices, etc.
*/

#ifndef ROL_PDEOPT_POISSONDATA_H
#define ROL_PDEOPT_POISSONDATA_H

#include "Teuchos_GlobalMPISession.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Version.hpp"
#include "MatrixMarket_Tpetra.hpp"

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

  Teuchos::RCP<const Teuchos::Comm<int> > commPtr_;
  int myRank_;
  int numProcs_;

  Teuchos::RCP<const Tpetra::Map<> >    myOverlapMap_;
  Teuchos::RCP<const Tpetra::Map<> >    myUniqueMap_;
  Teuchos::RCP<Tpetra::CrsGraph<> >     matGraph_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >    matA_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >    matB_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >    matM_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >    matR_;
  Teuchos::RCP<Tpetra::MultiVector<> >  vecUd_;

  shards::CellTopology cellType_;
  int spaceDim_;
  int numNodesPerCell_;
  int numCubPoints_;

  int totalNumCells_;
  int totalNumDofs_;
  int numCells_;

  Teuchos::RCP<Intrepid::FieldContainer<Real> > cubPoints_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cubWeights_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cellNodes_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cellJac_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cellJacInv_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cellJacDet_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cellWeightedMeasure_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > valReference_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradReference_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > valPhysical_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradPhysical_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > valPhysicalWeighted_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradPhysicalWeighted_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradgradMats_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > valvalMats_;

public:

  PoissonData(const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
              const Teuchos::RCP<Teuchos::ParameterList> &parlist,
              const Teuchos::RCP<std::ostream> &outStream) {

    /************************************/
    /*** Retrieve communication data. ***/
    /************************************/
    commPtr_  = comm;
    myRank_   = commPtr_->getRank();
    numProcs_ = commPtr_->getSize();
    *outStream << "Total number of processors: " << numProcs_ << std::endl;
    /************************************/
    /************************************/


    /****************************************************************************/
    /*** Initialize mesh / finite element fields / degree-of-freedom manager. ***/
    /****************************************************************************/

    // Mesh manager.
    meshMgr_ = Teuchos::rcp(new MeshManager_Rectangle<Real>(*parlist));
    printMeshData(*outStream);
    // Finite element fields.
    Teuchos::RCP<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> > > basisPtrQ1 =
      Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >);
    basisPtrs_.resize(1, Teuchos::null);
    basisPtrs_[0] = basisPtrQ1;
    // Degree-of-freedom manager.
    dofMgr_ = Teuchos::rcp(new DofManager<Real>(meshMgr_, basisPtrs_));
    // Retrieve total number of cells in the mesh.
    totalNumCells_ = meshMgr_->getNumCells();
    // Retrieve total number of degrees of freedom in the mesh.
    totalNumDofs_ = dofMgr_->getNumDofs();

    /****************************************************************************/
    /****************************************************************************/


    /****************************************************/
    /*** Build parallel communication infrastructure. ***/
    /****************************************************/

    // Partition the cells in the mesh.  We use a basic quasi-equinumerous partitioning,
    // where the remainder, if any, is assigned to the last processor.
    Teuchos::Array<int> myGlobIds_;
    Teuchos::Array<int> myCellIds_;
    Teuchos::Array<int> cellOffsets_(numProcs_, 0);
    int cellsPerProc = totalNumCells_ / numProcs_;
    numCells_ = cellsPerProc;
    int cellSplit = parlist->sublist("Geometry").get("Partition type", 1);
    switch(cellSplit) {
      case 0:
        if (myRank_ == 0) {  // remainder in the first
          numCells_ += totalNumCells_ % numProcs_;
        }
        for (int i=1; i<numProcs_; ++i) {
          cellOffsets_[i] = cellOffsets_[i-1] + cellsPerProc + (static_cast<int>(i==1))*(totalNumCells_ % numProcs_);
        }
        break;
      case 1:
        if (myRank_ == numProcs_-1) { // remainder in the last
          numCells_ += totalNumCells_ % numProcs_;
        }
        for (int i=1; i<numProcs_; ++i) {
          cellOffsets_[i] = cellOffsets_[i-1] + cellsPerProc;
        }
        break;
      case 2:
        if (myRank_ < (totalNumCells_%numProcs_)) { // spread remainder, starting from the first
          numCells_++;
        }
        for (int i=1; i<numProcs_; ++i) {
          cellOffsets_[i] = cellOffsets_[i-1] + cellsPerProc + (static_cast<int>(i-1<(totalNumCells_%numProcs_)));
        }
        break;
    }
    Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
    int numLocalDofs = cellDofs.dimension(1);
    *outStream << "Cell offsets across processors: " << cellOffsets_ << std::endl;
    for (int i=0; i<numCells_; ++i) {
      myCellIds_.push_back(cellOffsets_[myRank_]+i);
      for (int j=0; j<numLocalDofs; ++j) {
        myGlobIds_.push_back( cellDofs(cellOffsets_[myRank_]+i,j) );
      }
    }
    std::sort(myGlobIds_.begin(), myGlobIds_.end());
    myGlobIds_.erase( std::unique(myGlobIds_.begin(), myGlobIds_.end()), myGlobIds_.end() );

    // Build maps.
    myOverlapMap_ = Teuchos::rcp(new Tpetra::Map<>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                                                   myGlobIds_, 0, comm));
    //std::cout << std::endl << myOverlapMap_->getNodeElementList();
    /** One can also use the non-member function:
          myOverlapMap_ = Tpetra::createNonContigMap<int,int>(myGlobIds_, comm);
        to build the overlap map.
    **/
    myUniqueMap_ = Tpetra::createOneToOne<int,int>(myOverlapMap_);
    //std::cout << std::endl << myUniqueMap_->getNodeElementList() << std::endl;

    /****************************************************/
    /****************************************************/


    /****************************************************/
    /*** Set up local discretization data and arrays. ***/
    /****************************************************/

    // Retrieve some basic cell information.
    cellType_ = (basisPtrs_[0])->getBaseCellTopology();   // get the cell type from any basis
    spaceDim_ = cellType_.getDimension();                 // retrieve spatial dimension
    numNodesPerCell_ = cellType_.getNodeCount();          // retrieve number of nodes per cell

    // Cubature data.
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                                          // create cubature factory
    int cubDegree = 4;                                                                          // set cubature degree, e.g., 2
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
    valReference_         = Teuchos::rcp(new Intrepid::FieldContainer<Real>(lfs, numCubPoints_));  
    gradReference_        = Teuchos::rcp(new Intrepid::FieldContainer<Real>(lfs, numCubPoints_, spaceDim_));  
    valPhysical_          = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, lfs, numCubPoints_));
    gradPhysical_         = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, lfs, numCubPoints_, spaceDim_));
    valPhysicalWeighted_  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, lfs, numCubPoints_));
    gradPhysicalWeighted_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, lfs, numCubPoints_, spaceDim_));
    gradgradMats_         = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, lfs, lfs));
    valvalMats_           = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, lfs, lfs));

    // Geometric definition of the cells in the mesh, based on the cell-to-node map and the domain partition.
    Intrepid::FieldContainer<Real> &nodes = *meshMgr_->getNodes();
    Intrepid::FieldContainer<int>  &ctn   = *meshMgr_->getCellToNodeMap();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numNodesPerCell_; ++j) {
        for (int k=0; k<spaceDim_; ++k) {
          (*cellNodes_)(i, j, k) = nodes(ctn(myCellIds_[i],j), k);
        }
      }
    }

    /****************************************************/
    /****************************************************/


    /****************************************************************/
    /*** Assemble cellwise contributions to vectors and matrices. ***/
    /****************************************************************/

    cellCub->getCubature(*cubPoints_, *cubWeights_);                                         // retrieve cubature points and weights
    (*basisPtrs_[0]).getValues(*gradReference_, *cubPoints_, Intrepid::OPERATOR_GRAD);       // evaluate grad operator at cubature points
    (*basisPtrs_[0]).getValues(*valReference_, *cubPoints_, Intrepid::OPERATOR_VALUE);       // evaluate value operator at cubature points

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
    Intrepid::FunctionSpaceTools::integrate<Real>(*gradgradMats_,                            // compute local grad.grad (stiffness) matrices
                                                  *gradPhysical_,
                                                  *gradPhysicalWeighted_,
                                                  Intrepid::COMP_CPP);

    Intrepid::FunctionSpaceTools::HGRADtransformVALUE<Real>(*valPhysical_,                   // transform reference values into physical space
                                                            *valReference_);
    Intrepid::FunctionSpaceTools::multiplyMeasure<Real>(*valPhysicalWeighted_,               // multiply with weighted measure
                                                        *cellWeightedMeasure_,
                                                        *valPhysical_);
    Intrepid::FunctionSpaceTools::integrate<Real>(*valvalMats_,                              // compute local val.val (mass) matrices
                                                  *valPhysical_,
                                                  *valPhysicalWeighted_,
                                                  Intrepid::COMP_CPP);

    /****************************************************************/
    /****************************************************************/

    /****************************************/
    /*** Assemble global data structures. ***/
    /****************************************/

    // Assemble graph.
    matGraph_ = Teuchos::rcp(new Tpetra::CrsGraph<>(myUniqueMap_, 0));
    Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        matGraph_->insertGlobalIndices(cellDofs(myCellIds_[i],j), cellDofsArrayRCP(myCellIds_[i]*numLocalDofs, numLocalDofs));
      }
    }
    matGraph_->fillComplete();

    // Assemble matrices.
    matA_ = Tpetra::rcp(new Tpetra::CrsMatrix<>(matGraph_));
    int numLocalMatEntries = numLocalDofs*numLocalDofs;
    Teuchos::ArrayRCP<const double> gradgradArrayRCP = gradgradMats_->getData();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        matA_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                   cellDofsArrayRCP(myCellIds_[i]*numLocalDofs, numLocalDofs),
                                   gradgradArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
      }
    }
    matA_->fillComplete();

    matM_ = Tpetra::rcp(new Tpetra::CrsMatrix<>(matGraph_));
    Teuchos::ArrayRCP<const double> valvalArrayRCP = valvalMats_->getData();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        matM_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                   cellDofsArrayRCP(myCellIds_[i]*numLocalDofs, numLocalDofs),
                                   valvalArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
      }
    }
    matM_->fillComplete();

    matB_ = matM_;

    matR_ = matM_;

    // Assemble vectors.
    vecUd_ = Tpetra::rcp(new Tpetra::MultiVector<>(matA_->getRowMap(), 1, true));
    vecUd_->putScalar(1.0);

    /****************************************/
    /****************************************/

  }


  Teuchos::RCP<Tpetra::CrsMatrix<> > getMatA() const {
    return matA_;
  }


  Teuchos::RCP<Tpetra::CrsMatrix<> > getMatB() const {
    return matB_;
  }


  Teuchos::RCP<Tpetra::CrsMatrix<> > getMatM() const {
    return matM_;
  }


  Teuchos::RCP<Tpetra::CrsMatrix<> > getMatR() const {
    return matR_;
  }


  Teuchos::RCP<Tpetra::MultiVector<> > getVecUd() const {
    return vecUd_;
  }


  void printMeshData(std::ostream &outStream) const {
    Teuchos::RCP<Intrepid::FieldContainer<Real> > nodesPtr = meshMgr_->getNodes();
    Teuchos::RCP<Intrepid::FieldContainer<int> >  cellToNodeMapPtr = meshMgr_->getCellToNodeMap();
    Intrepid::FieldContainer<Real> &nodes = *nodesPtr;
    Intrepid::FieldContainer<int>   &cellToNodeMap = *cellToNodeMapPtr;
    outStream << "Number of nodes = " << meshMgr_->getNumNodes() << std::endl;
    outStream << "Number of cells = " << meshMgr_->getNumCells() << std::endl;
    outStream << "Number of edges = " << meshMgr_->getNumEdges() << std::endl;
    // Print mesh to file.
    if ((myRank_ == 0)) {
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
  }


  void outputTpetraData() const {
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> >   matWriter;
    Tpetra::MatrixMarket::Writer< Tpetra::MultiVector<> > vecWriter;
    matWriter.writeSparseFile("stiffness_mat", matA_);
    matWriter.writeSparseFile("mass_mat", matM_);
    vecWriter.writeDenseFile("vec", vecUd_);
  }


}; // class Poisson_Data

#endif
