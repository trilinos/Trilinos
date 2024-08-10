// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  sysbuilder.hpp
    \brief Builds block systems from individual matrices.
*/

#ifndef ROL_PDEOPT_SYSBUILDER_H
#define ROL_PDEOPT_SYSBUILDER_H

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Version.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include "MatrixMarket_Tpetra.hpp"

//// Global Timers.
#ifdef ROL_TIMERS
namespace ROL {
  namespace PDEOPT {

    ROL::Ptr<Teuchos::Time> AssembleBlockSystem = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble Block System");

  }
}
#endif


template<class Real>
class SysBuilder {

  typedef Tpetra::Map<>::local_ordinal_type   LO;
  typedef Tpetra::Map<>::global_ordinal_type  GO;
  typedef Tpetra::Map<>::node_type            NO;
  typedef Tpetra::MultiVector<Real,LO,GO,NO>  MV;
  typedef Teuchos::ArrayView<GO>::size_type   SZ;
  typedef Tpetra::global_size_t              GSZ;

private:



public:

  // Destructor.
  virtual ~SysBuilder() {}

  // Constructor.
  SysBuilder() {}   

  // Builds a KKT system matrix from individual matrices. 
  void buildMatrix(ROL::Ptr<Tpetra::CrsMatrix<> > &sysMat,
                   const ROL::Ptr<const Tpetra::CrsMatrix<> > &A,
                   const ROL::Ptr<const Tpetra::CrsMatrix<> > &B,
                   const ROL::Ptr<const Tpetra::CrsMatrix<> > &H11,
                   const ROL::Ptr<const Tpetra::CrsMatrix<> > &H12,
                   const ROL::Ptr<const Tpetra::CrsMatrix<> > &H21,
                   const ROL::Ptr<const Tpetra::CrsMatrix<> > &H22) {
    // Get individual maps.
    ROL::Ptr<const Tpetra::Map<> > mapSim = A->getDomainMap();
    ROL::Ptr<const Tpetra::Map<> > mapOpt = B->getDomainMap();
    ROL::Ptr<const Tpetra::Map<> > mapRes = A->getRangeMap();

    // Get individual indices.
    Teuchos::ArrayView<const GO> idxSim = mapSim->getLocalElementList();
    Teuchos::ArrayView<const GO> idxOpt = mapOpt->getLocalElementList();
    Teuchos::ArrayView<const GO> idxRes = mapRes->getLocalElementList();

    // Get communicator.
    ROL::Ptr<const Teuchos::Comm<int> > comm = mapSim->getComm();

    // Get an estimate of the number of entries per row.
    size_t totalNumEntries = A->getGlobalMaxNumRowEntries() + B->getGlobalMaxNumRowEntries();

    // Glue indices.
    Teuchos::Array<GO> idxSys;
    SZ offsetSim = 0, offsetOpt = idxSim.size(), offsetRes = idxSim.size() + idxOpt.size(); 
    for (SZ i=0; i<idxSim.size(); ++i) {
      idxSys.push_back(idxSim[i]+offsetSim);
    }
    for (SZ i=0; i<idxOpt.size(); ++i) {
      idxSys.push_back(idxOpt[i]+offsetOpt);
    }
    for (SZ i=0; i<idxRes.size(); ++i) {
      idxSys.push_back(idxRes[i]+offsetRes);
    }

    // Create block-system map.
    GSZ numSysIds  = idxSim.size() + idxOpt.size() + idxRes.size();
    ROL::Ptr<Tpetra::Map<> > sysMap = ROL::makePtr<Tpetra::Map<>>(numSysIds, idxSys, 0, comm);

    // Assemble global graph structure.
    ROL::Ptr<Tpetra::CrsGraph<> > sysGraph = ROL::makePtr<Tpetra::CrsGraph<>>(sysMap, totalNumEntries);

    // Create transposes of matrices A and B.
    Tpetra::RowMatrixTransposer<> transposerA(A);
    Tpetra::RowMatrixTransposer<> transposerB(B);
    ROL::Ptr<Tpetra::CrsMatrix<> > A_trans = transposerA.createTranspose();
    ROL::Ptr<Tpetra::CrsMatrix<> > B_trans = transposerB.createTranspose();

    // Insert indices of matrix H12 (Hessian_11), matrix H12 (Hessian_12) and matrix A (Jacobian_1) transposed.
    LO nRows = idxSim.size();
    size_t numEntries(0);
    typename Tpetra::CrsMatrix<>::nonconst_global_inds_host_view_type indices;
    typename Tpetra::CrsMatrix<>::nonconst_values_host_view_type values;
    for (LO i=0; i<nRows; ++i) {
      //
      if (H11 != ROL::nullPtr) {
        numEntries = H11->getNumEntriesInGlobalRow(idxSim[i]);
        Kokkos::resize(indices,numEntries);
        Kokkos::resize(values,numEntries);
        H11->getGlobalRowCopy(idxSim[i], indices, values, numEntries);
        sysGraph->insertGlobalIndices(i+offsetSim, indices.size(), indices.data());
      }
      //
      if (H12 != ROL::nullPtr) {
        numEntries = H12->getNumEntriesInGlobalRow(idxSim[i]);
        Kokkos::resize(indices,numEntries);
        Kokkos::resize(values,numEntries);
        H12->getGlobalRowCopy(idxSim[i], indices, values, numEntries);
        for (size_t j=0; j<numEntries; ++j) {
          indices(j) += offsetOpt;
        }
        sysGraph->insertGlobalIndices(i+offsetSim, indices.size(), indices.data());
      }
      //
      numEntries = A_trans->getNumEntriesInGlobalRow(idxSim[i]);
      Kokkos::resize(indices,numEntries);
      Kokkos::resize(values,numEntries);
      A_trans->getGlobalRowCopy(idxSim[i], indices, values, numEntries);
      for (size_t j=0; j<numEntries; ++j) {
        indices(j) += offsetRes;
      }
      sysGraph->insertGlobalIndices(i+offsetSim, indices.size(), indices.data());
    }
    // Insert indices of matrix H21 (Hessian_21), matrix H22 (Hessian_22) and matrix B (Jacobian_2) transposed.
    nRows = idxOpt.size();
    for (LO i=0; i<nRows; ++i) {
      //
      if (H21 != ROL::nullPtr) {
        numEntries = H21->getNumEntriesInGlobalRow(idxOpt[i]);
        Kokkos::resize(indices,numEntries);
        Kokkos::resize(values,numEntries);;
        H21->getGlobalRowCopy(idxOpt[i], indices, values, numEntries);
        sysGraph->insertGlobalIndices(i+offsetOpt, indices.size(), indices.data());
      }
      //
      if (H22 != ROL::nullPtr) {
        numEntries = H22->getNumEntriesInGlobalRow(idxOpt[i]);
        Kokkos::resize(indices,numEntries);
        Kokkos::resize(values,numEntries);
        H22->getGlobalRowCopy(idxOpt[i], indices, values, numEntries);
        for (size_t j=0; j<numEntries; ++j) {
          indices(j) += offsetOpt;
        }
        sysGraph->insertGlobalIndices(i+offsetOpt, indices.size(), indices.data());
      }
      //
      numEntries = B_trans->getNumEntriesInGlobalRow(idxOpt[i]);
      Kokkos::resize(indices,numEntries);
      Kokkos::resize(values,numEntries);
      B_trans->getGlobalRowCopy(idxOpt[i], indices, values, numEntries);
      for (size_t j=0; j<numEntries; ++j) {
        indices(j) += offsetRes;
      }
      sysGraph->insertGlobalIndices(i+offsetOpt, indices.size(), indices.data());
    }
    // Insert indices of matrix A (Jacobian_1) and matrix B (Jacobian_2).
    nRows = idxRes.size();
    for (LO i=0; i<nRows; ++i) {
      //
      numEntries = A->getNumEntriesInGlobalRow(idxRes[i]);
      Kokkos::resize(indices,numEntries);
      Kokkos::resize(values,numEntries);
      A->getGlobalRowCopy(idxRes[i], indices, values, numEntries);
      sysGraph->insertGlobalIndices(i+offsetRes, indices.size(), indices.data());
      //
      numEntries = B->getNumEntriesInGlobalRow(idxRes[i]);
      Kokkos::resize(indices,numEntries);
      Kokkos::resize(values,numEntries);
      B->getGlobalRowCopy(idxRes[i], indices, values, numEntries);
      for (size_t j=0; j<numEntries; ++j) {
        indices(j) += offsetOpt;
      }
      sysGraph->insertGlobalIndices(i+offsetRes, indices.size(), indices.data());
    }

    // Fill-complete the graph.
    sysGraph->fillComplete();

    // Create block-system matrix and fill with values.
    sysMat = ROL::makePtr<Tpetra::CrsMatrix<>>(sysGraph);
    sysMat->resumeFill(); sysMat->setAllToScalar(static_cast<Real>(0));

    // Insert values of matrix H12 (Hessian_11), matrix H12 (Hessian_12) and matrix A (Jacobian_1) transposed.
    nRows = idxSim.size();
    for (LO i=0; i<nRows; ++i) {
      //
      if (H11 != ROL::nullPtr) {
        numEntries = H11->getNumEntriesInGlobalRow(idxSim[i]);
        Kokkos::resize(indices,numEntries);
        Kokkos::resize(values,numEntries);
        H11->getGlobalRowCopy(idxSim[i], indices, values, numEntries);
        sysMat->replaceGlobalValues(i+offsetSim, indices, values);
      }
      //
      if (H12 != ROL::nullPtr) {
        numEntries = H12->getNumEntriesInGlobalRow(idxSim[i]);
        Kokkos::resize(indices,numEntries);
        Kokkos::resize(values,numEntries);
        H12->getGlobalRowCopy(idxSim[i], indices, values, numEntries);
        for (size_t j=0; j<numEntries; ++j) {
          indices(j) += offsetOpt;
        }
        sysMat->replaceGlobalValues(i+offsetSim, indices, values);
      }
      //
      numEntries = A_trans->getNumEntriesInGlobalRow(idxSim[i]);
      Kokkos::resize(indices,numEntries);
      Kokkos::resize(values,numEntries);
      A_trans->getGlobalRowCopy(idxSim[i], indices, values, numEntries);
      for (size_t j=0; j<numEntries; ++j) {
        indices(j) += offsetRes;
      }
      sysMat->replaceGlobalValues(i+offsetSim, indices, values);
    }
    // Insert values of matrix H21 (Hessian_21), matrix H22 (Hessian_22) and matrix B (Jacobian_2) transposed.
    nRows = idxOpt.size();
    for (LO i=0; i<nRows; ++i) {
      //
      if (H21 != ROL::nullPtr) {
        numEntries = H21->getNumEntriesInGlobalRow(idxOpt[i]);
        Kokkos::resize(indices,numEntries);
        Kokkos::resize(values,numEntries);
        H21->getGlobalRowCopy(idxOpt[i], indices, values, numEntries);
        sysMat->replaceGlobalValues(i+offsetOpt, indices, values);
      }
      //
      if (H22 != ROL::nullPtr) {
        numEntries = H22->getNumEntriesInGlobalRow(idxOpt[i]);
        Kokkos::resize(indices,numEntries);
        Kokkos::resize(values,numEntries);
        H22->getGlobalRowCopy(idxOpt[i], indices, values, numEntries);
        for (size_t j=0; j<numEntries; ++j) {
          indices(j) += offsetOpt;
        }
        sysMat->replaceGlobalValues(i+offsetOpt, indices, values);
      }
      //
      numEntries = B_trans->getNumEntriesInGlobalRow(idxOpt[i]);
      Kokkos::resize(indices,numEntries);
      Kokkos::resize(values,numEntries);
      B_trans->getGlobalRowCopy(idxOpt[i], indices, values, numEntries);
      for (size_t j=0; j<numEntries; ++j) {
        indices(j) += offsetRes;
      }
      sysMat->replaceGlobalValues(i+offsetOpt, indices, values);
    }
    // Insert values of matrix A (Jacobian_1) and matrix B (Jacobian_2).
    nRows = idxRes.size();
    for (LO i=0; i<nRows; ++i) {
      //
      numEntries = A->getNumEntriesInGlobalRow(idxRes[i]);
      Kokkos::resize(indices,numEntries);
      Kokkos::resize(values,numEntries);
      A->getGlobalRowCopy(idxRes[i], indices, values, numEntries);
      sysMat->replaceGlobalValues(i+offsetRes, indices, values);
      //
      numEntries = B->getNumEntriesInGlobalRow(idxRes[i]);
      Kokkos::resize(indices,numEntries);
      Kokkos::resize(values,numEntries);
      B->getGlobalRowCopy(idxRes[i], indices, values, numEntries);
      for (size_t j=0; j<numEntries; ++j) {
        indices(j) += offsetOpt;
      }
      sysMat->replaceGlobalValues(i+offsetRes, indices, values);
    }

    // Fill-complete the matrix.
    sysMat->fillComplete();

  }

  // Builds a block-system vector from individual vectors. 
  void buildVector(ROL::Ptr<MV> &sysVec,
                   const ROL::Ptr<const MV > &vecSim,
                   const ROL::Ptr<const MV > &vecOpt,
                   const ROL::Ptr<const MV > &vecRes) {
    // Get individual maps.
    ROL::Ptr<const Tpetra::Map<> > mapSim = vecSim->getMap();
    ROL::Ptr<const Tpetra::Map<> > mapOpt = vecOpt->getMap();
    ROL::Ptr<const Tpetra::Map<> > mapRes = vecRes->getMap();

    // Get individual indices.
    Teuchos::ArrayView<const GO> idxSim = mapSim->getLocalElementList();
    Teuchos::ArrayView<const GO> idxOpt = mapOpt->getLocalElementList();
    Teuchos::ArrayView<const GO> idxRes = mapRes->getLocalElementList();

    // Get communicator.
    ROL::Ptr<const Teuchos::Comm<int> > comm = mapSim->getComm();
    
    // Glue indices.
    Teuchos::Array<GO> idxSys;
    SZ offsetSim = 0, offsetOpt = idxSim.size(), offsetRes = idxSim.size() + idxOpt.size(); 
    for (SZ i=0; i<idxSim.size(); ++i) {
      idxSys.push_back(idxSim[i]+offsetSim);
    }
    for (SZ i=0; i<idxOpt.size(); ++i) {
      idxSys.push_back(idxOpt[i]+offsetOpt);
    }
    for (SZ i=0; i<idxRes.size(); ++i) {
      idxSys.push_back(idxRes[i]+offsetRes);
    }

    // Create block-system vector.
    size_t numVecs = vecSim->getNumVectors();
    GSZ numSysIds  = idxSim.size() + idxOpt.size() + idxRes.size();
    ROL::Ptr<Tpetra::Map<> > sysMap = ROL::makePtr<Tpetra::Map<>>(numSysIds, idxSys, 0, comm);
    sysVec = ROL::makePtr<MV>(sysMap, numVecs, true);

    // Fill block-system vector.
    for (size_t colId=0; colId<numVecs; ++colId) {
      Teuchos::ArrayRCP<const Real> dataSim = vecSim->getData(colId);
      for (size_t j=0; j<vecSim->getLocalLength(); ++j) {
        sysVec->replaceGlobalValue(idxSim[j]+offsetSim, colId, dataSim[j]);
      }
    }    
    for (size_t colId=0; colId<numVecs; ++colId) {
      Teuchos::ArrayRCP<const Real> dataOpt = vecOpt->getData(colId);
      for (size_t j=0; j<vecOpt->getLocalLength(); ++j) {
        sysVec->replaceGlobalValue(idxOpt[j]+offsetOpt, colId, dataOpt[j]);
      }
    }    
    for (size_t colId=0; colId<numVecs; ++colId) {
      Teuchos::ArrayRCP<const Real> dataRes = vecRes->getData(colId);
      for (size_t j=0; j<vecRes->getLocalLength(); ++j) {
        sysVec->replaceGlobalValue(idxRes[j]+offsetRes, colId, dataRes[j]);
      }
    }    
    
  }

}; // class SysBuilder

#endif

