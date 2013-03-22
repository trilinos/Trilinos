// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
// @HEADER

#ifndef TPETRA_CRSMATRIXSOLVEOP_DEF_HPP
#define TPETRA_CRSMATRIXSOLVEOP_DEF_HPP

#include "Tpetra_CrsMatrix.hpp"

#ifdef DOXYGEN_USE_ONLY
  #include "Tpetra_CrsMatrixSolveOp_decl.hpp"
#endif

/*! \file Tpetra_CrsMatrixSolveOp_def.hpp 

    The implementations for the members of Tpetra::CrsMatrixSolveOp and related non-member constructors.
 */

namespace Tpetra {

  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CrsMatrixSolveOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::CrsMatrixSolveOp(const Teuchos::RCP<const CrsMatrix<MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &A) 
  : matrix_(A) {
  }

  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CrsMatrixSolveOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::~CrsMatrixSolveOp() {
  }


  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void 
  CrsMatrixSolveOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::apply(
              const MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> & X,
                    MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> & Y,
                    Teuchos::ETransp mode, OpScalar alpha, OpScalar beta) const 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(!matrix_->isFillComplete(), std::runtime_error, 
        Teuchos::typeName(*this) << "::apply(): underlying matrix is not fill-complete.");
    TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
        Teuchos::typeName(*this) << "::apply(X,Y): X and Y must have the same number of vectors.");
    TEUCHOS_TEST_FOR_EXCEPTION(matrix_->isLowerTriangular() == false && matrix_->isUpperTriangular() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::apply() requires either upper or lower triangular structure in underlying matrix.");
    TEUCHOS_TEST_FOR_EXCEPTION( alpha != Teuchos::ScalarTraits<OpScalar>::one() || beta != Teuchos::ScalarTraits<OpScalar>::zero(), std::runtime_error,
        Teuchos::typeName(*this) << "::apply(): non-trivial alpha,beta not supported at this time.");
    if (mode == Teuchos::NO_TRANS) {
      applyNonTranspose(X,Y);
    }
    else {
      applyTranspose(X,Y);
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void 
  CrsMatrixSolveOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::applyNonTranspose(
      const MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> & X_in, 
            MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> & Y_in) const 
  {
    // Solve U X = Y  or  L X = Y
    // X belongs to domain map, while Y belongs to range map
    typedef Teuchos::ScalarTraits<OpScalar> ST;
    using Teuchos::null;

    const size_t numVectors = X_in.getNumVectors();
    Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > importer = matrix_->getGraph()->getImporter();
    Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > exporter = matrix_->getGraph()->getExporter();
    Teuchos::RCP<const MV> X;

    // it is okay if X and Y reference the same data, because we can perform a triangular solve in-situ.
    // however, we require that column access to each is strided.

    // set up import/export temporary multivectors
    if (importer != null) {
      if (importMV_ != null && importMV_->getNumVectors() != numVectors) importMV_ = null;
      if (importMV_ == null) {
        importMV_ = Teuchos::rcp( new MV(matrix_->getColMap(),numVectors) );
      }
    }
    if (exporter != null) {
      if (exportMV_ != null && exportMV_->getNumVectors() != numVectors) exportMV_ = null;
      if (exportMV_ == null) {
        exportMV_ = Teuchos::rcp( new MV(matrix_->getRowMap(),numVectors) );
      }
    }

    // solve(NO_TRANS): RangeMap -> DomainMap
    // lclMatSolve_: RowMap -> ColMap
    // importer: DomainMap -> ColMap
    // exporter: RowMap -> RangeMap
    // 
    // solve = reverse(exporter)  o   lclMatSolve_  o reverse(importer)
    //         RangeMap   ->    RowMap     ->     ColMap         ->    DomainMap
    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors
    if (exporter != null) {
      exportMV_->doImport(X_in, *exporter, INSERT);
      X = exportMV_;
    }
    else if (X_in.isConstantStride() == false) {
      // cannot handle non-constant stride right now
      // generate a copy of X_in
      X = Teuchos::rcp(new MV(X_in));
    }
    else {
      // just temporary, so this non-owning RCP is okay
      X = Teuchos::rcp( &X_in, false );
    }

    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    // We will compute solution into the to-be-exported MV
    if (importer != null) {
      matrix_->template localSolve<OpScalar,OpScalar>(*X,*importMV_,Teuchos::NO_TRANS);
      // Make sure target is zero: necessary because we are adding.
      Y_in.putScalar(ST::zero());  
      Y_in.doExport(*importMV_, *importer, ADD);
    }
    // otherwise, solve into Y
    else {
      // can't solve into non-strided multivector
      if (Y_in.isConstantStride() == false) {
        // generate a strided copy of Y
        MV Y(Y_in);
        matrix_->template localSolve<OpScalar,OpScalar>(*X,Y,Teuchos::NO_TRANS);
        Y_in = Y;
      }
      else {
        matrix_->template localSolve<OpScalar,OpScalar>(*X,Y_in,Teuchos::NO_TRANS);
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void 
  CrsMatrixSolveOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::applyTranspose(
               const MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> & X_in, 
                     MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> &Y_in) const 
  {
    typedef Teuchos::ScalarTraits<OpScalar> ST;
    using Teuchos::null;

    const size_t numVectors = X_in.getNumVectors();
    Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > importer = matrix_->getGraph()->getImporter();
    Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > exporter = matrix_->getGraph()->getExporter();
    Teuchos::RCP<const MV> X;

    // it is okay if X and Y reference the same data, because we can perform a triangular solve in-situ.
    // however, we require that column access to each is strided.

    // set up import/export temporary multivectors
    if (importer != null) {
      if (importMV_ != null && importMV_->getNumVectors() != numVectors) importMV_ = null;
      if (importMV_ == null) {
        importMV_ = Teuchos::rcp( new MV(matrix_->getColMap(),numVectors) );
      }
    }
    if (exporter != null) {
      if (exportMV_ != null && exportMV_->getNumVectors() != numVectors) exportMV_ = null;
      if (exportMV_ == null) {
        exportMV_ = Teuchos::rcp( new MV(matrix_->getRowMap(),numVectors) );
      }
    }

    // solve(TRANS): DomainMap -> RangeMap
    // lclMatSolve_(TRANS): ColMap -> RowMap
    // importer: DomainMap -> ColMap
    // exporter: RowMap -> RangeMap
    // 
    // solve = importer o   lclMatSolve_  o  exporter
    //         Domainmap -> ColMap     ->      RowMap -> RangeMap
    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if (importer != null) {
      importMV_->doImport(X_in,*importer,INSERT);
      X = importMV_;
    }
    else if (X_in.isConstantStride() == false) {
      // cannot handle non-constant stride right now
      // generate a copy of X_in
      X = Teuchos::rcp(new MV(X_in));
    }
    else {
      // just temporary, so this non-owning RCP is okay
      X = Teuchos::rcp( &X_in, false );
    }


    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    // We will compute solution into the to-be-exported MV; get a view
    if (exporter != null) {
      matrix_->template localSolve<OpScalar,OpScalar>(*X,*exportMV_,Teuchos::CONJ_TRANS);
      // Make sure target is zero: necessary because we are adding
      Y_in.putScalar(ST::zero());  
      Y_in.doExport(*importMV_, *importer, ADD);
    }
    // otherwise, solve into Y
    else {
      if (Y_in.isConstantStride() == false) {
        // generate a strided copy of Y
        MV Y(Y_in);
        matrix_->template localSolve<OpScalar,OpScalar>(*X,Y,Teuchos::CONJ_TRANS);
        Y_in = Y;
      }
      else {
        matrix_->template localSolve<OpScalar,OpScalar>(*X,Y_in,Teuchos::CONJ_TRANS);
      }
    }
  }

  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool 
  CrsMatrixSolveOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::hasTransposeApply() const {
    return true;
  }

  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsMatrixSolveOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getDomainMap() const {
    return matrix_->getRangeMap();
  }

  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsMatrixSolveOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getRangeMap() const {
    return matrix_->getDomainMap();
  }

} // end of namespace Tpetra

template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
Teuchos::RCP< Tpetra::CrsMatrixSolveOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> >
Tpetra::createCrsMatrixSolveOp(const Teuchos::RCP<const Tpetra::CrsMatrix<MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &A) {
  return Teuchos::rcp(new Tpetra::CrsMatrixSolveOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>(A) );
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

//! Explicit instantiation macro supporting the CrsMatrixSolveOp class. Instantiates the class, the non-member constructor, and the necessary CrsMatrix::solve() member.
#define TPETRA_CRSMATRIX_SOLVEOP_INSTANT(OPSCALAR,MATSCALAR,LO,GO,NODE) \
  \
  template class CrsMatrixSolveOp< OPSCALAR , MATSCALAR , LO , GO , NODE >; \
  \
  template Teuchos::RCP< Tpetra::CrsMatrixSolveOp<OPSCALAR,MATSCALAR,LO,GO,NODE> > \
  createCrsMatrixSolveOp(const Teuchos::RCP<const Tpetra::CrsMatrix<MATSCALAR,LO,GO,NODE> > &A); \
  \
  template void CrsMatrix<MATSCALAR,LO,GO,NODE>::localSolve<OPSCALAR,OPSCALAR>( \
        const MultiVector<OPSCALAR,LO,GO,NODE> &X, \
              MultiVector<OPSCALAR,LO,GO,NODE> &Y, \
              Teuchos::ETransp mode) const;

#define TPETRA_CRSMATRIX_SOLVEOP_INSTANT_SINGLE(SCALAR,LO,GO,NODE) \
        TPETRA_CRSMATRIX_SOLVEOP_INSTANT(SCALAR,SCALAR,LO,GO,NODE)

#endif // TPETRA_CRSMATRIXSOLVEOP_DEF_HPP
