//@HEADER
// ************************************************************************
// 
//               Tpetra: Templated Linear Algebra Services Package 
//                 Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef TPETRA_CRSMATRIXSOLVEOP_HPP
#define TPETRA_CRSMATRIXSOLVEOP_HPP

#include <Teuchos_RCP.hpp>
#include "Tpetra_CrsMatrix.hpp"

namespace Tpetra {

  //! \brief A class for wrapping a Tpetra::CrsMatrix solve in a Tpetra::Operator.
  template <class OpScalar, class MatScalar = OpScalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatVec = Kokkos::DefaultSparseMultiply<MatScalar,LocalOrdinal,Node>, class LocalMatSolve = Kokkos::DefaultSparseSolve<MatScalar,LocalOrdinal,Node> >
  class CrsMatrixSolveOp : public Operator<OpScalar,LocalOrdinal,GlobalOrdinal,Node> {
    public:
      //! @name Constructor/Destructor Methods
      //@{ 

      //! Constructor
      CrsMatrixSolveOp(const Teuchos::RCP<const CrsMatrix<MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> > &A);

      //! Destructor
      virtual ~CrsMatrixSolveOp();

      //@}

      //! @name Methods implementing Operator
      //@{ 

      //! Computes this matrix-vector multilication y = A x.
      //! This calls solve() on the underlying CrsMatrix object.
      void apply(const MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
                 Teuchos::ETransp mode = Teuchos::NO_TRANS, OpScalar alpha = Teuchos::ScalarTraits<OpScalar>::one(), OpScalar beta = Teuchos::ScalarTraits<OpScalar>::zero()) const;

      //! Indicates whether this operator supports inverting the adjoint operator.
      //! This is true.
      bool hasTransposeApply() const;

      //! \brief Returns the Map associated with the domain of this operator.
      //! This is the range map of the underlying CrsMatrix.
      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const;

      //! Returns the Map associated with the domain of this operator.
      //! This is the domain map of the underlying CrsMatrix.
      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const;

      //@}
    
    protected:
      typedef MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> MV;

      // underlying CrsMatrix
      const Teuchos::RCP<const CrsMatrix<MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> > matrix_;

      // multivectors used for import/export dest/source in apply()
      mutable Teuchos::RCP<MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> > importMV_, exportMV_;

      // private methods for transpose or non-transpose
      void applyNonTranspose(const MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> &Y) const;
      void applyTranspose(const MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> &Y) const;
  };


  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  CrsMatrixSolveOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::CrsMatrixSolveOp(const Teuchos::RCP<const CrsMatrix<MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> > &A) 
  : matrix_(A) {
  }

  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  CrsMatrixSolveOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::~CrsMatrixSolveOp() {
  }


  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void 
  CrsMatrixSolveOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::apply(
              const MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> & X,   
                    MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> & Y,
                    Teuchos::ETransp mode, OpScalar alpha, OpScalar beta) const 
  {
    TEST_FOR_EXCEPTION(!matrix_->isFillComplete(), std::runtime_error, 
        Teuchos::typeName(*this) << "::apply(): underlying matrix is not fill-complete.");
    TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
        Teuchos::typeName(*this) << "::apply(X,Y): X and Y must have the same number of vectors.");
    TEST_FOR_EXCEPTION(matrix_->isLowerTriangular() == false && matrix_->isUpperTriangular() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::apply() requires either upper or lower triangular structure in underlying matrix.");
    TEST_FOR_EXCEPTION( alpha != Teuchos::ScalarTraits<OpScalar>::one() || beta != Teuchos::ScalarTraits<OpScalar>::zero(), std::runtime_error,
        Teuchos::typeName(*this) << "::apply(): non-trivial alpha,beta not supported at this time.");
    if (mode == Teuchos::NO_TRANS) {
      applyNonTranspose(X,Y);
    }
    else {
      applyTranspose(X,Y);
    }
  }

  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void 
  CrsMatrixSolveOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::applyNonTranspose(
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
      matrix_->template solve<OpScalar,OpScalar>(*X,*importMV_,Teuchos::NO_TRANS);
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
        matrix_->template solve<OpScalar,OpScalar>(*X,Y,Teuchos::NO_TRANS);
        Y_in = Y;
      }
      else {
        matrix_->template solve<OpScalar,OpScalar>(*X,Y_in,Teuchos::NO_TRANS);
      }
    }
  }


  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void 
  CrsMatrixSolveOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::applyTranspose(
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
      matrix_->template solve<OpScalar,OpScalar>(*X,*exportMV_,Teuchos::CONJ_TRANS);
      // Make sure target is zero: necessary because we are adding
      Y_in.putScalar(ST::zero());  
      Y_in.doExport(*importMV_, *importer, ADD);
    }
    // otherwise, solve into Y
    else {
      if (Y_in.isConstantStride() == false) {
        // generate a strided copy of Y
        MV Y(Y_in);
        matrix_->template solve<OpScalar,OpScalar>(*X,Y,Teuchos::CONJ_TRANS);
        Y_in = Y;
      }
      else {
        matrix_->template solve<OpScalar,OpScalar>(*X,Y_in,Teuchos::CONJ_TRANS);
      }
    }
  }

  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool 
  CrsMatrixSolveOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::hasTransposeApply() const {
    return true;
  }

  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsMatrixSolveOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getDomainMap() const {
    return matrix_->getRangeMap();
  }

  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsMatrixSolveOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getRangeMap() const {
    return matrix_->getDomainMap();
  }

  //! Non-member function to create Tpetra::CrsMatrixSolveOp
  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  Teuchos::RCP< CrsMatrixSolveOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> >
  createCrsMatrixSolveOp(const Teuchos::RCP<const CrsMatrix<MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> > &A) {
    return rcp(new CrsMatrixSolveOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>(A) );
  }

} // end of namespace Tpetra

#endif // TPETRA_CRSMATRIXSOLVEOP_HPP
