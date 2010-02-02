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

#ifndef TPETRA_CRSMATRIXMULTIPLYOP_HPP
#define TPETRA_CRSMATRIXMULTIPLYOP_HPP

#include <Teuchos_RCP.hpp>
#include "Tpetra_CrsMatrix.hpp"

namespace Tpetra {

  //! \brief A class for wrapping a CrsMatrix multiply in a Operator.
  template <class OpScalar, class MatScalar = OpScalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatVec = Kokkos::DefaultSparseMultiply<MatScalar,LocalOrdinal,Node>, class LocalMatSolve = Kokkos::DefaultSparseSolve<MatScalar,LocalOrdinal,Node> >
  class CrsMatrixMultiplyOp : public Operator<OpScalar,LocalOrdinal,GlobalOrdinal,Node> {
    public:
      //! @name Constructor/Destructor Methods
      //@{ 

      //! Constructor
      CrsMatrixMultiplyOp(const Teuchos::RCP<const CrsMatrix<MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> > &A);

      //! Destructor
      virtual ~CrsMatrixMultiplyOp();

      //@}

      //! @name Methods implementing Operator
      //@{ 

      //! Computes this matrix-vector multilication Y = A X.
      //! This calls multiply<OpScalar,OpScalar>() on the underlying CrsMatrix object.
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
      void applyTranspose(const MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
                          OpScalar alpha, OpScalar beta) const;

      void applyNonTranspose(const MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
                             OpScalar alpha, OpScalar beta) const;
  };


  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  CrsMatrixMultiplyOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::CrsMatrixMultiplyOp(const Teuchos::RCP<const CrsMatrix<MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> > &A) 
  : matrix_(A) {
    // we don't require that A is fill complete; we will query for the importer/exporter at apply()-time
  }

  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  CrsMatrixMultiplyOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::~CrsMatrixMultiplyOp() {
  }

  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void 
  CrsMatrixMultiplyOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::apply(
              const MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> & X_in, MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> &Y_in,
              Teuchos::ETransp mode, OpScalar alpha, OpScalar beta) const 
  {
    TEST_FOR_EXCEPTION(!matrix_->isFillComplete(), std::runtime_error, 
        Teuchos::typeName(*this) << "::apply(): underlying matrix is not fill-complete.");
    TEST_FOR_EXCEPTION(X_in.getNumVectors() != Y_in.getNumVectors(), std::runtime_error,
        Teuchos::typeName(*this) << "::apply(X,Y): X and Y must have the same number of vectors.");
    TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<OpScalar>::isComplex && mode == Teuchos::TRANS, std::logic_error,
        Teuchos::typeName(*this) << "::apply() does not currently support transposed multiplications for complex scalar types.");
    if (mode == Teuchos::NO_TRANS) {
      applyNonTranspose(X_in, Y_in, alpha, beta);
    }
    else {
      applyTranspose(X_in, Y_in, alpha, beta);
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void 
  CrsMatrixMultiplyOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::applyNonTranspose(
      const MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> & X_in, 
            MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> & Y_in,
            OpScalar alpha, OpScalar beta) const 
  {
    typedef Teuchos::ScalarTraits<OpScalar> ST;
    using Teuchos::null;

    const int myImageID = Teuchos::rank(*matrix_->getComm());
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    if (myImageID == 0) {
      *out << "Entering CrsMatrixMultiplyOp::applyNonTranspose()" << std::endl
                << "Column Map: " << std::endl;
    }
    *out << matrix_->getColMap() << std::endl;
    if (myImageID == 0) {
      *out << "Initial input: " << std::endl;
    }
    X_in.describe(*out,Teuchos::VERB_EXTREME);
#endif

    const size_t numVectors = X_in.getNumVectors();
    // because of Views, it is difficult to determine if X and Y point to the same data. 
    // however, if they reference the exact same object, we will do the user the favor of copying X into new storage (with a warning)
    // we ony need to do this if we have trivial importers; otherwise, we don't actually apply the operator from X into Y
    const Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > importer = matrix_->getGraph()->getImporter();
    const Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > exporter = matrix_->getGraph()->getExporter();
    // access X indirectly, in case we need to create temporary storage
    Teuchos::RCP<const MV> X;

    // some parameters for below
    const bool Y_is_replicated = !Y_in.isDistributed(),
               Y_is_overwritten = (beta == ST::zero());
    if (Y_is_replicated && myImageID > 0) {
      beta = ST::zero();
    }

    // currently, cannot multiply from multivector of non-constant stride
    if (X_in.isConstantStride() == false && importer == null) {
      // generate a strided copy of X_in 
      X = Teuchos::rcp(new MV(X_in));
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) *out << "X is not constant stride, duplicating X results in a strided copy" << std::endl;
      X->describe(*out,Teuchos::VERB_EXTREME);
#endif
    }
    else {
      // just temporary, so this non-owning RCP is okay
      X = Teuchos::rcp(&X_in, false);
    }

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

    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if (importer != null) {
      importMV_->doImport(X_in, *importer, INSERT);
      // multiply out of importMV_
      X = importMV_;
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) {
        *out << "Performed import of X using importer..." << std::endl;
      }
      X->describe(*out,Teuchos::VERB_EXTREME);
#endif
    }

    // If we have a non-trivial exporter, we must export elements that are permuted or go to other processors
    // We will compute solution into the to-be-exported MV
    if (exporter != null) {
      // Do actual computation
      matrix_->template multiply<OpScalar,OpScalar>(*X, *exportMV_, Teuchos::NO_TRANS, alpha, ST::zero());
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) *out << "Export vector after multiply()..." << std::endl;
      exportMV_->describe(*out,Teuchos::VERB_EXTREME);
#endif
      if (Y_is_overwritten) Y_in.putScalar(ST::zero());
      else                  Y_in.scale(beta);
      Y_in.doExport(*exportMV_, *exporter, ADD);
    }
    // otherwise, multiply into Y
    else {
      // can't multiply in-situ; can't mutiply into non-strided multivector
      if (Y_in.isConstantStride() == false || X.getRawPtr() == &Y_in) {
        // generate a strided copy of Y 
        MV Y(Y_in);
        matrix_->template multiply<OpScalar,OpScalar>(*X, Y, Teuchos::NO_TRANS, alpha, beta);
        Y_in = Y;
      }
      else {
        matrix_->template multiply<OpScalar,OpScalar>(*X, Y_in, Teuchos::NO_TRANS, alpha, beta);
      }
    }
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
    if (myImageID == 0) *out << "Y_in vector after multiply/export..." << std::endl;
    Y_in.describe(*out,Teuchos::VERB_EXTREME);
#endif
    // Handle case of rangemap being a local replicated map: in this case, sum contributions from each processor
    if (Y_is_replicated) {
      Y_in.reduce();
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) *out << "Output vector is local; result after reduce()..." << std::endl;
      Y_in.describe(*out,Teuchos::VERB_EXTREME);
#endif
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void 
  CrsMatrixMultiplyOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::applyTranspose(
               const MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> & X_in, 
                     MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> & Y_in,
               OpScalar alpha, OpScalar beta) const 
  {
    typedef Teuchos::ScalarTraits<OpScalar> ST;
    using Teuchos::null;

    int myImageID = Teuchos::rank(*matrix_->getComm());
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    if (myImageID == 0) {
      *out << "Entering CrsMatrixMultiplyOp::applyTranspose()" << std::endl
                << "Column Map: " << std::endl;
    }
    *out << matrix_->getColMap() << std::endl;
    if (myImageID == 0) {
      *out << "Initial input: " << std::endl;
    }
    X_in.describe(*out,Teuchos::VERB_EXTREME);
#endif

    const size_t numVectors = X_in.getNumVectors();
    // because of Views, it is difficult to determine if X and Y point to the same data. 
    // however, if they reference the exact same object, we will do the user the favor of copying X into new storage (with a warning)
    // we ony need to do this if we have trivial importers; otherwise, we don't actually apply the operator from X into Y
    Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > importer = matrix_->getGraph()->getImporter();
    Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > exporter = matrix_->getGraph()->getExporter();
    // access X indirectly, in case we need to create temporary storage
    Teuchos::RCP<const MV> X;

    // some parameters for below
    const bool Y_is_replicated = !Y_in.isDistributed(),
               Y_is_overwritten = (beta == ST::zero());
    if (Y_is_replicated && myImageID > 0) {
      beta = ST::zero();
    }

    // currently, cannot multiply from multivector of non-constant stride
    if (X_in.isConstantStride() == false && importer==null) {
      // generate a strided copy of X_in 
      X = Teuchos::rcp(new MV(X_in));
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) *out << "X is not constant stride, duplicating X results in a strided copy" << std::endl;
      X->describe(*out,Teuchos::VERB_EXTREME);
#endif
    }
    else {
      // just temporary, so this non-owning RCP is okay
      X = Teuchos::rcp(&X_in, false);
    }

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

    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors
    if (exporter != null) {
      exportMV_->doImport(X_in,*exporter,INSERT);
      // multiply out of exportMV_
      X = exportMV_;
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) {
        *out << "Performed import of X using exporter..." << std::endl;
      }
      X->describe(*out,Teuchos::VERB_EXTREME);
#endif
    }

    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    // We will compute colutioni into the to-be-exported MV; get a view
    if (importer != null) {
      // Do actual computation
      matrix_->template multiply<OpScalar,OpScalar>(*X, *importMV_, Teuchos::CONJ_TRANS, alpha, ST::zero());
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) *out << "Import vector after multiply()..." << std::endl;
      importMV_->describe(*out,Teuchos::VERB_EXTREME);
#endif
      if (Y_is_overwritten) Y_in.putScalar(ST::zero());
      else                  Y_in.scale(beta);
      Y_in.doExport(*importMV_,*importer,ADD);
    }
    // otherwise, multiply into Y
    else {
      // can't multiply in-situ; can't multiply into non-strided multivector
      if (Y_in.isConstantStride() == false || X.getRawPtr() == &Y_in) {
        // generate a strided copy of Y
        MV Y(Y_in);
        matrix_->template multiply<OpScalar,OpScalar>(*X, Y, Teuchos::CONJ_TRANS, alpha, beta);
        Y_in = Y;
      }
      else {
        matrix_->template multiply<OpScalar,OpScalar>(*X, Y_in, Teuchos::CONJ_TRANS, alpha, beta);
      }
    }
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
    if (myImageID == 0) *out << "Y_in vector after multiply/export..." << std::endl;
    Y_in.describe(*out,Teuchos::VERB_EXTREME);
#endif
    // Handle case of rangemap being a local replicated map: in this case, sum contributions from each processor
    if (Y_is_replicated) {
      Y_in.reduce();
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) *out << "Output vector is local; result after reduce()..." << std::endl;
      Y_in.describe(*out,Teuchos::VERB_EXTREME);
#endif
    }
  }

  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool 
  CrsMatrixMultiplyOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::hasTransposeApply() const {
    return true;
  }

  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsMatrixMultiplyOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getDomainMap() const {
    return matrix_->getDomainMap();
  }

  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsMatrixMultiplyOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getRangeMap() const {
    return matrix_->getRangeMap();
  }

  //! Non-member function to create CrsMatrixMultiplyOp
  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  Teuchos::RCP< CrsMatrixMultiplyOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> >
  createCrsMatrixMultiplyOp(const Teuchos::RCP<const CrsMatrix<MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> > &A) {
    return rcp(new CrsMatrixMultiplyOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>(A) );
  }

} // end of namespace Tpetra

#endif // TPETRA_CRSMATRIXMULTIPLYOP_HPP
