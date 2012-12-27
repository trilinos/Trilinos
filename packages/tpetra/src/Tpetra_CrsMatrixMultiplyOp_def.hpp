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

#ifndef TPETRA_CRSMATRIXMULTIPLYOP_DEF_HPP
#define TPETRA_CRSMATRIXMULTIPLYOP_DEF_HPP

#include "Tpetra_CrsMatrix.hpp"

#ifdef DOXYGEN_USE_ONLY
  #include "Tpetra_CrsMatrixMultiplyOp_decl.hpp"
#endif

#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
  #include "Teuchos_VerboseObject.hpp"
#endif

/*! \file Tpetra_CrsMatrixMultiplyOp_def.hpp 

    The implementations for the members of Tpetra::CrsMatrixMultiplyOp and related non-member constructors.
 */

namespace Tpetra {

  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CrsMatrixMultiplyOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::CrsMatrixMultiplyOp(const Teuchos::RCP<const CrsMatrix<MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &A) 
  : matrix_(A) {
    // we don't require that A is fill complete; we will query for the importer/exporter at apply()-time
#ifdef HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_PROFILING
    importTimer_ = Teuchos::TimeMonitor::getNewTimer( "CrsMatrixMultiplyOp::import" );
    exportTimer_ = Teuchos::TimeMonitor::getNewTimer( "CrsMatrixMultiplyOp::export" );
#endif
  }

  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CrsMatrixMultiplyOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::~CrsMatrixMultiplyOp() {
  }

  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void 
  CrsMatrixMultiplyOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::apply(
              const MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> & X_in, 
                    MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> &Y_in,
                    Teuchos::ETransp mode, OpScalar alpha, OpScalar beta) const 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(!matrix_->isFillComplete(), std::runtime_error, 
        Teuchos::typeName(*this) << "::apply(): underlying matrix is not fill-complete.");
    TEUCHOS_TEST_FOR_EXCEPTION(X_in.getNumVectors() != Y_in.getNumVectors(), std::runtime_error,
        Teuchos::typeName(*this) << "::apply(X,Y): X and Y must have the same number of vectors.");
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<OpScalar>::isComplex && mode == Teuchos::TRANS, std::logic_error,
        Teuchos::typeName(*this) << "::apply() does not currently support transposed multiplications for complex scalar types.");
    if (mode == Teuchos::NO_TRANS) {
      applyNonTranspose(X_in, Y_in, alpha, beta);
    }
    else {
      applyTranspose(X_in, Y_in, alpha, beta);
    }
  }


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void 
  CrsMatrixMultiplyOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  gaussSeidel (const MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> &B,
	       MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> &X,
	       const MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> &D,
	       const OpScalar& dampingFactor,
	       const ESweepDirection direction) const
  {
    using Teuchos::null;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcpFromRef;
    typedef Teuchos::ScalarTraits<OpScalar> STS;
    typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
    typedef Export<LocalOrdinal, GlobalOrdinal, Node> export_type;
    typedef Import<LocalOrdinal, GlobalOrdinal, Node> import_type;
    typedef MultiVector<OpScalar, LocalOrdinal, GlobalOrdinal, Node> MV;

    // We don't need the Export object because this method assumes
    // that the row, domain, and range Maps are the same.  We do need
    // the Import object, if there is one, though.
    RCP<const import_type> importer = matrix_->getGraph()->getImporter();
    RCP<const export_type> exporter = matrix_->getGraph()->getExporter();
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! exporter.is_null (), 
      std::runtime_error,
      "gaussSeidel requires that the row, domain, and range Maps be the same.  "
      "This cannot be the case, because the matrix has a nontrivial Export object.");

#ifdef TEUCHOS_DEBUG
    {
      RCP<const map_type> domainMap = matrix_->getDomainMap ();
      RCP<const map_type> rangeMap = matrix_->getRangeMap ();
      RCP<const map_type> rowMap = matrix_->getGraph ()->getRowMap ();

      // The relation 'isSameAs' is transitive.  It's also a
      // collective, so we don't have to do a "shared" test for
      // exception (i.e., a global reduction on the test value).
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! domainMap->isSameAs (*rangeMap) || ! domainMap->isSameAs (*rowMap),
        std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidel requires that the row, domain, and "
	"range Maps of the matrix be the same.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! X.getMap ()->isSameAs (*domainMap),
        std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidel requires that the input / output "
	"multivector X be in the domain Map of the matrix.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! B.getMap ()->isSameAs (*rangeMap),
        std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidel requires that the input multivector B "
	"be in the range Map of the matrix.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! D.getMap ()->isSameAs (*rowMap),
        std::runtime_error,
        "Tpetra::CrsMatrix::gaussSeidel requires that the input multivector D "
	"be in the row Map of the matrix.");
    }
#endif // TEUCHOS_DEBUG

    RCP<const map_type> colMap = matrix_->getGraph ()->getColMap ();

    // The Gauss-Seidel / SOR kernel expects multivectors of constant
    // stride.  If either of them is not, we have to duplicate them.
    // Do this before making the column Map view of X, since the view
    // has to be the view of the multivector used in the kernel.
    // Issue a performance warning if we had to do this.  It's really
    // a problem with our kernel, not the user's problem, but it's
    // still a good idea to warn them.

    RCP<MV> X_in;
    RCP<MV> X_colMap;
    bool copiedInput;
    if (X.isConstantStride()) {
      X_in = rcpFromRef (X);
      copiedInput = false;
    }
    else {
      // Don't just make a copy of X; make a column Map multivector,
      // make a domain Map view of it, and copy X into the domain Map
      // view.
      X_colMap = rcp (new MV (colMap, X.getNumVectors ()));
      X_in = X_colMap->offsetViewNonConst (X.getMap (), 0);
      *X_in = X; // deep copy
      copiedInput = true;
      TPETRA_EFFICIENCY_WARNING(
        ! X.isConstantStride (), 
	std::runtime_error, 
	"gaussSeidel: The current implementation of the Gauss-Seidel kernel "
	"requires that X and B both have constant stride.  Since X does not "
	"have constant stride, we had to make a copy.  This is a limitation of "
	"the current implementation and not your fault, but we still report it "
	"as an efficiency warning for your information.");
    }

    RCP<const MV> B_in;
    if (B.isConstantStride()) {
      B_in = rcpFromRef (B);
    }
    else {
      B_in = rcp (new MV (B));
      TPETRA_EFFICIENCY_WARNING(
        ! B.isConstantStride (), 
	std::runtime_error, 
	"gaussSeidel: The current implementation of the Gauss-Seidel kernel "
	"requires that X and B both have constant stride.  Since B does not "
	"have constant stride, we had to make a copy.  This is a limitation of "
	"the current implementation and not your fault, but we still report it "
	"as an efficiency warning for your information.");
    }

    // After this block of code, X_colMap is a view of X, distributed
    // using the column Map of the matrix.
    if (importer.is_null ()) {
      // The domain and column Maps are the same; no need to Import.
      // Don't assign to X_colMap if we've already made it above.
      if (! copiedInput) {
	X_colMap = X_in;
      }
    } 
    else {
      // mfh 26 Dec 2012: Alas, I can't currently ask a MultiVector
      // whether it is a view of something else.  This would let me
      // avoid making a new view.  It would also let me check the
      // precondition that X be a domain Map view of a column Map
      // multivector.
      X_colMap = X_in->offsetViewNonConst (colMap, 0);
      X_colMap->doImport (*X_in, *importer, INSERT);
    }

    // Do local Gauss-Seidel.
    if (direction != Symmetric) {
      Kokkos::ESweepDirection localDirection;
      if (direction == Forward) {
	localDirection = Kokkos::Forward;
      } else if (direction == Backward) {
	localDirection = Kokkos::Backward;
      } else {
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "The 'direction' "
          "enum does not have any of its valid values: Forward, Backward, or "
          "Symmetric.");
      }
      matrix_->template localGaussSeidel<OpScalar,OpScalar> (*B_in, *X_colMap, D, dampingFactor, localDirection);
    }
    else { // direction == Symmetric
      matrix_->template localGaussSeidel<OpScalar,OpScalar> (*B_in, *X_colMap, D, dampingFactor, Kokkos::Forward);
      // Communicate again before the Backward sweep.
      X_colMap->doImport (*X_in, *importer, INSERT);
      matrix_->template localGaussSeidel<OpScalar,OpScalar> (*B_in, *X_colMap, D, dampingFactor, Kokkos::Backward);
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void 
  CrsMatrixMultiplyOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::applyNonTranspose(
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
      {
#ifdef HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_PROFILING
        Teuchos::TimeMonitor lcltimer(*importTimer_);
#endif
        importMV_->doImport(X_in, *importer, INSERT);
      }
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
      matrix_->template localMultiply<OpScalar,OpScalar>(*X, *exportMV_, Teuchos::NO_TRANS, alpha, ST::zero());
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) *out << "Export vector after localMultiply()..." << std::endl;
      exportMV_->describe(*out,Teuchos::VERB_EXTREME);
#endif
      if (Y_is_overwritten) Y_in.putScalar(ST::zero());
      else                  Y_in.scale(beta);
      {
#ifdef HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_PROFILING
        Teuchos::TimeMonitor lcltimer(*exportTimer_);
#endif
        Y_in.doExport(*exportMV_, *exporter, ADD);
      }
    }
    // otherwise, multiply into Y
    else {
      // can't multiply in-situ; can't mutiply into non-strided multivector
      if (Y_in.isConstantStride() == false || X.getRawPtr() == &Y_in) {
        // generate a strided copy of Y 
        MV Y(Y_in);
        matrix_->template localMultiply<OpScalar,OpScalar>(*X, Y, Teuchos::NO_TRANS, alpha, beta);
        Y_in = Y;
      }
      else {
        matrix_->template localMultiply<OpScalar,OpScalar>(*X, Y_in, Teuchos::NO_TRANS, alpha, beta);
      }
    }
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
    if (myImageID == 0) *out << "Y_in vector after local multiply/export..." << std::endl;
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
  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void 
  CrsMatrixMultiplyOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::applyTranspose(
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
      {
#ifdef HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_PROFILING
        Teuchos::TimeMonitor lcltimer(*importTimer_);
#endif
        exportMV_->doImport(X_in,*exporter,INSERT);
      }
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
    // We will compute solution into the to-be-exported MV; get a view
    if (importer != null) {
      // Do actual computation
      matrix_->template localMultiply<OpScalar,OpScalar>(*X, *importMV_, Teuchos::CONJ_TRANS, alpha, ST::zero());
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) *out << "Import vector after localMultiply()..." << std::endl;
      importMV_->describe(*out,Teuchos::VERB_EXTREME);
#endif
      if (Y_is_overwritten) Y_in.putScalar(ST::zero());
      else                  Y_in.scale(beta);
      //
      {
#ifdef HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_PROFILING
        Teuchos::TimeMonitor lcltimer(*importTimer_);
#endif
        Y_in.doExport(*importMV_,*importer,ADD);
      }
    }
    // otherwise, multiply into Y
    else {
      // can't multiply in-situ; can't multiply into non-strided multivector
      if (Y_in.isConstantStride() == false || X.getRawPtr() == &Y_in) {
        // generate a strided copy of Y
        MV Y(Y_in);
        matrix_->template localMultiply<OpScalar,OpScalar>(*X, Y, Teuchos::CONJ_TRANS, alpha, beta);
        Y_in = Y;
      }
      else {
        matrix_->template localMultiply<OpScalar,OpScalar>(*X, Y_in, Teuchos::CONJ_TRANS, alpha, beta);
      }
    }
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
    if (myImageID == 0) *out << "Y_in vector after local multiply/export..." << std::endl;
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

  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool 
  CrsMatrixMultiplyOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::hasTransposeApply() const {
    return true;
  }

  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsMatrixMultiplyOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getDomainMap() const {
    return matrix_->getDomainMap();
  }

  template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsMatrixMultiplyOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getRangeMap() const {
    return matrix_->getRangeMap();
  }

} // Tpetra namespace

template <class OpScalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
Teuchos::RCP< Tpetra::CrsMatrixMultiplyOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> >
Tpetra::createCrsMatrixMultiplyOp(const Teuchos::RCP<const Tpetra::CrsMatrix<MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &A) {
  return Teuchos::rcp(new Tpetra::CrsMatrixMultiplyOp<OpScalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>(A) );
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

//! Explicit instantiation macro supporting the CrsMatrixMultiplyOp class. Instantiates the class, the non-member constructor, and the necessary CrsMatrix::multiply() member.
#define TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(OPSCALAR,MATSCALAR,LO,GO,NODE) \
  \
  template class CrsMatrixMultiplyOp< OPSCALAR , MATSCALAR , LO , GO , NODE >; \
  \
  template Teuchos::RCP< Tpetra::CrsMatrixMultiplyOp<OPSCALAR,MATSCALAR,LO,GO,NODE> > \
  createCrsMatrixMultiplyOp(const Teuchos::RCP<const Tpetra::CrsMatrix<MATSCALAR,LO,GO,NODE> > &A); \
  \
  template void CrsMatrix<MATSCALAR,LO,GO,NODE>::localMultiply<OPSCALAR,OPSCALAR>( \
        const MultiVector<OPSCALAR,LO,GO,NODE> &X, \
              MultiVector<OPSCALAR,LO,GO,NODE> &Y, \
              Teuchos::ETransp mode,               \
              OPSCALAR alpha, OPSCALAR beta        \
              ) const; \
  template void CrsMatrix<MATSCALAR,LO,GO,NODE>::localGaussSeidel<OPSCALAR,OPSCALAR>( \
        const MultiVector<OPSCALAR,LO,GO,NODE> &B, \
              MultiVector<OPSCALAR,LO,GO,NODE> &X, \
        const MultiVector<OPSCALAR,LO,GO,NODE> &D, \
        const OPSCALAR& alpha, \
        const Kokkos::ESweepDirection direction) const;

#define TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT_SINGLE(SCALAR,LO,GO,NODE) \
        TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(SCALAR,SCALAR,LO,GO,NODE)

#endif // TPETRA_CRSMATRIXMULTIPLYOP_DEF_HPP
