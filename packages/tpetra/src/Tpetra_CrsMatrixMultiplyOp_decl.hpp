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

#ifndef TPETRA_CRSMATRIXMULTIPLYOP_DECL_HPP
#define TPETRA_CRSMATRIXMULTIPLYOP_DECL_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Operator.hpp"
#include <Teuchos_TimeMonitor.hpp>


/*! \file Tpetra_CrsMatrixMultiplyOp_decl.hpp 

    The declarations for the class Tpetra::CrsMatrixMultiplyOp and related non-member constructors.
 */

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS  
  // forward declaration
  template <class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class CrsMatrix;
#endif

  //! \brief A class for wrapping a CrsMatrix multiply in a Operator.
  template <class Scalar, class MatScalar = Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<MatScalar,LocalOrdinal,Node>::SparseOps >
  class CrsMatrixMultiplyOp : public Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
    public:
      //! @name Constructor/Destructor Methods
      //@{ 

      //! Constructor
      CrsMatrixMultiplyOp(const Teuchos::RCP<const CrsMatrix<MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &A);

      //! Destructor
      virtual ~CrsMatrixMultiplyOp();

      //@}

      //! @name Methods implementing Operator
      //@{ 

      //! Computes this matrix-vector multilication Y = A X.
      //! This calls multiply<Scalar,Scalar>() on the underlying CrsMatrix object.
      void apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
                 Teuchos::ETransp mode = Teuchos::NO_TRANS, Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(), Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;

      /// \brief "Hybrid" Jacobi + (Gauss-Seidel or SOR) on \f$B = A X\f$.
      ///
      /// "Hybrid" means Jacobi for interprocess communication, but
      /// Successive Over-Relaxation (SOR) or Gauss-Seidel for
      /// intraprocess computation.  Gauss-Seidel is a special case of
      /// SOR, where the damping factor is one.
      ///
      /// The Forward or Backward sweep directions have their usual
      /// SOR meaning within the process.  Interprocess communication
      /// occurs once before the sweep, as it would in Jacobi.
      ///
      /// The Symmetric sweep direction means first Forward, then
      /// Backward.  Before each sweep is an interprocess
      /// communication, as in Jacobi.  Thus, Symmetric results in two
      /// interprocess communication steps.
      ///
      /// \param B [in] Right-hand side(s).
      /// \param X [in/out] On input: initial guess(es).  On output:
      ///   result multivector(s).
      /// \param D [in] Inverse of diagonal entries of the matrix A.
      /// \param dampingFactor [in] SOR damping factor.  A damping
      ///   factor of one results in Gauss-Seidel.
      /// \param direction [in] Sweep direction: Forward, Backward, or
      ///   Symmetric.
      /// \param numSweeps [in] Number of sweeps.  We count each
      ///   Symmetric sweep (including both its Forward and its
      ///   Backward sweep) as one.
      void 
      gaussSeidel (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B,
		   MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
		   const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &D,
		   const Scalar& dampingFactor,
		   const ESweepDirection direction,
		   const int numSweeps) const;

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
      typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

      // underlying CrsMatrix
      const Teuchos::RCP<const CrsMatrix<MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > matrix_;

      // multivectors used for import/export dest/source in apply()
      mutable Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > importMV_, exportMV_;

#ifdef HAVE_KOKKOSCLASSIC_CUDA_NODE_MEMORY_PROFILING
    Teuchos::RCP<Teuchos::Time> importTimer_, exportTimer_;
#endif

      // private methods for transpose or non-transpose
      void applyTranspose(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
                          Scalar alpha, Scalar beta) const;

      void applyNonTranspose(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
                             Scalar alpha, Scalar beta) const;
  };

  /*! \brief Non-member function to create CrsMatrixMultiplyOp

      \relatesalso CrsMatrixMultiplyOp
   */
  template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::RCP< CrsMatrixMultiplyOp<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> >
  createCrsMatrixMultiplyOp(const Teuchos::RCP<const CrsMatrix<MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &A);

} // end of namespace Tpetra

#endif // TPETRA_CRSMATRIXMULTIPLYOP_DECL_HPP
