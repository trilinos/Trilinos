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

#ifdef HAVE_KOKKOS_CUDA_NODE_MEMORY_PROFILING
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
