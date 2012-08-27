/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_CONTAINER_HPP
#define IFPACK2_CONTAINER_HPP
#include "Ifpack2_ConfigDefs.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Describable.hpp"
#include <iostream>

//! Ifpack2::Container: a pure virtual class for creating and solving local linear problems.

/*!
Class Ifpack2::Container provides the abstract interfaces for 
containers. A "container" is an object that hosts all it is necessary
to create, populate, and solve local linear problems. The local
linear problem matrix, B, is a submatrix of the local components
of a distributed matrix, A. The idea of container is to
specify the rows of A that are contained in B, then extract B from A,
and compute all it is necessary to solve a linear system in B.
Then, set starting solution (if necessary) and right-hand side for B,
and solve the linear system in B.

<P>A container should be used in the following manner:
- Create an container object, specifying the number of rows of B.
- If necessary, set parameters for the solution using
  setParameters().
- Initialize the container by calling initialize().
- Specify the ID of the local rows of A that are contained in B,
  using ID().
- Prepare the linear system solver using compute().
- set X and/or Y elements using getX() and getY().
- Solve the linear system using apply().
- Get the components of the computed solution using getY().

The number of vectors can be set using setNumVectors(), and it
is defaulted to 1.

<P>Containers are currently used by class Ifpack2::BlockRelaxation.

<P>Ifpack2::Container is a pure virtual class.
Two concrete implementations are provided in classes
Ifpack2::SparseContainer (that stores matrices in sparse the format
Tpetra::CrsMatrix) and Ifpack2::DenseContainer
(for relatively small matrices, as matrices are stored as
Tpetra::SerialDenseMatrix's).

Note: The MatrixType and InverseType can have different Scalars, ordinals (and even nodes).
You can mix and match so long as implicit conversions are available.
The most obvious use case for this are:
1) MatrixGlobalOrdinal==long long and InverseGlobalOrdinal==short
2) MatrixScalar=float and InverseScalar=double 

\date Last update Aug-12.
  
*/

namespace Ifpack2 {

template<class MatrixType, class InverseType>
class Container : public Teuchos::Describable {

public:
  typedef typename MatrixType::scalar_type          MatrixScalar;
  typedef typename MatrixType::local_ordinal_type   MatrixLocalOrdinal;
  typedef typename MatrixType::global_ordinal_type  MatrixGlobalOrdinal;
  typedef typename MatrixType::node_type            MatrixNode;

  typedef typename InverseType::scalar_type         InverseScalar;
  typedef typename InverseType::local_ordinal_type  InverseLocalOrdinal;
  typedef typename InverseType::global_ordinal_type InverseGlobalOrdinal;
  typedef typename InverseType::node_type           InverseNode;

  //! Destructor.
  virtual ~Container() {};

  //! Returns the number of rows of the matrix and X/Y.
  virtual size_t getNumRows() const = 0;

  //! Returns the ID associated to local row i. 
  /*!
   * The set of (local) rows assigned to this container is defined
   * by calling ID(i) = j, where i (from 0 to NumRows()) indicates
   * the container-row, and j indicates the local row in the calling
   * process.
   *
   * This is usually used to recorder the local row ID (on calling process)
   * of the i-th row in the container.
   */
  virtual MatrixLocalOrdinal & ID(const size_t i) = 0;

  //! Initializes the container, by performing all operations that only require matrix structure.
  virtual void initialize() = 0;

  //! Finalizes the linear system matrix and prepares for the application of the inverse.
  virtual void compute(const Teuchos::RCP<const Tpetra::RowMatrix<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode> >& Matrix) = 0;

  //! Sets all necessary parameters.
  virtual void setParameters(const Teuchos::ParameterList& List) = 0;

  //! Returns \c true is the container has been successfully initialized.
  virtual bool isInitialized() const = 0;

  //! Returns \c true is the container has been successfully computed.
  virtual bool isComputed() const = 0;

  //! Computes Y = alpha * M^{-1} X + beta*Y
  /*! Here the X and Y are the size of the global problem the container was extracted from to begin with.
   */ 
  virtual void apply(const Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& X,
		     Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& Y,
		     Teuchos::ETransp mode=Teuchos::NO_TRANS,
		     MatrixScalar alpha=Teuchos::ScalarTraits< MatrixScalar >::one(),
		     MatrixScalar beta=Teuchos::ScalarTraits< MatrixScalar >::zero())=0;

  //! Computes Y = alpha * diag(D) * M^{-1} (diag(D) X) + beta*Y
  /*! Here D, X and Y are the size of the global problem the container was extracted from to begin with.  D has to be 
      a single Vector, while X and Y can be MultiVectors.  This function is designed to support techniques with overlap,
      such as Schwarz methods.
   */ 
  virtual void weightedApply(const Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& X,
			     Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& Y,
			     const Tpetra::Vector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& D,
			     Teuchos::ETransp mode=Teuchos::NO_TRANS,
			     MatrixScalar alpha=Teuchos::ScalarTraits< MatrixScalar >::one(),
			     MatrixScalar beta=Teuchos::ScalarTraits< MatrixScalar >::zero())=0;


  //! Prints out basic information about the container.
  virtual std::ostream& print(std::ostream& os) const = 0;

};

template <class MatrixType, class InverseType>
inline std::ostream& operator<<(std::ostream& os, const Ifpack2::Container<MatrixType,InverseType> & obj)
{
  return(obj.print(os));
}

}

#endif // IFPACK2_CONTAINER_HPP
