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

#ifndef TPETRA_ROWMATRIX_HPP
#define TPETRA_ROWMATRIX_HPP

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include "Tpetra_Operator.hpp"
#include "Tpetra_RowGraph.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"

namespace Tpetra 
{
  //! \brief A pure virtual interface for row-partitioned matrices.
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal=LocalOrdinal>
  class RowMatrix : public Operator<Scalar,LocalOrdinal,GlobalOrdinal>
  {
    public:
      //! @name Constructor/Destructor Methods
      //@{ 

      // !Destructor.
      virtual ~RowMatrix() {};

      //@}

      //! @name Matrix Query Methods
      //@{ 

      //! Returns \c true if fillComplete() has been called.
      virtual bool isFillComplete() const = 0;

      //! \brief If matrix indices are stored as local indices, this function returns true. Otherwise, it returns false.
      virtual bool indicesAreLocal() const = 0;

      //! \brief If matrix indices are stored as global indices, this function returns false. Otherwise, it returns true.
      virtual bool indicesAreGlobal() const = 0;

      //! \brief Indicates whether the matrix is lower triangular.
      virtual bool lowerTriangular() const = 0;

      //! \brief Indicates whether the matrix is upper triangular.
      virtual bool upperTriangular() const = 0;

      //! Returns the communicator.
      virtual Teuchos::RCP<const Teuchos::Comm<int> > getComm() const = 0;

      //! \brief Indicates whether the matrix has a well-defined column map. 
      /*! The column map does not exist until after fillComplete(), unless the matrix was constructed with one. */
      virtual bool hasColMap() const = 0;

      //! Returns the RowGraph associated with this matrix. 
      virtual const RowGraph<LocalOrdinal,GlobalOrdinal> &getGraph() const = 0;

      //! Returns the Map that describes the row distribution in this matrix.
      virtual const Map<LocalOrdinal,GlobalOrdinal> & getRowMap() const = 0;

      //! \brief Returns the Map that describes the column distribution in this matrix.
      /*! Cannot be called before fillComplete(), unless the matrix was constructed with a column map. */
      virtual const Map<LocalOrdinal,GlobalOrdinal> & getColMap() const = 0;

      //! Returns the number of global matrix rows. 
      virtual GlobalOrdinal numGlobalRows() const = 0;

      //! \brief Returns the number of global matrix columns. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      virtual GlobalOrdinal numGlobalCols() const = 0;

      //! Returns the number of matrix rows owned by the calling image. 
      virtual Teuchos_Ordinal numLocalRows() const = 0;

      //! \brief Returns the number of matrix columns needed by the calling image to apply the forward operator.
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      virtual Teuchos_Ordinal numLocalCols() const = 0;

      //! Returns the index base for global indices for this matrix. 
      virtual Teuchos_Ordinal getIndexBase() const = 0;

      //! \brief Returns the number of nonzero entries in the global matrix. 
      /*! Returns the number of global entries in the associated graph. */
      virtual GlobalOrdinal numGlobalEntries() const = 0;

      //! \brief Returns the number of nonzero entries in the calling image's portion of the matrix. 
      /*! Before fillComplete() is called, this could include duplicated entries. */
      virtual Teuchos_Ordinal numMyEntries() const = 0;

      //! \brief Returns the current number of nonzero entries on this node in the specified global row .
      /*! Throws exception std::runtime_error if the specified global row does not belong to this node. */
      virtual Teuchos_Ordinal numEntriesForGlobalRow(GlobalOrdinal globalRow) const = 0;

      //! Returns the current number of nonzero entries on this node in the specified local row.
      /*! Throws exception std::runtime_error if the specified local row is not valid for this node. */
      virtual Teuchos_Ordinal numEntriesForMyRow(LocalOrdinal localRow) const = 0;

      //! \brief Returns the number of global nonzero diagonal entries, based on global row/column index comparisons. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      virtual GlobalOrdinal numGlobalDiagonals() const = 0;

      //! \brief Returns the number of local nonzero diagonal entries, based on global row/column index comparisons. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      virtual Teuchos_Ordinal numMyDiagonals() const = 0;

      //! \brief Returns the maximum number of nonzero entries across all rows/columns on all images. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      virtual GlobalOrdinal globalMaxNumRowEntries() const = 0;

      //! \brief Returns the maximum number of nonzero entries across all rows/columns on this image. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      virtual Teuchos_Ordinal myMaxNumRowEntries() const = 0;

      //@}

      //! @name Data Access Methods
      // @{ 

      //! \brief Get a copy of the diagonal entries owned by this node, with local row idices.
      /*! Returns a distributed Vector object partitioned according to the matrix's row map, containing the 
          the zero and non-zero diagonals owned by this node. */
      virtual void getLocalDiagCopy(Vector<Scalar,LocalOrdinal,GlobalOrdinal> &diag) const = 0;

      //! Returns a copy of the specified (and locally owned) row, using local indices.
      /*! Before fillComplete(), the results will not include entries submitted to another node and may contain duplicated entries.
       * \pre hasColMap() == true
       */
      virtual void extractMyRowCopy(LocalOrdinal localRow, 
                                    const Teuchos::ArrayView<LocalOrdinal> &indices, 
                                    const Teuchos::ArrayView<Scalar> &values,
                                    Teuchos_Ordinal &numEntries) const  = 0;

      //! Returns a copy of the specified (and locally owned) row, using global indices.
      /*! Before fillComplete(), the results will not include entries submitted to another node and may contain duplicated entries. */
      virtual void extractGlobalRowCopy(GlobalOrdinal globalRow, 
                                        const Teuchos::ArrayView<GlobalOrdinal> &indices,
                                        const Teuchos::ArrayView<Scalar> &values,
                                        Teuchos_Ordinal &numEntries) const = 0;

      //@}

      //! @name I/O Methods
      //@{ 
      
      //! Prints the matrix on the specified stream. This is very verbose.
      virtual void print(std::ostream& os) const = 0;

      // @}

  }; // class RowMatrix

} // namespace Tpetra

#endif

