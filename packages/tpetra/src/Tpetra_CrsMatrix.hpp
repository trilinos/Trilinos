//@HEADER
// ************************************************************************
// 
//               Tpetra: Templated Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
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

#ifndef TPETRA_CRSMATRIX_HPP
#define TPETRA_CRSMATRIX_HPP

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_CompileTimeAssert.hpp>
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsGraph.hpp"

// TODO: Add static graph construction option to constructors, a la Epetra

namespace Tpetra {
  // struct for i,j,v triplets
  template <class Ordinal, class Scalar>
  struct CrsIJV {
    CrsIJV() {}
    CrsIJV(Ordinal row, Ordinal col, const Scalar &val) {
      i = row;
      j = col;
      v = val;
    }
    Ordinal i,j;
    Scalar  v;
  };
}

namespace Teuchos {
  // SerializationTraits specialization for CrsIJV, using DirectSerialization
  template<typename Ordinal, typename Scalar>
  class SerializationTraits<int,Tpetra::CrsIJV<Ordinal,Scalar> >
  : public DirectSerializationTraits<int,Tpetra::CrsIJV<Ordinal,Scalar> >
  {};
}

namespace std {
  template <class Ordinal, class Scalar>
  bool operator<(const Tpetra::CrsIJV<Ordinal,Scalar> &ijv1, const Tpetra::CrsIJV<Ordinal,Scalar> &ijv2) {
    return ijv1.i < ijv2.i;
  }
}

namespace Tpetra 
{
  //! Tpetra::CrsMatrix: A class for constructing and using sparse compressed index matrices and row access.
  /*!
   * This class allows the construction of sparse matrices with row-access. 
   * Method insertGlobalValues() can be used to set both locally
   * owned and non-local elements; the shipping of data is done with hardcoded
   * MPI calls when fillComplete() is called.
   *
   * The nonzero elements of  locally owned row can be accessed by method
   * extractMyRowCopy() or extractGlobalRowCopy(). The former returns the column
   * indices using local numbering, the latter using global numbering.
   *
   */
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal=LocalOrdinal>
  class CrsMatrix : public RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal>
  {
    public:
      //! @name Constructor/Destructor Methods
      //@{ 

      //! Constructor specifying the number of non-zeros for all rows.
      CrsMatrix(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, Teuchos_Ordinal numNNZ);

      //! Constructor specifying the number of non-zeros for each row.
      CrsMatrix(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, const Teuchos::ArrayView<Teuchos_Ordinal> &NNZPerRowToAlloc);

      // !Destructor.
      virtual ~CrsMatrix();

      //@}

      //! @name Methods implementing Operator
      //@{ 

      //! Returns the Map associated with the domain of this operator.
      const Map<LocalOrdinal,GlobalOrdinal> & getDomainMap() const;

      //! Returns the Map associated with the domain of this operator.
      const Map<LocalOrdinal,GlobalOrdinal> & getRangeMap() const;

      //! Computes the matrix-vector multilication y = A x.
      void apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>& X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

      //@}

      //! @name Methods implementing RowMatrix
      //@{ 

      //! Returns \c true if fillComplete() has been called.
      bool isFillComplete() const;

      //! Returns \c true if optimizeStorage() has been called.
      bool isStorageOptimized() const;

      //! Dictates that the graph is static, so that new entries cannot be added to the matrix. */
      bool isStaticGraph() const;

      //! \brief If matrix indices are stored as local indices, this function returns true. Otherwise, it returns false.
      bool indicesAreLocal() const;

      //! \brief If matrix indices are stored as global indices, this function returns false. Otherwise, it returns true.
      bool indicesAreGlobal() const;

      //! \brief Indicates whether the matrix is lower triangular.
      bool lowerTriangular() const;

      //! \brief Indicates whether the matrix is upper triangular.
      bool upperTriangular() const;

      //! Returns the communicator.
      Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

      //! \brief Indicates whether the matrix has a well-defined column map. 
      /*! The column map does not exist until after fillComplete(), unless the matrix was constructed with one. */
      bool hasColMap() const; 

      //! Returns the RowGraph associated with this matrix. 
      const RowGraph<LocalOrdinal,GlobalOrdinal> &getGraph() const;

      //! Returns the Map that describes the row distribution in this matrix.
      const Map<LocalOrdinal,GlobalOrdinal> & getRowMap() const;

      //! \brief Returns the Map that describes the column distribution in this matrix.
      /*! Throws exception if getColMap() == false. */
      const Map<LocalOrdinal,GlobalOrdinal> & getColMap() const;

      //! Returns the number of global matrix rows. 
      GlobalOrdinal numGlobalRows() const;

      //! \brief Returns the number of global matrix columns. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      GlobalOrdinal numGlobalCols() const;

      //! Returns the number of matrix rows owned by the calling image. 
      Teuchos_Ordinal numLocalRows() const;

      //! \brief Returns the number of matrix columns needed by the calling image to apply the forward operator.
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      Teuchos_Ordinal numLocalCols() const;

      //! Returns the index base for global indices for this matrix. 
      Teuchos_Ordinal getIndexBase() const;

      //! \brief Returns the number of nonzero entries in the global matrix. 
      /*! Returns the number of global entries in the associated graph. */
      GlobalOrdinal numGlobalEntries() const;

      //! \brief Returns the number of nonzero entries in the calling image's portion of the matrix. 
      /*! Before fillComplete() is called, this could include duplicated entries. */
      Teuchos_Ordinal numMyEntries() const;

      //! \brief Returns the current number of nonzero entries on this node in the specified global row .
      /*! Throws exception std::runtime_error if the specified global row does not belong to this node. */
      Teuchos_Ordinal numEntriesForGlobalRow(GlobalOrdinal globalRow) const;

      //! Returns the current number of nonzero entries on this node in the specified local row.
      /*! Throws exception std::runtime_error if the specified local row is not valid for this node. */
      Teuchos_Ordinal numEntriesForMyRow(LocalOrdinal localRow) const;

      //! \brief Returns the number of global nonzero diagonal entries, based on global row/column index comparisons. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      GlobalOrdinal numGlobalDiagonals() const;

      //! \brief Returns the number of local nonzero diagonal entries, based on global row/column index comparisons. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      Teuchos_Ordinal numMyDiagonals() const;

      //! \brief Returns the maximum number of nonzero entries across all rows/columns on all images. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      GlobalOrdinal globalMaxNumRowEntries() const;

      //! \brief Returns the maximum number of nonzero entries across all rows/columns on this image. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      Teuchos_Ordinal myMaxNumRowEntries() const;

      //@}

      //! @name Data Entry Methods
      //@{ 

      //! Submit multiple entries, using global IDs.
      /*! All index values must be in the global space. */
      void insertGlobalValues(GlobalOrdinal globalRow, 
                         const Teuchos::ArrayView<const GlobalOrdinal> &cols,
                         const Teuchos::ArrayView<const Scalar>        &vals);


      //! Submit a single entry, using global IDs.
      /*! All index values must be in the global space. */
      void insertGlobalValue(GlobalOrdinal globalRow, GlobalOrdinal globalCol, Scalar value);

      //! Set all matrix entries equal to scalarThis.
      void setAllToScalar(const Scalar &alpha);

      //! Scale the current values of a matrix, this = alpha*this. 
      void scale(const Scalar &alpha);

      /*! \brief Signal that data entry is complete, specifying domain and range maps. 
          Off-node entries are distributed, repeated entries are summed, and global indices are transformed to local indices.
          If \c OptimizeStorage is true, then optimizeStorage() is called as well.
       */ 
      void fillComplete(const Map<LocalOrdinal,GlobalOrdinal> &domainMap, const Map<LocalOrdinal,GlobalOrdinal> &rangeMap, bool OptimizeStorage=false);

      /*! \brief Signal that data entry is complete. 
          Off-node entries are distributed, repeated entries are summed, and global indices are transformed to local indices.
          If \c OptimizeStorage is true, then optimizeStorage() is called as well.
          \note This method calls fillComplete( getRowMap(), getRowMap() ).
       */
      void fillComplete(bool OptimizeStorage=false);

      //! \brief Re-allocate the data into contiguous storage.
      void optimizeStorage();

      // @}

      //! @name Data Access Methods
      // @{ 

      //! \brief Get a copy of the diagonal entries owned by this node, with local row idices.
      /*! Returns a distributed Vector object partitioned according to the matrix's row map, containing the 
          the zero and non-zero diagonals owned by this node. */
      void getLocalDiagCopy(Vector<Scalar,LocalOrdinal,GlobalOrdinal> &diag) const;

      //! Returns a copy of the specified (and locally owned) row, using local indices.
      /*! Before fillComplete(), the results will not include entries submitted to another node and may contain duplicated entries.
       * \pre hasColMap() == true
       */
      void extractMyRowCopy(Teuchos_Ordinal localRow, 
                            const Teuchos::ArrayView<LocalOrdinal> &indices, 
                            const Teuchos::ArrayView<Scalar> &values,
                            Teuchos_Ordinal &numEntries) const;

      //! Returns a copy of the specified (and locally owned) row, using global indices.
      /*! Before fillComplete(), the results will not include entries submitted to another node and may contain duplicated entries. */
      void extractGlobalRowCopy(GlobalOrdinal globalRow, 
                                const Teuchos::ArrayView<GlobalOrdinal> &indices,
                                const Teuchos::ArrayView<Scalar> &values,
                                Teuchos_Ordinal &numEntries) const;

      //@}

      //! @name I/O Methods
      //@{ 
      
      //! Prints the matrix on the specified stream. This is very verbose.
      void print(std::ostream& os) const;

      // @}

    protected:
      // copy constructor disabled
      CrsMatrix(const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal> &Source);
      // operator= disabled
      CrsMatrix& operator=(const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal> &rhs);
      // useful typedefs
      typedef Teuchos::OrdinalTraits<LocalOrdinal>    LOT;
      typedef Teuchos::OrdinalTraits<GlobalOrdinal>   GOT;
      typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> MV;

      void allocateValues();
      void sortEntries();
      void mergeRedundantEntries();

      // Performs importing of off-processor elements and adds them to the locally owned elements.
      void globalAssemble();
      // multiplication routines
      void GeneralMV (typename MV::const_pointer x       , typename MV::pointer y       ) const;
      void GeneralMM (typename MV::const_double_pointer X, typename MV::double_pointer Y, Teuchos_Ordinal numVectors) const;
      void GeneralMhV(typename MV::const_pointer x       , typename MV::pointer y       ) const;
      void GeneralMhM(typename MV::const_double_pointer X, typename MV::double_pointer Y, Teuchos_Ordinal numVectors) const;

      CrsGraph<LocalOrdinal,GlobalOrdinal> graph_;
      bool fillComplete_;

      // values
      // before optimizeStorage()
      Teuchos::Array<Teuchos::ArrayRCP<Scalar> > values_;
      // after optimizeStorage()
      Teuchos::ArrayRCP<Scalar> contigValues_;
      /* valuesPtrs_[j] is always the begin() iterator from an ArrayView of the source or contigValues_
         of the appropriate length. in a debug build, it is an ArrayRCP, which does bounds checking. in an optimized
         build, it is a C pointer. valuesPtrs_ is allocated to numLocalRows()+1; the span of the jth row begins with
         valuesPtrs_[j] and ends before valuesPtrs_[j+1] */
      Teuchos::Array<typename Teuchos::ArrayRCP<Scalar>::iterator> valuesPtrs_;
      // multivectors used for import/export dest/source in apply()
      mutable Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> > importMV_, exportMV_;

      // a map between a (non-local) row and a list of (col,val)
      // TODO: this functionality will be divided between CrsGraph and CrsMatrix after the former comes into existence
      std::map<GlobalOrdinal, std::list<std::pair<GlobalOrdinal,Scalar> > > nonlocals_;

  }; // class CrsMatrix



  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::CrsMatrix(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, Teuchos_Ordinal maxNNZPerRow)
  : graph_(rowMap,maxNNZPerRow)
  , fillComplete_(false)
  , values_(rowMap.getNumMyEntries())
  , contigValues_(Teuchos::null)
  , valuesPtrs_(0)
  , importMV_(Teuchos::null)
  , exportMV_(Teuchos::null)
  {
    TEST_FOR_EXCEPTION(maxNNZPerRow < 1 && maxNNZPerRow != 0, std::runtime_error,
        Teuchos::typeName(*this) << "::CrsMatrix(rowMap,maxNNZPerRow): maxNNZPerRow must be non-negative.");
    allocateValues();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::CrsMatrix(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, const Teuchos::ArrayView<Teuchos_Ordinal> &NNZPerRowToAlloc)
  : graph_(rowMap,NNZPerRowToAlloc)
  , fillComplete_(false)
  , values_(rowMap.getNumMyEntries())
  , contigValues_(Teuchos::null)
  , valuesPtrs_(0)
  , importMV_(Teuchos::null)
  , exportMV_(Teuchos::null)
  {
    TEST_FOR_EXCEPTION(NNZPerRowToAlloc.size() != rowMap.getNumMyEntries(), std::runtime_error,
        Teuchos::typeName(*this) << "::CrsMatrix(rowMap,NNZPerRowToAlloc): NNZPerRowToAlloc must have as many entries as specified by rowMap for this node.");
    allocateValues();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::allocateValues() 
  {
    Teuchos_Ordinal numlocal = getRowMap().getNumMyEntries();
    for (Teuchos_Ordinal i=0; i<numlocal; ++i) {
      Teuchos_Ordinal nnzi = graph_.numAllocatedEntriesForMyRow(i);
      if (nnzi > 0) {
        values_[i] = Teuchos::arcp<Scalar>(nnzi);
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::~CrsMatrix()
  {}


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::isFillComplete() const
  { return fillComplete_; }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::isStorageOptimized() const
  { return (contigValues_ != Teuchos::null); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::indicesAreLocal() const
  { return graph_.indicesAreLocal(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::indicesAreGlobal() const
  { return graph_.indicesAreGlobal(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::hasColMap() const
  { return graph_.hasColMap(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP<const Teuchos::Comm<int> > 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getComm() const
  { return graph_.getComm(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::numGlobalEntries() const
  { return graph_.numGlobalEntries(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::numMyEntries() const
  { return graph_.numMyEntries(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::numGlobalRows() const
  { return graph_.numGlobalRows(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::numGlobalCols() const
  { return graph_.numGlobalCols(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::numLocalRows() const
  { return graph_.numLocalRows(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::numLocalCols() const
  { return graph_.numLocalCols(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::numGlobalDiagonals() const
  { return graph_.numGlobalDiagonals(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::numMyDiagonals() const
  { return graph_.numMyDiagonals(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::numEntriesForGlobalRow(GlobalOrdinal globalRow) const
  { return graph_.numEntriesForGlobalRow(globalRow); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::numEntriesForMyRow(LocalOrdinal localRow) const
  { return graph_.numEntriesForMyRow(localRow); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::globalMaxNumRowEntries() const
  { return graph_.globalMaxNumRowEntries(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::myMaxNumRowEntries() const
  { return graph_.myMaxNumRowEntries(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getIndexBase() const
  { return getRowMap().getIndexBase(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  const Map<LocalOrdinal,GlobalOrdinal> & 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getRowMap() const 
  { return graph_.getRowMap(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  const Map<LocalOrdinal,GlobalOrdinal> & 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getColMap() const 
  { return graph_.getColMap(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  const Map<LocalOrdinal,GlobalOrdinal> & 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getDomainMap() const
  { return graph_.getDomainMap(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  const Map<LocalOrdinal,GlobalOrdinal> & 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getRangeMap() const
  { return graph_.getRangeMap(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::insertGlobalValues(GlobalOrdinal globalRow, 
                         const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                         const Teuchos::ArrayView<const Scalar>  &values)
  {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(isFillComplete() == true || isStorageOptimized() == true, std::runtime_error,
      Teuchos::typeName(*this) << "::insertGlobalValues(): fillComplete() has already been called.");
    /* this version can always be called:
       a) either indices are still global, or
       b) indices are local, in which case we have a column map (either from construction or from 
          our graph), so that we can translate the given global indices */
    TEST_FOR_EXCEPTION(values.size() != indices.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::insertGlobalValues(): values.size() must equal indices.size().");
    typename Teuchos::ArrayView<const GlobalOrdinal>::iterator ind = indices.begin();
    typename Teuchos::ArrayView<const Scalar       >::iterator val =  values.begin();
    if (getRowMap().isMyGlobalIndex(globalRow)) {
      Teuchos_Ordinal myRow = getRowMap().getLocalIndex(globalRow),
                     rowNNZ = numEntriesForMyRow(myRow),
                      toAdd = indices.size(),
                   rowAlloc = graph_.numAllocatedEntriesForMyRow(myRow);
      if (rowNNZ+toAdd > rowAlloc) {
#       if defined(THROW_TPETRA_EFFICIENCY_WARNINGS) || defined(PRINT_TPETRA_EFFICIENCY_WARNINGS)
          std::string err = Teuchos::typeName(*this) + "::insertGlobalValues(): Pre-allocated space has been exceeded, requiring new allocation. To improve efficiency, suggest larger allocation.";
#         if defined(THROW_TPETRA_EFFICIENCY_WARNINGS)
            TEST_FOR_EXCEPTION(true, std::runtime_error, err);
#         else
            std::cerr << err << std::endl;
#         endif
#       endif
        // assumption: the number allocated for each row is stored in the graph.
        // increase the allocation, copy old entries to new storage
        rowAlloc = rowNNZ+toAdd;
        ArrayRCP<Scalar> newVals = Teuchos::arcp<Scalar>(Teuchos::as<Teuchos_Ordinal>(rowAlloc));
        std::copy(values_[myRow].begin(), values_[myRow].begin()+rowNNZ, newVals.begin());
        values_[myRow] = newVals;
      }
      // insert indices and values
      graph_.insertGlobalIndices(globalRow,indices);
      std::copy( values.begin(), values.end(),  values_[myRow].begin()+rowNNZ);
#ifdef HAVE_TPETRA_DEBUG
      // the assumption is that graph_.numAllocatedEntriesForMyRow is the allocated size here
      TEST_FOR_EXCEPTION( rowAlloc != graph_.numAllocatedEntriesForMyRow(myRow) 
                          || rowNNZ+toAdd != numEntriesForMyRow(myRow), std::logic_error,
          Teuchos::typeName(*this) << "::insertGlobalValues(): Internal logic error. Please contact Tpetra team.");
#endif
    }
    else {
      for (; val != values.end(); ++val, ++ind) {
        nonlocals_[globalRow].push_back(std::make_pair(*ind, *val));
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::insertGlobalValue(GlobalOrdinal globalRow, GlobalOrdinal globalCol, Scalar value)
  { insertGlobalValues(globalRow,Teuchos::arrayView(&globalCol,1),Teuchos::arrayView(&value,1)); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::setAllToScalar(const Scalar &alpha)
  { 
    TEST_FOR_EXCEPT(isStorageOptimized());
    for (Teuchos_Ordinal r=0; r<numLocalRows(); ++r) {
      std::fill(values_[r].begin(), values_[r].begin()+numEntriesForMyRow(r), alpha);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::scale(const Scalar &alpha)
  { 
    TEST_FOR_EXCEPT(isStorageOptimized());
    for (Teuchos_Ordinal r=0; r<numLocalRows(); ++r) {
      typename Teuchos::ArrayRCP<Scalar>::iterator val, vend;
      val = values_[r].begin();
      vend = values_[r].begin()+numEntriesForMyRow(r);
      while (val != vend) {
        (*val++) *= alpha;
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::extractMyRowCopy(Teuchos_Ordinal myRow, 
                                                                     const Teuchos::ArrayView<LocalOrdinal> &indices, 
                                                                     const Teuchos::ArrayView<Scalar>       &values,
                                                                     Teuchos_Ordinal &numEntries) const
  {
    numEntries = numEntriesForMyRow(myRow);
    TEST_FOR_EXCEPTION(getRowMap().isMyLocalIndex(myRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::extractMyRowCopy(myRow,...): specified row (==" << myRow << ") is not valid on this node.");
    TEST_FOR_EXCEPTION(indices.size() < numEntries || values.size() < numEntries, std::runtime_error, 
        Teuchos::typeName(*this) << "::extractMyRowCopy(myRow,indices,values): size of indices,values must be sufficient to store the specified row.");
#ifdef HAVE_TPETRA_DEBUG
    Teuchos_Ordinal nnzagain;
    graph_.extractMyRowCopy(myRow,indices,nnzagain);
    TEST_FOR_EXCEPTION(nnzagain != numEntries, std::logic_error, 
        Teuchos::typeName(*this) << "::extractMyRowCopy(): Internal logic error. Please contact Tpetra team.");
#else
    graph_.extractMyRowCopy(myRow,indices,numEntries);
#endif
    std::copy( values_[myRow].begin(), values_[myRow].begin()+numEntries, values.begin() );
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::extractGlobalRowCopy(GlobalOrdinal globalRow, 
                                                                      const Teuchos::ArrayView<GlobalOrdinal> &indices,
                                                                      const Teuchos::ArrayView<Scalar>  &values,
                                                                      Teuchos_Ordinal &numEntries) const
  {
    // Only locally owned rows can be queried, otherwise complain
    Teuchos_Ordinal myRow = getRowMap().getLocalIndex(globalRow);
    TEST_FOR_EXCEPTION(myRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::extractGlobalRowCopy(globalRow,...): globalRow does not belong to this node.");
    numEntries = graph_.numAllocatedEntriesForMyRow(myRow);
    TEST_FOR_EXCEPTION(
        indices.size() < numEntries || values.size() < numEntries, std::runtime_error, 
        Teuchos::typeName(*this) << "::extractGlobalRowCopy(globalRow,indices,values): size of indices,values must be sufficient to store the specified row.");
#ifdef HAVE_TPETRA_DEBUG
    Teuchos_Ordinal nnzagain;
    graph_.extractGlobalRowCopy(globalRow,indices,nnzagain);
    TEST_FOR_EXCEPTION(nnzagain != numEntries, std::logic_error, 
        Teuchos::typeName(*this) << "::extractMyRowCopy(): Internal logic error. Please contact Tpetra team.");
#else
    graph_.extractGlobalRowCopy(globalRow,indices,numEntries);
#endif
    std::copy( values_[myRow].begin(), values_[myRow].begin()+numEntries, values.begin() );
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &X, 
                                                                 MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &Y, 
                                                                 Teuchos::ETransp mode) const
  {
    // TODO: add support for alpha,beta term coefficients: Y = alpha*A*X + beta*Y
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::null;
    using Teuchos::ArrayView;
    TEST_FOR_EXCEPTION(!isFillComplete(), std::runtime_error, 
        Teuchos::typeName(*this) << ": cannot call apply() until fillComplete() has been called.");
    TEST_FOR_EXCEPTION(X.numVectors() != Y.numVectors(), std::runtime_error,
        Teuchos::typeName(*this) << "::apply(X,Y): X and Y must have the same number of vectors.");

#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
    int myImageID = Teuchos::rank(*getComm());
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    if (myImageID == 0) {
      *out << "Entering CrsMatrix::apply()" << std::endl
                << "Column Map: " << std::endl;
    }
    *out << this->getColMap() << std::endl;
    if (myImageID == 0) {
      *out << "Initial input: " << std::endl;
    }
    X.print(*out); X.printValues(*out);
#   endif

    Teuchos_Ordinal numVectors = X.numVectors();
    // because of Views, it is difficult to determine if X and Y point to the same data. 
    // however, if they reference the exact same object, we will do the user the favor of copying X into new storage (with a warning)
    // we ony need to do this if we have trivial importers; otherwise, we don't actually apply the operator from X into Y
    Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal> > importer = graph_.getImporter();
    Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal> > exporter = graph_.getExporter();
    Teuchos::RCP<const MV> Xcopy;
    typename MV::const_double_pointer Xdata = X.extractConstView2D();
    typename MV::double_pointer       Ydata = Y.extractView2D();
    if (&X==&Y && importer==null && exporter==null) {
#     if defined(THROW_TPETRA_EFFICIENCY_WARNINGS) || defined(PRINT_TPETRA_EFFICIENCY_WARNINGS)
      std::string err = Teuchos::typeName(*this) + "::apply(X,Y): If X and Y are the same, it necessitates a temporary copy of X, which is inefficient.";
#     if defined(THROW_TPETRA_EFFICIENCY_WARNINGS)
      TEST_FOR_EXCEPTION(true, std::runtime_error, err);
#     else
      std::cerr << err << std::endl;
#     endif
#     endif
      // generate a copy of X 
      Xcopy = Teuchos::rcp(new MV(X));
      Xdata = Xcopy->extractConstView2D();
#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) *out << "X and Y are co-located, duplicating X results in a stride copy" << std::endl;
      *out << this->getColMap() << std::endl;
      Xcopy->print(*out); Xcopy->printValues(*out);
#   endif
    }
    if (importer != null) {
      if (importMV_ != null && importMV_->numVectors() != numVectors) importMV_ = null;
      if (importMV_ == null) {
        importMV_ = Teuchos::rcp( new MV(getColMap(),numVectors) );
      }
    }
    if (exporter != null) {
      if (exportMV_ != null && exportMV_->numVectors() != numVectors) exportMV_ = null;
      if (exportMV_ == null) {
        exportMV_ = Teuchos::rcp( new MV(getRowMap(),numVectors) );
      }
    }
    if (mode == Teuchos::NO_TRANS) {
      // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
      if (importer != null) {
        importMV_->doImport(X, *importer, INSERT);
        Xdata = importMV_->extractConstView2D();
#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
        if (myImageID == 0) {
          *out << "Performed import of X using importer..." << std::endl;
        }
        importMV_->print(*out); importMV_->printValues(*out);
#   endif
      }
      // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
      // We will compute solution into the to-be-exported MV; get a view
      if (exporter != null) {
        Ydata = exportMV_->extractView2D();
      }
      // Do actual computation
      if (numVectors==1) {
        GeneralMV(Xdata[0],Ydata[0]);
      }
      else {
        GeneralMM(Xdata,Ydata,numVectors);
      }
#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) *out << "Matrix-MV product..." << std::endl;
      if (exportMV_ != null) {
        exportMV_->print(*out); exportMV_->printValues(*out);
      } else {
        Y.print(*out); Y.printValues(*out);
      }
#   endif
      // do the export
      if (exporter != null) {
        Y.putScalar(ST::zero());  // Make sure target is zero: necessary because we are adding. may need adjusting for alpha,beta apply()
        Y.doExport(*exportMV_, *exporter, ADD); // Fill Y with Values from export vector
#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
        if (myImageID == 0) *out << "Output vector after export() using exporter..." << std::endl;
        Y.print(*out); Y.printValues(*out);
#   endif
      }
      // Handle case of rangemap being a local replicated map: in this case, sum contributions from each processor
      if (Y.isDistributed() == false) {
        Y.reduce();
#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
        if (myImageID == 0) *out << "Output vector is local; result after reduce()..." << std::endl;
        Y.print(*out); Y.printValues(*out);
#   endif
      }
    }
    else {
      // mode == CONJ_TRANS or TRANS
      TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<Scalar>::isComplex && mode == Teuchos::TRANS, std::logic_error,
          Teuchos::typeName(*this) << "::apply() does not currently support transposed multiplications for complex scalar types.");
      // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors
      if (exporter != null) {
        exportMV_->doImport(X,*exporter,INSERT);
        Xdata = exportMV_->extractConstView2D();
#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
        if (myImageID == 0) {
          *out << "Performed import of X using exporter..." << std::endl;
        }
        exportMV_->print(*out); exportMV_->printValues(*out);
#   endif
      }
      // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
      // We will compute colutioni into the to-be-exported MV; get a view
      if (importer != null) {
        Ydata = importMV_->extractView2D();
      }
      // Do the actual transposed multiplication
      if (numVectors==1) {
        GeneralMhV(Xdata[0],Ydata[0]);
      }
      else {
        GeneralMhM(Xdata,Ydata,numVectors);
      }
#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) *out << "Matrix-MV product..." << std::endl;
      if (importMV_ != null) {
        importMV_->print(*out); importMV_->printValues(*out);
      } else {
        Y.print(*out); Y.printValues(*out);
      }
#   endif
      if (importer != null) {
        Y.putScalar(ST::zero()); // Make sure target is zero: necessary because we are adding. may need adjusting for alpha,beta apply()
        Y.doExport(*importMV_,*importer,ADD);
#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
        if (myImageID == 0) *out << "Output vector after export() using importer..." << std::endl;
        Y.print(*out); Y.printValues(*out);
#   endif
      }
      // Handle case of rangemap being a local replicated map: in this case, sum contributions from each processor
      if (Y.isDistributed() == false) {
        Y.reduce();
#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
        if (myImageID == 0) *out << "Output vector is local; result after reduce()..." << std::endl;
        Y.print(*out); Y.printValues(*out);
#   endif
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::print(std::ostream& os) const 
  {
    // support this operation whether fillComplete() or not
    using std::endl;
    int myImageID = Teuchos::rank(*getComm());
    if (myImageID == 0)
    {
      os << "Tpetra::CrsMatrix, label = " << this->getObjectLabel() << endl;
      os << "Number of global rows    = " << numGlobalRows() << endl;
      if (isFillComplete())
      {
        os << "Number of global columns    = " << numGlobalCols() << endl;
        os << "Status = fill complete" << endl;
        os << "Number of global nonzeros   = " << numGlobalEntries() << endl;
        os << "Global max nonzeros per row = " << globalMaxNumRowEntries() << endl;
      }
      else
      {
        os << "Status = fill not complete" << endl;
      }
    }
    if (isFillComplete())
    {
      for (int pid=0; pid < Teuchos::size(*getComm()); ++pid)
      {
        if (pid == myImageID)
        {
          Teuchos::Array<GlobalOrdinal> indices(myMaxNumRowEntries());
          Teuchos::Array<Scalar>         values(myMaxNumRowEntries());
          Teuchos_Ordinal rowSize;
          os << "% Number of rows on image " << myImageID << " = " << numLocalRows() << endl;
          for (Teuchos_Ordinal i=0; i < numLocalRows(); ++i)
          {
            GlobalOrdinal globalRow = getRowMap().getGlobalIndex(i);
            extractGlobalRowCopy(globalRow, indices(), values(), rowSize);
            if (rowSize > Teuchos::OrdinalTraits<Teuchos_Ordinal>::zero()) {
              for (Teuchos_Ordinal j=0; j < rowSize; ++j) {
                os << "Matrix(" << globalRow << ", " << indices[j] << ") = " << values[j] << endl;
              }
            }
          }
        }
        Teuchos::barrier(*getComm());
        Teuchos::barrier(*getComm());
        Teuchos::barrier(*getComm());
      }
    }
  }



  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::fillComplete(const Map<LocalOrdinal,GlobalOrdinal> &domainMap, const Map<LocalOrdinal,GlobalOrdinal> &rangeMap, bool OptimizeStorage)
  {
    if (getComm()->getSize() > 1) {
      globalAssemble();
    }
    else {
      TEST_FOR_EXCEPTION(nonlocals_.size() > 0, std::runtime_error,
          Teuchos::typeName(*this) << "::fillComplete(): cannot have non-local entries on a serial run. Invalid entry was submitted to the CrsMatrix.");
    }

    if (graph_.isFillComplete()) {
      if (!isStaticGraph() && !isFillComplete()) {
        return;
      }
    }

    if (!isStaticGraph()) {
      graph_.makeIndicesLocal(domainMap,rangeMap);
    }

    sortEntries();
    mergeRedundantEntries();

    if (!isStaticGraph()) {
      graph_.fillComplete(domainMap,rangeMap);
    }
  
    fillComplete_ = true;

    if (OptimizeStorage) optimizeStorage();
  }



  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::fillComplete(bool OptimizeStorage)
  {
    fillComplete(getRowMap(),getRowMap(),OptimizeStorage);
  }



  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::optimizeStorage()
  {
    TEST_FOR_EXCEPT(true);
  }



  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::globalAssemble()
  {
    using Teuchos::arcp;
    using Teuchos::OrdinalTraits;
    using Teuchos::Array;
    using Teuchos::SerialDenseMatrix;
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    using Teuchos::rcp;
    using std::list;
    using std::pair;
    using std::make_pair;
    typedef typename std::map<GlobalOrdinal,std::list<pair<GlobalOrdinal,Scalar> > >::const_iterator NLITER;
    int numImages = getComm()->getSize();
    int myImageID = getComm()->getRank();
    // Determine if any nodes have global entries to share
    Teuchos_Ordinal MyNonlocals = nonlocals_.size(), 
                    MaxGlobalNonlocals;
    Teuchos::reduceAll<Teuchos_Ordinal>(*getComm(),Teuchos::REDUCE_MAX,MyNonlocals,&MaxGlobalNonlocals);
    if (MaxGlobalNonlocals == 0) return;  // no entries to share

    // compute a list of NLRs from nonlocals_ and use it to compute:
    //      IdsAndRows: a vector of (id,row) pairs
    //          NLR2Id: a map from NLR to the Id that owns it
    // globalNeighbors: a global graph of connectivity between images: globalNeighbors(i,j) indicates that j sends to i
    //         sendIDs: a list of all images I send to
    //         recvIDs: a list of all images I receive from (constructed later)
    Array<pair<int,GlobalOrdinal> > IdsAndRows;
    std::map<GlobalOrdinal,int> NLR2Id;
    SerialDenseMatrix<int,char> globalNeighbors;
    Array<int> sendIDs, recvIDs;
    {
      // nonlocals_ contains the entries we are holding for all non-local rows
      // we want a list of the rows for which we have data
      Array<GlobalOrdinal> NLRs;
      std::set<GlobalOrdinal> setOfRows;
      for (NLITER iter = nonlocals_.begin(); iter != nonlocals_.end(); ++iter)
      {
        setOfRows.insert(iter->first);
      }
      // copy the elements in the set into an Array
      NLRs.resize(setOfRows.size());
      std::copy(setOfRows.begin(), setOfRows.end(), NLRs.begin());

      // get a list of ImageIDs for the non-local rows (NLRs)
      Array<int> NLRIds(NLRs.size());
      {
        bool invalidGIDs = getRowMap().getRemoteIndexList(NLRs(),NLRIds());
        char lclerror = ( invalidGIDs ? 1 : 0 );
        char gblerror;
        Teuchos::reduceAll(*getComm(),Teuchos::REDUCE_MAX,lclerror,&gblerror);
        TEST_FOR_EXCEPTION(gblerror, std::runtime_error,
            Teuchos::typeName(*this) << "::globalAssemble(): non-local entries correspond to invalid rows.");
      }

      // build up a list of neighbors, as well as a map between NLRs and Ids
      // localNeighbors[i] != 0 iff I have data to send to image i
      // put NLRs,Ids into an array of pairs
      IdsAndRows.reserve(NLRs.size());
      Array<char> localNeighbors(numImages,0);
      typename Array<GlobalOrdinal>::const_iterator nlr;
      typename Array<int>::const_iterator id;
      for (nlr = NLRs.begin(), id = NLRIds.begin();
           nlr != NLRs.end(); ++nlr, ++id) 
      {
        NLR2Id[*nlr] = *id;
        localNeighbors[*id] = 1;
        IdsAndRows.push_back(make_pair<int,GlobalOrdinal>(*id,*nlr));
      }
      for (int j=0; j<numImages; ++j)
      {
        if (localNeighbors[j]) {
          sendIDs.push_back(j);
        }
      }
      // sort IdsAndRows, by Ids first, then rows
      std::sort(IdsAndRows.begin(),IdsAndRows.end());
      // gather from other nodes to form the full graph
      globalNeighbors.shapeUninitialized(numImages,numImages);
      Teuchos::gatherAll(*getComm(),numImages,localNeighbors.getRawPtr(),numImages*numImages,globalNeighbors.values());
      // globalNeighbors at this point contains (on all images) the
      // connectivity between the images. 
      // globalNeighbors(i,j) != 0 means that j sends to i/that i receives from j
    }

    ////////////////////////////////////////////////////////////////////////////////////// 
    // FIGURE OUT WHO IS SENDING TO WHOM AND HOW MUCH
    // DO THIS IN THE PROCESS OF PACKING ALL OUTGOING DATA ACCORDING TO DESTINATION ID
    ////////////////////////////////////////////////////////////////////////////////////// 

    // loop over all columns to know from which images I can expect to receive something
    for (int j=0; j<numImages; ++j)
    {
      if (globalNeighbors(myImageID,j)) {
        recvIDs.push_back(j);
      }
    }
    Teuchos_Ordinal numRecvs = recvIDs.size();

    // we know how many we're sending to already
    // form a contiguous list of all data to be sent
    // track the number of entries for each ID
    Array<CrsIJV<GlobalOrdinal,Scalar> > IJVSendBuffer;
    Array<Teuchos_Ordinal> sendSizes(sendIDs.size(), 0);
    Teuchos_Ordinal numSends = 0;
    for (typename Array<pair<int,GlobalOrdinal> >::const_iterator IdAndRow = IdsAndRows.begin();
         IdAndRow != IdsAndRows.end(); ++IdAndRow) 
    {
      int            id = IdAndRow->first;
      GlobalOrdinal row = IdAndRow->second;
      // have we advanced to a new send?
      if (sendIDs[numSends] != id) {
        numSends++;
        TEST_FOR_EXCEPTION(sendIDs[numSends] != id, std::logic_error, Teuchos::typeName(*this) << "::globalAssemble(): internal logic error. Contact Tpetra team.");
      }
      // copy data for row into contiguous storage
      for (typename list<pair<GlobalOrdinal,Scalar> >::const_iterator jv = nonlocals_[row].begin(); jv != nonlocals_[row].end(); ++jv)
      {
        IJVSendBuffer.push_back( CrsIJV<GlobalOrdinal,Scalar>(row,jv->first,jv->second) );
        sendSizes[numSends]++;
      }
    }
    if (IdsAndRows.size() > 0) {
      numSends++; // one last increment, to make it a count instead of an index
    }
    TEST_FOR_EXCEPTION(Teuchos::as<typename Array<int>::size_type>(numSends) != sendIDs.size(), std::logic_error, Teuchos::typeName(*this) << "::globalAssemble(): internal logic error. Contact Tpetra team.");

    // don't need this data anymore
    nonlocals_.clear();

    ////////////////////////////////////////////////////////////////////////////////////// 
    // TRANSMIT SIZE INFO BETWEEN SENDERS AND RECEIVERS
    ////////////////////////////////////////////////////////////////////////////////////// 
    // perform non-blocking sends: send sizes to our recipients
    Array<Teuchos::RCP<Teuchos::CommRequest> > sendRequests;
    for (int s=0; s < numSends ; ++s)
    {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      sendRequests.push_back( Teuchos::isend<int,int>(*getComm(),rcp(&sendSizes[s],false),sendIDs[s]) );
    }
    // perform non-blocking receives: receive sizes from our senders
    Array<Teuchos::RCP<Teuchos::CommRequest> > recvRequests;
    Array<Teuchos_Ordinal> recvSizes(numRecvs);
    for (int r=0; r < numRecvs; ++r) {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      recvRequests.push_back( Teuchos::ireceive(*getComm(),rcp(&recvSizes[r],false),recvIDs[r]) );
    }
    // wait on all 
    if (!sendRequests.empty()) {
      Teuchos::waitAll(*getComm(),sendRequests());
    }
    if (!recvRequests.empty()) {
      Teuchos::waitAll(*getComm(),recvRequests());
    }
    Teuchos::barrier(*getComm());
    sendRequests.clear();
    recvRequests.clear();

    ////////////////////////////////////////////////////////////////////////////////////
    // NOW SEND/RECEIVE ALL ROW DATA
    ////////////////////////////////////////////////////////////////////////////////////
    // from the size info, build the ArrayViews into IJVSendBuffer
    Array<ArrayView<CrsIJV<GlobalOrdinal,Scalar> > > sendBuffers(numSends,Teuchos::null);
    {
      Teuchos_Ordinal cur = 0;
      for (Teuchos_Ordinal s=0; s<numSends; ++s) {
        sendBuffers[s] = IJVSendBuffer(cur,sendSizes[s]);
        cur += sendSizes[s];
      }
    }
    // perform non-blocking sends
    for (Teuchos_Ordinal s=0; s < numSends ; ++s)
    {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      ArrayRCP<CrsIJV<GlobalOrdinal,Scalar> > tmparcp = arcp(sendBuffers[s].getRawPtr(),0,sendBuffers[s].size(),false);
      sendRequests.push_back( Teuchos::isend<int,CrsIJV<GlobalOrdinal,Scalar> >(*getComm(),tmparcp,sendIDs[s]) );
    }
    // calculate amount of storage needed for receives
    // setup pointers for the receives as well
    Teuchos_Ordinal totalRecvSize = std::accumulate(recvSizes.begin(),recvSizes.end(),0);
    Array<CrsIJV<GlobalOrdinal,Scalar> > IJVRecvBuffer(totalRecvSize);
    // from the size info, build the ArrayViews into IJVRecvBuffer
    Array<ArrayView<CrsIJV<GlobalOrdinal,Scalar> > > recvBuffers(numRecvs,Teuchos::null);
    {
      Teuchos_Ordinal cur = 0;
      for (Teuchos_Ordinal r=0; r<numRecvs; ++r) {
        recvBuffers[r] = IJVRecvBuffer(cur,recvSizes[r]);
        cur += recvSizes[r];
      }
    }
    // perform non-blocking recvs
    for (Teuchos_Ordinal r=0; r < numRecvs ; ++r)
    {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      ArrayRCP<CrsIJV<GlobalOrdinal,Scalar> > tmparcp = arcp(recvBuffers[r].getRawPtr(),0,recvBuffers[r].size(),false);
      recvRequests.push_back( Teuchos::ireceive(*getComm(),tmparcp,recvIDs[r]) );
    }
    // perform waits
    if (!sendRequests.empty()) {
      Teuchos::waitAll(*getComm(),sendRequests());
    }
    if (!recvRequests.empty()) {
      Teuchos::waitAll(*getComm(),recvRequests());
    }
    Teuchos::barrier(*getComm());
    sendRequests.clear();
    recvRequests.clear();


    ////////////////////////////////////////////////////////////////////////////////////
    // NOW PROCESS THE RECEIVED ROW DATA
    ////////////////////////////////////////////////////////////////////////////////////
    // TODO: instead of adding one entry at a time, add one row at a time.
    //       this requires resorting; they arrived sorted by sending node, so that entries could be non-contiguous if we received
    //       multiple entries for a particular row from different processors.
    //       it also requires restoring the data, which may make it not worth the trouble.
    for (typename Array<CrsIJV<GlobalOrdinal,Scalar> >::const_iterator ijv = IJVRecvBuffer.begin(); ijv != IJVRecvBuffer.end(); ++ijv)
    {
      insertGlobalValue(ijv->i, ijv->j, ijv->v);
    }

    // WHEW! THAT WAS TIRING!
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::GeneralMV(typename MV::const_pointer x, typename MV::pointer y) const
  {
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    ArrayView<const LocalOrdinal> cinds;
    typename ArrayView<const LocalOrdinal>::iterator cind;
    typename ArrayRCP<const Scalar>::const_iterator aval;
    for (Teuchos_Ordinal r=0; r < getRowMap().getNumMyEntries(); ++r) {
      graph_.extractMyRowConstView(r,cinds);
      Scalar sum = ST::zero();
      cind = cinds.begin();
      aval = values_[r].begin(); 
      while (cind != cinds.end()) 
      {
        sum += (*aval++) * x[*cind++];
      }
      y[r] = sum;
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::GeneralMM (typename MV::const_double_pointer X, typename MV::double_pointer Y, Teuchos_Ordinal numVectors) const
  {
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    Teuchos_Ordinal nlrs = numLocalRows();
    ArrayView<const LocalOrdinal> cinds;
    typename ArrayView<const LocalOrdinal>::iterator cind;
    typename ArrayRCP<const Scalar>::const_iterator aval;
    for (Teuchos_Ordinal r=0; r < nlrs; ++r) {
      graph_.extractMyRowConstView(r,cinds);
      for (Teuchos_Ordinal j=0; j<numVectors; ++j) {
        typename MV::pointer       yvals = Y[j];
        typename MV::const_pointer xvals = X[j]; 
        cind = cinds.begin(); 
        aval = values_[r].begin(); 
        Scalar sum = ST::zero();
        while (cind != cinds.end()) {
          sum += (*aval++) * xvals[*cind++];
        }
        yvals[r] = sum;
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::GeneralMhV(typename MV::const_pointer x, typename MV::pointer y) const
  { 
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    Teuchos_Ordinal nlrs = numLocalRows(),
                    nlcs = numLocalCols();
    // Initialize y for transpose multiply
    std::fill( y, y+nlcs, ST::zero() );
    // apply conjugate transpose of matrix to x
    // use column triad formulation, accumulating into y each column of A^H (each row of A) times each entry in x
    ArrayView<const LocalOrdinal> cinds;
    typename ArrayView<const LocalOrdinal>::iterator cind;
    typename ArrayRCP<const Scalar>::iterator aval;
    for (Teuchos_Ordinal r=0; r < nlrs; ++r) {
      graph_.extractMyRowConstView(r,cinds);
      // loop over entries in this column of A^H (this row of A)
      cind = cinds.begin();
      aval = values_[r].begin();
      while (cind != cinds.end()) {
        y[*cind++] += ST::conjugate(*aval++) * x[r];
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::GeneralMhM(typename MV::const_double_pointer X, typename MV::double_pointer Y, Teuchos_Ordinal numVectors) const
  {
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    Teuchos_Ordinal nlrs = numLocalRows(),
                    nlcs = numLocalCols();
    // Initialize Y for transpose multiply
    for (Teuchos_Ordinal j=0; j<numVectors; ++j) {
      std::fill( Y[j], Y[j]+nlcs, ST::zero() );
    }
    // apply conjugate transpose of matrix to X
    // use outer-product formulation, hitting Y with a rank-1 update comprised of each column of A^H (each row of A) and each row of X
    ArrayView<const LocalOrdinal> cinds;
    typename ArrayView<const LocalOrdinal>::iterator cind;
    typename ArrayRCP<const Scalar>::iterator aval;
    for (Teuchos_Ordinal r=0; r < nlrs; ++r) {
      graph_.extractMyRowConstView(r,cinds);
      // loop over numvectors
      for (Teuchos_Ordinal j=0; j<numVectors; ++j) {
        typename MV::pointer       yvals = Y[j];
        typename MV::const_pointer xvals = X[j]; 
        // loop over entries in this column of A^H (this row of A)
        cind = cinds.begin();
        aval = values_[r].begin();
        while (cind != cinds.end()) {
          yvals[*cind++] += ST::conjugate(*aval++) * xvals[r];
        }
      }
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getLocalDiagCopy(Vector<Scalar,LocalOrdinal,GlobalOrdinal> &dvec) const
  {
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
        Teuchos::typeName(*this) << ": cannot call getLocalDiagCopy() until fillComplete() has been called.");
    TEST_FOR_EXCEPTION(!dvec.getMap().isSameAs(getRowMap()), std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalDiagCopy(dvec): dvec must have the same map as the CrsMatrix.");
#ifdef HAVE_TPETRA_DEBUG
    int numDiagFound = 0;
#endif
    Teuchos::ArrayView<Scalar> values;
    dvec.extractView1D(values);
    Teuchos_Ordinal nlrs = numLocalRows();
    typename Teuchos::ArrayRCP<const Scalar>::iterator v;
    typename Teuchos::ArrayView<Scalar>::iterator ov;
    Teuchos::ArrayView<const LocalOrdinal> cinds;
    typename Teuchos::ArrayView<const LocalOrdinal>::iterator i;
    ov = values.begin();
    for (Teuchos_Ordinal r=0; r < nlrs; ++r) {
      *ov = Teuchos::ScalarTraits<Scalar>::zero();
      GlobalOrdinal rgid = getRowMap().getGlobalIndex(r);
      if (getColMap().isMyGlobalIndex(rgid)) {
        LocalOrdinal rlid = getColMap().getLocalIndex(rgid);
        graph_.extractMyRowConstView(r,cinds);
        i = cinds.begin();
        v = values_[r].begin();
        while (i != cinds.end()) {
          if (*i == rlid) {
            *ov = *v;
#ifdef HAVE_TPETRA_DEBUG
            ++numDiagFound;
#endif
            break;
          }
          else if (*i > rlid) {
            // indices are sorted; we have passed the diagonal, so quit
            break;
          }
          v++;
          i++;
        }
      }
      ++ov;
    }
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(numDiagFound != graph_.numMyDiagonals(), std::logic_error, 
        "CrsMatrix::getLocalDiagCopy(): logic error. Please contact Tpetra team.");
#endif
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::isStaticGraph() const
  { return false; }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  const RowGraph<LocalOrdinal,GlobalOrdinal> &
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getGraph() const
  { return graph_; }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::sortEntries()
  {
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(graph_.indicesAreGlobal() == true, std::logic_error,
        Teuchos::typeName(*this) << "::sortEntries(): sortEntries() must be called after indices are transformed to local.\n"
        << "Likely internal logic error. Please contact Tpetra team.");
    if (graph_.indicesAreSorted()) return;
    if (graph_.indicesAreAllocated()) {
      const Teuchos_Ordinal nlrs = numLocalRows();
      for (Teuchos_Ordinal r=0; r < nlrs; ++r)
      {
        ArrayView<LocalOrdinal> inds;
        typename ArrayRCP<Scalar>::iterator vals = values_[r].begin();
        graph_.extractMyRowView(r,inds);
        sort2(inds.begin(), inds.end(), vals);
      }
    }
    graph_.indicesAreSorted(true);
    return;
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::lowerTriangular() const
  { return graph_.lowerTriangular(); }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::upperTriangular() const
  { return graph_.upperTriangular(); }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::mergeRedundantEntries() 
  {
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(graph_.indicesAreSorted() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::mergeRedundantEntries() cannot be called before indices are sorted.\n"
        << "Likely interal logic error. Please contact Tpetra team.");
    for (Teuchos_Ordinal r=0; r<numLocalRows(); ++r) 
    {
      Teuchos_Ordinal rnnz = numEntriesForMyRow(r);
      if (rnnz > 1) {
        ArrayView<const LocalOrdinal> inds;
        graph_.extractMyRowConstView(r,inds);
        typename ArrayRCP<Scalar>::iterator vals = values_[r].begin();
        Teuchos_Ordinal curEntry = 0;
        Scalar curValue = vals[0];
        for (Teuchos_Ordinal k=1; k<rnnz; ++k) {
          if (inds[k] == inds[k-1]) {
            curValue += vals[k];
          }
          else {
            vals[curEntry++] = curValue;
            curValue = vals[k];
          }
        }
        vals[curEntry] = curValue;
      }
    }
    graph_.removeRedundantIndices();
  }

} // namespace Tpetra

#endif
