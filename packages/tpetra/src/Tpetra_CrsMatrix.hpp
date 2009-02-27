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
#include "Tpetra_Operator.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"

// Assumption: sizeof(GlobalOrdinal) >= sizeof(LocalOrdinal)
// Assumption: max(GlobalOrdinal) >= max(Teuchos_Ordinal) >= max(LocalOrdinal)
// TODO: Add CompileTime tests of these assumptions, add these assumptions to documentation

// TODO: Add static graph construction option to constructors, a là Epetra

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
   * Method submitEntries() can be used to set both locally
   * owned and non-local elements; the shipping of data is done with hardcoded
   * MPI calls when fillComplete() is called.
   *
   * The nonzero elements of  locally owned row can be accessed by method
   * getMyRowCopy() or getGlobalRowCopy(). The former returns the column
   * indices using local numbering, the latter using global numbering.
   *
   */
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal=LocalOrdinal>
  class CrsMatrix : public Operator<Scalar,LocalOrdinal,GlobalOrdinal>
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
      //! @name Query Methods
      //@{ 
      
      //! Returns \c true if fillComplete() has been called.
      inline bool isFillComplete() const;

      //! Returns \c true if optimizeStorage() has been called.
      inline bool isStorageOptimized() const;

      //! \brief If matrix indices have been transformed to local, this function returns true. Otherwise, it returns false.
      /*! Indices are transformed to local after calling fillComplete(). Furthermore, indices are local if the 
          matrix was constructed with a completed CrsGraph. */
      inline bool indicesAreLocal() const;

      //! \brief If matrix indices have been transformed to local, this function returns false. Otherwise, it returns true.
      /*! Indices are global after calling fillComplee() or if the matrix was constructed with a completed CrsGraph. */
      inline bool indicesAreGlobal() const;

      //! Returns the communicator.
      inline Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

      //! Returns the number of global matrix rows. 
      inline GlobalOrdinal getNumGlobalRows() const;

      //! \brief Returns the number of global matrix columns. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      inline GlobalOrdinal getNumGlobalCols() const;

      //! Returns the number of matrix rows owned by the calling image. 
      inline Teuchos_Ordinal getNumLocalRows() const;

      //! Returns the index base for global indices for this matrix. 
      inline Teuchos_Ordinal getIndexBase() const;

      //! Returns the Map that describes the row distribution in this matrix.
      inline const Map<LocalOrdinal,GlobalOrdinal> & getRowMap() const;

      //! \brief Returns the Map that describes the column distribution in this matrix.
      /*! Cannot be called before fillComplete(), unless the matrix was constructed with a column map. */
      inline const Map<LocalOrdinal,GlobalOrdinal> & getColMap() const;

      //! \brief Indicates whether the matrix has a well-defined column map. 
      /*! The column map does not exist until after fillComplete(), unless the matrix was constructed with one. */
      inline bool hasColMap() const; 

      //! \brief Returns the number of matrix columns needed by the calling image to apply the forward operator.
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      inline Teuchos_Ordinal getNumLocalCols() const;

      //! \brief Returns the number of nonzero entries in the global matrix. 
      /*! Returns the number of global entries in the associated graph. */
      inline GlobalOrdinal getGlobalNNZ() const;

      //! \brief Returns the number of nonzero entries in the calling image's portion of the matrix. 
      /*! Before fillComplete() is called, this could include duplicated entries. */
      inline Teuchos_Ordinal getLocalNNZ() const;

      //! \brief Returns the current number of nonzero entries on this node in the specified global row .
      /*! Throws exception std::runtime_error if the specified global row does not belong to this node. */
      inline Teuchos_Ordinal getRowNNZGlobal(GlobalOrdinal globalRow) const;

      //! Returns the current number of nonzero entries on this node in the specified local row.
      /*! Throws exception std::runtime_error if the specified local row is not valid for this node. */
      inline Teuchos_Ordinal getRowNNZLocal(LocalOrdinal localRow) const;

      //! \brief Returns the number of global nonzero diagonal entries, based on global row/column index comparisons. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      inline GlobalOrdinal getNumGlobalNZDiags() const;

      //! \brief Returns the number of local nonzero diagonal entries, based on global row/column index comparisons. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      inline Teuchos_Ordinal getNumLocalNZDiags() const;

      //! \brief Returns the maximum number of nonzero entries across all rows/columns on all images. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      inline GlobalOrdinal getGlobalMaxRowNNZ() const;

      //! \brief Returns the maximum number of nonzero entries across all rows/columns on this image. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      inline Teuchos_Ordinal getLocalMaxRowNNZ() const;

      //@}

      //! @name Methods implementing Tpetra::Operator
      //@{ 

      //! Returns the Map associated with the domain of this operator.
      const Map<LocalOrdinal,GlobalOrdinal> & getDomainMap() const;

      //! Returns the Map associated with the domain of this operator.
      const Map<LocalOrdinal,GlobalOrdinal> & getRangeMap() const;

      //! Computes the matrix-vector multilication y = A x.
      void apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>& X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

      //@}

      //! @name Data Entry Methods
      //@{ 

      //! Submit multiple entries, using global IDs.
      /*! All index values must be in the global space. */
      void submitEntries(GlobalOrdinal globalRow, 
                         const Teuchos::ArrayView<const GlobalOrdinal> &cols,
                         const Teuchos::ArrayView<const Scalar>        &vals);

      //! \brief Submit a single entry, using global IDs.
      /*! The index values must be inthe global space. When entering multiple entries in a single row, 
          this method is less efficient than submitEntries(). */
      void submitEntry(GlobalOrdinal globalRow, GlobalOrdinal globalCol, Scalar value);

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
      Teuchos::RCP<Vector<Scalar,LocalOrdinal,GlobalOrdinal> > getLocalDiagCopy() const;

      //! Returns a copy of the specified (and locally owned) row, using local indices.
      /*! Before fillComplete(), the results will not include entries submitted to another node and may contain duplicated entries.
       * \pre haveColMap() == true
       */
      void getLocalRowCopy(LocalOrdinal localRow, 
                           const Teuchos::ArrayView<LocalOrdinal> &indices, 
                           const Teuchos::ArrayView<Scalar> &values) const;

      //! Returns a copy of the specified (and locally owned) row, using global indices.
      /*! Before fillComplete(), the results will not include entries submitted to another node and may contain duplicated entries. */
      void getGlobalRowCopy(GlobalOrdinal globalRow, 
                            const Teuchos::ArrayView<GlobalOrdinal> &indices,
                            const Teuchos::ArrayView<Scalar> &values) const;

      //@}

      //! @name I/O Methods
      //@{ 
      
      //! Prints the matrix on the specified stream. This is very verbose.
      void print(std::ostream& os) const;

      // @}

    private:
      // copy constructor disabled
      CrsMatrix(const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal> &Source);
      // operator= disabled
      CrsMatrix& operator=(const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal> &rhs);
      // useful typedefs
      typedef Teuchos::OrdinalTraits<LocalOrdinal>    LOT;
      typedef Teuchos::OrdinalTraits<GlobalOrdinal>   GOT;
      typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> MV;

      // Performs importing of off-processor elements and adds them to the locally owned elements.
      void globalAssemble();
      // multiplication routines
      void GeneralMV (typename MV::const_pointer x       , typename MV::pointer y       ) const;
      void GeneralMM (typename MV::const_double_pointer X, typename MV::double_pointer Y, Teuchos_Ordinal numVectors) const;
      void GeneralMhV(typename MV::const_pointer x       , typename MV::pointer y       ) const;
      void GeneralMhM(typename MV::const_double_pointer X, typename MV::double_pointer Y, Teuchos_Ordinal numVectors) const;

      // Graph Info: to be moved to CrsGraph
      Teuchos::RCP<const Teuchos::Comm<int> > comm_;
      Map<LocalOrdinal,GlobalOrdinal> rowMap_, colMap_;
      Map<LocalOrdinal,GlobalOrdinal> rangeMap_, domainMap_;
      GlobalOrdinal     globalNNZ_, numGlobalDiags_, globalMaxNumEntries_;
      Teuchos_Ordinal    localNNZ_,  numLocalDiags_,  localMaxNumEntries_;
      // structures used before optimizeStorage()
      Teuchos::Array<Teuchos::ArrayRCP<GlobalOrdinal> > colGInds_;                      // empty after fillComplete()
      Teuchos::Array<Teuchos::ArrayRCP<LocalOrdinal> >  colLInds_;                      // empty before fillComplete(), after optimizeStorage()
      Teuchos::Array<LocalOrdinal> rowNNZ_, rowAlloc_;                                  // empty after OptimizeStorage()
      // structure used after optimizeStorage()
      Teuchos::ArrayRCP<LocalOrdinal> contigColInds_;                                   // allocated during optimizeStorage()
      /* colIndsPtrs_[j] is always the begin() iterator from an ArrayView of the source or contigColInds_
         of the appropriate length. in a debug build, it is an ArrayRCP, which does bounds checking. in an optimized
         build, it is a C pointer. colIndsPtrs_ is allocated to numLocalRows()+1; the span of the jth row begins with
         colIndsPtrs_[j] and ends before colIndsPtrs_[j+1] */
      Teuchos::Array<typename Teuchos::ArrayRCP<LocalOrdinal>::iterator> colIndsPtrs_;
      // these are RCPs because they are optional
      // importer is needed if DomainMap is not sameas ColumnMap
      Teuchos::RCP<Import<LocalOrdinal,GlobalOrdinal> > importer_;
      // exporter is needed if RowMap is not sameas RangeMap
      Teuchos::RCP<Export<LocalOrdinal,GlobalOrdinal> > exporter_;

      // values
      // before optimizeStorage()
      Teuchos::Array<Teuchos::ArrayRCP<Scalar> > values_;
      // after optimizeStorage()
      Teuchos::ArrayRCP<Scalar> contigValues_;
      /* valuesPtrs_[j] is always the begin() iterator from an ArrayView of the source or contigColInds_
         of the appropriate length. in a debug build, it is an ArrayRCP, which does bounds checking. in an optimized
         build, it is a C pointer. valuesPtrs_ is allocated to numLocalRows()+1; the span of the jth row begins with
         valuesPtrs_[j] and ends before valuesPtrs_[j+1] */
      Teuchos::Array<typename Teuchos::ArrayRCP<Scalar>::iterator> valuesPtrs_;
      // multivectors used for import/export dest/source in apply()
      mutable Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> > importMV_, exportMV_;

      // a map between a (non-local) row and a list of (col,val)
      // TODO: this functionality will be divided between CrsGraph and CrsMatrix after the former comes into existence
      std::map<GlobalOrdinal, std::list<std::pair<GlobalOrdinal,Scalar> > > nonlocals_;

      bool fillComplete_;

  }; // class CrsMatrix



  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::CrsMatrix(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, Teuchos_Ordinal maxNNZPerRow)
  // CrsGraph members
  : comm_(rowMap.getComm())
  , rowMap_(rowMap)
  , colMap_(rowMap)     //
  , rangeMap_(rowMap)   //  these have to be set to something; we'll set them appropriately later
  , domainMap_(rowMap)  // 
  , globalNNZ_(GOT::invalid())
  , numGlobalDiags_(GOT::zero())
  , globalMaxNumEntries_(GOT::zero())
  , localNNZ_(0)
  , numLocalDiags_(0)
  , localMaxNumEntries_(0)
  , colGInds_(rowMap.getNumMyEntries())
  , colLInds_(0)
  , rowNNZ_(rowMap.getNumMyEntries(),0)
  , rowAlloc_(rowMap.getNumMyEntries(),maxNNZPerRow)
  , contigColInds_(Teuchos::null)
  , colIndsPtrs_(0)
  , importer_(Teuchos::null)
  , exporter_(Teuchos::null)
  // CrsMatrix members
  , values_(rowMap.getNumMyEntries())
  , contigValues_(Teuchos::null)
  , valuesPtrs_(0)
  , importMV_(Teuchos::null)
  , exportMV_(Teuchos::null)
  , fillComplete_(false)
  {
    if (maxNNZPerRow > 0) {
      for (Teuchos_Ordinal i=0; i<rowMap.getNumMyEntries(); ++i) {
        colGInds_[i] = Teuchos::arcp<GlobalOrdinal>(rowAlloc_[i]);
        values_[i]   = Teuchos::arcp<Scalar       >(rowAlloc_[i]);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::CrsMatrix(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, const Teuchos::ArrayView<Teuchos_Ordinal> &NNZPerRowToAlloc)
  // CrsGraph members
  : comm_(rowMap.getComm())
  , rowMap_(rowMap)
  , colMap_(rowMap)     //
  , rangeMap_(rowMap)   //  these have to be set to something; we'll set them appropriately later
  , domainMap_(rowMap)  // 
  , globalNNZ_(GOT::invalid())
  , numGlobalDiags_(GOT::zero())
  , globalMaxNumEntries_(GOT::zero())
  , localNNZ_(0)
  , numLocalDiags_(0)
  , localMaxNumEntries_(0)
  , colGInds_(rowMap.getNumMyEntries())
  , colLInds_(0)
  , rowNNZ_(rowMap.getNumMyEntries(),0)
  , rowAlloc_(rowMap.getNumMyEntries(),0)
  , contigColInds_(Teuchos::null)
  , colIndsPtrs_(0)
  , importer_(Teuchos::null)
  , exporter_(Teuchos::null)
  // CrsMatrix members
  , values_(rowMap.getNumMyEntries())
  , contigValues_(Teuchos::null)
  , valuesPtrs_(0)
  , importMV_(Teuchos::null)
  , exportMV_(Teuchos::null)
  , fillComplete_(false)
  {
    TEST_FOR_EXCEPTION(NNZPerRowToAlloc.size() != rowMap.getNumMyEntries(), std::runtime_error,
        Teuchos::typeName(*this) << "::CrsMatrix(rowMap,NNZPerRowToAlloc): NNZPerRowToAlloc must have as many entries as specified by rowMap for this node.");
    typename Teuchos::ArrayView<Teuchos_Ordinal>::iterator nnz    = NNZPerRowToAlloc.begin(),
                                                           nnzend = NNZPerRowToAlloc.end();
    typename Teuchos::Array<Teuchos::ArrayRCP<LocalOrdinal> >::iterator rinds  = colGInds_.begin();
    typename Teuchos::Array<Teuchos::ArrayRCP<Scalar> >::iterator       rvals  = values_.begin();
    typename Teuchos::Array<Teuchos_Ordinal>::iterator                  ralloc = rowAlloc_.begin();
    for (; nnz != nnzend; ++nnz, ++rinds, ++rvals, ++ralloc) {
      if (*nnz > 0) {
        (*rinds) = Teuchos::arcp<GlobalOrdinal>(*nnz);
        (*rvals) = Teuchos::arcp<Scalar       >(*nnz);
        (*ralloc) = *nnz;
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
  { return !isFillComplete(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::indicesAreGlobal() const
  { return isFillComplete(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::hasColMap() const
  { return isFillComplete(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP<const Teuchos::Comm<int> > 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getComm() const
  { return comm_; }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getGlobalNNZ() const
  {
    return globalNNZ_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getLocalNNZ() const
  {
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
      Teuchos::typeName(*this) << ":CrsMatrix: cannot call getLocalNNZ() until fillComplete() has been called.");
    return localNNZ_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getNumGlobalRows() const
  { return rowMap_.getNumGlobalEntries(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getNumGlobalCols() const
  {
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
      Teuchos::typeName(*this) << ": cannot call getNumGlobalCols() until fillComplete() has been called.");
    return colMap_.getNumGlobalEntries();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getNumLocalRows() const
  { return rowMap_.getNumMyEntries(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getNumLocalCols() const
  {
    TEST_FOR_EXCEPTION(hasColMap() == false, std::runtime_error,
      Teuchos::typeName(*this) << ": cannot call getNumMyCols() without a valid column map, typically after fillComplete() has been called.");
    return colMap_.getNumMyEntries();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getNumGlobalNZDiags() const
  {
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
      Teuchos::typeName(*this) << ": cannot call getNumGlobalNZDiags() until fillComplete() has been called.");
    return numGlobalDiags_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getNumLocalNZDiags() const
  {
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
      Teuchos::typeName(*this) << ": cannot call getNumLocalNZDiags() until fillComplete() has been called.");
    return numLocalDiags_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getRowNNZGlobal(GlobalOrdinal globalRow) const
  {
    LocalOrdinal myRow = rowMap_.getLocalIndex(globalRow);
    TEST_FOR_EXCEPTION(myRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::getRowNNZGlobal(globalRow): globalRow does not below to this node.");
    return (!isStorageOptimized() ? rowNNZ_[myRow] : colIndsPtrs_[myRow+1]-colIndsPtrs_[myRow]);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getRowNNZLocal(LocalOrdinal localRow) const
  {
    TEST_FOR_EXCEPTION(!rowMap_.isMyLocalIndex(localRow), std::runtime_error,
        Teuchos::typeName(*this) << "::getRowNNZLocal(localRow): localRow not valid for this node.");
    return (!isStorageOptimized() ? rowNNZ_[localRow] : colIndsPtrs_[localRow+1]-colIndsPtrs_[localRow]);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getGlobalMaxRowNNZ() const
  {
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
      Teuchos::typeName(*this) << ": cannot call getGlobalMaxRowNNZ() until fillComplete() has been called.");
    return globalMaxNumEntries_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getLocalMaxRowNNZ() const
  {
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
      Teuchos::typeName(*this) << ": cannot call getLocalMaxRowNNZ() until fillComplete() has been called.");
    return localMaxNumEntries_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getIndexBase() const
  { return rowMap_.getIndexBase(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  const Map<LocalOrdinal,GlobalOrdinal> & 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getRowMap() const 
  { return rowMap_; }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  const Map<LocalOrdinal,GlobalOrdinal> & 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getColMap() const 
  { 
    TEST_FOR_EXCEPTION(hasColMap() == false, std::runtime_error,
        Teuchos::typeName(*this) << ": matrix does not currently have a valid column map.");
    return colMap_; 
  }
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  const Map<LocalOrdinal,GlobalOrdinal> & 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getDomainMap() const
  {
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getDomainMap(): matrix does not currently have range or domain maps.");
    return domainMap_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  const Map<LocalOrdinal,GlobalOrdinal> & 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getRangeMap() const
  {
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getRangeMap(): matrix does not currently have range or domain maps.");
    return rangeMap_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::submitEntries(GlobalOrdinal globalRow, 
                         const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                         const Teuchos::ArrayView<const Scalar>  &values)
  {
    /* this method adds a new entry to the graph, which requires that the graph was not given and the fillComplete() has not been called */
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(isFillComplete() == true || isStorageOptimized() == true, std::runtime_error,
      Teuchos::typeName(*this) << "::submitEntries(): fillComplete() has already been called.");
    /* this version can always be called:
       a) either indices are still global, or
       b) indices are local, in which case we have a column map (either from construction or from 
          our graph), so that we can translate the given global indices */
    TEST_FOR_EXCEPTION(values.size() != indices.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::submitEntries(): values.size() must equal indices.size().");
    typename Teuchos::ArrayView<const GlobalOrdinal>::iterator ind = indices.begin();
    typename Teuchos::ArrayView<const Scalar       >::iterator val =  values.begin();
    if (rowMap_.isMyGlobalIndex(globalRow)) {
      Teuchos_Ordinal myRow = rowMap_.getLocalIndex(globalRow),
                     rowNNZ = rowNNZ_[myRow], 
                      toAdd = indices.size();
      if (rowNNZ+toAdd > rowAlloc_[myRow]) {
#       if defined(THROW_TPETRA_EFFICIENCY_WARNINGS) || defined(PRINT_TPETRA_EFFICIENCY_WARNINGS)
          std::string err = Teuchos::typeName(*this) + "::submitEtnry(): Pre-allocated space has been exceeded, requiring new allocation. To improve efficiency, suggest larger allocation.";
#         if defined(THROW_TPETRA_EFFICIENCY_WARNINGS)
            TEST_FOR_EXCEPTION(true, std::runtime_error, err);
#         else
            std::cerr << err << std::endl;
#         endif
#       endif
        // increase allocation to necessary amount, copy previous data, reassign ArrayRCPs
        rowAlloc_[myRow] = rowNNZ+toAdd;
        ArrayRCP<GlobalOrdinal> newInds = Teuchos::arcp<GlobalOrdinal>(Teuchos::as<Teuchos_Ordinal>(rowAlloc_[myRow]));
        ArrayRCP<Scalar>        newVals = Teuchos::arcp<Scalar       >(Teuchos::as<Teuchos_Ordinal>(rowAlloc_[myRow]));
        std::copy(colGInds_[myRow].begin(),colGInds_[myRow].begin()+rowNNZ,newInds.begin());
        std::copy(  values_[myRow].begin(),  values_[myRow].begin()+rowNNZ,newVals.begin());
        colGInds_[myRow] = newInds;
        values_[myRow] = newVals;
      }
      std::copy(indices.begin(),indices.end(),colGInds_[myRow].begin()+rowNNZ);
      std::copy( values.begin(), values.end(),  values_[myRow].begin()+rowNNZ);
      rowNNZ_[myRow] += toAdd;
    }
    else {
      for (; val != values.end(); ++val, ++ind) {
        nonlocals_[globalRow].push_back(std::make_pair(*ind, *val));
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::submitEntry(GlobalOrdinal globalRow, GlobalOrdinal globalCol, Scalar value)
  { submitEntries(globalRow,Teuchos::arrayView(&globalCol,1),Teuchos::arrayView(&value,1)); }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::setAllToScalar(const Scalar &alpha)
  { 
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    if (isStorageOptimized() == true) {
      // contiguous storage: the easy (and fast) way
      std::fill(contigValues_.begin(), contigValues_.end(), alpha);
    }
    else {
      // non-contiguous storage: the hard (and slow) way
      typename Array<ArrayRCP<Scalar> >::const_iterator rowvals;
      typename Array<LocalOrdinal>::const_iterator rownnz;
      for (rowvals = values_.begin(), rownnz = rowNNZ_.begin(); rowvals != values_.end(); ++rowvals, ++rownnz) 
      {
        std::fill(rowvals->begin(),rowvals->begin()+(*rownnz),alpha);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::scale(const Scalar &alpha)
  { 
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    if (isStorageOptimized() == true) {
      // contiguous storage: the easy (and fast) way
      for (typename Teuchos::ArrayRCP<Scalar>::iterator val  = contigValues_.begin(); 
                                                        val != contigValues_.end(); ++val) 
      {
        (*val) *= alpha;
      }
    }
    else {
      // non-contiguous storage: the hard (and slow) way
      typename Array<ArrayRCP<Scalar> >::const_iterator rowvals;
      typename Array<Teuchos_Ordinal>::const_iterator rownnz;
      for (rowvals = values_.begin(), rownnz = rowNNZ_.begin(); rowvals != values_.end(); ++rowvals, ++rownnz) 
      {
        for (typename Teuchos::ArrayRCP<Scalar>::iterator val = rowvals->begin(); val != rowvals->begin() + (*rownnz); ++val)
        {
          (*val) *= alpha;
        }
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getLocalRowCopy(LocalOrdinal myRow, 
                                                                     const Teuchos::ArrayView<LocalOrdinal> &indices, 
                                                                     const Teuchos::ArrayView<Scalar>       &values) const
  {
    TEST_FOR_EXCEPTION(hasColMap() == false, std::runtime_error,
        Teuchos::typeName(*this) << ": cannot call getLocalRowCopy() without a valid column map.");
    Teuchos_Ordinal rownnz = (!isStorageOptimized() ? rowNNZ_[myRow] : colIndsPtrs_[myRow+1]-colIndsPtrs_[myRow]);
    TEST_FOR_EXCEPTION(indices.size() != rownnz || values.size() != rownnz, std::runtime_error, 
        Teuchos::typeName(*this) << "::getLocalRowCopy(indices,values): size of indices,values must be sufficient to store the values for the specified row.");
    if (indicesAreLocal()) {
      if (isStorageOptimized()) {
        std::copy(colIndsPtrs_[myRow],colIndsPtrs_[myRow+1],indices.begin());
        std::copy( valuesPtrs_[myRow], valuesPtrs_[myRow+1], values.begin());
      }
      else {
        std::copy(colLInds_[myRow].begin(),colLInds_[myRow].begin()+rownnz,indices.begin());
        std::copy(  values_[myRow].begin(),  values_[myRow].begin()+rownnz, values.begin());
      }
    }
    else {
      TEST_FOR_EXCEPTION(isStorageOptimized() == true, std::logic_error,
          Teuchos::typeName(*this) << "::getLocalRowCopy(): internal logic error. Contact Tpetra team.");
      typename Teuchos::ArrayView<LocalOrdinal> il = indices.begin();
      typename Teuchos::ArrayRCP<LocalOrdinal>  gi = colGInds_[myRow].begin(), gend = gi+rownnz;
      for (; gi != gend; ++il, ++gi) {
        (*il) = colMap_.getLocalIndex(*gi);
      }
      std::copy(  values_[myRow].begin(),  values_[myRow].begin()+rownnz, values.begin());
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getGlobalRowCopy(GlobalOrdinal globalRow, 
                                                                      const Teuchos::ArrayView<GlobalOrdinal> &indices,
                                                                      const Teuchos::ArrayView<Scalar>  &values) const
  {
    // Only locally owned rows can be queried, otherwise complain
    Teuchos_Ordinal myRow = rowMap_.getLocalIndex(globalRow);
    TEST_FOR_EXCEPTION(!rowMap_.isMyGlobalIndex(globalRow), std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowCOpy(globalRow): globalRow does not belong to this node.");
    Teuchos_Ordinal rownnz = Teuchos::as<Teuchos_Ordinal>(rowNNZ_[myRow]);
    TEST_FOR_EXCEPTION(
        indices.size() != rownnz || values.size() != rownnz, std::runtime_error, 
        Teuchos::typeName(*this) << "::getGlobalRowCopy(indices,values): size of indices,values must be sufficient to store the values for the specified row.");
    std::copy(values_[myRow].begin(),values_[myRow].begin()+rownnz,values.begin());
    if (isFillComplete())
    {
      // translate local IDs back to global IDs
      typename Teuchos::ArrayView<GlobalOrdinal>::iterator  gind = indices.begin();
      typename Teuchos::ArrayRCP<LocalOrdinal>::const_iterator lind = colLInds_[myRow].begin(),
                                                               lend = colLInds_[myRow].begin() + rownnz;
      for (; lind != lend; ++lind, ++gind) 
      {
        (*gind) = colMap_.getGlobalIndex(*lind);
      }
    }
    else
    {
      std::copy(colGInds_[myRow].begin(),colGInds_[myRow].begin()+rownnz,indices.begin());
    }
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
    int myImageID = Teuchos::rank(*comm_);
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
    Teuchos::RCP<const MV> Xcopy;
    typename MV::const_double_pointer Xdata = X.extractConstView2D();
    typename MV::double_pointer       Ydata = Y.extractView2D();
    if (&X==&Y && importer_==null && exporter_==null) {
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
    if (importer_ != null) {
      if (importMV_ != null && importMV_->numVectors() != numVectors) importMV_ = null;
      if (importMV_ == null) {
        importMV_ = Teuchos::rcp( new MV(colMap_,numVectors) );
      }
    }
    if (exporter_ != null) {
      if (exportMV_ != null && exportMV_->numVectors() != numVectors) exportMV_ = null;
      if (exportMV_ == null) {
        exportMV_ = Teuchos::rcp( new MV(rowMap_,numVectors) );
      }
    }
    if (mode == Teuchos::NO_TRANS) {
      // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
      if (importer_ != null) {
        importMV_->doImport(X, *importer_, INSERT);
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
      if (exporter_ != null) {
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
      if (exporter_ != null) {
        Y.putScalar(ST::zero());  // Make sure target is zero: necessary because we are adding. may need adjusting for alpha,beta apply()
        Y.doExport(*exportMV_, *exporter_, ADD); // Fill Y with Values from export vector
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
      if (exporter_ != null) {
        exportMV_->doImport(X,*exporter_,INSERT);
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
      if (importer_ != null) {
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
      if (importer_ != null) {
        Y.putScalar(ST::zero()); // Make sure target is zero: necessary because we are adding. may need adjusting for alpha,beta apply()
        Y.doExport(*importMV_,*importer_,ADD);
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
    int myImageID = Teuchos::rank(*comm_);
    if (myImageID == 0)
    {
      os << "Tpetra::CrsMatrix, label = " << this->getObjectLabel() << endl;
      os << "Number of global rows    = " << getNumGlobalRows() << endl;
      if (isFillComplete())
      {
        os << "Number of global columns    = " << getNumGlobalCols() << endl;
        os << "Status = fill complete" << endl;
        os << "Number of global nonzeros   = " << getGlobalNNZ() << endl;
        os << "Global max nonzeros per row = " << getGlobalMaxRowNNZ() << endl;
      }
      else
      {
        os << "Status = fill not complete" << endl;
      }
    }
    if (isFillComplete())
    {
      for (int pid=0; pid < Teuchos::size(*comm_); ++pid)
      {
        if (pid == myImageID)
        {
          Teuchos::Array<GlobalOrdinal> indices(getLocalMaxRowNNZ());
          Teuchos::Array<Scalar>         values(getLocalMaxRowNNZ());
          os << "% Number of rows on image " << myImageID << " = " << getNumLocalRows() << endl;
          for (Teuchos_Ordinal i=0; i < getNumLocalRows(); ++i)
          {
            GlobalOrdinal globalRow = rowMap_.getGlobalIndex(i);
            Teuchos_Ordinal rowSize = getNumRowEntries(globalRow);
            if (rowSize > Teuchos::OrdinalTraits<Teuchos_Ordinal>::zero()) {
              getGlobalRowCopy(globalRow, indices(0,rowSize), values(0,rowSize));
              for (Teuchos_Ordinal j=0; j < rowSize; ++j) {
                os << "Matrix(" << globalRow << ", " << indices[j] << ") = " << values[j] << endl;
              }
            }
          }
        }
        Teuchos::barrier(*comm_);
        Teuchos::barrier(*comm_);
        Teuchos::barrier(*comm_);
      }
    }
  }



  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::fillComplete(const Map<LocalOrdinal,GlobalOrdinal> &domainMap, const Map<LocalOrdinal,GlobalOrdinal> &rangeMap, bool OptimizeStorage)
  {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(isFillComplete() == true, std::runtime_error,
      Teuchos::typeName(*this) << "::fillComplete(): fillComplete() has already been called.");

    typename Teuchos::ArrayView<const GlobalOrdinal> myGlobalEntries = rowMap_.getMyGlobalEntries();
    domainMap_ = domainMap;
    rangeMap_  = rangeMap;
    Teuchos_Ordinal numLocalRows = rowMap_.getNumMyEntries();

    // =============================== //
    // Part 0: send off-image elements //
    // =============================== //
    if (comm_->getSize() > 1) {
      globalAssemble();
    }
    else {
      TEST_FOR_EXCEPTION(nonlocals_.size() > 0, std::runtime_error,
          Teuchos::typeName(*this) << "::fillComplete(): cannot have non-local entries on a serial run. Invalid entry was submitted to the CrsMatrix.");
    }

    // =============================== //
    // Part I: remove repeated indices //
    // =============================== //
    //
    // Load all matrix entries in a hash table, then re-fill
    // the row with the last inserted value.
    // This also has the effect of sorting the indices.
    for (Teuchos_Ordinal r = 0; r < numLocalRows; ++r)
    {
      // put them in the map from colGInds_,values_
      std::map<GlobalOrdinal,Scalar> row;
      {
        typename ArrayRCP<GlobalOrdinal>::const_iterator cind = colGInds_[r].begin(),
                                                         cend = colGInds_[r].begin() + rowNNZ_[r];
        typename ArrayRCP<Scalar>::const_iterator         val = values_[r].begin();
        for (; cind != cend; ++cind, ++val)
        {
          typename std::map<GlobalOrdinal,Scalar>::iterator loc = row.find(*cind);
          if (loc == row.end()) {
            // insert as is; unfortunately, this is necessary for Scalar who don't initialize to zero (e.g., dd_real)
            row[*cind] = *val;
          }
          else {
            (*loc).second += *val;
          }
        }
      }
      // get them out of the map, back to colGInds_,values_
      Teuchos_Ordinal count = 0;
      typename std::map<GlobalOrdinal,Scalar>::const_iterator miter = row.begin(),
                                                               mend = row.end();
      typename Teuchos::ArrayRCP<GlobalOrdinal>::const_iterator citer = colGInds_[r].begin();
      typename Teuchos::ArrayRCP<Scalar       >::const_iterator viter = values_[r].begin();
      for (; miter != mend; ++miter, ++citer, ++viter)
      {
        (*citer) = miter->first;
        (*viter) = miter->second;
        ++count;
      }
      rowNNZ_[r] = count;
    }

    // =============================== //
    // Part II: build the column space //
    // =============================== //
    // construct a list of columns for which we have a non-zero
    TEST_FOR_EXCEPTION(hasColMap() == true,std::logic_error,"Incomplete feature.");
    {
      std::set<GlobalOrdinal> nnzcols;
      for (Teuchos_Ordinal r=0; r < numLocalRows; ++r)
      {
        typename ArrayRCP<GlobalOrdinal>::const_iterator cind = colGInds_[r].begin(),
                                                         cend = colGInds_[r].begin() + rowNNZ_[r];
        for (; cind != cend; ++cind) {
          nnzcols.insert(*cind);
        }
      }
      Array<GlobalOrdinal> myColumns(nnzcols.size());
      std::copy(nnzcols.begin(),nnzcols.end(),myColumns.begin());
      colMap_ = Map<LocalOrdinal,GlobalOrdinal>(GOT::invalid(), myColumns(), rowMap_.getIndexBase(), comm_);
    }

    // create import, export
    if (!domainMap_.isSameAs(colMap_)) {
      importer_ = Teuchos::rcp( new Import<LocalOrdinal,GlobalOrdinal>(domainMap_,colMap_) );
    }
    else {
      importer_ = Teuchos::null;
    }
    if (!rangeMap_.isSameAs(rowMap_)) {
      exporter_ = Teuchos::rcp( new Export<LocalOrdinal,GlobalOrdinal>(rowMap_,rangeMap_) );
    }
    else {
      exporter_ = Teuchos::null;
    }

    // ============================== //
    // Part IV: move to local indices //
    // Reuse the space allocated for  //
    // the global indicies.           //
    // ============================== //
    colLInds_.resize(numLocalRows);
    for (Teuchos_Ordinal r=0; r < numLocalRows; ++r)
    {
      colLInds_[r] = Teuchos::arcp_reinterpret_cast<LocalOrdinal>(colGInds_[r]);
      typename ArrayRCP<GlobalOrdinal>::const_iterator cindG = colGInds_[r].begin(),
                                                       cendG = colGInds_[r].begin() + rowNNZ_[r];
      typename ArrayRCP<LocalOrdinal>::iterator        cindL = colLInds_[r].begin();
      for (; cindG != cendG; ++cindG, ++cindL)
      {
        (*cindL) = colMap_.getLocalIndex(*cindG);
      }
      colGInds_[r] = Teuchos::null;   // data is owned by colLInds_[r] now
    }

    // ================================================ //
    // Part V: compute some local and global quantities //
    // ================================================ //
    for (Teuchos_Ordinal r=0; r < numLocalRows; ++r)
    {
      GlobalOrdinal rgid = myGlobalEntries[r];
      Teuchos_Ordinal numEntries = rowNNZ_[r];
      localNNZ_ += numEntries;
      if ( colMap_.isMyGlobalIndex(rgid) ) {
        if ( std::find(colLInds_[r].begin(),colLInds_[r].begin()+numEntries,colMap_.getLocalIndex(rgid)) != colLInds_[r].begin()+numEntries ) {
          ++numLocalDiags_;
        }
      }
      localMaxNumEntries_ = std::max(localMaxNumEntries_,numEntries);
    }
    {
      GlobalOrdinal lcl[2], gbl[2];
      lcl[0] = localNNZ_; lcl[1] = numLocalDiags_;
      Teuchos::reduceAll<int,GlobalOrdinal>(*comm_,Teuchos::REDUCE_SUM,2,lcl,gbl);
      globalNNZ_ = gbl[0]; numGlobalDiags_ = gbl[1];
    }
    Teuchos::reduceAll<int,GlobalOrdinal>(*comm_,Teuchos::REDUCE_MAX,localMaxNumEntries_,&globalMaxNumEntries_);

    // mark transformation as successfully completed
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
    int numImages = comm_->getSize();
    int myImageID = comm_->getRank();
    // Determine if any nodes have global entries to share
    Teuchos_Ordinal MyNonlocals = nonlocals_.size(), 
                    MaxGlobalNonlocals;
    Teuchos::reduceAll<Teuchos_Ordinal>(*comm_,Teuchos::REDUCE_MAX,MyNonlocals,&MaxGlobalNonlocals);
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
        bool invalidGIDs = rowMap_.getRemoteIndexList(NLRs(),NLRIds());
        char lclerror = ( invalidGIDs ? 1 : 0 );
        char gblerror;
        Teuchos::reduceAll(*comm_,Teuchos::REDUCE_MAX,lclerror,&gblerror);
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
      Teuchos::gatherAll(*comm_,numImages,localNeighbors.getRawPtr(),numImages*numImages,globalNeighbors.values());
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
      sendRequests.push_back( Teuchos::isend<int,int>(*comm_,rcp(&sendSizes[s],false),sendIDs[s]) );
    }
    // perform non-blocking receives: receive sizes from our senders
    Array<Teuchos::RCP<Teuchos::CommRequest> > recvRequests;
    Array<Teuchos_Ordinal> recvSizes(numRecvs);
    for (int r=0; r < numRecvs; ++r) {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      recvRequests.push_back( Teuchos::ireceive(*comm_,rcp(&recvSizes[r],false),recvIDs[r]) );
    }
    // wait on all 
    if (!sendRequests.empty()) {
      Teuchos::waitAll(*comm_,sendRequests());
    }
    if (!recvRequests.empty()) {
      Teuchos::waitAll(*comm_,recvRequests());
    }
    Teuchos::barrier(*comm_);
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
      sendRequests.push_back( Teuchos::isend<int,CrsIJV<GlobalOrdinal,Scalar> >(*comm_,tmparcp,sendIDs[s]) );
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
      recvRequests.push_back( Teuchos::ireceive(*comm_,tmparcp,recvIDs[r]) );
    }
    // perform waits
    if (!sendRequests.empty()) {
      Teuchos::waitAll(*comm_,sendRequests());
    }
    if (!recvRequests.empty()) {
      Teuchos::waitAll(*comm_,recvRequests());
    }
    Teuchos::barrier(*comm_);
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
      submitEntries(ijv->i, Teuchos::arrayView(&ijv->j,1), Teuchos::arrayView(&ijv->v,1));
    }

    // WHEW! THAT WAS TIRING!
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::GeneralMV(typename MV::const_pointer x, typename MV::pointer y) const
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    typename Teuchos::ArrayRCP<LocalOrdinal>::const_iterator cind, cend;
    typename Teuchos::ArrayRCP<Scalar      >::const_iterator aval;
    for (Teuchos_Ordinal r=0; r < rowMap_.getNumMyEntries(); ++r) {
      cind = colLInds_[r].begin(); 
      cend = colLInds_[r].begin() + rowNNZ_[r];
      aval =  values_[r].begin(); 
      Scalar sum = ST::zero();
      for (; cind != cend; ++cind, ++aval) {
        sum += (*aval) * x[*cind];
      }
      y[r] = sum;
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::GeneralMM (typename MV::const_double_pointer X, typename MV::double_pointer Y, Teuchos_Ordinal numVectors) const
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    Teuchos_Ordinal numLocalRows = rowMap_.getNumMyEntries();
    typename Teuchos::ArrayRCP<LocalOrdinal>::const_iterator cind, cend;
    typename Teuchos::ArrayRCP<Scalar      >::const_iterator aval;
    for (Teuchos_Ordinal r=0; r < numLocalRows; ++r) {
      for (Teuchos_Ordinal j=0; j<numVectors; ++j) {
        typename MV::pointer       yvals = Y[j];
        typename MV::const_pointer xvals = X[j]; 
        cind = colLInds_[r].begin(); 
        cend = colLInds_[r].begin() + rowNNZ_[r];
        aval =  values_[r].begin(); 
        Scalar sum = ST::zero();
        for (; cind != cend; ++cind, ++aval) {
          sum += (*aval) * xvals[*cind];
        }
        yvals[r] = sum;
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::GeneralMhV(typename MV::const_pointer x, typename MV::pointer y) const
  { 
    typedef Teuchos::ScalarTraits<Scalar> ST;
    Teuchos_Ordinal numLocalRows = getNumLocalRows(),
                    numLocalCols = getNumLocalCols();
    // Initialize y for transpose multiply
    std::fill( y, y+numLocalCols, ST::zero() );
    // apply conjugate transpose of matrix to x
    // use column triad formulation, accumulating into y each column of A^H (each row of A) times each entry in x
    typename Teuchos::ArrayRCP<LocalOrdinal>::const_iterator cind, cend;
    typename Teuchos::ArrayRCP<Scalar      >::const_iterator aval;
    for (Teuchos_Ordinal r=0; r < numLocalRows; ++r) {
      // loop over entries in this column of A^H (this row of A)
      cind = colLInds_[r].begin();
      cend = colLInds_[r].begin() + rowNNZ_[r];
      aval =  values_[r].begin();
      for (; cind != cend; ++cind, ++aval) {
        y[*cind] += ST::conjugate(*aval) * x[r];
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::GeneralMhM(typename MV::const_double_pointer X, typename MV::double_pointer Y, Teuchos_Ordinal numVectors) const
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    Teuchos_Ordinal numLocalRows = getNumLocalRows(),
                    numLocalCols = getNumLocalCols();
    // Initialize Y for transpose multiply
    for (Teuchos_Ordinal j=0; j<numVectors; ++j) {
      std::fill( Y[j], Y[j]+numLocalCols, ST::zero() );
    }
    // apply conjugate transpose of matrix to X
    // use outer-product formulation, hitting Y with a rank-1 update comprised of each column of A^H (each row of A) and each row of X
    typename Teuchos::ArrayRCP<LocalOrdinal>::const_iterator cind, cend;
    typename Teuchos::ArrayRCP<Scalar      >::const_iterator aval;
    for (Teuchos_Ordinal r=0; r < numLocalRows; ++r) {
      // loop over numvectors
      for (Teuchos_Ordinal j=0; j<numVectors; ++j) {
        typename MV::pointer       yvals = Y[j];
        typename MV::const_pointer xvals = X[j]; 
        // loop over entries in this column of A^H (this row of A)
        cind = colLInds_[r].begin();
        cend = colLInds_[r].begin() + rowNNZ_[r];
        aval =  values_[r].begin();
        for (; cind != cend; ++cind, ++aval) {
          yvals[*cind] += ST::conjugate(*aval) * xvals[r];
        }
      }
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP<Vector<Scalar,LocalOrdinal,GlobalOrdinal> > 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getLocalDiagCopy() const
  {
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
        Teuchos::typeName(*this) << ": cannot call getLocalDiagCopy() until fillComplete() has been called.");
    Teuchos::RCP<Vector<Scalar,LocalOrdinal,GlobalOrdinal> > dvec = 
      Teuchos::rcp(new Vector<Scalar,LocalOrdinal,GlobalOrdinal>(rowMap_,false));
#ifdef TPETRA_DEBUG
    int numDiagFound = 0;
#endif
    Teuchos::ArrayView<Scalar> values;
    dvec->extractView1D(values);
    Teuchos_Ordinal numLocalRows = rowMap_.getNumMyEntries();
    typename Teuchos::ArrayRCP<Scalar>::const_iterator v;
    typename Teuchos::ArrayView<Scalar>::iterator ov;
    ov = values.begin();
    for (Teuchos_Ordinal r=0; r < numLocalRows; ++r) {
      *ov = Teuchos::ScalarTraits<Scalar>::zero();
      GlobalOrdinal rgid = rowMap_.getGlobalIndex(r);
      if (colMap_.isMyGlobalIndex(rgid)) {
        LocalOrdinal rlid = colMap_.getLocalIndex(rgid);
        typename Teuchos::ArrayRCP<LocalOrdinal>::const_iterator i = colLInds_[r].begin(),
                                                              iend = colLInds_[r].begin() + rowNNZ_[r];
        for (v = values_[r].begin(); i != iend; ++i, ++v) {
          if (*i == rlid) {
            *ov = *v;
#ifdef TPETRA_DEBUG
            ++numDiagFound;
#endif
            break;
          }
        }
      }
      ++ov;
    }
#ifdef TPETRA_DEBUG
    TEST_FOR_EXCEPTION(numDiagFound != numLocalDiags_, std::logic_error, "CrsMatrix::getLocalDiagCopy(): logic error. Please contact Tpetra team.");
#endif
    return dvec;
  }


} // namespace Tpetra

#endif
