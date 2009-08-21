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

#ifndef TPETRA_CRSMATRIX_HPP
#define TPETRA_CRSMATRIX_HPP

// TODO: add support for alpha,beta term coefficients: Y = alpha*A*X + beta*Y
// TODO: row-wise insertion of entries in globalAssemble()

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_getRawPtr.hpp>
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
  //! \brief A class for constructing and using sparse compressed index matrices with row access.
  /*!
   * This class allows the construction of sparse matrices with row-access. 
   * Method insertGlobalValues() can be used to set both locally
   * owned and non-local elements; the shipping of data is done with hardcoded
   * MPI calls when fillComplete() is called.
   *
   * The nonzero elements of  locally owned row can be accessed by method
   * getLocalRowCopy() or getGlobalRowCopy(). The former returns the column
   * indices using local numbering, the latter using global numbering.
   *
   * This class is templated on \c Scalar, \c LocalOrdinal and \c GlobalOrdinal. 
   * The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
   * type, if omitted, defaults to the \c LocalOrdinal type.
   * The class utilizes CrsGraph object which has the same local and global ordinal types.
   *
   */
  template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class CrsMatrix : public RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>
  {
    public:
      //! @name Constructor/Destructor Methods
      //@{ 

      //! Constructor specifying the number of non-zeros for all rows.
      CrsMatrix(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, size_t numNNZ, bool staticProfile = false);

      //! Constructor specifying the number of non-zeros for each row.
      CrsMatrix(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, const Teuchos::ArrayView<size_t> &NNZPerRowToAlloc, bool staticProfile = false);

      //! Constructor specifying a column map and the number of non-zeros for all rows.
      CrsMatrix(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, const Map<LocalOrdinal,GlobalOrdinal> &colMap, size_t numNNZ, bool staticProfile = false);

      //! Constructor specifying a column map and the number of non-zeros for each row.
      CrsMatrix(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, const Map<LocalOrdinal,GlobalOrdinal> &colMap, const Teuchos::ArrayView<size_t> &NNZPerRowToAlloc, bool staticProfile = false);

      //! Constructor specifying a pre-constructed graph.
      CrsMatrix(const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> &graph);

      // !Destructor.
      virtual ~CrsMatrix();

      //@}

      //! @name Methods implementing Operator
      //@{ 

      //! \brief Returns the Map associated with the domain of this operator.
      //! This will be equal to the row map until fillComplete() is called.
      const Map<LocalOrdinal,GlobalOrdinal> & getDomainMap() const;

      //! Returns the Map associated with the domain of this operator.
      //! This will be equal to the row map until fillComplete() is called.
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
      bool isLocallyIndexed() const;

      //! \brief If matrix indices are stored as global indices, this function returns false. Otherwise, it returns true.
      bool isGloballyIndexed() const;

      //! \brief Indicates whether the matrix is lower triangular.
      bool isLowerTriangular() const;

      //! \brief Indicates whether the matrix is upper triangular.
      bool isUpperTriangular() const;

      //! Returns the communicator.
      Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

      //! Returns the underlying node.
      virtual Node& getNode() const = 0;

      //! \brief Indicates whether the matrix has a well-defined column map. 
      /*! The column map does not exist until after fillComplete(), unless the matrix was constructed with one. */
      bool hasColMap() const; 

      //! Returns the RowGraph associated with this matrix. 
      const RowGraph<LocalOrdinal,GlobalOrdinal> &getGraph() const;

      //! Returns the Map that describes the row distribution in this matrix.
      const Map<LocalOrdinal,GlobalOrdinal> & getRowMap() const;

      //! \brief Returns the Map that describes the column distribution in this matrix.
      //! This will be equal to the row map if hasColMap() == false.
      const Map<LocalOrdinal,GlobalOrdinal> & getColMap() const;

      //! Returns the number of global matrix rows. 
      GlobalOrdinal getNumGlobalRows() const;

      //! \brief Returns the number of global matrix columns. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      GlobalOrdinal getNumGlobalCols() const;

      //! Returns the number of matrix rows owned by the calling image. 
      size_t getNumLocalRows() const;

      //! \brief Returns the number of matrix columns needed by the calling image to apply the forward operator.
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      size_t getNumLocalCols() const;

      //! Returns the index base for global indices for this matrix. 
      GlobalOrdinal getIndexBase() const;

      //! \brief Returns the number of nonzero entries in the global matrix. 
      /*! Returns the number of global entries in the associated graph. */
      global_size_t getNumGlobalEntries() const;

      //! \brief Returns the number of nonzero entries in the calling image's portion of the matrix. 
      /*! Before fillComplete() is called, this could include duplicated entries. */
      size_t getNumLocalEntries() const;

      //! \brief Returns the current number of nonzero entries on this node in the specified global row .
      /*! Throws exception std::runtime_error if the specified global row does not belong to this node. */
      size_t getNumEntriesForGlobalRow(GlobalOrdinal globalRow) const;

      //! Returns the current number of nonzero entries on this node in the specified local row.
      /*! Throws exception std::runtime_error if the specified local row is not valid for this node. */
      size_t getNumEntriesForLocalRow(LocalOrdinal localRow) const;

      //! \brief Returns the number of global nonzero diagonal entries, based on global row/column index comparisons. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      GlobalOrdinal getNumGlobalDiags() const;

      //! \brief Returns the number of local nonzero diagonal entries, based on global row/column index comparisons. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      size_t getNumLocalDiags() const;

      //! \brief Returns the maximum number of nonzero entries across all rows/columns on all images. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      GlobalOrdinal getGlobalMaxNumRowEntries() const;

      //! \brief Returns the maximum number of nonzero entries across all rows/columns on this image. 
      /*! May not be called before fillComplete(), unless the matrix was constructed with a column map. */
      size_t getLocalMaxNumRowEntries() const;

      //@}

      //! @name Data Entry Methods
      //@{ 

      //! Submit matrix entries, using local IDs.
      /*! All index values must be in the local space. */
      void insertMyValues(LocalOrdinal localRow, 
                         const Teuchos::ArrayView<const LocalOrdinal> &cols,
                         const Teuchos::ArrayView<const Scalar>       &vals);

      //! Submit matrix entries, using global IDs.
      /*! All index values must be in the global space. */
      void insertGlobalValues(GlobalOrdinal globalRow, 
                         const Teuchos::ArrayView<const GlobalOrdinal> &cols,
                         const Teuchos::ArrayView<const Scalar>        &vals);

      //! Replace matrix entries, using global IDs.
      /*! All index values must be in the global space. 
         If (globalRow,cols[i]) corresponds to an entry 
         that is duplicated in the matrix (likely because it was inserted more than once and fillComplete() 
         has not been called), the behavior of this function is not defined. */
      void replaceGlobalValues(GlobalOrdinal globalRow, 
                         const Teuchos::ArrayView<const GlobalOrdinal> &cols,
                         const Teuchos::ArrayView<const Scalar>        &vals);

      //! Sum into multiple entries, using global IDs.
      /*! All index values must be in the global space. */
      void sumIntoGlobalValues(GlobalOrdinal globalRow, 
                         const Teuchos::ArrayView<const GlobalOrdinal> &cols,
                         const Teuchos::ArrayView<const Scalar>        &vals);


      //! Set all matrix entries equal to scalarThis.
      void setAllToScalar(const Scalar &alpha);

      //! Scale the current values of a matrix, this = alpha*this. 
      void scale(const Scalar &alpha);

      //! \brief Communicate non-local contributions to other nodes.
      void globalAssemble();

      /*! \brief Signal that data entry is complete, specifying domain and range maps. 
          Off-node entries are distributed, repeated entries are summed, and global indices are transformed to local indices.
          If \c OptimizeStorage is true, then optimizeStorage() is called as well.
       */ 
      void fillComplete(const Map<LocalOrdinal,GlobalOrdinal> &domainMap, const Map<LocalOrdinal,GlobalOrdinal> &rangeMap, bool OptimizeStorage=true);

      /*! \brief Signal that data entry is complete. 
          Off-node entries are distributed, repeated entries are summed, and global indices are transformed to local indices.
          If \c OptimizeStorage is true, then optimizeStorage() is called as well.
          \note This method calls fillComplete( getRowMap(), getRowMap() ).
       */
      void fillComplete(bool OptimizeStorage=true);

      //! \brief Re-allocate the data into contiguous storage.
      void optimizeStorage();

      // @}

      //! @name Extraction Methods
      // @{ 

      //! \brief Get a copy of the diagonal entries owned by this node, with local row idices.
      /*! Returns a distributed Vector object partitioned according to the matrix's row map, containing the 
          the zero and non-zero diagonals owned by this node. */
      void getLocalDiagCopy(Vector<Scalar,LocalOrdinal,GlobalOrdinal> &diag) const;

      //! Returns a copy of the specified (and locally owned) row, using local indices.
      /*! Before fillComplete(), the results will not include entries submitted to another node and may contain duplicated entries.
       * \pre hasColMap() == true
       */
      void getLocalRowCopy(LocalOrdinal localRow, 
                            const Teuchos::ArrayView<LocalOrdinal> &indices, 
                            const Teuchos::ArrayView<Scalar> &values,
                            size_t &numEntries) const;

      //! Returns a copy of the specified (and locally owned) row, using global indices.
      /*! Before fillComplete(), the results will not include entries submitted to another node and may contain duplicated entries. */
      void getGlobalRowCopy(GlobalOrdinal globalRow, 
                                const Teuchos::ArrayView<GlobalOrdinal> &indices,
                                const Teuchos::ArrayView<Scalar> &values,
                                size_t &numEntries) const;

      //! Get a non-persisting view of the elements in a specified global row of the matrix.
      /*!
        \param GlobalRow - (In) Global row from which to retrieve matrix entries.
        \param Indices - (Out) Indices for the global row.
        \param Values - (Out) Values for the global row.

         Note: If \c GlobalRow does not belong to this node, then \c indices and \c values are set to <tt>Teuchos::null</tt>.

        \pre isLocallyIndexed()==false
       */
      void extractGlobalRowView(GlobalOrdinal GlobalRow, 
                                Teuchos::ArrayView<GlobalOrdinal> &indices, 
                                Teuchos::ArrayView<Scalar> &values);

      //! Get a view of the elements in a specified local row of the graph.
      /*!
        \param LocalRow - (In) Local row from which to retrieve matrix entries.
        \param Indices - (Out) Indices for the local row.
        \param Values - (Out) Values for the local row.

         Note: If \c LocalRow is not valid for this node, then \c indices and \c values are set to <tt>Teuchos::null</tt>.

        \pre isLocallyIndexed()==true
       */
      void extractMyRowView(LocalOrdinal LocalRow, 
                                Teuchos::ArrayView<LocalOrdinal> &indices, 
                                Teuchos::ArrayView<Scalar> &values);

      //! Get a non-persisting view of the elements in a specified global row of the graph.
      /*!
        \param GlobalRow - (In) Global row from which to retrieve matrix entries.
        \param Indices - (Out) Indices for the global row.
        \param Values - (Out) Values for the global row.

         Note: If \c GlobalRow does not belong to this node, then \c indices and \c values are set to <tt>Teuchos::null</tt>.

        \pre isLocallyIndexed()==false
       */
      void getGlobalRowView(GlobalOrdinal GlobalRow, 
                                     Teuchos::ArrayView<const GlobalOrdinal> &indices,
                                     Teuchos::ArrayView<const Scalar> &values) const;

      //! Get a view of the elements in a specified local row of the graph.
      /*!
        \param LocalRow - (In) Local row from which to retrieve matrix entries.
        \param Indices - (Out) Indices for the local row.
        \param Values - (Out) Values for the local row.

         Note: If \c LocalRow is not valid for this node, then \c indices and \c values are set to <tt>Teuchos::null</tt>.

        \pre isLocallyIndexed()==true
       */
      void getLocalRowView(LocalOrdinal LocalRow, 
                           Teuchos::ArrayView<const LocalOrdinal> &indices,
                           Teuchos::ArrayView<const Scalar> &values) const;

      //! Returns the CrsGraph associated with this matrix. 
      const CrsGraph<LocalOrdinal,GlobalOrdinal> &getCrsGraph() const;

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

      // multiplication routines
      void GeneralMV (typename MV::const_pointer x       , typename MV::pointer y       ) const;
      void GeneralMM (typename MV::const_double_pointer X, typename MV::double_pointer Y, size_t numVectors) const;
      void GeneralMhV(typename MV::const_pointer x       , typename MV::pointer y       ) const;
      void GeneralMhM(typename MV::const_double_pointer X, typename MV::double_pointer Y, size_t numVectors) const;
      inline typename Teuchos::ArrayRCP<const Scalar>::iterator getVptr(size_t row) const;
      inline typename Teuchos::ArrayRCP<Scalar>::iterator getVptr(size_t row);

      CrsGraph<LocalOrdinal,GlobalOrdinal> graph_;
      bool staticGraph_;
      bool constructedWithFilledGraph_;
      bool fillComplete_;
      bool storageOptimized_;

      //
      // Unoptimized structure
      // These are allocated if storage is not optimized or allocation is not static
      //
      Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > values_;

      //
      // Optimized structure
      // Structure used if allocation is static or after optimizeStorage()
      //
      Teuchos::ArrayRCP<Scalar> contigValues_;
      /* valuesPtrs[j] is the begin() iterator from an ArrayView of 
         contigValues_ corresponding to the proper row, of the appropriate length.
         In a debug build, it is an ArrayRCP, which does bounds checking. in an optimized
         build, it is a C pointer. valuesPtrs is allocated to getNumLocalRows()+1; the span of the jth row begins with
         valuesPtrs[j] and ends before valuesPtrs[j+1] */
      Teuchos::ArrayRCP<typename Teuchos::ArrayRCP<Scalar>::iterator> valuesPtrs_;

      // multivectors used for import/export dest/source in apply()
      mutable Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> > importMV_, exportMV_;

      // a map between a (non-local) row and a list of (col,val)
      std::map<GlobalOrdinal, std::list<std::pair<GlobalOrdinal,Scalar> > > nonlocals_;

  }; // class CrsMatrix



  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::CrsMatrix(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, size_t maxNNZPerRow, bool staticProfile)
  : graph_(rowMap,maxNNZPerRow,staticProfile)
  , staticGraph_(false)
  , constructedWithFilledGraph_(false)
  , fillComplete_(false)
  , storageOptimized_(false)
  , contigValues_(Teuchos::null)
  , importMV_(Teuchos::null)
  , exportMV_(Teuchos::null)
  {
    TEST_FOR_EXCEPTION(maxNNZPerRow < 1 && maxNNZPerRow != 0, std::runtime_error,
        Teuchos::typeName(*this) << "::CrsMatrix(rowMap,maxNNZPerRow): maxNNZPerRow must be non-negative.");
    allocateValues();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::CrsMatrix(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, const Teuchos::ArrayView<size_t> &NNZPerRowToAlloc, bool staticProfile)
  : graph_(rowMap,NNZPerRowToAlloc,staticProfile)
  , staticGraph_(false)
  , constructedWithFilledGraph_(false)
  , fillComplete_(false)
  , storageOptimized_(false)
  , contigValues_(Teuchos::null)
  , importMV_(Teuchos::null)
  , exportMV_(Teuchos::null)
  {
    TEST_FOR_EXCEPTION(NNZPerRowToAlloc.size() != rowMap.getNumMyEntries(), std::runtime_error,
        Teuchos::typeName(*this) << "::CrsMatrix(rowMap,NNZPerRowToAlloc): NNZPerRowToAlloc must have as many entries as specified by rowMap for this node.");
    allocateValues();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::CrsMatrix(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, const Map<LocalOrdinal,GlobalOrdinal> &colMap, size_t maxNNZPerRow, bool staticProfile)
  : graph_(rowMap,colMap,maxNNZPerRow,staticProfile)
  , staticGraph_(false)
  , constructedWithFilledGraph_(false)
  , fillComplete_(false)
  , storageOptimized_(false)
  , contigValues_(Teuchos::null)
  , importMV_(Teuchos::null)
  , exportMV_(Teuchos::null)
  {
    TEST_FOR_EXCEPTION(maxNNZPerRow < 1 && maxNNZPerRow != 0, std::runtime_error,
        Teuchos::typeName(*this) << "::CrsMatrix(rowMap,maxNNZPerRow): maxNNZPerRow must be non-negative.");
    allocateValues();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::CrsMatrix(const Map<LocalOrdinal,GlobalOrdinal> &rowMap, const Map<LocalOrdinal,GlobalOrdinal> &colMap, const Teuchos::ArrayView<size_t> &NNZPerRowToAlloc, bool staticProfile)
  : graph_(rowMap,colMap,NNZPerRowToAlloc,staticProfile)
  , staticGraph_(false)
  , constructedWithFilledGraph_(false)
  , fillComplete_(false)
  , storageOptimized_(false)
  , contigValues_(Teuchos::null)
  , importMV_(Teuchos::null)
  , exportMV_(Teuchos::null)
  {
    TEST_FOR_EXCEPTION(NNZPerRowToAlloc.size() != rowMap.getNumMyEntries(), std::runtime_error,
        Teuchos::typeName(*this) << "::CrsMatrix(rowMap,NNZPerRowToAlloc): NNZPerRowToAlloc must have as many entries as specified by rowMap for this node.");
    allocateValues();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::CrsMatrix(const CrsGraph<LocalOrdinal,GlobalOrdinal> &graph)
  : graph_(graph)
  , staticGraph_(true)
  , fillComplete_(false)
  , storageOptimized_(false)
  , contigValues_(Teuchos::null)
  , importMV_(Teuchos::null)
  , exportMV_(Teuchos::null)
  {
    // we won't prohibit the case where the graph is not yet filled, but we will check below to ensure that the
    // graph isn't filled between now and when fillComplete() is called on this matrix
    constructedWithFilledGraph_ = graph_.isFillComplete();
    allocateValues();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::allocateValues() 
  {
    size_t numlocal = getRowMap().getNumMyEntries();
    if (numlocal > 0) {
      if (graph_.isStaticProfile()) {
        const size_t nta = graph_.totalAllocation();
        size_t sofar = 0;
        valuesPtrs_ = Teuchos::arcp<typename Teuchos::ArrayRCP<Scalar>::iterator>(numlocal+1);
        if (nta) {
          contigValues_ = Teuchos::arcp<Scalar>(nta);
          for (size_t r=0; r<numlocal; ++r) {
            size_t ntarow = graph_.numAllocatedEntriesForMyRow(r);
            valuesPtrs_[r] = contigValues_.persistingView(sofar,ntarow).begin();
            sofar += ntarow;
          }
          valuesPtrs_[numlocal] = contigValues_.end();
#ifdef HAVE_TPETRA_DEBUG
          TEST_FOR_EXCEPTION(sofar != nta, std::logic_error,
              Teuchos::typeName(*this) << "::allocateValues(): Internal logic error. Please contact Tpetra team.");
#endif
        }
        else {
          std::fill(valuesPtrs_.begin(),valuesPtrs_.end(),
                    Teuchos::NullIteratorTraits<typename Teuchos::ArrayRCP<Scalar>::iterator>::getNull());
        }
      }
      else {
        values_ = Teuchos::arcp< Teuchos::ArrayRCP<Scalar> >(numlocal);
        for (size_t r=0; r<numlocal; ++r) {
          size_t ntarow = graph_.numAllocatedEntriesForMyRow(r);
          if (ntarow > 0) {
            values_[r] = Teuchos::arcp<Scalar>(ntarow);
          }
        }
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
  { return storageOptimized_; }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::isLocallyIndexed() const
  { return graph_.isLocallyIndexed(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::isGloballyIndexed() const
  { return graph_.isGloballyIndexed(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::hasColMap() const
  { return graph_.hasColMap(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP<const Teuchos::Comm<int> > 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getComm() const
  { return graph_.getComm(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  global_size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getNumGlobalEntries() const
  { return graph_.getNumGlobalEntries(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getNumLocalEntries() const
  { return graph_.getNumLocalEntries(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getNumGlobalRows() const
  { return graph_.getNumGlobalRows(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getNumGlobalCols() const
  { return graph_.getNumGlobalCols(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getNumLocalRows() const
  { return graph_.getNumLocalRows(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getNumLocalCols() const
  { return graph_.getNumLocalCols(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getNumGlobalDiags() const
  { return graph_.getNumGlobalDiags(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getNumLocalDiags() const
  { return graph_.getNumLocalDiags(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getNumEntriesForGlobalRow(GlobalOrdinal globalRow) const
  { return graph_.getNumEntriesForGlobalRow(globalRow); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getNumEntriesForLocalRow(LocalOrdinal localRow) const
  { 
    using Teuchos::OrdinalTraits;
    if (!getRowMap().isNodeLocalElement(localRow)) return OrdinalTraits<size_t>::invalid();
    return graph_.RNNZ(localRow);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getGlobalMaxNumRowEntries() const
  { return graph_.getGlobalMaxNumRowEntries(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getLocalMaxNumRowEntries() const
  { return graph_.getLocalMaxNumRowEntries(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  GlobalOrdinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getIndexBase() const
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
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::insertMyValues(LocalOrdinal localRow, 
                         const Teuchos::ArrayView<const LocalOrdinal> &indices,
                         const Teuchos::ArrayView<const Scalar>       &values) 
  {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(isStorageOptimized() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertMyValues(): cannot insert new values after optimizeStorage() has been called.");
    TEST_FOR_EXCEPTION(graph_.isGloballyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertMyValues(): graph indices are global; use insertGlobalValues().");
    TEST_FOR_EXCEPTION(hasColMap() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::insertMyValues(): cannot insert local indices without a column map; ");
    TEST_FOR_EXCEPTION(isStaticGraph() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertMyValues(): matrix was constructed with static graph; cannot insert new entries.");
    TEST_FOR_EXCEPTION(values.size() != indices.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::insertMyValues(): values.size() must equal indices.size().");
    TEST_FOR_EXCEPTION(getRowMap().isNodeLocalElement(localRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::insertMyValues(): row does not belong to this node.");
    Teuchos::Array<LocalOrdinal> finds;
    Teuchos::Array<Scalar>       fvals;
    // use column map to filter the entries:
    const Tpetra::Map<LocalOrdinal,GlobalOrdinal> &cmap = getColMap();
    for (size_t i=0; i<indices.size(); ++i) {
      if (cmap.isNodeLocalElement(indices[i])) {
        finds.push_back(indices[i]);
        fvals.push_back(values[i]);
      }
    }
    Teuchos::ArrayView<const LocalOrdinal> findices = finds();
    Teuchos::ArrayView<const Scalar      > fvalues  = fvals();
    size_t rowNNZ = getNumEntriesForLocalRow(localRow),
                     toAdd = findices.size(),
                  rowAlloc = graph_.numAllocatedEntriesForMyRow(localRow);
    if (rowNNZ+toAdd > rowAlloc) {
      TEST_FOR_EXCEPTION(graph_.isStaticProfile() == true, std::runtime_error,
          Teuchos::typeName(*this) << "::insertMyValues(): new indices exceed statically allocated graph structure.");
      TPETRA_EFFICIENCY_WARNING(true,std::runtime_error,
          "::insertMyValues(): Pre-allocated space has been exceeded, requiring new allocation. To improve efficiency, suggest larger allocation.");
#ifdef HAVE_TPETRA_DEBUG
      // assumption: the number currently allocated in values_ is that stored in the graph.
      TEST_FOR_EXCEPTION( rowAlloc != values_[localRow].size(), std::logic_error,
          Teuchos::typeName(*this) << "::insertMyValues(): Internal logic error or unsupported use case. Please contact Tpetra team.");
#endif
      // increase the allocation, copy old entries to new storage
      rowAlloc = rowNNZ+toAdd;
      ArrayRCP<Scalar> newVals = Teuchos::arcp<Scalar>(rowAlloc);
      std::copy(values_[localRow].begin(), values_[localRow].begin()+rowNNZ, newVals.begin());
      values_[localRow] = newVals;
    }
    // insert indices and values
    graph_.insertMyIndices(localRow,findices);  // this will enlarge allocation and append new indices
    std::copy( fvalues.begin(), fvalues.end(),  getVptr(localRow)+rowNNZ);
#ifdef HAVE_TPETRA_DEBUG
    // the assumption is that graph_.numAllocatedEntriesForMyRow is the allocated size here
    TEST_FOR_EXCEPTION( rowAlloc != graph_.numAllocatedEntriesForMyRow(localRow) 
                        || rowNNZ+toAdd != getNumEntriesForLocalRow(localRow), std::logic_error,
                        Teuchos::typeName(*this) << "::insertMyValues(): Internal logic error or unsupported use case. Please contact Tpetra team.");
#endif
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::insertGlobalValues(GlobalOrdinal globalRow, 
                         const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                         const Teuchos::ArrayView<const Scalar>  &values)
  {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(graph_.isLocallyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertGlobalValues(): graph indices are local; use insertMyValues().");
    TEST_FOR_EXCEPTION(isStaticGraph() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertGlobalValues(): matrix was constructed with static graph. Cannot insert new entries.");
    TEST_FOR_EXCEPTION(values.size() != indices.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::insertGlobalValues(): values.size() must equal indices.size().");

    LocalOrdinal myRow = getRowMap().getLocalIndex(globalRow);

    // if we have a column map, use it to filter the entries.
    // only filter if this is our row.
    Teuchos::Array<GlobalOrdinal> finds;
    Teuchos::Array<Scalar>        fvals;
    Teuchos::ArrayView<const GlobalOrdinal> findices = indices;
    Teuchos::ArrayView<const Scalar       > fvalues  = values;
    if (hasColMap() && myRow != LOT::invalid()) {
      const Tpetra::Map<LocalOrdinal,GlobalOrdinal> &cmap = getColMap();
      for (size_t i=0; i<indices.size(); ++i) {
        if (cmap.isNodeGlobalElement(indices[i])) {
          finds.push_back(indices[i]);
          fvals.push_back(values[i]);
        }
      }
      findices = finds();
      fvalues  = fvals();
    }

    // add the new indices and values
    if (myRow != LOT::invalid()) {
      size_t rowNNZ = getNumEntriesForLocalRow(myRow),
                      toAdd = findices.size(),
                   rowAlloc = graph_.numAllocatedEntriesForMyRow(myRow);
      if (rowNNZ+toAdd > rowAlloc) {
        TEST_FOR_EXCEPTION(graph_.isStaticProfile() == true, std::runtime_error,
            Teuchos::typeName(*this) << "::insertGlobalValues(): new indices exceed statically allocated graph structure.");
        TPETRA_EFFICIENCY_WARNING(true,std::runtime_error,
            "::insertGlobalValues(): Pre-allocated space has been exceeded, requiring new allocation. To improve efficiency, suggest larger allocation.");
#ifdef HAVE_TPETRA_DEBUG
        // assumption: the number currently allocated in values_ is that stored in the graph.
        TEST_FOR_EXCEPTION( rowAlloc != values_[myRow].size(), std::logic_error,
            Teuchos::typeName(*this) << "::insertGlobalValues(): Internal logic error or unsupported use case. Please contact Tpetra team.");
#endif
        // increase the allocation, copy old entries to new storage
        rowAlloc = rowNNZ+toAdd;
        ArrayRCP<Scalar> newVals = Teuchos::arcp<Scalar>(rowAlloc);
        std::copy(values_[myRow].begin(), values_[myRow].begin()+rowNNZ, newVals.begin());
        values_[myRow] = newVals;
      }
      // insert indices and values
      graph_.insertGlobalIndices(globalRow,findices);  // this will enlarge allocation and append new indices
      std::copy( fvalues.begin(), fvalues.end(),  getVptr(myRow)+rowNNZ);
#ifdef HAVE_TPETRA_DEBUG
      // the assumption is that graph_.numAllocatedEntriesForMyRow is the allocated size here as well
      TEST_FOR_EXCEPTION( rowAlloc != graph_.numAllocatedEntriesForMyRow(myRow) 
                          || rowNNZ+toAdd != getNumEntriesForLocalRow(myRow), std::logic_error,
          Teuchos::typeName(*this) << "::insertGlobalValues(): Internal logic error or unsupported use case. Please contact Tpetra team.");
#endif
    }
    else {
      typename Teuchos::ArrayView<const GlobalOrdinal>::iterator ind = findices.begin();
      typename Teuchos::ArrayView<const Scalar       >::iterator val =  fvalues.begin();
      for (; val != fvalues.end(); ++val, ++ind) {
        nonlocals_[globalRow].push_back(std::make_pair(*ind, *val));
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::replaceGlobalValues(GlobalOrdinal globalRow, 
                         const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                         const Teuchos::ArrayView<const Scalar>  &values)
  {
    // find the values for the specified indices
    // if the row is not ours, throw an exception
    // ignore values not in the matrix (indices not found)
    // operate whether indices are local or global
    const size_t TOINV = Teuchos::OrdinalTraits<size_t>::invalid();
    TEST_FOR_EXCEPTION(values.size() != indices.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::replaceGlobalValues(): values.size() must equal indices.size().");
    typename Teuchos::ArrayView<const GlobalOrdinal>::iterator ind = indices.begin();
    typename Teuchos::ArrayView<const        Scalar>::iterator val = values.begin();
    LocalOrdinal lrow = getRowMap().getLocalIndex(globalRow);
    TEST_FOR_EXCEPTION(lrow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::replaceGlobalValues(): specified global row does not belong to this processor.");
    typename Teuchos::ArrayRCP<Scalar>::iterator vptr = getVptr(lrow);
    if (isLocallyIndexed() == true) {
      while (ind != indices.end()) {
        LocalOrdinal lind = getColMap().getLocalIndex(*ind);
        size_t loc = graph_.findMyIndex(lrow,lind);
        if (loc != TOINV) {
          vptr[loc] = (*val);
        }
        ++ind;
        ++val;
      }
    }
    else if (isGloballyIndexed() == true) {
      while (ind != indices.end()) {
        size_t loc = graph_.findGlobalIndex(lrow,*ind);
        if (loc != TOINV) {
          vptr[loc] = (*val);
        }
        ++ind;
        ++val;
      }
    }
    //else {
      // indices are not allocated; nothing to do
    //}
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::sumIntoGlobalValues(GlobalOrdinal globalRow, 
                         const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                         const Teuchos::ArrayView<const Scalar>  &values)
  {
    // find the values for the specified indices
    // if the row is not ours, throw an exception
    // ignore values not in the matrix (indices not found)
    // operate whether indices are local or global
    const size_t TOINV = Teuchos::OrdinalTraits<size_t>::invalid();
    TEST_FOR_EXCEPTION(values.size() != indices.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::sumIntoGlobalValues(): values.size() must equal indices.size().");
    typename Teuchos::ArrayView<const GlobalOrdinal>::iterator ind = indices.begin();
    typename Teuchos::ArrayView<const        Scalar>::iterator val = values.begin();
    LocalOrdinal lrow = getRowMap().getLocalIndex(globalRow);
    TEST_FOR_EXCEPTION(lrow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::sumIntoGlobalValues(): specified global row does not belong to this processor.");
    typename Teuchos::ArrayRCP<Scalar>::iterator vptr = getVptr(lrow);
    if (isLocallyIndexed() == true) {
      while (ind != indices.end()) {
        LocalOrdinal lind = getColMap().getLocalIndex(*ind);
        size_t loc = graph_.findMyIndex(lrow,lind);
        if (loc != TOINV) {
          vptr[loc] += (*val);
        }
        ++ind;
        ++val;
      }
    }
    else if (isGloballyIndexed() == true) {
      while (ind != indices.end()) {
        size_t loc = graph_.findGlobalIndex(lrow,*ind);
        if (loc != TOINV) {
          vptr[loc] += (*val);
        }
        ++ind;
        ++val;
      }
    }
    //else {
      // indices are not allocated; nothing to do
    //}
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::scale(const Scalar &alpha)
  { 
    TEST_FOR_EXCEPT(isStorageOptimized());
    for (size_t r=0; r<getNumLocalRows(); ++r) {
      typename Teuchos::ArrayRCP<Scalar>::iterator val, vend;
      val = getVptr(r);
      vend = val+getNumEntriesForLocalRow(r);
      while (val != vend) {
        (*val++) *= alpha;
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getLocalRowCopy(LocalOrdinal LocalRow, 
                                                                     const Teuchos::ArrayView<LocalOrdinal> &indices, 
                                                                     const Teuchos::ArrayView<Scalar>       &values,
                                                                     size_t &numEntries) const
  {
    numEntries = getNumEntriesForLocalRow(LocalRow);
    TEST_FOR_EXCEPTION(getRowMap().isNodeLocalElement(LocalRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowCopy(LocalRow,...): specified row (==" << LocalRow << ") is not valid on this node.");
    TEST_FOR_EXCEPTION(indices.size() < numEntries || values.size() < numEntries, std::runtime_error, 
        Teuchos::typeName(*this) << "::getLocalRowCopy(LocalRow,indices,values): size of indices,values must be sufficient to store the specified row.");
#ifdef HAVE_TPETRA_DEBUG
    size_t nnzagain;
    graph_.getLocalRowCopy(LocalRow,indices,nnzagain);
    TEST_FOR_EXCEPTION(nnzagain != numEntries, std::logic_error, 
        Teuchos::typeName(*this) << "::getLocalRowCopy(): Internal logic error. Please contact Tpetra team.");
#else
    graph_.getLocalRowCopy(LocalRow,indices,numEntries);
#endif
    typename Teuchos::ArrayRCP<const Scalar>::iterator vptr = getVptr(LocalRow);
    std::copy( vptr, vptr+numEntries, values.begin() );
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getGlobalRowCopy(GlobalOrdinal globalRow, 
                                                                      const Teuchos::ArrayView<GlobalOrdinal> &indices,
                                                                      const Teuchos::ArrayView<Scalar>  &values,
                                                                      size_t &numEntries) const
  {
    // Only locally owned rows can be queried, otherwise complain
    size_t myRow = getRowMap().getLocalIndex(globalRow);
    TEST_FOR_EXCEPTION(myRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowCopy(globalRow,...): globalRow does not belong to this node.");
    numEntries = getNumEntriesForLocalRow(myRow);
    TEST_FOR_EXCEPTION(
        indices.size() < numEntries || values.size() < numEntries, std::runtime_error, 
        Teuchos::typeName(*this) << "::getGlobalRowCopy(globalRow,indices,values): size of indices,values must be sufficient to store the specified row.");
#ifdef HAVE_TPETRA_DEBUG
    size_t nnzagain;
    graph_.getGlobalRowCopy(globalRow,indices,nnzagain);
    TEST_FOR_EXCEPTION(nnzagain != numEntries, std::logic_error, 
        Teuchos::typeName(*this) << "::getLocalRowCopy(): Internal logic error. Please contact Tpetra team.");
#else
    graph_.getGlobalRowCopy(globalRow,indices,numEntries);
#endif
    typename Teuchos::ArrayRCP<const Scalar>::iterator vptr = getVptr(myRow);
    std::copy( vptr, vptr+numEntries, values.begin() );
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::extractGlobalRowView(GlobalOrdinal GlobalRow, 
                                Teuchos::ArrayView<GlobalOrdinal> &indices, 
                                Teuchos::ArrayView<Scalar> &values)
  {
    TEST_FOR_EXCEPTION(isLocallyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::extractGlobalRowView(): global indices do not exist; call extractMyRowView().");
    size_t lrow = getRowMap().getLocalIndex(GlobalRow);
    TEST_FOR_EXCEPTION(lrow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::extractGlobalRowView(GlobalRow,...): GlobalRow (== " << GlobalRow << ") does not belong to this node.");
    size_t rnnz = getNumEntriesForLocalRow(lrow);
    graph_.extractGlobalRowView(GlobalRow,indices);
    if (rnnz == 0) {
      values = Teuchos::ArrayView<Scalar>(Teuchos::null);
    }
    else {
      if (isStorageOptimized() == true || graph_.isStaticProfile() == true) {
        values = Teuchos::arrayView<Scalar>(Teuchos::getRawPtr(valuesPtrs_[lrow]),rnnz);
      }
      else {
        values = values_[lrow](0,rnnz);
      }
    }
    return;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::extractMyRowView(LocalOrdinal LocalRow, 
                                Teuchos::ArrayView<LocalOrdinal> &indices, 
                                Teuchos::ArrayView<Scalar> &values)
  {
    TEST_FOR_EXCEPTION(isGloballyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::extractMyRowView(): local indices do not exist; call extractGlobalRowVie().");
    TEST_FOR_EXCEPTION(getRowMap().isNodeLocalElement(LocalRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::extractMyRowView(LocalRow,...): LocalRow (== " << LocalRow << ") is not valid on this node.");
    size_t rnnz = getNumEntriesForLocalRow(LocalRow);
    graph_.extractMyRowView(LocalRow,indices);
    if (rnnz == 0) {
      values = Teuchos::ArrayView<Scalar>(Teuchos::null);
    }
    else {
      if (isStorageOptimized() == true || graph_.isStaticProfile() == true) {
        values = Teuchos::arrayView<Scalar>(Teuchos::getRawPtr(valuesPtrs_[LocalRow]),rnnz);
      }
      else {
        values = values_[LocalRow](0,rnnz);
      }
    }
    return;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getGlobalRowView(GlobalOrdinal GlobalRow, 
                                     Teuchos::ArrayView<const GlobalOrdinal> &indices,
                                     Teuchos::ArrayView<const Scalar> &values) const
  {
    TEST_FOR_EXCEPTION(isLocallyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowView(): global indices do not exist; call getLocalRowView().");
    size_t lrow = getRowMap().getLocalIndex(GlobalRow);
    TEST_FOR_EXCEPTION(lrow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowView(GlobalRow,...): GlobalRow (== " << GlobalRow << ") does not belong to this node.");
    size_t rnnz = getNumEntriesForLocalRow(lrow);
    graph_.getGlobalRowView(GlobalRow,indices);
    if (rnnz == 0) {
      values = Teuchos::ArrayView<Scalar>(Teuchos::null);
    }
    else {
      if (isStorageOptimized() == true || graph_.isStaticProfile() == true) {
        values = Teuchos::arrayView<Scalar>(Teuchos::getRawPtr(valuesPtrs_[lrow]),rnnz);
      }
      else {
        values = values_[lrow](0,rnnz);
      }
    }
    return;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getLocalRowView(LocalOrdinal LocalRow, 
                                 Teuchos::ArrayView<const LocalOrdinal> &indices,
                                 Teuchos::ArrayView<const Scalar> &values) const
  {
    TEST_FOR_EXCEPTION(isGloballyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowView(): local indices do not exist; call getGlobalRowView().");
    TEST_FOR_EXCEPTION(getRowMap().isNodeLocalElement(LocalRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowView(LocalRow,...): LocalRow (== " << LocalRow << ") is not valid on this node.");
    size_t rnnz = getNumEntriesForLocalRow(LocalRow);
    graph_.getLocalRowView(LocalRow,indices);
    if (rnnz == 0) {
      values = Teuchos::ArrayView<Scalar>(Teuchos::null);
    }
    else {
      if (isStorageOptimized() == true || graph_.isStaticProfile() == true) {
        values = Teuchos::arrayView<Scalar>(Teuchos::getRawPtr(valuesPtrs_[LocalRow]),rnnz);
      }
      else {
        values = values_[LocalRow](0,rnnz);
      }
    }
    return;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &X, 
                                                                 MultiVector<Scalar,LocalOrdinal,GlobalOrdinal> &Y, 
                                                                 Teuchos::ETransp mode) const
  {
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

    const size_t numVectors = X.numVectors();
    // because of Views, it is difficult to determine if X and Y point to the same data. 
    // however, if they reference the exact same object, we will do the user the favor of copying X into new storage (with a warning)
    // we ony need to do this if we have trivial importers; otherwise, we don't actually apply the operator from X into Y
    Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal> > importer = graph_.getImporter();
    Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal> > exporter = graph_.getExporter();
    Teuchos::RCP<const MV> Xcopy;
    typename MV::const_double_pointer Xdata = X.extractConstView2D();
    typename MV::double_pointer       Ydata = Y.extractView2D();
    if (&X==&Y && importer==null && exporter==null) {
      TPETRA_EFFICIENCY_WARNING(true,std::runtime_error,
          "::apply(X,Y): If X and Y are the same, it necessitates a temporary copy of X, which is inefficient.");
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
      os << "Number of global rows    = " << getNumGlobalRows() << endl;
      if (isFillComplete())
      {
        os << "Number of global columns    = " << getNumGlobalCols() << endl;
        os << "Status = fill complete" << endl;
        os << "Number of global nonzeros   = " << getNumGlobalEntries() << endl;
        os << "Global max nonzeros per row = " << getGlobalMaxNumRowEntries() << endl;
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
          Teuchos::Array<GlobalOrdinal> indices(getLocalMaxNumRowEntries());
          Teuchos::Array<Scalar>         values(getLocalMaxNumRowEntries());
          size_t rowSize;
          os << "% Number of rows on image " << myImageID << " = " << getNumLocalRows() << endl;
          for (size_t i=0; i < getNumLocalRows(); ++i)
          {
            GlobalOrdinal globalRow = getRowMap().getGlobalIndex(i);
            getGlobalRowCopy(globalRow, indices(), values(), rowSize);
            if (rowSize > Teuchos::OrdinalTraits<size_t>::zero()) {
              for (size_t j=0; j < rowSize; ++j) {
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
      // if the graph is filled, it should be because 
      // - it was given at construction, already filled
      // - it was filled on a previous call to fillComplete()
      // if neither, it means that the graph has been filled by the user.
      // this means that a graph passed at construction has been filled in the interim.
      // in this case, indices have been sorted and merged, so that we are no longer able
      // to sort or merge the associated values. nothing we can do here.
      if ((constructedWithFilledGraph_ || isFillComplete()) == false) {
        TPETRA_ABUSE_WARNING(true,std::runtime_error,
            "::fillComplete(): fillComplete() has been called on graph since matrix was constructed.");
        return;
      }
    }

    if (!isStaticGraph()) {
      graph_.makeIndicesLocal(domainMap,rangeMap);
    }

    sortEntries();
    mergeRedundantEntries();

    if (!isStaticGraph()) {
      graph_.fillComplete(domainMap,rangeMap,OptimizeStorage);
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
    // optimizeStorage will perform two functions:
    // 1) create a single allocation of memory
    // 2) pack data in that allocation
    // if isStaticProfile() == true, then 1) has already been done
    if (isStorageOptimized() == true) return;
    TEST_FOR_EXCEPTION(isFillComplete() == false || graph_.indicesAreSorted() == false || graph_.noRedundancies() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::optimizeStorage(): fillComplete() must be called before optimizeStorage().");
    
    // 1) allocate single memory block
    const size_t nlrs = getNumLocalRows();
    if (nlrs > 0) {
      if (graph_.isStaticProfile() == false) {
        // need to allocate storage, create pointers, copy data, and delete old data
        const size_t nE = getNumLocalEntries();
        valuesPtrs_ = Teuchos::arcp<typename Teuchos::ArrayRCP<Scalar>::iterator>(nlrs+1);
        if (nE > 0) {
          contigValues_ = Teuchos::arcp<Scalar>(nE);
          size_t sofar = 0;
          for (size_t r=0; r<nlrs; ++r) {
            size_t rne = graph_.RNNZ(r);
            valuesPtrs_[r] = contigValues_.persistingView(sofar,rne).begin();
            std::copy(values_[r].begin(),values_[r].begin()+rne,valuesPtrs_[r]);
            values_[r] = Teuchos::null;
            sofar += rne;
          }
          valuesPtrs_[nlrs] = contigValues_.end();
        }
        else {
          std::fill(valuesPtrs_.begin(),valuesPtrs_.end(),
                    Teuchos::NullIteratorTraits<typename Teuchos::ArrayRCP<Scalar>::iterator>::getNull());
        }
        values_ = Teuchos::null;
      }
      else {
        // storage is already allocated and pointers are set; just need to pack
        // need to change pointers and move data
        // remember, these aren't just pointers, but also bounds-checked 
        const size_t nE = getNumLocalEntries();
        if (nE > 0) {
          size_t sofar = 0;
          typename Teuchos::ArrayRCP<Scalar>::iterator newptr, oldptr;
          for (size_t r=0; r<nlrs; ++r) {
            size_t rne = graph_.RNNZ(r);
            newptr = contigValues_.persistingView(sofar,rne).begin();
            oldptr = valuesPtrs_[r];
            if (newptr != oldptr) {
              for (size_t j=0; j<rne; ++j) {
                newptr[j] = oldptr[j];
              }
              valuesPtrs_[r] = newptr;
            }
            sofar += rne;
          }
          valuesPtrs_[nlrs] = contigValues_.persistingView(sofar,0).begin();
        }
      }
    }
    storageOptimized_ = true;
#ifdef HAVE_TPETRA_DEBUG
    // check that we only have the data structures that we need
    TEST_FOR_EXCEPTION( values_ != Teuchos::null, std::logic_error,
        Teuchos::typeName(*this) << "::optimizeStorage(): Internal logic error. Please contact Tpetra team.");
#endif
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
    using Teuchos::tuple;
    typedef typename std::map<GlobalOrdinal,std::list<pair<GlobalOrdinal,Scalar> > >::const_iterator NLITER;
    int numImages = getComm()->getSize();
    int myImageID = getComm()->getRank();
    // Determine if any nodes have global entries to share
    size_t MyNonlocals = nonlocals_.size(), 
                    MaxGlobalNonlocals;
    Teuchos::reduceAll<size_t>(*getComm(),Teuchos::REDUCE_MAX,MyNonlocals,&MaxGlobalNonlocals);
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
    size_t numRecvs = recvIDs.size();

    // we know how many we're sending to already
    // form a contiguous list of all data to be sent
    // track the number of entries for each ID
    Array<CrsIJV<GlobalOrdinal,Scalar> > IJVSendBuffer;
    Array<size_t> sendSizes(sendIDs.size(), 0);
    size_t numSends = 0;
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
    Array<size_t> recvSizes(numRecvs);
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
      size_t cur = 0;
      for (size_t s=0; s<numSends; ++s) {
        sendBuffers[s] = IJVSendBuffer(cur,sendSizes[s]);
        cur += sendSizes[s];
      }
    }
    // perform non-blocking sends
    for (size_t s=0; s < numSends ; ++s)
    {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      ArrayRCP<CrsIJV<GlobalOrdinal,Scalar> > tmparcp = arcp(sendBuffers[s].getRawPtr(),0,sendBuffers[s].size(),false);
      sendRequests.push_back( Teuchos::isend<int,CrsIJV<GlobalOrdinal,Scalar> >(*getComm(),tmparcp,sendIDs[s]) );
    }
    // calculate amount of storage needed for receives
    // setup pointers for the receives as well
    size_t totalRecvSize = std::accumulate(recvSizes.begin(),recvSizes.end(),0);
    Array<CrsIJV<GlobalOrdinal,Scalar> > IJVRecvBuffer(totalRecvSize);
    // from the size info, build the ArrayViews into IJVRecvBuffer
    Array<ArrayView<CrsIJV<GlobalOrdinal,Scalar> > > recvBuffers(numRecvs,Teuchos::null);
    {
      size_t cur = 0;
      for (size_t r=0; r<numRecvs; ++r) {
        recvBuffers[r] = IJVRecvBuffer(cur,recvSizes[r]);
        cur += recvSizes[r];
      }
    }
    // perform non-blocking recvs
    for (size_t r=0; r < numRecvs ; ++r)
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
      try {
        insertGlobalValues(ijv->i, tuple(ijv->j), tuple(ijv->v));
      }
      catch (std::runtime_error &e) {
        std::ostringstream omsg;
        omsg << e.what() << std::endl
          << "caught in globalAssemble() in " << __FILE__ << ":" << __LINE__ << std::endl ;
        throw std::runtime_error(omsg.str());
      }
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
    typename ArrayRCP<const Scalar>::iterator aval;
    for (size_t r=0; r < getRowMap().getNumMyEntries(); ++r) {
      graph_.getLocalRowView(r,cinds);
      Scalar sum = ST::zero();
      cind = cinds.begin();
      aval = getVptr(r);
      while (cind != cinds.end()) 
      {
        sum += (*aval++) * x[*cind++];
      }
      y[r] = sum;
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::GeneralMM (typename MV::const_double_pointer X, typename MV::double_pointer Y, size_t numVectors) const
  {
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    size_t nlrs = getNumLocalRows();
    ArrayView<const LocalOrdinal> cinds;
    typename ArrayView<const LocalOrdinal>::iterator cind;
    typename ArrayRCP<const Scalar>::iterator aval;
    for (size_t r=0; r < nlrs; ++r) {
      graph_.getLocalRowView(r,cinds);
      for (size_t j=0; j<numVectors; ++j) {
        typename MV::pointer       yvals = Y[j];
        typename MV::const_pointer xvals = X[j]; 
        cind = cinds.begin(); 
        aval = getVptr(r);
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
    size_t nlrs = getNumLocalRows(),
                    nlcs = getNumLocalCols();
    // Initialize y for transpose multiply
    std::fill( y, y+nlcs, ST::zero() );
    // apply conjugate transpose of matrix to x
    // use column triad formulation, accumulating into y each column of A^H (each row of A) times each entry in x
    ArrayView<const LocalOrdinal> cinds;
    typename ArrayView<const LocalOrdinal>::iterator cind;
    typename ArrayRCP<const Scalar>::iterator aval;
    for (size_t r=0; r < nlrs; ++r) {
      graph_.getLocalRowView(r,cinds);
      // loop over entries in this column of A^H (this row of A)
      cind = cinds.begin();
      aval = getVptr(r);
      while (cind != cinds.end()) {
        y[*cind++] += ST::conjugate(*aval++) * x[r];
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::GeneralMhM(typename MV::const_double_pointer X, typename MV::double_pointer Y, size_t numVectors) const
  {
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    size_t nlrs = getNumLocalRows(),
                    nlcs = getNumLocalCols();
    // Initialize Y for transpose multiply
    for (size_t j=0; j<numVectors; ++j) {
      std::fill( Y[j], Y[j]+nlcs, ST::zero() );
    }
    // apply conjugate transpose of matrix to X
    // use outer-product formulation, hitting Y with a rank-1 update comprised of each column of A^H (each row of A) and each row of X
    ArrayView<const LocalOrdinal> cinds;
    typename ArrayView<const LocalOrdinal>::iterator cind;
    typename ArrayRCP<const Scalar>::iterator aval;
    for (size_t r=0; r < nlrs; ++r) {
      graph_.getLocalRowView(r,cinds);
      // loop over numvectors
      for (size_t j=0; j<numVectors; ++j) {
        typename MV::pointer       yvals = Y[j];
        typename MV::const_pointer xvals = X[j]; 
        // loop over entries in this column of A^H (this row of A)
        cind = cinds.begin();
        aval = getVptr(r);
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
    size_t nlrs = getNumLocalRows();
    typename Teuchos::ArrayRCP<const Scalar>::iterator v;
    typename Teuchos::ArrayView<Scalar>::iterator ov;
    Teuchos::ArrayView<const LocalOrdinal> cinds;
    typename Teuchos::ArrayView<const LocalOrdinal>::iterator i;
    ov = values.begin();
    for (size_t r=0; r < nlrs; ++r) {
      *ov = Teuchos::ScalarTraits<Scalar>::zero();
      GlobalOrdinal rgid = getRowMap().getGlobalIndex(r);
      if (getColMap().isNodeGlobalElement(rgid)) {
        LocalOrdinal rlid = getColMap().getLocalIndex(rgid);
        graph_.getLocalRowView(r,cinds);
        i = cinds.begin();
        v = getVptr(r);
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
    TEST_FOR_EXCEPTION(numDiagFound != graph_.getNumLocalDiags(), std::logic_error, 
        "CrsMatrix::getLocalDiagCopy(): logic error. Please contact Tpetra team.");
#endif
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::isStaticGraph() const
  { return staticGraph_; }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  const RowGraph<LocalOrdinal,GlobalOrdinal> &
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getGraph() const
  { return graph_; }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  const CrsGraph<LocalOrdinal,GlobalOrdinal> &
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getCrsGraph() const
  { return graph_; }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::sortEntries()
  {
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(graph_.isGloballyIndexed() == true, std::logic_error,
        Teuchos::typeName(*this) << "::sortEntries(): sortEntries() must be called after indices are transformed to local.\n"
        << "Likely internal logic error. Please contact Tpetra team.");
    if (graph_.indicesAreSorted()) return;
    if (graph_.indicesAreAllocated()) {
      const size_t nlrs = getNumLocalRows();
      for (size_t r=0; r < nlrs; ++r)
      {
        ArrayView<LocalOrdinal> inds;
        typename ArrayRCP<Scalar>::iterator vals = getVptr(r);
        graph_.extractMyRowView(r,inds);
        sort2(inds.begin(), inds.end(), vals);
      }
    }
    graph_.indicesAreSorted(true);  // we just sorted them
    return;
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::isLowerTriangular() const
  { return graph_.isLowerTriangular(); }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::isUpperTriangular() const
  { return graph_.isUpperTriangular(); }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::mergeRedundantEntries() 
  {
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(graph_.indicesAreSorted() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::mergeRedundantEntries() cannot be called before indices are sorted.\n"
        << "Likely interal logic error. Please contact Tpetra team.");
    if (graph_.noRedundancies()) return;
    for (size_t r=0; r<getNumLocalRows(); ++r) 
    {
      size_t rnnz = getNumEntriesForLocalRow(r);
      if (rnnz > 1) {
        ArrayView<const LocalOrdinal> inds;
        graph_.getLocalRowView(r,inds);
        typename ArrayRCP<Scalar>::iterator vals = getVptr(r);
        size_t curEntry = 0;
        Scalar curValue = vals[0];
        for (size_t k=1; k<rnnz; ++k) {
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

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  typename Teuchos::ArrayRCP<const Scalar>::iterator 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getVptr(size_t row) const
  {
    if (graph_.isStaticProfile() || isStorageOptimized()) {
      return valuesPtrs_[row];
    }
    else {
      return values_[row].begin();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  typename Teuchos::ArrayRCP<Scalar>::iterator 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal>::getVptr(size_t row)
  {
    if (graph_.isStaticProfile() || isStorageOptimized()) {
      return valuesPtrs_[row];
    }
    else {
      return values_[row].begin();
    }
  }

} // namespace Tpetra

#endif
