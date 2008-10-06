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
  class SerializationTraits<Ordinal,Tpetra::CrsIJV<Ordinal,Scalar> >
  : public DirectSerializationTraits<Ordinal,Tpetra::CrsIJV<Ordinal,Scalar> >
  {};
}

namespace Tpetra 
{
  //! Tpetra::CrsMatrix: A class for constructing and using sparse compressed index matrices and row access.
  /*!
   * This class allows the construction of sparse matrices with row-access. 
   * Methods submitEntry() and submitEntries() can be used to set both locally
   * owned and non-local elements; the shipping of data is done with hardcoded
   * MPI calls when fillComplete() is called.
   *
   * The nonzero elements of  locally owned row can be accessed by method
   * getMyRowCopy() or getGlobalRowCopy(). The former returns the column
   * indices using local numbering, the latter using global numbering.
   *
   */
  template<class Ordinal, class Scalar>
  class CrsMatrix : public Operator<Ordinal,Scalar>
  {
    public:
      //! @name Constructor/Destructor Methods
      //@{ 

      //! Constructor specifying the domain map only.
      CrsMatrix(const Map<Ordinal> &rowMap);

      // !Destructor.
      virtual ~CrsMatrix();

      //@}
      //! @name Query Methods
      //@{ 
      
      //! Returns \c true if the matrix has already been fill completed.
      inline bool isFillCompleted() const;

      //! Returns the communicator.
      Teuchos::RCP<const Teuchos::Comm<Ordinal> > getComm() const;

      //! Returns the number of nonzero entries in the global matrix. 
      inline Ordinal getNumGlobalNonzeros() const;

      //! Returns the number of nonzero entries in the calling image's portion of the matrix. 
      inline Ordinal getNumMyNonzeros() const;

      //! Returns the number of global matrix rows. 
      inline Ordinal getNumGlobalRows() const;

      //! Returns the number of global matrix columns. 
      inline Ordinal getNumGlobalCols() const;

      //! Returns the number of matrix rows owned by the calling image. 
      inline Ordinal getNumMyRows() const;

      //! Returns the number of matrix columns referenced by the calling image. 
      inline Ordinal getNumMyCols() const;

      //! Returns the number of global nonzero diagonal entries, based on global row/column index comparisons. 
      inline Ordinal getNumGlobalDiagonals() const;

      //! Returns the number of local nonzero diagonal entries, based on global row/column index comparisons. 
      inline Ordinal getNumMyDiagonals() const;

      //! Returns the maximum number of nonzero entries across all rows/columns on all images. 
      inline Ordinal getGlobalMaxNumEntries() const;

      //! Returns the maximum number of nonzero entries across all rows/columns on this image. 
      inline Ordinal getMyMaxNumEntries() const;

      //! Returns the index base for global indices for this matrix. 
      inline Ordinal getIndexBase() const;
      
      //! Returns the Map that describes the row distribution in this matrix.
      const Map<Ordinal> & getRowMap() const;

      //! Returns the Map that describes the column distribution in this matrix.
      const Map<Ordinal> & getColMap() const;

      //@}

      //! @name Methods implementing Tpetra::Operator
      //@{ 
      
      //! Returns the Map associated with the domain of this operator.
      const Map<Ordinal> & getDomainMap() const;

      //! Returns the Map associated with the domain of this operator.
      const Map<Ordinal> & getRangeMap() const;

      //! Computes the matrix-vector multilication y = A x.
      void apply(const MultiVector<Ordinal,Scalar>& X, MultiVector<Ordinal,Scalar> &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

      //@}

      //! @name Construction Methods
      //@{ 

      //! Signals that data entry is complete. Matrix data is converted into a more optimized form.
      void fillComplete();

      //! Submits one local or nonlocal entry to the matrix using global IDs.
      void submitEntry(Ordinal globalRow, Ordinal globalCol,
                       const Scalar &value);

      //! Submit multiple entries, using global IDs.
      /*! All index values must be in the global space. Behavoir is defined by the CombineMode passed in. */
      void submitEntries(Ordinal globalRow, 
                         const Teuchos::ArrayView<const Ordinal> &cols,
                         const Teuchos::ArrayView<const Scalar> &vals);

      //! Set all matrix entries equal to scalarThis.
      void setAllToScalar(const Scalar &alpha);

      //! Scale the current values of a matrix, this = alpha*this. 
      void scale(const Scalar &alpha);

      // @}

      //! @name Data Access Methods
      // @{ 

      //! Returns the current number of nonzero entries in specified global index on this image. 
      inline Ordinal getNumRowEntries(Ordinal globalRow) const;

      //! Returns a copy of the specified local row, column indices are local.
      void getMyRowCopy(Ordinal myRow, const Teuchos::ArrayView<Ordinal> &indices, 
                                       const Teuchos::ArrayView<Scalar> &values) const;

      //! Returns a copy of the specified (and locally owned) global row, column indices are global.
      void getGlobalRowCopy(Ordinal globalRow, 
                            const Teuchos::ArrayView<Ordinal> &indices,
                            const Teuchos::ArrayView<Scalar> &values) const;

      //@}

      //! @name I/O Methods
      //@{ 
      
      //! Prints the matrix on the specified stream. This is very verbose.
      void print(std::ostream& os) const;

      // @}

    private:
      // copy constructor disabled
      CrsMatrix(const CrsMatrix<Ordinal,Scalar> &Source);
      // operator= disabled
      CrsMatrix& operator=(const CrsMatrix<Ordinal, Scalar> &rhs);
      // useful typedefs
      typedef Teuchos::OrdinalTraits<Ordinal> OT;

      // Performs importing of off-processor elements and adds them to the locally owned elements.
      void globalAssemble();
      // multiplication routines
      void GeneralMV(const Teuchos::ArrayView<const Scalar> &X, const Teuchos::ArrayView<Scalar> &Y) const;
      void GeneralMM(const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> > &X, const Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> > &Y) const;

      const Teuchos::RCP<const Teuchos::Comm<Ordinal> > comm_;
      const Map<Ordinal>& rowMap_;
      Teuchos::RCP<const Map<Ordinal> > colMap_;

      Ordinal numMyRows_, numGlobalRows_;

      // column indices
      Teuchos::Array<Teuchos::Array<Ordinal> > colinds_;
      // values
      Teuchos::Array<Teuchos::Array<Scalar> > values_;

      // TODO: consider using a contiguous storage
      // costs are: insertions are more difficult
      mutable Teuchos::RCP<MultiVector<Ordinal,Scalar> > importMV_, exportMV_;

      // importer is needed if DomainMap is not sameas ColumnMap
      // FINISH: currently, domain map == range map == row map which is typically != column map
      //         ergo, we will usually have an importer
      Teuchos::RCP<Import<Ordinal> > importer_;
      // exporter is needed if RowMap is not sameas DomainMap
      // FINISH: currently, domain map == range map == row map, so that we never have an exporter
      Teuchos::RCP<Export<Ordinal> > exporter_;

      // a map between a (non-local) row and a list of (col,val)
      // TODO: switch this to a hash-based map (instead of a sort-based map) as soon as one is available
      std::map<Ordinal, std::list<std::pair<Ordinal,Scalar> > > nonlocals_;

      bool fillCompleted_;

      Ordinal globalNNZ_;
      Ordinal myNNZ_;

      Ordinal numGlobalDiags_;
      Ordinal numMyDiags_;

      Ordinal globalMaxNumEntries_;
      Ordinal myMaxNumEntries_;

  }; // class CrsMatrix



  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  template<class Ordinal, class Scalar>
  CrsMatrix<Ordinal,Scalar>::CrsMatrix(const Map<Ordinal> &rowMap) 
  : comm_(rowMap.getComm())
  , rowMap_(rowMap)
  , numMyRows_(rowMap.getNumMyEntries())
  , numGlobalRows_(rowMap.getNumGlobalEntries())
  , colinds_(numMyRows_)
  ,  values_(numMyRows_)
  {
    fillCompleted_       = false;
    globalNNZ_           = OT::zero();
    myNNZ_               = OT::zero();
    numGlobalDiags_      = OT::zero();
    numMyDiags_          = OT::zero();
    globalMaxNumEntries_ = OT::zero();
    myMaxNumEntries_     = OT::zero();
  }


  template<class Ordinal, class Scalar>
  CrsMatrix<Ordinal,Scalar>::~CrsMatrix()
  {}


  template<class Ordinal, class Scalar>
  bool CrsMatrix<Ordinal,Scalar>::isFillCompleted() const
  { return fillCompleted_; }


  template<class Ordinal, class Scalar>
  Teuchos::RCP<const Teuchos::Comm<Ordinal> > 
  CrsMatrix<Ordinal,Scalar>::getComm() const
  { return comm_; }


  template<class Ordinal, class Scalar>
  Ordinal CrsMatrix<Ordinal,Scalar>::getNumGlobalNonzeros() const
  {
    TEST_FOR_EXCEPTION(isFillCompleted() == false, std::runtime_error,
      "Tpetra::CrsMatrix: cannot call getNumGlobalNonzeros() until fillComplete() has been called.");
    return globalNNZ_;
  }


  template<class Ordinal, class Scalar>
  Ordinal CrsMatrix<Ordinal,Scalar>::getNumMyNonzeros() const
  {
    TEST_FOR_EXCEPTION(isFillCompleted() == false, std::runtime_error,
      "Tpetra::CrsMatrix: cannot call getNumMyNonzeros() until fillComplete() has been called.");
    return myNNZ_;
  }


  template<class Ordinal, class Scalar>
  Ordinal CrsMatrix<Ordinal,Scalar>::getNumGlobalRows() const
  { return rowMap_.getNumGlobalEntries(); }


  template<class Ordinal, class Scalar>
  Ordinal CrsMatrix<Ordinal,Scalar>::getNumGlobalCols() const
  {
    TEST_FOR_EXCEPTION(isFillCompleted() == false, std::runtime_error,
      "Tpetra::CrsMatrix: cannot call getNumGlobalCols() until fillComplete() has been called.");
    return colMap_->getNumGlobalEntries();
  }


  template<class Ordinal, class Scalar>
  Ordinal CrsMatrix<Ordinal,Scalar>::getNumMyRows() const
  { return rowMap_.getNumMyEntries(); }


  template<class Ordinal, class Scalar>
  Ordinal CrsMatrix<Ordinal,Scalar>::getNumMyCols() const
  {
    TEST_FOR_EXCEPTION(isFillCompleted() == false, std::runtime_error,
      "Tpetra::CrsMatrix: cannot call getNumMyCols() until fillComplete() has been called.");
    return colMap_->getNumMyEntries();
  }


  template<class Ordinal, class Scalar>
  Ordinal CrsMatrix<Ordinal,Scalar>::getNumGlobalDiagonals() const
  {
    TEST_FOR_EXCEPTION(isFillCompleted() == false, std::runtime_error,
      "Tpetra::CrsMatrix: cannot call getNumGlobalDiagonals() until fillComplete() has been called.");
    return numGlobalDiags_;
  }


  template<class Ordinal, class Scalar>
  Ordinal CrsMatrix<Ordinal,Scalar>::getNumMyDiagonals() const
  {
    TEST_FOR_EXCEPTION(isFillCompleted() == false, std::runtime_error,
      "Tpetra::CrsMatrix: cannot call getNumMyDiagonals() until fillComplete() has been called.");
    return numMyDiags_;
  }


  template<class Ordinal, class Scalar>
  Ordinal CrsMatrix<Ordinal,Scalar>::getNumRowEntries(Ordinal globalRow) const
  {
    TEST_FOR_EXCEPTION(!rowMap_.isMyGlobalIndex(globalRow), std::runtime_error,
        "Tpetra::CrsMatrix::getNumRowEntries(globalRow): globalRow does not belong to this node.");
    return values_[rowMap_.getLocalIndex(globalRow)].size();
  }


  template<class Ordinal, class Scalar>
  Ordinal CrsMatrix<Ordinal,Scalar>::getGlobalMaxNumEntries() const
  {
    TEST_FOR_EXCEPTION(isFillCompleted() == false, std::runtime_error,
      "Tpetra::CrsMatrix: cannot call getGlobalMaxNumEntries() until fillComplete() has been called.");
    return globalMaxNumEntries_;
  }


  template<class Ordinal, class Scalar>
  Ordinal CrsMatrix<Ordinal,Scalar>::getMyMaxNumEntries() const
  {
    TEST_FOR_EXCEPTION(isFillCompleted() == false, std::runtime_error,
      "Tpetra::CrsMatrix: cannot call getMyMaxNumEntries() until fillComplete() has been called.");
    return myMaxNumEntries_;
  }


  template<class Ordinal, class Scalar>
  Ordinal CrsMatrix<Ordinal,Scalar>::getIndexBase() const
  { return rowMap_.getIndexBase(); }


  template<class Ordinal, class Scalar>
  const Map<Ordinal> & CrsMatrix<Ordinal,Scalar>::getRowMap() const 
  { return rowMap_; }


  template<class Ordinal, class Scalar>
  const Map<Ordinal> & CrsMatrix<Ordinal,Scalar>::getColMap() const 
  { 
    TEST_FOR_EXCEPTION(isFillCompleted() == false, std::runtime_error,
        "Tpetra::CrsMatrix: cannot call getColMap() until fillComplete() has been called.");
    return *colMap_; 
  }

  template<class Ordinal, class Scalar>
  const Map<Ordinal> & CrsMatrix<Ordinal,Scalar>::getDomainMap() const
  {
    return getColMap();
  }

  template<class Ordinal, class Scalar>
  const Map<Ordinal> & CrsMatrix<Ordinal,Scalar>::getRangeMap() const
  {
    return getRowMap();
  }

  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::submitEntry(Ordinal globalRow, Ordinal globalCol, const Scalar &value)
  {
    TEST_FOR_EXCEPTION(isFillCompleted() == true, std::runtime_error,
      "Tpetra::CrsMatrix::submitEntry(): fillComplete() has already been called.");
    if (rowMap_.isMyGlobalIndex(globalRow)) {
      Ordinal myRow = rowMap_.getLocalIndex(globalRow);
      colinds_[myRow].push_back(globalCol);
      values_[myRow].push_back(value);
    }
    else {
      nonlocals_[globalRow].push_back(std::make_pair(globalCol, value));
    }
  }


  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::submitEntries(Ordinal globalRow, 
                         const Teuchos::ArrayView<const Ordinal> &indices,
                         const Teuchos::ArrayView<const Scalar>  &values)
  {
    TEST_FOR_EXCEPTION(values.size() != indices.size(), std::runtime_error,
        "Tpetra::CrsMatrix::submitEntries(): values.size() must equal indices.size().");
    typename  Teuchos::ArrayView<const Scalar>::iterator val =  values.begin();
    typename Teuchos::ArrayView<const Ordinal>::iterator ind = indices.begin();
    for (; val != values.end(); ++val, ++ind) 
    {
      submitEntry(globalRow, *ind, *val);
    }
  }


  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::setAllToScalar(const Scalar &alpha)
  { 
    using Teuchos::Array;
    TEST_FOR_EXCEPTION(isFillCompleted() == false, std::runtime_error,
        "Tpetra::CrsMatrix: cannot call setAllToScalar() until fillComplete() has been called.");
    for (typename Array<Array<Scalar> >::const_iterator row = values_.begin(); row != values_.end(); ++row) 
    {
      std::fill(row->begin(),row->end(),alpha);
    }
  }


  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::scale(const Scalar &alpha)
  { 
    using Teuchos::Array;
    TEST_FOR_EXCEPTION(isFillCompleted() == false, std::runtime_error,
        "Tpetra::CrsMatrix: cannot call scale() until fillComplete() has been called.");
    for (typename Array<Array<Scalar> >::const_iterator row = values_.begin(); row != values_.end(); ++row) 
    {
      for (typename Array<Scalar>::const_iterator val = row->begin(); val != row->end(); ++val)
      {
        *val *= alpha;
      }
    }
  }


  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::getMyRowCopy(Ordinal myRow, 
                                               const Teuchos::ArrayView<Ordinal> &indices, 
                                               const Teuchos::ArrayView<Scalar>  &values) const
  {
    TEST_FOR_EXCEPTION(isFillCompleted() == false, std::runtime_error,
        "Tpetra::CrsMatrix: cannot call getMyRowCopy() until fillComplete() has been called.");
    TEST_FOR_EXCEPTION(indices.size() != colinds_[myRow].size() || values.size() != values_[myRow].size(), std::runtime_error, 
        "Tpetra::CrsMatrix::getMyRowCopy(indices,values): size of indices,values must be sufficient to store the values for the specified row.");
    std::copy( values_[myRow].begin(), values_[myRow].end(), values.begin());
    std::copy(colinds_[myRow].begin(),colinds_[myRow].end(),indices.begin());
  }


  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::getGlobalRowCopy(Ordinal globalRow, 
                                                   const Teuchos::ArrayView<Ordinal> &indices,
                                                   const Teuchos::ArrayView<Scalar>  &values) const
  {
    // TODO: is it possible to even call this one? we can't preallocate for the output buffer without knowing the size, 
    //       which (above) seems to require fillCompleted()
    // 
    // Only locally owned rows can be queried, otherwise complain
    TEST_FOR_EXCEPTION(!rowMap_.isMyGlobalIndex(globalRow), std::runtime_error,
        "Tpetra::CrsMatrix::getGlobalRowCOpy(globalRow): globalRow does not belong to this node.");
    Ordinal myRow = rowMap_.getLocalIndex(globalRow);
    TEST_FOR_EXCEPTION(
        Teuchos::as<typename Teuchos::Array<Ordinal>::size_type>(indices.size()) != colinds_[myRow].size() 
          || Teuchos::as<typename Teuchos::Array<Scalar>::size_type>(values.size()) != values_[myRow].size(), std::runtime_error, 
        "Tpetra::CrsMatrix::getGlobalRowCopy(indices,values): size of indices,values must be sufficient to store the values for the specified row.");
    std::copy(values_[myRow].begin(),values_[myRow].end(),values.begin());
    if (isFillCompleted())
    {
      // translate local IDs back to global IDs
      typename Teuchos::ArrayView<Ordinal>::iterator gind = indices.begin();
      for (typename Teuchos::Array<Ordinal>::const_iterator lind = colinds_[myRow].begin();
           lind != colinds_[myRow].end(); ++lind, ++gind) 
      {
        (*gind) = colMap_->getGlobalIndex(*lind);
      }
    }
    else
    {
      std::copy(colinds_[myRow].begin(),colinds_[myRow].end(),indices.begin());
    }
  }


  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::apply(const MultiVector<Ordinal,Scalar> &X, MultiVector<Ordinal,Scalar> &Y, Teuchos::ETransp mode) const
  {
    // TODO: add support for alpha,beta term coefficients: Y = alpha*A*X + beta*Y
    using Teuchos::null;
    using Teuchos::ArrayView;
    TEST_FOR_EXCEPTION(!isFillCompleted(), std::runtime_error, 
        "Tpetra::CrsMatrix: cannot call apply() until fillComplete() has been called.");
    TEST_FOR_EXCEPTION(X.numVectors() != Y.numVectors(), std::runtime_error,
        "Tpetra::CrsMatrix::apply(X,Y): X and Y must have the same number of vectors.");

    // FINISH: in debug mode, we need to establish that X and Y are compatible with (sameas?) DomainMap and RangeMap (respectively)
    // FINISH: or maybe not; I guess this should happen below in the calls to import(), export()

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

    Ordinal numVectors = X.numVectors();
    // because of Views, it is difficult to determine if X and Y point to the same data. 
    // however, if they reference the exact same object, we will do the user the favor of copying X into new storage (with a warning)
    // we ony need to do this if we have trivial importers; otherwise, we don't actually apply the operator from X into Y
    Teuchos::RCP<const MultiVector<Ordinal,Scalar> > Xcopy;
    ArrayView<const ArrayView<const Scalar> > Xdata = X.extractConstView();
    ArrayView<const ArrayView<      Scalar> > Ydata = Y.extractView();
    if (&X==&Y && importer_==null && exporter_==null) {
#     if defined(THROW_TPETRA_EFFICIENCY_WARNINGS) || defined(PRINT_TPETRA_EFFICIENCY_WARNINGS)
      std::string err;
      err += "Tpetra::CrsMatrix<" + Teuchos::TypeNameTraits<Ordinal>::name() + "," + Teuchos::TypeNameTraits<Ordinal>::name()
        + ">::apply(X,Y): If X and Y are the same, it necessitates a temporary copy of X, which is inefficient.";
#     if defined(THROW_TPETRA_EFFICIENCY_WARNINGS)
      TEST_FOR_EXCEPTION(true, std::runtime_error, err);
#     else
      std::cerr << err << std::endl;
#     endif
#     endif
      // generate a copy of X 
      Xcopy = Teuchos::rcp(new MultiVector<Ordinal,Scalar>(X));
      Xdata = Xcopy->extractConstView();
#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) *out << "X and Y are co-located, duplicating X results in a stride copy" << std::endl;
      *out << this->getColMap() << std::endl;
      Xcopy->print(*out); Xcopy->printValues(*out);
#   endif
    }
    if (importer_ != null) {
      if (importMV_ != null && importMV_->numVectors() != numVectors) importMV_ = null;
      if (importMV_ == null) {
        importMV_ = Teuchos::rcp( new MultiVector<Ordinal,Scalar>(*colMap_,numVectors) );
      }
    }
    if (exporter_ != null) {
      if (exportMV_ != null && exportMV_->numVectors() != numVectors) exportMV_ = null;
      if (exportMV_ == null) {
        exportMV_ = Teuchos::rcp( new MultiVector<Ordinal,Scalar>(rowMap_,numVectors) );
      }
    }
    // only support NO_TRANS currently
    TEST_FOR_EXCEPTION(mode != Teuchos::NO_TRANS, std::logic_error,
        "Tpetra::CrsMatrix::apply() does not currently support transposed multiplications.");
    if (mode == Teuchos::NO_TRANS) {
      // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
      if (importer_ != null) {
        importMV_->doImport(X, *importer_, INSERT);
        Xdata = importMV_->extractConstView();
#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
        if (myImageID == 0) {
          *out << "Performed import of X..." << std::endl;
        }
        importMV_->print(*out); importMV_->printValues(*out);
#   endif
      }
      // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
      // We will compute solution into the to-be-exported MV; get a view
      if (exporter_ != null) {
        Ydata = exportMV_->extractView();
      }
      // Do actual computation
      if (numVectors==Teuchos::OrdinalTraits<Ordinal>::one()) {
        GeneralMV(Xdata[0],Ydata[0]);
      }
      else {
        GeneralMM(Xdata,Ydata);
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
        Y.putScalar(0.0);  // Make sure target is zero: necessary because we are adding. may need adjusting for alpha,beta apply()
        Y.doExport(*exportMV_, *exporter_, ADD); // Fill Y with Values from export vector
#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
        if (myImageID == 0) *out << "Output vector after export()..." << std::endl;
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


  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::print(std::ostream& os) const 
  {
    using std::endl;
    int myImageID = Teuchos::rank(*comm_);
    if (myImageID == 0)
    {
      os << "Tpetra::CrsMatrix, label = " << this->label() << endl;
      os << "Number of global rows    = " << getNumGlobalRows() << endl;
      if (isFillCompleted())
      {
        os << "Number of global columns    = " << getNumGlobalCols() << endl;
        os << "Status = fillCompleted" << endl;
        os << "Number of global nonzeros   = " << getNumGlobalNonzeros() << endl;
        os << "Global max nonzeros per row = " << getGlobalMaxNumEntries() << endl;
      }
      else
      {
        os << "Status = not fillCompleted" << endl;
      }
    }
    if (isFillCompleted())
    {
      for (int pid=0; pid < Teuchos::size(*comm_); ++pid)
      {
        if (pid == myImageID)
        {
          Teuchos::Array<Ordinal> indices(getMyMaxNumEntries());
          Teuchos::Array<Scalar>   values(getMyMaxNumEntries());
          os << "% Number of rows on image " << myImageID << " = " << getNumMyRows() << endl;
          for (Ordinal i=OT::zero(); i < getNumMyRows(); ++i)
          {
            Ordinal globalRow = rowMap_.getGlobalIndex(i);
            Ordinal rowSize = getNumRowEntries(globalRow);
            if (rowSize > Teuchos::OrdinalTraits<Ordinal>::zero()) {
              getGlobalRowCopy(globalRow, indices(0,rowSize), values(0,rowSize));
              for (Ordinal j=OT::zero(); j < rowSize; ++j) {
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
  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::fillComplete()
  {
    using Teuchos::Array;
    TEST_FOR_EXCEPTION(isFillCompleted() == true, std::runtime_error,
      "Tpetra::CrsMatrix::fillComplete(): fillComplete() has already been called.");

    Teuchos::ArrayView<const Ordinal> myGlobalEntries = rowMap_.getMyGlobalEntries();

    // =============================== //
    // Part 0: send off-image elements //
    // =============================== //
    if (comm_->getSize() > 1) {
      globalAssemble();
    }
    else {
      TEST_FOR_EXCEPTION(nonlocals_.size() > 0, std::runtime_error,
          "Tpetra::CrsMatrix::fillComplete(): cannot have non-local entries on a serial run. Invalid entry was submitted to the CrsMatrix.");
    }

    // =============================== //
    // Part I: remove repeated indices //
    // =============================== //
    //
    // I load all matrix entries in a hash table, then I re-fill
    // the row with the last inserted value.
    // This also has the effect of sorting the indices.
    for (Ordinal r = OT::zero(); r < numMyRows_ ; ++r)
    {
      // put them in the map from colinds_,values_
      std::map<Ordinal,Scalar> row;
      {
        typename Array<Ordinal>::const_iterator cind = colinds_[r].begin();
        typename  Array<Scalar>::const_iterator val  = values_[r].begin();
        for (; cind != colinds_[r].end(); ++cind, ++val)
        {
          row[*cind] += *val;
        }
      }
      // get them out of the map, back to colinds_,values_
      typename std::map<Ordinal,Scalar>::size_type count = 0;
      for (typename std::map<Ordinal,Scalar>::iterator iter = row.begin(); 
           iter != row.end(); ++iter)
      {
        colinds_[r][count] = iter->first;
        values_[r][count] = iter->second;
        ++count;
      }
      colinds_[r].resize(count);
      values_[r].resize(count);
    }

    // =============================== //
    // Part II: build the column space //
    // =============================== //
    // construct a list of columns for which we have a non-zero, but are not locally owned by our row map
    std::set<Ordinal> nnzcols;
    for (Ordinal r=OT::zero(); r < numMyRows_ ; ++r)
    {
      for (typename Array<Ordinal>::const_iterator cind = colinds_[r].begin();
           cind != colinds_[r].end(); ++cind)
      {
        if (!rowMap_.isMyGlobalIndex(*cind))
        {
          nnzcols.insert(*cind);
        }
      }
    }
    Array<Ordinal> myPaddedGlobalEntries(myGlobalEntries.begin(),myGlobalEntries.end());
    for (typename std::set<Ordinal>::iterator iter = nnzcols.begin(); iter != nnzcols.end() ; ++iter)
    {
      myPaddedGlobalEntries.push_back(*iter);
    }
    colMap_ = Teuchos::rcp(new Map<Ordinal>(Teuchos::as<Ordinal>(-1), myPaddedGlobalEntries(),
                                            rowMap_.getIndexBase(), *rowMap_.getPlatform()) );
    importer_ = Teuchos::rcp( new Import<Ordinal>(rowMap_,*colMap_) );

    // ============================== //
    // Part IV: move to local indices //
    // ============================== //
    for (Ordinal r=OT::zero(); r < numMyRows_ ; ++r)
    {
      for (typename Array<Ordinal>::iterator cind = colinds_[r].begin();
           cind != colinds_[r].end(); ++cind) 
      {
        (*cind) = colMap_->getLocalIndex(*cind);
      }
    }

    // ================================================ //
    // Part V: compute some local and global quantities //
    // ================================================ //
    for (Ordinal r=OT::zero(); r < numMyRows_; ++r)
    {
      Ordinal numEntries = colinds_[r].size();
      myNNZ_ += numEntries;
      if ( std::find(colinds_[r].begin(),colinds_[r].end(),r) != colinds_[r].end() )
      {
        ++numMyDiags_;
      }
      myMaxNumEntries_ = std::max(myMaxNumEntries_,numEntries);
    }
    // TODO: consider eliminating one of these comms in favor of a pair
    Teuchos::reduceAll<Ordinal>(*comm_,Teuchos::REDUCE_SUM,myNNZ_,&globalNNZ_);
    Teuchos::reduceAll<Ordinal>(*comm_,Teuchos::REDUCE_SUM,numMyDiags_,&numGlobalDiags_);
    Teuchos::reduceAll<Ordinal>(*comm_,Teuchos::REDUCE_MAX,myMaxNumEntries_,&globalMaxNumEntries_);

    // mark transformation as successfully completed
    fillCompleted_ = true;
  }



  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::globalAssemble()
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
    typedef typename std::map<Ordinal,std::list<pair<Ordinal,Scalar> > >::const_iterator NLITER;
    int numImages = comm_->getSize();
    int myImageID = comm_->getRank();
    // Determine if any nodes have global entries to share
    Ordinal MyNonlocals = nonlocals_.size(),
            MaxGlobalNonlocals;
    Teuchos::reduceAll<Ordinal>(*comm_,Teuchos::REDUCE_MAX,MyNonlocals,&MaxGlobalNonlocals);
    if (MaxGlobalNonlocals == OT::zero()) return;  // no entries to share

    // compute a list of NLRs from nonlocals_ and use it to compute:
    //      IdsAndRows: a vector of (id,row) pairs
    //          NLR2Id: a map from NLR to the Id that owns it
    // globalNeighbors: a global graph of connectivity between images: globalNeighbors(i,j) indicates that j sends to i
    //         sendIDs: a list of all images I send to
    //         recvIDs: a list of all images I receive from (constructed later)
    Array<pair<Ordinal,Ordinal> > IdsAndRows;
    std::map<Ordinal,Ordinal> NLR2Id;
    SerialDenseMatrix<Ordinal,char> globalNeighbors;
    Array<Ordinal> sendIDs, recvIDs;
    {
      // nonlocals_ contains the entries we are holding for all non-local rows
      // we want a list of the rows for which we have data
      Array<Ordinal> NLRs;
      std::set<Ordinal> setOfRows;
      for (NLITER iter = nonlocals_.begin(); iter != nonlocals_.end(); ++iter)
      {
        setOfRows.insert(iter->first);
      }
      // copy the elements in the set into an Array
      NLRs.resize(setOfRows.size());
      std::copy(setOfRows.begin(), setOfRows.end(), NLRs.begin());

      // get a list of ImageIDs for the non-local rows (NLRs)
      Array<Ordinal> NLRIds(NLRs.size());
      {
        bool invalidGIDs = rowMap_.getRemoteIndexList(NLRs(),NLRIds());
        char lclerror = ( invalidGIDs ? OrdinalTraits<char>::one() : OrdinalTraits<char>::zero() );
        char gblerror;
        Teuchos::reduceAll(*comm_,Teuchos::REDUCE_MAX,lclerror,&gblerror);
        TEST_FOR_EXCEPTION(gblerror, std::runtime_error,
            "Tpetra::CrsMatrix::globalAssemble(): non-local entries correspond to invalid rows.");
      }

      // build up a list of neighbors, as well as a map between NLRs and Ids
      // localNeighbors[i] != 0 iff I have data to send to image i
      // put NLRs,Ids into an array of pairs
      IdsAndRows.reserve(NLRs.size());
      Array<char> localNeighbors(numImages,OrdinalTraits<char>::zero());
      for (typename Array<Ordinal>::const_iterator nlr = NLRs.begin(), id = NLRIds.begin();
          nlr != NLRs.end(); ++nlr, ++id) 
      {
        NLR2Id[*nlr] = *id;
        localNeighbors[*id] = OrdinalTraits<char>::one();
        IdsAndRows.push_back(make_pair(*id,*nlr));
      }
      for (Ordinal j=OT::zero(); j<numImages; ++j)
      {
        if (localNeighbors[j]) {
          sendIDs.push_back(j);
        }
      }
      // sort IdsAndRows, by Ids first, then rows
      std::sort(IdsAndRows.begin(),IdsAndRows.end());
      // gather from other nodes to form the full graph
      globalNeighbors.shapeUninitialized(numImages,numImages);
      Teuchos::gatherAll<Ordinal>(*comm_,numImages,localNeighbors.getRawPtr(),numImages*numImages,globalNeighbors.values());
      // globalNeighbors at this point contains (on all images) the
      // connectivity between the images. 
      // globalNeighbors(i,j) != 0 means that j sends to i/that i receives from j
    }

    ////////////////////////////////////////////////////////////////////////////////////// 
    // FIGURE OUT WHO IS SENDING TO WHOM AND HOW MUCH
    // DO THIS IN THE PROCESS OF PACKING ALL OUTGOING DATA ACCORDING TO DESTINATION ID
    ////////////////////////////////////////////////////////////////////////////////////// 

    // loop over all columns to know from which images I can expect to receive something
    for (Ordinal j=OT::zero(); j<numImages; ++j)
    {
      if (globalNeighbors(myImageID,j)) {
        recvIDs.push_back(j);
      }
    }
    Ordinal numRecvs = recvIDs.size();

    // we know how many we're sending to already
    // form a contiguous list of all data to be sent
    // track the number of entries for each ID
    Array<CrsIJV<Ordinal,Scalar> > IJVSendBuffer;
    Array<Ordinal> sendSizes(sendIDs.size(), 0);
    Ordinal numSends = 0;
    for (typename Array<pair<Ordinal,Ordinal> >::const_iterator IdAndRow = IdsAndRows.begin();
         IdAndRow != IdsAndRows.end(); ++IdAndRow) 
    {
      Ordinal  id = IdAndRow->first;
      Ordinal row = IdAndRow->second;
      // have we advanced to a new send?
      if (sendIDs[numSends] != id) {
        numSends++;
        TEST_FOR_EXCEPTION(sendIDs[numSends] != id, std::logic_error, "Tpetra::CrsMatrix::globalAssemble(): internal logic error. Contact Tpetra team.");
      }
      // copy data for row into contiguous storage
      for (typename list<pair<Ordinal,Scalar> >::const_iterator jv = nonlocals_[row].begin(); jv != nonlocals_[row].end(); ++jv)
      {
        IJVSendBuffer.push_back( CrsIJV<Ordinal,Scalar>(row,jv->first,jv->second) );
        sendSizes[numSends]++;
      }
    }
    if (IdsAndRows.size() > 0) {
      numSends++; // one last increment, to make it a count instead of an index
    }
    TEST_FOR_EXCEPTION(Teuchos::as<typename Array<Ordinal>::size_type>(numSends) != sendIDs.size(), std::logic_error, "Tpetra::CrsMatrix::globalAssemble(): internal logic error. Contact Tpetra team.");

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
      sendRequests.push_back( Teuchos::isend<Ordinal,Ordinal>(*comm_,rcp(&sendSizes[s],false),sendIDs[s]) );
    }
    // perform non-blocking receives: receive sizes from our senders
    Array<Teuchos::RCP<Teuchos::CommRequest> > recvRequests;
    Array<Ordinal> recvSizes(numRecvs);
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
    Array<ArrayView<CrsIJV<Ordinal,Scalar> > > sendBuffers(numSends,Teuchos::null);
    {
      Ordinal cur = 0;
      for (int s=0; s<numSends; ++s) {
        sendBuffers[s] = IJVSendBuffer(cur,sendSizes[s]);
        cur += sendSizes[s];
      }
    }
    // perform non-blocking sends
    for (int s=0; s < numSends ; ++s)
    {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      ArrayRCP<CrsIJV<Ordinal,Scalar> > tmparcp = arcp(sendBuffers[s].getRawPtr(),OT::zero(),sendBuffers[s].size(),false);
      sendRequests.push_back( Teuchos::isend<Ordinal,CrsIJV<Ordinal,Scalar> >(*comm_,tmparcp,sendIDs[s]) );
    }
    // calculate amount of storage needed for receives
    // setup pointers for the receives as well
    Ordinal totalRecvSize = std::accumulate(recvSizes.begin(),recvSizes.end(),OT::zero());
    Array<CrsIJV<Ordinal,Scalar> > IJVRecvBuffer(totalRecvSize);
    // from the size info, build the ArrayViews into IJVRecvBuffer
    Array<ArrayView<CrsIJV<Ordinal,Scalar> > > recvBuffers(numRecvs,Teuchos::null);
    {
      Ordinal cur = 0;
      for (int r=0; r<numRecvs; ++r) {
        recvBuffers[r] = IJVRecvBuffer(cur,recvSizes[r]);
        cur += recvSizes[r];
      }
    }
    // perform non-blocking recvs
    for (int r=0; r < numRecvs ; ++r)
    {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      ArrayRCP<CrsIJV<Ordinal,Scalar> > tmparcp = arcp(recvBuffers[r].getRawPtr(),OT::zero(),recvBuffers[r].size(),false);
      recvRequests.push_back( Teuchos::ireceive<Ordinal,CrsIJV<Ordinal,Scalar> >(*comm_,tmparcp,recvIDs[r]) );
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
    for (typename Array<CrsIJV<Ordinal,Scalar> >::const_iterator ijv = IJVRecvBuffer.begin(); ijv != IJVRecvBuffer.end(); ++ijv)
    {
      submitEntry(ijv->i, ijv->j, ijv->v);
    }

    // WHEW! THAT WAS TIRING!
  }


  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::GeneralMV(const Teuchos::ArrayView<const Scalar> &x, const Teuchos::ArrayView<Scalar> &y) const
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    typename Teuchos::Array<Ordinal>::const_iterator col;
    typename Teuchos::Array<Scalar >::const_iterator val;
    for (Ordinal r=0; r < numMyRows_; ++r) {
      col = colinds_[r].begin(); 
      val =  values_[r].begin(); 
      Scalar sum = ST::zero();
      for (; col != colinds_[r].end(); ++col, ++val) {
        sum += (*val) * x[*col];
      }
      y[r] = sum;
    }
  }


  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::GeneralMM(const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> > &X, const Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> > &Y) const
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    Ordinal numVectors = X.size();
    typename Teuchos::Array<Ordinal>::const_iterator col;
    typename Teuchos::Array<Scalar >::const_iterator val;
    for (Ordinal r=0; r < numMyRows_; ++r) {
      for (int j=0; j<numVectors; ++j) {
        const Teuchos::ArrayView<const Scalar> &xvals = X[j]; 
        col = colinds_[r].begin(); 
        val = values_[r].begin(); 
        Scalar sum = ST::zero();
        for (; col != colinds_[r].end(); ++col, ++val) {
          sum += (*val) * xvals[*col];
        }
        Y[j][r] = sum;
      }
    }
  }


} // namespace Tpetra

#endif
