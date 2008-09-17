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

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include "Tpetra_Operator.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"

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
      CrsMatrix(const Teuchos::RCP<const Teuchos::Comm<Ordinal> > & comm, Map<Ordinal> &rowMap);

      // !Destructor.
      virtual ~CrsMatrix();

      //@}
      //! @name Query Methods
      //@{ 
      
      //! Returns \c true if the matrix has already been fill completed.
      bool isFillCompleted() const;

      //! Returns the communicator.
      Teuchos::RCP<const Teuchos::Comm<Ordinal> > getComm() const;

      //! Returns the number of nonzero entries in the global matrix. 
      Ordinal getNumGlobalNonzeros() const;

      //! Returns the number of nonzero entries in the calling image's portion of the matrix. 
      Ordinal getNumMyNonzeros() const;

      //! Returns the number of global matrix rows. 
      Ordinal getNumGlobalRows() const;

      //! Returns the number of global matrix columns. 
      Ordinal getNumGlobalCols() const;

      //! Returns the number of matrix rows owned by the calling image. 
      Ordinal getNumMyRows() const;

      //! Returns the number of matrix columns owned by the calling image. 
      Ordinal getNumMyCols() const;

      //! Returns the number of global nonzero diagonal entries, based on global row/column index comparisons. 
      Ordinal getNumGlobalDiagonals() const;

      //! Returns the number of local nonzero diagonal entries, based on global row/column index comparisons. 
      Ordinal getNumMyDiagonals() const;

      //! Returns the current number of nonzero entries in specified global index on this image. 
      Ordinal getNumRowEntries(Ordinal globalRow) const;

      //! Returns the maximum number of nonzero entries across all rows/columns on all images. 
      Ordinal getGlobalMaxNumEntries() const;

      //! Returns the maximum number of nonzero entries across all rows/columns on this image. 
      Ordinal getMyMaxNumEntries() const;

      //! Returns the index base for global indices for this matrix. 
      Ordinal getIndexBase() const;
      
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
      void apply(const MultiVector<Ordinal,Scalar>& X, MultiVector<Ordinal,Scalar> &Y,
                 Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

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
                         const Teuchos::ArrayView<const Scalar> &values, 
                         const Teuchos::ArrayView<const Ordinal> &indices);

      // @}
      //! @name Computational Methods
      // @{ 

      //! Returns a copy of the specified local row, column indices are local.
      void getMyRowCopy(Ordinal localRow, const Teuchos::ArrayView<Ordinal> &indices, 
                                          const Teuchos::ArrayView<Scalar> &values) const;

      //! Returns a copy of the specified (and locally owned) global row, column indices are global.
      void getGlobalRowCopy(Ordinal globalRow, 
                            const Teuchos::ArrayView<Ordinal> &indices,
                            const Teuchos::ArrayView<Scalar> &values) const;

      //@}
      //! @name I/O Methods
      //@{ 
      
      //! Prints the matrix on the specified stream.
      void print(std::ostream& os) const;

      //! Set all matrix entries equal to scalarThis.
      void setAllToScalar(Scalar scalarThis);

      //! Scale the current values of a matrix, this = scalarThis*this. 
      void scale(Scalar scalarThis);

      //! Basic print, for debugging purposes only.
      void printValues();

      // @}

    private:
      //! copy constructor.
      CrsMatrix(const CrsMatrix<Ordinal,Scalar> &Source);

      CrsMatrix& operator=(const CrsMatrix<Ordinal, Scalar> &rhs);

      // struct for i,j,v
      struct CrsIJV {
        CrsIJV(Ordinal row, Ordinal col, const Scalar &val) {
          i = row;
          j = col;
          v = val;
        }
        Ordinal i,j;
        Scalar  v;
      };

      //! Performs importing of off-processor elements and adds them to the locally owned elements.
      void globalAssemble();

      const Teuchos::RCP<const Teuchos::Comm<Ordinal> > comm_;
      const Map<Ordinal>& rowMap_;
      Teuchos::RCP<Map<Ordinal> > colMap_;

      Ordinal numMyRows_, numGlobalRows_;
      const std::vector<Ordinal>& myGlobalEntries_;

      std::vector<std::vector<Ordinal> > indices_;
      std::vector<std::vector<Scalar> > values_;

      MultiVector<Ordinal, Scalar>* paddedMV_;
      Import<Ordinal>* importer_;

      // a map between a (non-local) row and a list of (col,val)
      // TODO: switch this to a hash-based map (instead of a sort-based map) as soon as one is available
      std::map<Ordinal, std::list<std::pair<Ordinal,Scalar> > > nonlocals_;

      bool fillCompleted_;

      Ordinal numGlobalNonzeros_;
      Ordinal numMyNonzeros_;

      Ordinal numGlobalDiagonals_;
      Ordinal numMyDiagonals_;

      Ordinal globalMaxNumEntries_;
      Ordinal myMaxNumEntries_;

  }; // class CrsMatrix



  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  template<class Ordinal, class Scalar>
  CrsMatrix<Ordinal,Scalar>::CrsMatrix(
      const Teuchos::RCP<const Teuchos::Comm<Ordinal> > & comm, 
            Map<Ordinal> &rowMap) 
  : comm_(comm)
  , rowMap_(rowMap)
  , numMyRows_(rowMap.getNumMyEntries())
  , numGlobalRows_(rowMap.getNumGlobalEntries())
  , myGlobalEntries_(rowMap.getMyGlobalEntries())
  {
    const Ordinal ZERO = Teuchos::OrdinalTraits<Ordinal>::zero();
    indices_.resize(numMyRows_);
    values_.resize(numMyRows_);
    fillCompleted_ = false;
    numGlobalNonzeros_ = ZERO;
    numMyNonzeros_     = ZERO;
    numGlobalDiagonals_ = ZERO;
    numMyDiagonals_     = ZERO;
    globalMaxNumEntries_ = ZERO;
    myMaxNumEntries_     = ZERO;
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
    return numGlobalNonzeros_;
  }


  template<class Ordinal, class Scalar>
  Ordinal CrsMatrix<Ordinal,Scalar>::getNumMyNonzeros() const
  {
    TEST_FOR_EXCEPTION(isFillCompleted() == false, std::runtime_error,
      "Tpetra::CrsMatrix: cannot call getNumMyNonzeros() until fillComplete() has been called.");
    return numMyNonzeros_;
  }


  template<class Ordinal, class Scalar>
  Ordinal CrsMatrix<Ordinal,Scalar>::getNumGlobalRows() const
  { return rowMap_.getNumGlobalEntries(); }


  template<class Ordinal, class Scalar>
  Ordinal CrsMatrix<Ordinal,Scalar>::getNumGlobalCols() const
  {
    TEST_FOR_EXCEPTION(isFillCompleted() == false, std::runtime_error,
      "Tpetra::CrsMatrix: cannot call getNumGlobalCols() until fillComplete() has been called.");
    return colMap_->getGlobalEntries();
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
    return numGlobalDiagonals_;
  }


  template<class Ordinal, class Scalar>
  Ordinal CrsMatrix<Ordinal,Scalar>::getNumMyDiagonals() const
  {
    TEST_FOR_EXCEPTION(isFillCompleted() == false, std::runtime_error,
      "Tpetra::CrsMatrix: cannot call getNumMyDiagonals() until fillComplete() has been called.");
    return numMyDiagonals_;
  }


  template<class Ordinal, class Scalar>
  Ordinal CrsMatrix<Ordinal,Scalar>::getNumRowEntries(Ordinal globalRow) const
  {
    TEST_FOR_EXCEPTION(isFillCompleted() == false, std::runtime_error,
      "Tpetra::CrsMatrix: cannot call getNumRowEntries() until fillComplete() has been called.");
    return indices_[rowMap_.getLID(globalRow)].size();
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
  { return colMap_; }


  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::fillComplete()
  {
    TEST_FOR_EXCEPTION(isFillCompleted() == true, std::runtime_error,
      "Tpetra::CrsMatrix::fillComplete(): fillComplete() has already been called.");

    // =============================== //
    // Part 0: send off-image elements //
    // =============================== //

    if (getComm().getSize() > 1) globalAssemble();

    // =============================== //
    // Part I: remove repeated indices //
    // =============================== //

    // I load all matrix entries in a hash table, then I re-fill
    // the row with the last inserted value.
    for (Ordinal i = ordinalZero() ; i < numMyRows_ ; ++i)
    {
      std::map<Ordinal, Scalar> singleRow;

      for (Ordinal j = ordinalZero() ; j < indices_[i].size() ; ++j)
      {
        singleRow[indices_[i][j]] += values_[i][j];
      }

      Ordinal count = ordinalZero();

      for (typename std::map<Ordinal,Scalar>::iterator iter = singleRow.begin() ; 
          iter != singleRow.end() ; ++iter)
      {
        indices_[i][count] = iter->first;
        values_[i][count] = iter->second;
        ++count;
      }

      indices_[i].resize(count);
      values_[i].resize(count);
    }

    // =============================== //
    // Part II: build the column space //
    // =============================== //

    // I have to find the list of non-locally owned columns

    std::map<Ordinal, bool> container; // replace with a hash table

    for (Ordinal i = ordinalZero() ; i < numMyRows_ ; ++i)
    {
      for (Ordinal j = ordinalZero() ; j < indices_[i].size() ; ++j)
      {
        Ordinal what = indices_[i][j];
        if (rowMap_.isMyGID(what)) continue;
        else
          container[what] = true;
      }
    }

    vector<Ordinal> MyPaddedGlobalEntries(myGlobalEntries_);

    for (typename std::map<Ordinal, bool>::iterator iter = container.begin() ; 
        iter != container.end() ; ++iter)
    {
      MyPaddedGlobalEntries.push_back(iter->first);
    }

    // now I can build the column space

    colMap_ = new ElementSpace<Ordinal>(-1, MyPaddedGlobalEntries.size(),
        MyPaddedGlobalEntries, rowMap_.getIndexBase(), rowMap_.platform());

    VectorColSpace_ = new VectorSpace<Ordinal, Scalar>(*colMap_, rowMap_.platform());

    paddedMV_ = new Vector<Ordinal, Scalar>(*VectorColSpace_);
    importer_ = new Import<Ordinal>(rowMap_, *colMap_);

    // ============================== //
    // Part IV: move to local indices //
    // ============================== //

    for (Ordinal i = ordinalZero() ; i < numMyRows_ ; ++i)
    {
      for (Ordinal j = ordinalZero() ; j < indices_[i].size() ; ++j)
      {
        indices_[i][j] = colMap_->getLID(indices_[i][j]);
      }
    }

    // ================================================ //
    // Part V: compute some local and global quantities //
    // ================================================ //

    for (Ordinal i = ordinalZero() ; i < numMyRows_ ; ++i)
    {
      Ordinal NumEntries = indices_[i].size();

      numMyNonzeros_ += NumEntries;

      for (Ordinal j = ordinalZero() ; j < NumEntries; ++j)
      {
        Ordinal col = indices_[i][j];
        if (col == i)
        {
          ++numMyDiagonals_;
          break;
        }
      }

      if (myMaxNumEntries_ < NumEntries) myMaxNumEntries_ = NumEntries;
    }

    // workaround, should be fixed in Comm
#ifdef HAVE_MPI_THIS_IS_BROKEN // FINISH
    MPI_Comm MpiCommunicator = (dynamic_cast<const MpiComm<Ordinal, Scalar>&>(getComm())).getMpiComm();
    MPI_Allreduce((void*)&numMyNonzeros_, (void*)&numGlobalNonzeros_, MpiTraits<Ordinal>::count(1),
        MpiTraits<Ordinal>::datatype(), MPI_SUM, MpiCommunicator);

    MPI_Allreduce((void*)&numMyDiagonals_, (void*)&numGlobalDiagonals_, MpiTraits<Ordinal>::count(1),
        MpiTraits<Ordinal>::datatype(), MPI_SUM, MpiCommunicator);

    MPI_Allreduce((void*)&myMaxNumEntries_, (void*)&globalMaxNumEntries_, MpiTraits<Ordinal>::count(1),
        MpiTraits<Ordinal>::datatype(), MPI_MAX, MpiCommunicator);
#else
    numGlobalNonzeros_ = numMyNonzeros_;
    numGlobalDiagonals_ = numMyDiagonals_;
    globalMaxNumEntries_ = myMaxNumEntries_;
#endif

    // mark transformation as successfully completed

    fillCompleted_ = true;
  }


  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::submitEntry(Ordinal globalRow, Ordinal globalCol,
                                              const Scalar &value)
  {
    TEST_FOR_EXCEPTION(isFillCompleted() == true, std::runtime_error,
      "Tpetra::CrsMatrix::submitEntry(): fillComplete() has already been called.");
    if (rowMap_.isMyGID(globalRow)) {
      Ordinal myRow = rowMap_.getLID(globalRow);
      indices_[myRow].push_back(globalCol);
      values_[myRow].push_back(value);
    }
    else {
      nonlocals_[globalRow].push_back(make_pair(globalCol, value));
    }
  }


  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::submitEntries(Ordinal globalRow, 
                         const Teuchos::ArrayView<const Scalar> &values, 
                         const Teuchos::ArrayView<const Ordinal> &indices)
  {
    TEST_FOR_EXCEPTION(values.size() != indices.size(), std::runtime_error,
        "Tpetra::CrsMatrix::submitEntries(): values.size() must equal indices.size().");
    typename Teuchos::ArrayView<Scalar>::iterator val = values.begin();
    typename Teuchos::ArrayView<Ordinal>::iterator ind = indices.begin();
    for (val != values.end(); ++val, ++ind) 
    {
      submitEntry(globalRow, *ind, *val);
    }
  }


  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::getMyRowCopy(Ordinal localRow, 
                                               const Teuchos::ArrayView<Ordinal> &indices, 
                                               const Teuchos::ArrayView<Scalar> &values) const;
  {
    TEST_FOR_EXCEPTION(isFillCompleted() == false, std::runtime_error,
        "Tpetra::CrsMatrix: cannot call getMyRowCopy() until fillComplete() has been called.");
    TEST_FOR_EXCEPTION(indices.size() < indices_[localRow].size() || values.size() < indices_[localRow].size(), std::runtime_error, 
        "Tpetra::CrsMatrix::getMyRowCopy(indices,values): size of indices,values must be sufficient to store the values for the specified row.");

    Ordinal NumEntries = indices_[MyRow].size();

    for (Ordinal i = ordinalZero() ; i < NumEntries ; ++i)
    {
      Indices[i] = indices_[MyRow][i];
      Values[i] = values_[MyRow][i];
    }
  }


  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::getGlobalRowCopy(Ordinal globalRow, 
                                                   const Teuchos::ArrayView<Ordinal> &indices,
                                                   const Teuchos::ArrayView<Scalar> &values) const
  {
    // Only locally owned rows can be queried, otherwise complain
    TEST_FOR_EXCEPTION(!rowMap_.isMyGID(globalRow), std::runtime_error,
        "Tpetra::CrsMatrix::getGlobalRowCOpy(globalRow): globalRow does not belong to this node.");

    Ordinal myRow = rowMap_.getLID(GlobalRow);

    Ordinal NumEntries = indices_[myRow].size();
    TEST_FOR_EXCEPTION(indices.size() < NumEntries || values.size() < NumEntries, std::runtime_error, 
        "Tpetra::CrsMatrix::getGlobalRowCopy(indices,values): size of indices,values must be sufficient to store the values for the specified row.");

    if (isFillCompleted())
    {
      for (Ordinal i = ordinalZero() ; i < NumEntries ; ++i)
      {
        Indices[i] = colMap_->getGID(indices_[myRow][i]);
        Values[i] = values_[myRow][i];
      }
    }
    else
    {
      for (Ordinal i = ordinalZero() ; i < NumEntries ; ++i)
      {
        Indices[i] = indices_[MyRow][i];
        Values[i] = values_[MyRow][i];
      }
    }
  }


  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::apply(const MultiVector<Ordinal, Scalar> &X, MultiVector<Ordinal,Scalar> &Y,
                                        Teuchos::ETransp mode = Teuchos::NO_TRANS) const
  {
    assert (Mode == AsIs);
    y.setAllToScalar(scalarZero());
    paddedMV_->doImport(x, *importer_, Insert);
    for (Ordinal i = ordinalZero() ; i < numMyRows_ ; ++i)
    {
      for (Ordinal j = ordinalZero() ; j < indices_[i].size() ; ++j)
      {
        Ordinal col = indices_[i][j];
        Scalar val = values_[i][j];
        y[i] += val * (*paddedMV_)[col];
      }
    }
  }


  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::print(std::ostream& os) const 
  {
    int MyImageID = rowMap_.comm().getRank();

    if (MyImageID == 0)
    {
      os << "Tpetra::CrsMatrix, label = " << label() << endl;
      os << "Number of global rows    = " << rowMap_.getNumGlobalEntries() << endl;
      if (isFillCompleted())
      {
        os << "Number of global columns = " << colMap_->getNumGlobalEntries() << endl;
        os << "Status = FillCompleted" << endl;
        os << "MyMaxNumEntries = " << getMyMaxNumEntries() << endl;
        os << "GlobalMaxNumEntries = " << getGlobalMaxNumEntries() << endl;
      }
      else
      {
        os << "Status = not FillCompleted" << endl;
      }
    }

    if (isFillCompleted())
    {
      for (int pid = 0 ; pid < rowMap_.comm().getSize() ; ++pid)
      {
        if (pid == MyImageID)
        {
          vector<Ordinal> Indices(getMyMaxNumEntries()); // FIXME
          vector<Scalar>  Values(getMyMaxNumEntries());
          Ordinal NumEntries;

          os << "% Number of rows on image " << MyImageID << " = " << rowMap_.getNumMyEntries() << endl;
          for (Ordinal i = ordinalZero() ; i < rowMap_.getNumMyEntries() ; ++i)
          {
            Ordinal GlobalRow = rowMap_.getGID(i);

            getGlobalRowCopy(GlobalRow, NumEntries, Indices, Values);
            for (Ordinal j = ordinalZero() ; j < NumEntries ; ++j)
              os << "Matrix(" << GlobalRow << ", " << Indices[j] << ") = " << Values[j] << ";" << endl;
          }
        }
        rowMap_.comm().barrier();
      }
    }
  }


  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::setAllToScalar(Scalar scalarThis)
  { TEST_FOR_EXCEPT(true); }


  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::scale(Scalar scalarThis)
  { TEST_FOR_EXCEPT(true); }


  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::printValues()
  {
    // this prints out the structure as they are
    for (int i = 0 ; i < indices_.size() ; ++i)
    {
      for (int j = 0 ; j < indices_[i].size() ; ++j)
      {
        cout << "local row " << i << ", col " << indices_[i][j] << ", val " << values_[i][j] << endl;
      }
    }
  }


  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::globalAssemble()
  {
    using Teuchos::OrdinalTraits;
    using Teuchos::Array;
    using Teuchos::SerialDenseMatrix;
    using Teuchos::ArrayView;
    using std::pair;
    using std::make_pair;
    typedef OrdinalTraits<Ordinal> OT;
    typedef typename std::map<Ordinal,std::list<pair<Ordinal,Scalar> > >::const_iterator NLITER;
    int numImages = comm_->getSize();
    int myImageID = comm_->getRank();
    // Determine if any nodes have global entries to share
    Ordinal MyNonlocals = nonlocals_.size(),
            MaxGlobalNonlocals;
    Teuchos::reduceAll<Ordinal>(*comm_,TEUCHOS::REDUCE_MAX,MyNonlocals,MaxGlobalNonlocals);
    if (MaxGlobalNonlocals == OT::OrdinalTraits<Ordinal>::zero()) return;  // no entries to share

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
      rowMap_.getRemoteIndexList(NLRs(),NLRIds());
      // build up a list of neighbors, as well as a map between NLRs and Ids
      // localNeighbors[i] != 0 iff I have data to send to image i
      // put NLRs,Ids into an array of pairs
      IdsAndRows.reserve(NLRs.size());
      Array<char> localNeighbors(numImages,OrdinalTraits<char>::zero());
      for (Array<Ordinal>::const_iterator nlr = NLRs.begin(), id = NLRIds.begin();
          nlr != NLRs.end(); ++nlr, ++id) 
      {
        NLR2Id[*nlr] = *id;
        localNeighbors[*id] = OrdinalTraits<char>::one();
        IdsAndRows.push_back(make_pair(*id,*nlr));
      }
      for (Ordinal j=OT::zero(); j<numImages; ++j) {
      {
        if (localNeighbors[j]) {
          sendIDs.push_back(j);
        }
      }
      // sort IdsAndRows, by Ids first, then rows
      std::sort(IdsAndRows.begin(),IdsAndRows.end());
      // gather from other nodes to form the full graph
      globalNeighbors.shapeUninitialized(numImages,numImages);
      Teuchos::gatherAll(*comm_,numImages,localNeighbors.getRawPtr(),numImages,globalNeighbors.values());
      // globalNeighbors at this point contains (on all images) the
      // connectivity between the images. 
      // globalNeighbors(i,j) != 0 means that j send to i/that i receives from j
    }

    ////////////////////////////////////////////////////////////////////////////////////// 
    // FIGURE OUT WHO IS SENDING TO WHOM AND HOW MUCH
    // DO THIS IN THE PROCESS OF PACKING ALL OUTGOING DATA ACCORDING TO DESTINATION ID
    ////////////////////////////////////////////////////////////////////////////////////// 

    // loop over all columns to know from which images I can expect to receive something
    for (Ordinal j=OT::zero(); j<numImages; ++j)
    {
      if (globalNeighbors(i,j)) {
        recvIDs.push_back(j);
      }
    }

    // we know how many we're sending to already
    // form a contiguous list of all data to be sent
    // track the number of entries for each ID
    Array<CrsIJV> IJVcontig;
    Array<Ordinal> sendSizes(sendIDs.size(), 0);
    typename Array<Ordinal>::size_type numSends = 0;
    for (typename Array<pair<Ordinal,Ordinal> >::const_iterator IdAndRow = IdsAndRows.begin();
         IdAndRow != IdsAndRows.end(); ++IdAndRow) 
    {
      Ordinal  id = IdAndRow->first;
      Ordinal row = IdAndRow->second;
      // have we advanced to a new send?
      if (sendIDs[numSends] != id) {
        numSends++;
      }
      // HERE
      // copy data for row into contiguous storage
      // for () {}
      // keep track of the amount, assign it to sendSizes[numSends]
    }
    numSends++; // one last increment, to make it a count instead of an index
    TEST_FOR_EXCEPTION(numSends != sendIDs.size(), std::logic_error, "Tpetra::CrsMatrix::globalAssemble(): internal logic error. Contact Tpetra team.");

    // from the size info, build the ArrayViews into IJVcontig

    // now I pack what has to be sent to the other images
    std::map<Ordinal,vector<Ordinal> > sendRows;
    std::map<Ordinal,vector<Ordinal> > sendCols;
    std::map<Ordinal,vector<Scalar> >  sendVals;
    for (NLITER iter = nonlocals_.begin(); iter != nonlocals_.end(); ++iter)
    {
      Ordinal row   = iter->first;
      int     image = NLR2Id[row];

      for (Ordinal i = ordinalZero() ; i < iter->second.size() ; ++i)
      {
        Ordinal col   = iter->second[i].first;
        Scalar  val   = iter->second[i].second;

        sendRows[image].push_back(row);
        sendCols[image].push_back(col);
        sendVals[image].push_back(val);
      }
    }

    int MyImageID = rowMap_.comm().getRank();

    vector<MPI_Request> send_requests(numImages * 3);
    vector<MPI_Status>  send_status(numImages * 3);

    Ordinal send_count = 0;


    for (int j = 0 ; j < numImages ; ++j)
    {
      int what = globalNeighbors[j + numImages * rowMap_.comm().getRank()];
      if (what > 0)
      {
        sendImages.push_back(j);
        send_sizes[send_count] = MpiTraits<Ordinal>::count(sendRows[j].size());

        MPI_Isend(&(send_sizes[send_count]), MpiTraits<Ordinal>::count(1), MpiTraits<Ordinal>::datatype(), 
            j, 23, MpiCommunicator, &(send_requests[send_count]));
        ++send_count;
      }
    }

    // Now receive the actual sizes
    vector<MPI_Request> recv_requests(numImages * 3);
    vector<MPI_Status>  recv_status(numImages * 3);

    vector<Ordinal> recv_sizes(numImages);
    vector<Ordinal> recv_images(numImages);

    Ordinal recv_count = 0;
    for (int j = 0 ; j < numImages ; ++j)
    {
      int what = globalNeighbors[j * numImages + rowMap_.comm().getRank()];
      if (what > 0)
      {
        recv_images[recv_count] = j;
        MPI_Irecv(&(recv_sizes[recv_count]), MpiTraits<Ordinal>::count(1), MpiTraits<Ordinal>::datatype(), 
            j, 23, MpiCommunicator, &(recv_requests[recv_count]));
        ++recv_count;
      }
    }

    MPI_Waitall(send_count, &(send_requests[0]), &(send_status[0]));
    MPI_Waitall(recv_count, &(recv_requests[0]), &(recv_status[0]));

    MPI_Barrier(MpiCommunicator);

    ////////////////////////////////////////////////////////////////////////////////////
    // NOW RECEIVE THE DATA BASED ON THE INFO FROM ABOVE
    ////////////////////////////////////////////////////////////////////////////////////

    std::map<Ordinal, vector<Ordinal> > recvRows;
    std::map<Ordinal, vector<Ordinal> > recvCols;
    std::map<Ordinal, vector<Scalar> >  recvVals;

    std::map<Ordinal, Ordinal> xxx;

    for (int i = 0 ; i < recv_count ; ++i)
    {
      Ordinal image = recv_images[i];
      recvRows[image].resize(recv_sizes[i]);
      recvCols[image].resize(recv_sizes[i]);
      recvVals[image].resize(recv_sizes[i]);
      xxx[image] = recv_sizes[i];
    }

    // At this point I know:
    // - I have to receive from `recv_count' images;
    // - image `i' will send recv_count[i] things, split in
    //   two vectors of Ordinal and a vector of Scalar.
    // First I start sending, then receiving

    send_count = 0;
    for (int j = 0 ; j < numImages ; ++j)
    {
      int what = globalNeighbors[j + numImages * rowMap_.comm().getRank()];
      if (what > 0)
      {
        // want to send to image `j', first Rows, then Cols, then Vals
        int osize = MpiTraits<Ordinal>::count(sendRows[j].size());
        int ssize = MpiTraits<Scalar>::count(sendRows[j].size());

        MPI_Isend(&(sendRows[j][0]), osize, MpiTraits<Ordinal>::datatype(), j, 32, MpiCommunicator, &(send_requests[send_count]));
        ++send_count;

        MPI_Isend(&(sendCols[j][0]), osize, MpiTraits<Ordinal>::datatype(), j, 33, MpiCommunicator, &(send_requests[send_count]));
        ++send_count;

        MPI_Isend(&(sendVals[j][0]), ssize, MpiTraits<Scalar>::datatype(), j, 34, MpiCommunicator, &(send_requests[send_count]));
        ++send_count;
      }
    }

    recv_count = 0;
    for (int j = 0 ; j < numImages ; ++j)
    {
      int what = globalNeighbors[j * numImages + rowMap_.comm().getRank()];
      if (what > 0)
      {
        int osize = MpiTraits<Ordinal>::count(xxx[j]);
        int ssize = MpiTraits<Scalar>::count(xxx[j]);

        // I want to receive from image `j', first Rows, then Cols, then Vals.
        MPI_Irecv(&(recvRows[j][0]), osize, MpiTraits<Ordinal>::datatype(), j, 32, MpiCommunicator, &(recv_requests[recv_count]));
        ++recv_count;

        MPI_Irecv(&(recvCols[j][0]), osize, MpiTraits<Ordinal>::datatype(), j, 33, MpiCommunicator, &(recv_requests[recv_count]));
        ++recv_count;

        MPI_Irecv(&(recvVals[j][0]), ssize, MpiTraits<Scalar>::datatype(), j, 34, MpiCommunicator, &(recv_requests[recv_count]));
        ++recv_count;
      }
    }

    MPI_Waitall(send_count, &(send_requests[0]), &(send_status[0]));
    MPI_Waitall(recv_count, &(recv_requests[0]), &(recv_status[0]));

    MPI_Barrier(MpiCommunicator);

    // now I add all the received elements to the list of local elements.

    for (int i = 0 ; i < recv_images.size() ; ++i)
    {
      int image = recv_images[i];
      for (int j = 0 ; j < recv_sizes[i] ; ++j)
      {
        submitEntry(Tpetra::Add, recvRows[image][j], recvCols[image][j], recvVals[image][j]);
      }
    }

    // don't need this data anymore
    nonLocals_.clear();
  }

} // namespace Tpetra

#endif
