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
  class CrsMatrix : public Teuchos::Object, public Teuchos::CompObject
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
      const Teuchos::Comm<Ordinal>& getComm() const;

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
      Ordinal getNumGlobalDiagonals() const

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
      void submitEntry(const Ordinal &globalRow, const Ordinal &globalCol,
                       const Scalar &value);

      //! Submit multiple entries, using global IDs.
      /*! All index values must be in the global space. Behavoir is defined by the CombineMode passed in. */
      void submitEntries(const Ordinal &globalRow, 
                         const Teuchos::ArrayView<const Scalar> &values, 
                         const Teuchos::ArrayView<const Ordinal> &indices);

      // @}
      //! @name Computational Methods
      // @{ 

      //! Returns a copy of the specified local row, column indices are local.
      void getMyRowCopy(const Ordinal &localRow, const Teuchos::ArrayView<Ordinal> &indices, 
                                                 const Teuchos::ArrayView<Scalar> &values) const;

      //! Returns a copy of the specified (and locally owned) global row, column indices are global.
      void getGlobalRowCopy(const Ordinal &globalRow, 
                            const Teuchos::ArrayView<Ordinal> &indices,
                            const Teuchos::ArrayView<Scalar> &values) const;

      //@}
      //! @name I/O Methods
      //@{ 
      
      //! Prints the matrix on the specified stream.
      virtual void print(std::ostream& os) const;

      //! Set all matrix entries equal to scalarThis.
      void setAllToScalar(Scalar scalarThis);

      //! Scale the current values of a matrix, this = scalarThis*this. 
      void scale(Scalar scalarThis);

      //! Basic print, for debugging purposes only.
      void rawPrint();

      // @}

    private:
      //! copy constructor.
      CrsMatrix(const CrsMatrix<Ordinal,Scalar> &Source);

      CrsMatrix& operator=(const CrsMatrix<Ordinal, Scalar> &rhs);

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

      std::map<Ordinal, std::vector<pair<Ordinal, Scalar> > > nonlocals_;

      bool fillCompleted_;

      Ordinal numGlobalNonzeros_;
      Ordinal numMyNonzeros_;

      Ordinal numGlobalDiagonals_;
      Ordinal numMyDiagonals_;

      Ordinal globalMaxNumEntries_;
      Ordinal myMaxNumEntries_;


  }; // class CrsMatrix

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
  {
    return fillCompleted_;
  }

  template<class Ordinal, class Scalar>
  Teuchos::RCP<const Teuchos::Comm<Ordinal> > 
  CrsMatrix<Ordinal,Scalar>::getComm() const
  {
    return comm_;
  }

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
  {
    return rowMap_.getNumGlobalEntries();
  }

  template<class Ordinal, class Scalar>
  Ordinal CrsMatrix<Ordinal,Scalar>::getNumGlobalCols() const
  {
    TEST_FOR_EXCEPTION(isFillCompleted() == false, std::runtime_error,
      "Tpetra::CrsMatrix: cannot call getNumGlobalCols() until fillComplete() has been called.");
    return colMap_->getGlobalEntries();
  }

  template<class Ordinal, class Scalar>
  Ordinal CrsMatrix<Ordinal,Scalar>::getNumMyRows() const
  {
    return rowMap_.getNumMyEntries();
  }

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
  {
    return rowMap_.getIndexBase();
  }

  template<class Ordinal, class Scalar>
  const Map<Ordinal> & CrsMatrix<Ordinal,Scalar>::getRowMap() const 
  {
    return rowMap_;
  }

  template<class Ordinal, class Scalar>
  const Map<Ordinal> & CrsMatrix<Ordinal,Scalar>::getColMap() const 
  {
    return colMap_;
  }

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
  void CrsMatrix<Ordinal,Scalar>::submitEntry(const Ordinal &globalRow, const Ordinal &globalCol,
                                              const Scalar &value)
  {
    if (fillCompleted_)
      throw(-1);
    if (rowMap_.isMyGID(globalRow))
    {
      Ordinal myRow = rowMap_.getLID(globalRow);
      indices_[myRow].push_back(globalCol);
      values_[myRow].push_back(value);
    }
    else
    {
      nonlocals_[globalRow].push_back(pair<Ordinal, Scalar>(globalCol, value));
    }
  }

  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::submitEntries(const Ordinal &globalRow, 
                         const Teuchos::ArrayView<Scalar> &values, 
                         const Teuchos::ArrayView<Ordinal> &indices)
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
  void CrsMatrix<Ordinal,Scalar>::getMyRowCopy(const Ordinal &localRow, 
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
  void CrsMatrix<Ordinal,Scalar>::getGlobalRowCopy(const Ordinal &globalRow, 
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
  void CrsMatrix<Ordinal,Scalar>::print(ostream& os) const 
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
  {
    throw(-2); // not yet implemented
  }

  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::scale(Scalar scalarThis)
  {
    throw(-2); // not yet implemented
  }

  template<class Ordinal, class Scalar>
  void CrsMatrix<Ordinal,Scalar>::rawPrint()
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
#ifdef HAVE_MPI_THIS_IS_BROKEN // FINISH
    MPI_Comm MpiCommunicator;
    try
    {
      MpiCommunicator = (dynamic_cast<const MpiComm<Ordinal, Scalar>&>(getComm())).getMpiComm();
    }
    catch(std::bad_cast bc) 
    {
      cerr << "Bad cast" << endl;
      throw(-1);
    }

    // First I want to check that we actually need to do this; it may be
    // that the user has only inserted locally owned elements.

    Ordinal MyNonlocals = nonlocals_.size(), GlobalNonlocals;

    MPI_Allreduce((void*)&MyNonlocals, (void*)&GlobalNonlocals, MpiTraits<Ordinal>::count(1),
        MpiTraits<Ordinal>::datatype(), MPI_MAX, MpiCommunicator);

    if (GlobalNonlocals == ordinalZero()) return;

    // Ok, so we need to do the hard work.

    int NumImages = getComm().getSize();

    std::map<Ordinal, Ordinal> containter_map;

    // this is a list of non-locally owned rows, in a map (should become a
    // hash some day for faster access)
    for (typename std::map<Ordinal, vector<pair<Ordinal, Scalar> > >::iterator iter = nonlocals_.begin() ; 
        iter != nonlocals_.end() ; ++iter)
    {
      containter_map[iter->first] = ordinalOne();
    }

    // convert the map to a vector so that I can use get getRemoteIDList()
    vector<Ordinal> container_vector;

    for (typename std::map<Ordinal, Ordinal>::iterator iter = containter_map.begin() ;
        iter != containter_map.end() ; ++iter)
    {
      container_vector.push_back(iter->first);
    }

    vector<int> image_vector(container_vector.size());

    rowMap_.getRemoteIDList (container_vector, image_vector);

    std::map<Ordinal, int> image_map;

    for (Ordinal i = ordinalZero() ; i < image_vector.size() ; ++i)
    {
      image_map[container_vector[i]] = image_vector[i];
    }

    vector<int> local_neighbors(rowMap_.comm().getSize());
    for (int i = 0 ; i < local_neighbors.size() ; ++i) local_neighbors[i] = 0;

    for (int i = 0 ; i < image_vector.size() ; ++i)
    {
      local_neighbors[image_vector[i]] = 1;
    }

    vector<int> global_neighbors(NumImages * NumImages);

    rowMap_.comm().gatherAll(&local_neighbors[0], &global_neighbors[0], NumImages);

    // `global_neighbors' at this point contains (on all images) the
    // connectivity between the images. On the row `i', a nonzero on col `j' means
    // that image i will send something to image j. On the column `j', a
    // nonzero on row `i' means that image j will receive something from
    // image i.

    // now I loop over all columns to know which image is supposed to send
    // me something
    vector<int> recvImages;

    for (int j = 0 ; j < NumImages ; ++j)
    {
      int what = global_neighbors[j * NumImages + rowMap_.comm().getRank()];
      if (what > 0)
      {
        recvImages.push_back(j);
      }
    }

    // do the same but with send
    vector<int> sendImages;

    // now I pack what has to be sent to the other images
    std::map<Ordinal, vector<Ordinal> > sendRows;
    std::map<Ordinal, vector<Ordinal> > sendCols;
    std::map<Ordinal, vector<Scalar> >  sendVals;

    for (typename std::map<Ordinal, vector<pair<Ordinal, Scalar> > >::iterator iter = nonlocals_.begin() ; 
        iter != nonlocals_.end() ; ++iter)
    {
      Ordinal row   = iter->first;
      int image = image_map[row];

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

    vector<MPI_Request> send_requests(NumImages * 3);
    vector<MPI_Status>  send_status(NumImages * 3);

    Ordinal send_count = 0;

    vector<Ordinal> send_sizes(NumImages); // because Isend is not buffered

    for (int j = 0 ; j < NumImages ; ++j)
    {
      int what = global_neighbors[j + NumImages * rowMap_.comm().getRank()];
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
    vector<MPI_Request> recv_requests(NumImages * 3);
    vector<MPI_Status>  recv_status(NumImages * 3);

    vector<Ordinal> recv_sizes(NumImages);
    vector<Ordinal> recv_images(NumImages);

    Ordinal recv_count = 0;
    for (int j = 0 ; j < NumImages ; ++j)
    {
      int what = global_neighbors[j * NumImages + rowMap_.comm().getRank()];
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
    for (int j = 0 ; j < NumImages ; ++j)
    {
      int what = global_neighbors[j + NumImages * rowMap_.comm().getRank()];
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
    for (int j = 0 ; j < NumImages ; ++j)
    {
      int what = global_neighbors[j * NumImages + rowMap_.comm().getRank()];
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
#endif
  }

} // namespace Tpetra

#endif
