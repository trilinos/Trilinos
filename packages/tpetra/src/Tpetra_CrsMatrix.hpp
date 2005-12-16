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

#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_VectorSpace.hpp"
#ifdef TPETRA_MPI
#include "Tpetra_MpiPlatform.hpp"
#include "Tpetra_MpiComm.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_SerialComm.hpp"
#endif
#include "Teuchos_ScalarTraits.hpp"
#include "Tpetra_Object.hpp"
#include "Tpetra_Import.hpp"

namespace Tpetra 
{
  enum ApplyMode 
  {
    AsIs,      // multiply with the matrix as it is stored
    Transpose, // multiply with its transpose
    Hermitian  // multiply with its Hermitian
  };

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
   * \date Last updated on 16-Dec-05
   *
   * \warning Class Tpetra::CisMatrix should be preferred to Tpetra::CrsMatrix. This class has been little tested.
   */
  template<class OrdinalType, class ScalarType>
  class CrsMatrix : public Object, public Teuchos::CompObject
  {
    public:
      //@{ \name Constructor/Destructor Methods

      //! Constructor specifying the primary distribution only.
      CrsMatrix(const Comm<OrdinalType, ScalarType>& Comm, 
                VectorSpace<OrdinalType, ScalarType>& VectorRowSpace) :
        Comm_(Comm),
        VectorRowSpace_(VectorRowSpace),
        RowSpace_(VectorRowSpace.elementSpace()),
        NumMyRows_(RowSpace_.getNumMyElements()),
        NumGlobalRows_(RowSpace_.getNumGlobalElements()),
        MyGlobalElements_(RowSpace_.getMyGlobalElements())
      {
        Indices_.resize(NumMyRows_);
        Values_.resize(NumMyRows_);

        FillCompleted_ = false;

        NumGlobalNonzeros_ = ordinalOne();
        NumMyNonzeros_     = ordinalOne();

        NumGlobalDiagonals_ = ordinalOne();
        NumMyDiagonals_     = ordinalOne();

        GlobalMaxNumEntries_ = ordinalOne();
        MyMaxNumEntries_     = ordinalOne();
      }

      // !Destructor.
      virtual ~CrsMatrix()
      {
        // do-nothing here
      }

      //@}
      //@{ \name Query Methods
      
      //! Returns \c true if the matrix has already been fill completed.
      bool isFillCompleted() const
      {
        return(FillCompleted_);
      }

      //! Returns the communicator.
      inline const Comm<OrdinalType, ScalarType>& getComm() const
      {
        return(Comm_);
      }

      //! Returns the number of nonzero entries in the global matrix. 
      inline OrdinalType getNumGlobalNonzeros () const
      {
        if(!isFillCompleted())
          throw(-1);
        
        return(NumGlobalNonzeros_);
      }

      //! Returns the number of nonzero entries in the calling image's portion of the matrix. 
      inline OrdinalType getNumMyNonzeros () const
      {
        if(!isFillCompleted())
          throw(-1);
        
        return(NumMyNonzeros_);
      }

      //! Returns the number of global matrix rows. 
      inline OrdinalType getNumGlobalRows () const
      {
        return(RowSpace_.getGlobalElements());
      }

      //! Returns the number of global matrix columns. 
      inline OrdinalType getNumGlobalCols () const
      {
        if(!isFillCompleted())
          throw(-1);
        
        return(ColSpace_->getGlobalElements());
      }

      //! Returns the number of matrix rows owned by the calling image. 
      inline OrdinalType getNumMyRows () const
      {
        return(RowSpace_.getMyElements());
      }

      //! Returns the number of matrix columns owned by the calling image. 
      inline OrdinalType getNumMyCols () const
      {
        if(!isFillCompleted())
          throw(-1);
        
        return(ColSpace_->getMyElements());
      }

      //! Returns the number of global nonzero diagonal entries, based on global row/column index comparisons. 
      inline OrdinalType getNumGlobalDiagonals () const
      {
        if(!isFillCompleted())
          throw(-1);
        
        return(NumGlobalDiagonals_);
      }

      //! Returns the number of local nonzero diagonal entries, based on global row/column index comparisons. 
      inline OrdinalType getNumMyDiagonals () const
      {
        if(!isFillCompleted())
          throw(-1);
        
        return(NumMyDiagonals_);
      }

      //! Returns the current number of nonzero entries in specified global index on this image. 
      inline OrdinalType getNumEntries (OrdinalType index) const
      {
        if(!isFillCompleted())
          throw(-1);
        
        return(Indices_[index].size());
      }

      //! Returns the maximum number of nonzero entries across all rows/columns on all images. 
      inline OrdinalType getGlobalMaxNumEntries () const
      {
        if(!isFillCompleted())
          throw(-1);
        
        return(GlobalMaxNumEntries_);
      }

      //! Returns the maximum number of nonzero entries across all rows/columns on this image. 
      inline OrdinalType getMyMaxNumEntries () const
      {
        if(!isFillCompleted())
          throw(-1);
        
        return(MyMaxNumEntries_);
      }

      //! Returns the index base for global indices for this matrix. 
      inline OrdinalType getIndexBase () const
      {
        return(RowSpace_.getIndexBase());
      }
      
      //! Returns true if this matrix is row-oriented, and false if this matrix is column-oriented.
      bool isRowOriented() const 
      {
        return(true);
      }

      //! Returns the VectorSpace that describes the primary distribution in this matrix.
      VectorSpace<OrdinalType, ScalarType> const& getPrimaryDist() const 
      {
        return(VectorRowSpace_);
      }

      VectorSpace<OrdinalType, ScalarType> const& getSecondaryDist() const 
      {
        if (isFillCompleted() == false)
          throw(-1); // not available at this point
        return(*VectorColSpace_);
      }

      VectorSpace<OrdinalType, ScalarType> const& getRowDist() const 
      {
        return(getPrimaryDist());
      }

      VectorSpace<OrdinalType, ScalarType> const& getColumnDist() const 
      {
        return(getSecondaryDist());
      }

      //! Returns the VectorSpace associated with the domain of this matrix.
      VectorSpace<OrdinalType, ScalarType> const& getDomainDist() const 
      {
        return(getSecondaryDist());
      }

      //! Returns the VectorSpace associated with the domain of this matrix.
      VectorSpace<OrdinalType, ScalarType> const& getRangeDist() const 
      {
        return(getPrimaryDist());
      }






      //@}
      //@{ \name Construction Methods
      
      //! Signals that data entry is complete. Matrix data is converted into a more optimized form.
      void fillComplete()
      {
        if (isFillCompleted())
          throw(-1);

        // =============================== //
        // Part 0: send off-image elements //
        // =============================== //

        if (getComm().getNumImages() > 1) globalAssemble();

        // =============================== //
        // Part I: remove repeated indices //
        // =============================== //
        
        // I load all matrix entries in a hash table, then I re-fill
        // the row with the last inserted value.
        for (OrdinalType i = ordinalZero() ; i < NumMyRows_ ; ++i)
        {
          std::map<OrdinalType, ScalarType> singleRow;

          for (OrdinalType j = ordinalZero() ; j < Indices_[i].size() ; ++j)
          {
            singleRow[Indices_[i][j]] += Values_[i][j];
          }

          OrdinalType count = ordinalZero();

          for (typename std::map<OrdinalType,ScalarType>::iterator iter = singleRow.begin() ; 
               iter != singleRow.end() ; ++iter)
          {
            Indices_[i][count] = iter->first;
            Values_[i][count] = iter->second;
            ++count;
          }

          Indices_[i].resize(count);
          Values_[i].resize(count);
        }

        // =============================== //
        // Part II: build the column space //
        // =============================== //
        
        // I have to find the list of non-locally owned columns

        map<OrdinalType, bool> container; // replace with a hash table

        for (OrdinalType i = ordinalZero() ; i < NumMyRows_ ; ++i)
        {
          for (OrdinalType j = ordinalZero() ; j < Indices_[i].size() ; ++j)
          {
            OrdinalType what = Indices_[i][j];
            if (RowSpace_.isMyGID(what)) continue;
            else
              container[what] = true;
          }
        }

        vector<OrdinalType> MyPaddedGlobalElements(MyGlobalElements_);

        for (typename std::map<OrdinalType, bool>::iterator iter = container.begin() ; 
             iter != container.end() ; ++iter)
        {
          MyPaddedGlobalElements.push_back(iter->first);
        }

        // now I can build the column space

        ColSpace_ = new ElementSpace<OrdinalType>(-1, MyPaddedGlobalElements.size(),
                                                  MyPaddedGlobalElements, RowSpace_.getIndexBase(), RowSpace_.platform());

        VectorColSpace_ = new VectorSpace<OrdinalType, ScalarType>(*ColSpace_, VectorRowSpace_.platform());

        PaddedVector_ = new Vector<OrdinalType, ScalarType>(*VectorColSpace_);
        Importer_ = new Import<OrdinalType>(RowSpace_, *ColSpace_);

        // ============================== //
        // Part IV: move to local indices //
        // ============================== //
        
        for (OrdinalType i = ordinalZero() ; i < NumMyRows_ ; ++i)
        {
          for (OrdinalType j = ordinalZero() ; j < Indices_[i].size() ; ++j)
          {
            Indices_[i][j] = ColSpace_->getLID(Indices_[i][j]);
          }
        }

        // ================================================ //
        // Part V: compute some local and global quantities //
        // ================================================ //
        
        for (OrdinalType i = ordinalZero() ; i < NumMyRows_ ; ++i)
        {
          OrdinalType NumEntries = Indices_[i].size();

          NumMyNonzeros_ += NumEntries;

          for (OrdinalType j = ordinalZero() ; j < NumEntries; ++j)
          {
            OrdinalType col = Indices_[i][j];
            if (col == i)
            {
              ++NumMyDiagonals_;
              break;
            }
          }
          
          if (MyMaxNumEntries_ < NumEntries) MyMaxNumEntries_ = NumEntries;
        }

        // workaround, should be fixed in Comm
#ifdef HAVE_MPI
        MPI_Comm MpiCommunicator = (dynamic_cast<const MpiComm<OrdinalType, ScalarType>&>(getComm())).getMpiComm();
        MPI_Allreduce((void*)&NumMyNonzeros_, (void*)&NumGlobalNonzeros_, MpiTraits<OrdinalType>::count(1),
                      MpiTraits<OrdinalType>::datatype(), MPI_SUM, MpiCommunicator);

        MPI_Allreduce((void*)&NumMyDiagonals_, (void*)&NumGlobalDiagonals_, MpiTraits<OrdinalType>::count(1),
                      MpiTraits<OrdinalType>::datatype(), MPI_SUM, MpiCommunicator);

        MPI_Allreduce((void*)&MyMaxNumEntries_, (void*)&GlobalMaxNumEntries_, MpiTraits<OrdinalType>::count(1),
                      MpiTraits<OrdinalType>::datatype(), MPI_MAX, MpiCommunicator);
#else
        NumGlobalNonzeros_ = NumMyNonzeros_;
        NumGlobalDiagonals_ = NumMyDiagonals_;
        GlobalMaxNumEntries_ = MyMaxNumEntries_;
#endif

        // mark transformation as successfully completed

        FillCompleted_ = true;
      }

      //! Signals that data entry is complete. Matrix data is converted into a more optimized form.
      void fillComplete(VectorSpace<OrdinalType, ScalarType> const& domainSpace, 
                        VectorSpace<OrdinalType, ScalarType> const& rangeSpace) 
      {
        throw(-1); // not implemented
      }

      //! Submits one local or nonlocal entry to the matrix using global IDs.
      inline void submitEntry(const CombineMode CM, const OrdinalType& GlobalRow, const OrdinalType& GlobalCol,
                              const ScalarType& Value)
      {
        if (CM != Tpetra::Add)
          throw(-1); // other modes not implemented yet.

        if (FillCompleted_)
          throw(-1);

        if (RowSpace_.isMyGID(GlobalRow))
        {
          OrdinalType MyRow = RowSpace_.getLID(GlobalRow);

          Indices_[MyRow].push_back(GlobalCol);
          Values_[MyRow].push_back(Value);
        }
        else
        {
          nonlocals_[GlobalRow].push_back(pair<OrdinalType, ScalarType>(GlobalCol, Value));
        }
      }

      //! Submit multiple entries, using global IDs.
      /*! All index values must be in the global space. Behavoir is defined by the CombineMode passed in. */
      void submitEntries(const CombineMode CM, OrdinalType const myRowOrColumn, 
                         OrdinalType const numEntries, ScalarType const* values, 
                         OrdinalType const* indices) 
      {
        for (OrdinalType i = ordinalZero() ; i < NumEntries ; ++i)
          submitEntry(CM, myRowOrColumn, indices[i], values[i]);
      }

      // @}
      // @{ \name Computational Methods
      
      //! Returns a copy of the specified local row, column indices are local.
      inline void getMyRowCopy(const OrdinalType MyRow, OrdinalType& NumEntries,
                               vector<OrdinalType>& Indices, vector<ScalarType>& Values) const
      {
        if (FillCompleted_)
          throw(-1);

        NumEntries = Indices_[MyRow].size();

        if (Indices.size() < NumEntries || Values.size() < NumEntries)
          throw(-1);

        for (OrdinalType i = ordinalZero() ; i < NumEntries ; ++i)
        {
          Indices[i] = Indices_[MyRow][i];
          Values[i] = Values_[MyRow][i];
        }
      }

      //! Returns a copy of the specified (and locally owned) global row, column indices are global.
      inline void getGlobalRowCopy(const OrdinalType GlobalRow, OrdinalType& NumEntries,
                                   vector<OrdinalType>& Indices, vector<ScalarType>& Values) const
      {
        // Only locally owned rows can be queried, otherwise complain
        if (!RowSpace_.isMyGID(GlobalRow))
          throw(-1);

        OrdinalType MyRow = RowSpace_.getLID(GlobalRow);

        NumEntries = Indices_[MyRow].size();

        if (Indices.size() < NumEntries || Values.size() < NumEntries)
          throw(-1);

        if (isFillCompleted())
        {
          for (OrdinalType i = ordinalZero() ; i < NumEntries ; ++i)
          {
            Indices[i] = ColSpace_->getGID(Indices_[MyRow][i]);
            Values[i] = Values_[MyRow][i];
          }
        }
        else
        {
          for (OrdinalType i = ordinalZero() ; i < NumEntries ; ++i)
          {
            Indices[i] = Indices_[MyRow][i];
            Values[i] = Values_[MyRow][i];
          }
        }
      }

      //! Computes the matrix-vector multilication y = A x.
      void apply(const Vector<OrdinalType, ScalarType>& x, Vector<OrdinalType, ScalarType> y,
                 ApplyMode Mode = AsIs) const
      {
        assert (Mode == AsIs);

        y.setAllToScalar(scalarZero());

        PaddedVector_->doImport(x, *Importer_, Insert);

        for (OrdinalType i = ordinalZero() ; i < NumMyRows_ ; ++i)
        {
          for (OrdinalType j = ordinalZero() ; j < Indices_[i].size() ; ++j)
          {
            OrdinalType col = Indices_[i][j];
            ScalarType val = Values_[i][j];
            y[i] += val * (*PaddedVector_)[col];
          }
        }
      }

      //@}
      //@{ \name I/O Methods
      
      //! Prints the matrix on the specified stream.
      virtual void print(ostream& os) const 
      {
        int MyImageID = RowSpace_.comm().getMyImageID();

        if (MyImageID == 0)
        {
          os << "Tpetra::CrsMatrix, label = " << label() << endl;
          os << "Number of global rows    = " << RowSpace_.getNumGlobalElements() << endl;
          if (isFillCompleted())
          {
            os << "Number of global columns = " << ColSpace_->getNumGlobalElements() << endl;
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
          for (int pid = 0 ; pid < RowSpace_.comm().getNumImages() ; ++pid)
          {
            if (pid == MyImageID)
            {
              vector<OrdinalType> Indices(getMyMaxNumEntries()); // FIXME
              vector<ScalarType>  Values(getMyMaxNumEntries());
              OrdinalType NumEntries;

              os << "% Number of rows on image " << MyImageID << " = " << RowSpace_.getNumMyElements() << endl;
              for (OrdinalType i = ordinalZero() ; i < RowSpace_.getNumMyElements() ; ++i)
              {
                OrdinalType GlobalRow = RowSpace_.getGID(i);

                getGlobalRowCopy(GlobalRow, NumEntries, Indices, Values);
                for (OrdinalType j = ordinalZero() ; j < NumEntries ; ++j)
                  os << "Matrix(" << GlobalRow << ", " << Indices[j] << ") = " << Values[j] << ";" << endl;
              }
            }
            RowSpace_.comm().barrier();
          }
        }
      }

      //! Set all matrix entries equal to scalarThis.
      inline void setAllToScalar (ScalarType scalarThis)
      {
        throw(-2); // not yet implemented
      }

      //! Scale the current values of a matrix, this = scalarThis*this. 
      inline void scale (ScalarType scalarThis)
      {
        throw(-2); // not yet implemented
      }


      //! Basic print, for debugging purposes only.
      void rawPrint()
      {
        // this prints out the structure as they are
        for (int i = 0 ; i < Indices_.size() ; ++i)
        {
          for (int j = 0 ; j < Indices_[i].size() ; ++j)
          {
            cout << "local row " << i << ", col " << Indices_[i][j] << ", val " << Values_[i][j] << endl;
          }
        }
      }
      // @}

    private:
      //! copy constructor.
      CrsMatrix(CrsMatrix<OrdinalType, ScalarType> const& Source)
      {
        throw(-1); // not implemented
      }

      CrsMatrix& operator=(const CrsMatrix<OrdinalType, ScalarType>& rhs)
      {
        throw(-1); // not implemented
      }

      //! Returns `1' for the ordinal type.
      inline OrdinalType ordinalOne() const
      {
        return(Teuchos::ScalarTraits<OrdinalType>::one());
      }

      //! Returns `0' for the ordinal type.
      inline OrdinalType ordinalZero() const
      {
        return(Teuchos::ScalarTraits<OrdinalType>::zero());
      }

      //! Returns `-1' for the ordinal type.
      inline OrdinalType ordinalMinusOne() const
      {
        return(ordinalZero() - ordinalOne());
      }

      //! Returns `1' for the scalar type.
      inline ScalarType scalarOne() const
      {
        return(Teuchos::ScalarTraits<ScalarType>::one());
      }

      //! Returns `0' for the scalar type.
      inline ScalarType scalarZero() const
      {
        return(Teuchos::ScalarTraits<ScalarType>::zero());
      }

      //! Returns `-1' for the scalar type.
      inline ScalarType scalarMinusOne() const
      {
        return(scalarZero() - scalarOne());
      }

      const vector<vector<OrdinalType> >& getIndices() const
      {
        return(Indices_);
      }

      const vector<vector<ScalarType> >& getValues() const
      {
        return(Values_);
      }

      //! Performs importing of off-processor elements and adds them to the locally owned elements.
      void globalAssemble()
      {
#ifdef HAVE_MPI
        MPI_Comm MpiCommunicator;
        try
        {
          MpiCommunicator = (dynamic_cast<const MpiComm<OrdinalType, ScalarType>&>(getComm())).getMpiComm();
        }
        catch(std::bad_cast bc) 
        {
          cerr << "Bad cast" << endl;
          throw(-1);
        }

        // First I want to check that we actually need to do this; it may be
        // that the user has only inserted locally owned elements.
        
        OrdinalType MyNonlocals = nonlocals_.size(), GlobalNonlocals;

        MPI_Allreduce((void*)&MyNonlocals, (void*)&GlobalNonlocals, MpiTraits<OrdinalType>::count(1),
                      MpiTraits<OrdinalType>::datatype(), MPI_MAX, MpiCommunicator);

        if (GlobalNonlocals == ordinalZero()) return;

        // Ok, so we need to do the hard work.
        
        int NumImages = getComm().getNumImages();

        map<OrdinalType, OrdinalType> containter_map;

        // this is a list of non-locally owned rows, in a map (should become a
        // hash some day for faster access)
        for (typename map<OrdinalType, vector<pair<OrdinalType, ScalarType> > >::iterator iter = nonlocals_.begin() ; 
             iter != nonlocals_.end() ; ++iter)
        {
          containter_map[iter->first] = ordinalOne();
        }

        // convert the map to a vector so that I can use get getRemoteIDList()
        vector<OrdinalType> container_vector;

        for (typename map<OrdinalType, OrdinalType>::iterator iter = containter_map.begin() ;
             iter != containter_map.end() ; ++iter)
        {
          container_vector.push_back(iter->first);
        }

        vector<int> image_vector(container_vector.size());

        RowSpace_.getRemoteIDList (container_vector, image_vector);

        map<OrdinalType, int> image_map;

        for (OrdinalType i = ordinalZero() ; i < image_vector.size() ; ++i)
        {
          image_map[container_vector[i]] = image_vector[i];
        }

        vector<int> local_neighbors(RowSpace_.comm().getNumImages());
        for (int i = 0 ; i < local_neighbors.size() ; ++i) local_neighbors[i] = 0;

        for (int i = 0 ; i < image_vector.size() ; ++i)
        {
          local_neighbors[image_vector[i]] = 1;
        }

        vector<int> global_neighbors(NumImages * NumImages);

        RowSpace_.comm().gatherAll(&local_neighbors[0], &global_neighbors[0], NumImages);

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
          int what = global_neighbors[j * NumImages + RowSpace_.comm().getMyImageID()];
          if (what > 0)
          {
            recvImages.push_back(j);
          }
        }

        // do the same but with send
        vector<int> sendImages;

        // now I pack what has to be sent to the other images
        map<OrdinalType, vector<OrdinalType> > sendRows;
        map<OrdinalType, vector<OrdinalType> > sendCols;
        map<OrdinalType, vector<ScalarType> >  sendVals;

        for (typename map<OrdinalType, vector<pair<OrdinalType, ScalarType> > >::iterator iter = nonlocals_.begin() ; 
             iter != nonlocals_.end() ; ++iter)
        {
          OrdinalType row   = iter->first;
          int image = image_map[row];

          for (OrdinalType i = ordinalZero() ; i < iter->second.size() ; ++i)
          {
            OrdinalType col   = iter->second[i].first;
            ScalarType  val   = iter->second[i].second;

            sendRows[image].push_back(row);
            sendCols[image].push_back(col);
            sendVals[image].push_back(val);
          }
        }

        int MyImageID = RowSpace_.comm().getMyImageID();

        vector<MPI_Request> send_requests(NumImages * 3);
        vector<MPI_Status>  send_status(NumImages * 3);

        OrdinalType send_count = 0;

        vector<OrdinalType> send_sizes(NumImages); // because Isend is not buffered

        for (int j = 0 ; j < NumImages ; ++j)
        {
          int what = global_neighbors[j + NumImages * RowSpace_.comm().getMyImageID()];
          if (what > 0)
          {
            sendImages.push_back(j);
            send_sizes[send_count] = MpiTraits<OrdinalType>::count(sendRows[j].size());

            MPI_Isend(&(send_sizes[send_count]), MpiTraits<OrdinalType>::count(1), MpiTraits<OrdinalType>::datatype(), 
                      j, 23, MpiCommunicator, &(send_requests[send_count]));
            ++send_count;
          }
        }

        // Now receive the actual sizes
        vector<MPI_Request> recv_requests(NumImages * 3);
        vector<MPI_Status>  recv_status(NumImages * 3);

        vector<OrdinalType> recv_sizes(NumImages);
        vector<OrdinalType> recv_images(NumImages);

        OrdinalType recv_count = 0;
        for (int j = 0 ; j < NumImages ; ++j)
        {
          int what = global_neighbors[j * NumImages + RowSpace_.comm().getMyImageID()];
          if (what > 0)
          {
            recv_images[recv_count] = j;
            MPI_Irecv(&(recv_sizes[recv_count]), MpiTraits<OrdinalType>::count(1), MpiTraits<OrdinalType>::datatype(), 
                      j, 23, MpiCommunicator, &(recv_requests[recv_count]));
            ++recv_count;
          }
        }

        MPI_Waitall(send_count, &(send_requests[0]), &(send_status[0]));
        MPI_Waitall(recv_count, &(recv_requests[0]), &(recv_status[0]));

        MPI_Barrier(MpiCommunicator);

        map<OrdinalType, vector<OrdinalType> > recvRows;
        map<OrdinalType, vector<OrdinalType> > recvCols;
        map<OrdinalType, vector<ScalarType> >  recvVals;

        map<OrdinalType, OrdinalType> xxx;

        for (int i = 0 ; i < recv_count ; ++i)
        {
          OrdinalType image = recv_images[i];
          recvRows[image].resize(recv_sizes[i]);
          recvCols[image].resize(recv_sizes[i]);
          recvVals[image].resize(recv_sizes[i]);

          xxx[image] = recv_sizes[i];
        }

        // At this point I know:
        // - I have to receive from `recv_count' images;
        // - image `i' will send recv_count[i] things, split in
        //   two vectors of OrdinalType and a vector of ScalarType.
        // First I start sending, then receiving

        send_count = 0;
        for (int j = 0 ; j < NumImages ; ++j)
        {
          int what = global_neighbors[j + NumImages * RowSpace_.comm().getMyImageID()];
          if (what > 0)
          {
            // want to send to image `j', first Rows, then Cols, then Vals
            int osize = MpiTraits<OrdinalType>::count(sendRows[j].size());
            int ssize = MpiTraits<ScalarType>::count(sendRows[j].size());
            
            MPI_Isend(&(sendRows[j][0]), osize, MpiTraits<OrdinalType>::datatype(), j, 32, MpiCommunicator, &(send_requests[send_count]));
            ++send_count;

            MPI_Isend(&(sendCols[j][0]), osize, MpiTraits<OrdinalType>::datatype(), j, 33, MpiCommunicator, &(send_requests[send_count]));
            ++send_count;

            MPI_Isend(&(sendVals[j][0]), ssize, MpiTraits<ScalarType>::datatype(), j, 34, MpiCommunicator, &(send_requests[send_count]));
            ++send_count;
          }
        }

        recv_count = 0;
        for (int j = 0 ; j < NumImages ; ++j)
        {
          int what = global_neighbors[j * NumImages + RowSpace_.comm().getMyImageID()];
          if (what > 0)
          {
            int osize = MpiTraits<OrdinalType>::count(xxx[j]);
            int ssize = MpiTraits<ScalarType>::count(xxx[j]);

            // I want to receive from image `j', first Rows, then Cols, then Vals.
            MPI_Irecv(&(recvRows[j][0]), osize, MpiTraits<OrdinalType>::datatype(), j, 32, MpiCommunicator, &(recv_requests[recv_count]));
            ++recv_count;

            MPI_Irecv(&(recvCols[j][0]), osize, MpiTraits<OrdinalType>::datatype(), j, 33, MpiCommunicator, &(recv_requests[recv_count]));
            ++recv_count;

            MPI_Irecv(&(recvVals[j][0]), ssize, MpiTraits<ScalarType>::datatype(), j, 34, MpiCommunicator, &(recv_requests[recv_count]));
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

      const Comm<OrdinalType, ScalarType>& Comm_;
      const VectorSpace<OrdinalType, ScalarType>& VectorRowSpace_;
      const ElementSpace<OrdinalType>& RowSpace_;

      OrdinalType NumMyRows_, NumGlobalRows_;
      const vector<OrdinalType>& MyGlobalElements_;

      std::vector<vector<OrdinalType> > Indices_;
      std::vector<vector<ScalarType> > Values_;

      ElementSpace<OrdinalType>* ColSpace_;
      VectorSpace<OrdinalType, ScalarType>* VectorColSpace_;

      Vector<OrdinalType, ScalarType>* PaddedVector_;
      Import<OrdinalType>* Importer_;

      map<OrdinalType, vector<pair<OrdinalType, ScalarType> > > nonlocals_;

      bool FillCompleted_;

      OrdinalType NumGlobalNonzeros_;
      OrdinalType NumMyNonzeros_;

      OrdinalType NumGlobalDiagonals_;
      OrdinalType NumMyDiagonals_;

      OrdinalType GlobalMaxNumEntries_;
      OrdinalType MyMaxNumEntries_;


  }; // class CrsMatrix

} // namespace Tpetra
