/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef EPETRAEXT_DISTARRAY_H
#define EPETRAEXT_DISTARRAY_H

#include "EpetraExt_ConfigDefs.h"
#include "EpetraExt_Exception.h"
#include "Epetra_Map.h"
#include "Epetra_DistObject.h"

namespace EpetraExt 
{
/*! 
\brief DistArray<T>: A class to store row-oriented multivectors of type T

Class DistArray allows the construction and usage of multivectors. These 
vectors contain element of type T, and the storage is row-oriented, and not
column-oriented as in class Epetra_MultiVector. As such, this class should
be used as a container for data, on which no BLAS-like operations are
performed. 

DistArray objects are indentified by an \c Epetra_Map and a \c RowSize. The 
map specifies the distribution of the elements across the processors and
therefore the number of local elements, while the RowSize gives the
total number of data assigned to each node. RowSize is constant for all
elements.

DistArray is derived from Epetra_DistObject, and it can therefore be
redistributed using Import/Export instructions.

The typical usage of this class is as follows:
\code
EpetraExt::DistArray<double> COO(VertexMap, NumDimensions);

// set the value of the j-th dimension of the i-th local node i:
COO(i, j) = 1.24
\endcode

\author Marzio Sala, ETHZ/D-INFK.

\date Last updated on Mar-06.
*/
  template<class T>
  class DistArray : public Epetra_DistObject
  {
    public:
      // @{ \name Constructors and Destructors

      //! Constructor for a given \c Map and \c RowSize.
      DistArray(const Epetra_Map& Map, const int RowSize) :
        Epetra_DistObject(Map)
      {
        // only Epetra_Map's with constant element size of 1 are allowed
        if (Map.MaxElementSize() != 1) 
          throw(Exception(__FILE__, __LINE__,
                          "Map.MaxElementSize() != 1"));
        if (!Map.ConstantElementSize()) 
          throw(Exception(__FILE__, __LINE__,
                          "Map.ConstantElementSize() != true"));

        MyLength_     = Map.NumMyElements();
        GlobalLength_ = Map.NumGlobalElements();
        RowSize_      = RowSize;
        count_        = 0;
        values_.resize(MyLength_ * RowSize_);
      }

      // @}
      // @{ \name Query methods

      //! Returns the length of the locally owned array.
      inline int MyLength() const
      {
        return(MyLength_);
      }

      //! Returns the global length of the array.
      inline int GlobalLength() const
      {
        return(GlobalLength_);
      }

      //! Returns the row size, that is, the amount of data associated with each element.
      inline int RowSize() const
      {
        return(RowSize_);
      }

      //! Returns a reference to the \c ID component of the \c LEID local element.
      inline T& operator()(const int LEID, const int ID)
      {
        assert (ID <= RowSize_);
        return(values_[LEID * RowSize_ + ID]);
      }

      inline T& operator()(const int GEID, const int ID, const bool isLocal)
      {
        int LEID = Map().LID(GEID);
        assert (LEID != -1);
        assert (ID <= RowSize_);
        return(values_[LEID * RowSize_ + ID]);
      }

      //! Prints the array on the specified stream.
      void Print(std::ostream& os) const
      {
        os << "DistArray object, label   = " << this->Label() << std::endl;
        os << "Number of local elements  = " << Map().NumMyElements() << std::endl;
        os << "Number of global elements = " << Map().NumGlobalElements() << std::endl;
        os << std::endl;

        for (int iproc=0; iproc < Comm().NumProc(); iproc++) 
        {
          if (iproc == 0)
          {
            os << "GEID\t";
            for (int i = 0; i < RowSize(); ++i) os << "V\t";
            os << std::endl;
          }

          if (Comm().MyPID() == iproc) 
          {
            for (int i = 0; i < Map().NumMyElements(); ++i)
            {
              os << Map().GID(i) << '\t';
              for (int j = 0; j < RowSize_; ++j)
                os << values_[i * RowSize_ + j] << '\t';
              os << std::endl;
            }
          }
        }
        os << std::endl;
      }

      int NextGID()
      {
        ++count_;
        if (count_ < Map().NumMyElements())
          return(Map().GID(count_));
        else
          return(-1);
      }

      int FirstGID()
      {
        count_ = 0;
        return(Map().GID(0));
      }

      //! Extracts a view of the array.
      const std::vector<T>& ExtractView() const
      {
        return(values_);
      }

      //! Returns a pointer to the internally stored data (non-const version).
      T* Values() 
      {
        return(&values_[0]);
      }

      //! Returns a pointer to the internally stored data (const version).
      const T* Values() const
      {
        return(&values_[0]);
      }

      // @}
    private:
      // @{ \name Epetra_DistObject methods

      virtual int CheckSizes(const Epetra_SrcDistObject& Source)
      {
        return(0);
      }

      virtual int CopyAndPermute(const Epetra_SrcDistObject& Source,
                                 int NumSameIDs,
                                 int NumPermuteIDs,
                                 int * PermuteToLIDs,
                                 int * PermuteFromLIDs,
                                 const Epetra_OffsetIndex * Indexor)
      {
        const DistArray& S = dynamic_cast<const DistArray&>(Source);
        const std::vector<T>& From = S.ExtractView();

        std::vector<T>& To = values_;

        //int * ToFirstPointInElementList = 0;
        //int * FromFirstPointInElementList = 0;
        //int * FromElementSizeList = 0;

        int j;

        int NumSameEntries;

        NumSameEntries = NumSameIDs;

        // Short circuit for the case where the source and target std::vector is the same.
        if (To==From) NumSameEntries = 0;

        // Do copy first
        if (NumSameIDs>0)
          if (To!=From) {
            for (j=0; j<NumSameEntries * RowSize_; j++)
            {
              To[j] = From[j];
            }
          }

        // Do local permutation next
        if (NumPermuteIDs>0) {

          for (j=0; j<NumPermuteIDs * RowSize_; j++) 
            To[PermuteToLIDs[j]] = From[PermuteFromLIDs[j]];
          // constant element size case
        }

        return(0);
      }

      virtual int PackAndPrepare(const Epetra_SrcDistObject& Source,
                                 int NumExportIDs,
                                 int* ExportLIDs,
                                 int& LenExports,
                                 char*& Exports,
                                 int& SizeOfPacket,
                                 int* Sizes,
                                 bool & VarSizes,
                                 Epetra_Distributor& Distor)
      {
        const DistArray& S = dynamic_cast<const DistArray&>(Source);
        const std::vector<T>& From = S.ExtractView();

        std::vector<T> To = values_;

        //int * FromFirstPointInElementList = 0;
        //int * FromElementSizeList = 0;

        SizeOfPacket = RowSize_ * sizeof(T); 

        if(NumExportIDs*SizeOfPacket>LenExports) {
          if (LenExports>0) delete [] Exports;
          LenExports = NumExportIDs*SizeOfPacket;
          Exports = new char[LenExports];
        }

        T* ptr;

        if (NumExportIDs>0) {
          ptr = (T*) Exports;

          // Point entry case
          for (int j=0; j<NumExportIDs; j++) 
            for (int k = 0; k < RowSize_ ; ++k)
              *ptr++ = From[ExportLIDs[j] * RowSize_ + k];
        }

        return(0);
      }

      virtual int UnpackAndCombine(const Epetra_SrcDistObject& Source,
                                   int NumImportIDs,
                                   int* ImportLIDs,
                                   int LenImports,
                                   char* Imports,
                                   int& SizeOfPacket,
                                   Epetra_Distributor& Distor,
                                   Epetra_CombineMode CombineMode,
                                   const Epetra_OffsetIndex * Indexor)
      {
        int j;

        if (CombineMode != Insert)
          EPETRA_CHK_ERR(-1); //Unsupported CombinedMode, will default to Zero

        std::cout << NumImportIDs << std::endl;
        if (NumImportIDs<=0) return(0);

        T* To = &values_[0];
        //int * ToFirstPointInElementList = 0;
        //int * ToElementSizeList = 0;

        T* ptr;
        // Unpack it...

        ptr = (T*) Imports;

        for (j=0; j<NumImportIDs; j++) 
          for (int k = 0; k < RowSize_ ; ++k)
            To[ImportLIDs[j] * RowSize_ + k] = *ptr++;

        return(0);
      }

      // @}
      // @{ \name Private data
      
      //! Container of local data
      std::vector<T> values_;
      //! Length of the locally owned array.
      int MyLength_;
      //! Length of the global array.
      int GlobalLength_;
      //! Amount of data associated with each element.
      int RowSize_;
      int count_;
      int last_;
      // @}
  };

} // namespace EpetraExt

#endif
