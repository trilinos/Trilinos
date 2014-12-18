// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
// @HEADER



#ifndef TPETRAUTILS_MATRIXGENERATOR_HPP
#define TPETRAUTILS_MATRIXGENERATOR_HPP

/// \file TpetraUtils_MatrixGenerator.hpp
/// \brief MiniFE test problem generator.
///
/// \warning If you don't know what this file is for, don't use it!

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_ComputeGatherMap.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>
#include <stdexcept>

namespace Tpetra {
  namespace Utils {

    template<class SparseMatrixType>
    class MatrixGenerator {
    public:
      //! This class' template parameter; a specialization of CrsMatrix.
      typedef SparseMatrixType sparse_matrix_type;

      /// Type of the entries of the sparse matrix.
      /// The first template parameter of CrsMatrix and MultiVector.
      typedef typename SparseMatrixType::scalar_type scalar_type;
      /// Type of the local indices of the sparse matrix.
      /// The second template parameter of CrsMatrix and MultiVector.
      typedef typename SparseMatrixType::local_ordinal_type local_ordinal_type;
      /// Type of the global indices of the sparse matrix.
      ///
      /// The third template parameter of CrsMatrix and MultiVector.
      /// This is also the type of indices as read from the Matrix
      /// Market file.  Indices of the sparse matrix are read in as
      /// global ordinals, since Matrix Market files represent the
      /// whole matrix and don't have a notion of distribution.
      typedef typename SparseMatrixType::global_ordinal_type
        global_ordinal_type;
      /// The Kokkos Node type.
      /// The fourth template parameter of CrsMatrix and MultiVector.
      typedef typename SparseMatrixType::node_type node_type;

      //! The MultiVector specialization associated with SparseMatrixType.
      typedef MultiVector<scalar_type,
                          local_ordinal_type,
                          global_ordinal_type,
                          node_type> multivector_type;

      //! The Vector specialization associated with SparseMatrixType.
      typedef Vector<scalar_type,
                     local_ordinal_type,
                     global_ordinal_type,
                     node_type> vector_type;

      typedef Teuchos::Comm<int> comm_type;
      typedef Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

    private:
      /// \typedef size_type
      /// \brief Handy typedef for entries of arrays such as rowPtr.
      typedef typename Teuchos::ArrayRCP<global_ordinal_type>::size_type size_type;

      static Teuchos::RCP<const map_type>
      makeRangeMap (const Teuchos::RCP<const comm_type>& pComm,
                    const Teuchos::RCP<node_type>& pNode,
                    const global_ordinal_type numRows)
      {
        using Teuchos::rcp;
        // A conventional, uniformly partitioned, contiguous map.
        return rcp (new map_type (static_cast<global_size_t> (numRows),
                                  static_cast<global_ordinal_type> (0),
                                  pComm, GloballyDistributed, pNode));
      }

      /// \brief Compute initial row map, or verify an existing one.
      ///
      /// The typical case when reading a sparse matrix from a file is
      /// for the reader itself to create a new row map, in particular
      /// a standard uniform contiguous one-to-one row map.  However,
      /// we also give the option to use an existing row map, if you
      /// are already using a particular distribution for (say) vector
      /// data and don't want to stop using it.  In the latter case
      /// (pRowMap is not null), we validate the communicator and node
      /// of the existing row map that you pass in.  In either case,
      /// you need to know the (global) number of rows in the matrix.
      ///
      /// \param pRowMap [in] If non-null, test pRowMap for validity,
      ///   and return it if valid.  Otherwise, if pRowMap is null,
      ///   initialize and return a (uniform contiguous one-to-one)
      ///   row map.  "Validity" here means that the map's
      ///   communicator and node are the same objects (pointerwise)
      ///   as the corresponding arguments.  (Note that the global
      ///   number of elements may not be the same as the number of
      ///   rows; a row map is not required to be one-to-one.)  The
      ///   typical case is to pass in null here, which is why we call
      ///   this routine "makeRowMap".
      /// \param pComm [in] Global communicator.
      /// \param pNode [in] Kokkos Node object.
      /// \param numRows [in] Global number of rows in the matrix.  If
      ///   pRowMap is nonnull, used only for error checking.
      ///
      /// \return If pRowMap is null, a new row map, otherwise pRowMap.
      static Teuchos::RCP<const map_type>
      makeRowMap (const Teuchos::RCP<const map_type>& pRowMap,
                  const Teuchos::RCP<const comm_type>& pComm,
                  const Teuchos::RCP<node_type>& pNode,
                  const global_ordinal_type numRows)
      {
        using Teuchos::rcp;
        // If the caller didn't provide a map, return a conventional,
        // uniformly partitioned, contiguous map.
        if (pRowMap.is_null()) {
          return rcp (new map_type (static_cast<global_size_t> (numRows),
                                    static_cast<global_ordinal_type> (0),
                                    pComm, GloballyDistributed, pNode));
        } else {
          TEUCHOS_TEST_FOR_EXCEPTION(
            ! pRowMap->isDistributed () && pComm->getSize () > 1,
            std::invalid_argument, "The specified row Map is not distributed, "
            "but the given communicator includes more than one process (in "
            "fact, there are " << pComm->getSize () << " processes).");
          // FIXME (mfh 25 Jun 2014) The communicators really don't
          // have to be identical; they only have to have the same
          // number of processes in the same order.
          TEUCHOS_TEST_FOR_EXCEPTION(
            pRowMap->getComm () != pComm, std::invalid_argument,
            "The specified row map's communicator (pRowMap->getComm()) is "
            "different than the given separately supplied communicator pComm.");
          // FIXME (mfh 25 Jun 2014) This really doesn't matter either.
          TEUCHOS_TEST_FOR_EXCEPTION(
            pRowMap->getNode () != pNode, std::invalid_argument,
            "The specified row map's node (pRowMap->getNode()) is different "
            "than the given separately supplied node pNode.");
          return pRowMap;
        }
      }

      /// \brief Compute domain map.
      ///
      /// Domain maps must always be one-to-one.  We will use this map
      /// when we call fillComplete() on the CrsMatrix that the reader
      /// constructs.
      ///
      /// \param pRangeMap [in] Valid range map of the matrix,
      ///   as returned by \c makeRangeMap().
      /// \param numRows [in] Global number of rows in the matrix.
      /// \param numCols [in] Global number of columns in the matrix.
      ///
      /// \return The domain map.  If numRows == numCols, this is
      ///   identical to the range map, otherwise we make a new map
      ///   for the domain.
      static Teuchos::RCP<const map_type>
      makeDomainMap (const Teuchos::RCP<const map_type>& pRangeMap,
                     const global_ordinal_type numRows,
                     const global_ordinal_type numCols)
      {
        // Abbreviations so that the map creation call isn't too long.
        typedef local_ordinal_type LO;
        typedef global_ordinal_type GO;
        typedef node_type NT;

        if (numRows == numCols) {
          return pRangeMap;
        } else {
          return createUniformContigMapWithNode<LO,GO,NT> (numCols,
                                                           pRangeMap->getComm (),
                                                           pRangeMap->getNode ());
        }
      }



      /// \brief Given my proc's data, return the completed sparse matrix.
      ///
      /// Each proc inserts its data into the sparse matrix, and then,
      /// if callFillComplete is true, all procs call fillComplete().
      /// (For whatever reason, you might not be done with the matrix
      /// yet, so you might want to call fillComplete() yourself.
      /// CrsMatrix::fillResume() doesn't currently work as you might
      /// expect when storage optimization is enabled; it fixes the
      /// graph of the matrix, so that you can't add new entries.)
      ///
      /// Column indices are zero-based on input.  This method will
      /// change them in place to match the index base of the input
      /// row Map (\c pRowMap).
      static Teuchos::RCP<sparse_matrix_type>
      makeMatrix (Teuchos::ArrayRCP<size_t>& myNumEntriesPerRow,
                  Teuchos::ArrayRCP<size_t>& myRowPtr,
                  Teuchos::ArrayRCP<global_ordinal_type>& myColInd,
                  Teuchos::ArrayRCP<scalar_type>& myValues,
                  const Teuchos::RCP<const map_type>& pRowMap,
                  const Teuchos::RCP<const map_type>& pRangeMap,
                  const Teuchos::RCP<const map_type>& pDomainMap,
                  const bool callFillComplete = true)
      {
        using Teuchos::ArrayView;
        using Teuchos::RCP;
        using Teuchos::rcp;
        using std::cerr;
        using std::endl;
        // Typedef to make certain type declarations shorter.
        typedef global_ordinal_type GO;

        // The row pointer array always has at least one entry, even
        // if the matrix has zero rows.  myNumEntriesPerRow, myColInd,
        // and myValues would all be empty arrays in that degenerate
        // case, but the row and domain maps would still be nonnull
        // (though they would be trivial maps).
        TEUCHOS_TEST_FOR_EXCEPTION(myRowPtr.is_null(), std::logic_error,
          "makeMatrix: myRowPtr array is null.  "
          "Please report this bug to the Tpetra developers.");
        TEUCHOS_TEST_FOR_EXCEPTION(pDomainMap.is_null(), std::logic_error,
          "makeMatrix: domain map is null.  "
          "Please report this bug to the Tpetra developers.");
        TEUCHOS_TEST_FOR_EXCEPTION(pRangeMap.is_null(), std::logic_error,
          "makeMatrix: range map is null.  "
          "Please report this bug to the Tpetra developers.");
        TEUCHOS_TEST_FOR_EXCEPTION(pRowMap.is_null(), std::logic_error,
          "makeMatrix: row map is null.  "
          "Please report this bug to the Tpetra developers.");

        // Construct the CrsMatrix, using the row map, with the
        // constructor specifying the number of nonzeros for each row.
        // Create with DynamicProfile, so that the fillComplete() can
        // do first-touch reallocation (a NUMA (Non-Uniform Memory
        // Access) optimization on multicore CPUs).
        RCP<sparse_matrix_type> A =
          rcp (new sparse_matrix_type (pRowMap, myNumEntriesPerRow,
                                       DynamicProfile));

        // List of the global indices of my rows.
        // They may or may not be contiguous.
        ArrayView<const GO> myRows = pRowMap->getNodeElementList ();
        const size_type myNumRows = myRows.size ();

        // Add this processor's matrix entries to the CrsMatrix.
        const GO indexBase = pRowMap->getIndexBase ();
        for (size_type i = 0; i < myNumRows; ++i) {
          const size_type myCurPos = myRowPtr[i];
          const local_ordinal_type curNumEntries = myNumEntriesPerRow[i];
          ArrayView<GO> curColInd = myColInd.view (myCurPos, curNumEntries);
          ArrayView<scalar_type> curValues = myValues.view (myCurPos, curNumEntries);

          // Modify the column indices in place to have the right index base.
          for (size_type k = 0; k < curNumEntries; ++k) {
            curColInd[k] += indexBase;
          }
          // Avoid constructing empty views of ArrayRCP objects.
          if (curNumEntries > 0) {
            A->insertGlobalValues (myRows[i], curColInd, curValues);
          }
        }
        // We've entered in all our matrix entries, so we can delete
        // the original data.  This will save memory when we call
        // fillComplete(), so that we never keep more than two copies
        // of the matrix's data in memory at once.
        myNumEntriesPerRow = null;
        myRowPtr = null;
        myColInd = null;
        myValues = null;

        if (callFillComplete) {
          A->fillComplete (pDomainMap, pRangeMap);
        }
        return A;
      }

      /// \brief Variant of makeMatrix() that takes parameters for
      ///   CrsMatrix's constructor and for fillComplete().
      ///
      /// Each process inserts its data into the sparse matrix, and
      /// then all processes call fillComplete().
      static Teuchos::RCP<sparse_matrix_type>
      makeMatrix (Teuchos::ArrayRCP<size_t>& myNumEntriesPerRow,
                  Teuchos::ArrayRCP<size_t>& myRowPtr,
                  Teuchos::ArrayRCP<global_ordinal_type>& myColInd,
                  Teuchos::ArrayRCP<scalar_type>& myValues,
                  const Teuchos::RCP<const map_type>& pRowMap,
                  const Teuchos::RCP<const map_type>& pRangeMap,
                  const Teuchos::RCP<const map_type>& pDomainMap,
                  const Teuchos::RCP<Teuchos::ParameterList>& constructorParams,
                  const Teuchos::RCP<Teuchos::ParameterList>& fillCompleteParams)
      {
        using Teuchos::ArrayView;
        using Teuchos::RCP;
        using Teuchos::rcp;
        using std::cerr;
        using std::endl;
        // Typedef to make certain type declarations shorter.
        typedef global_ordinal_type GO;

        // The row pointer array always has at least one entry, even
        // if the matrix has zero rows.  myNumEntriesPerRow, myColInd,
        // and myValues would all be empty arrays in that degenerate
        // case, but the row and domain maps would still be nonnull
        // (though they would be trivial maps).
        TEUCHOS_TEST_FOR_EXCEPTION(
          myRowPtr.is_null(), std::logic_error,
          "makeMatrix: myRowPtr array is null.  "
          "Please report this bug to the Tpetra developers.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          pDomainMap.is_null(), std::logic_error,
          "makeMatrix: domain map is null.  "
          "Please report this bug to the Tpetra developers.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          pRangeMap.is_null(), std::logic_error,
          "makeMatrix: range map is null.  "
          "Please report this bug to the Tpetra developers.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          pRowMap.is_null(), std::logic_error,
          "makeMatrix: row map is null.  "
          "Please report this bug to the Tpetra developers.");

        // Construct the CrsMatrix, using the row map, with the
        // constructor specifying the number of nonzeros for each row.
        // Create with DynamicProfile, so that the fillComplete() can
        // do first-touch reallocation (a NUMA (Non-Uniform Memory
        // Access) optimization on multicore CPUs).
        RCP<sparse_matrix_type> A =
          rcp (new sparse_matrix_type (pRowMap, myNumEntriesPerRow,
                                       DynamicProfile, constructorParams));

        // List of the global indices of my rows.
        // They may or may not be contiguous.
        ArrayView<const GO> myRows = pRowMap->getNodeElementList();
        const size_type myNumRows = myRows.size();

        // Add this processor's matrix entries to the CrsMatrix.
        const GO indexBase = pRowMap->getIndexBase ();
        for (size_type i = 0; i < myNumRows; ++i) {
          const size_type myCurPos = myRowPtr[i];
          const local_ordinal_type curNumEntries = myNumEntriesPerRow[i];
          ArrayView<GO> curColInd = myColInd.view (myCurPos, curNumEntries);
          ArrayView<scalar_type> curValues = myValues.view (myCurPos, curNumEntries);

          // Modify the column indices in place to have the right index base.
          for (size_type k = 0; k < curNumEntries; ++k) {
            curColInd[k] += indexBase;
          }
          if (curNumEntries > 0) {
            A->insertGlobalValues (myRows[i], curColInd, curValues);
          }
        }
        // We've entered in all our matrix entries, so we can delete
        // the original data.  This will save memory when we call
        // fillComplete(), so that we never keep more than two copies
        // of the matrix's data in memory at once.
        myNumEntriesPerRow = null;
        myRowPtr = null;
        myColInd = null;
        myValues = null;

        A->fillComplete (pDomainMap, pRangeMap, fillCompleteParams);
        return A;
      }

      /// \brief Variant of makeMatrix() that takes an optional column Map.
      ///
      /// This method computes \c colMap only if it is null on input,
      /// and if \c callFillComplete is true.
      static Teuchos::RCP<sparse_matrix_type>
      makeMatrix (Teuchos::ArrayRCP<size_t>& myNumEntriesPerRow,
                  Teuchos::ArrayRCP<size_t>& myRowPtr,
                  Teuchos::ArrayRCP<global_ordinal_type>& myColInd,
                  Teuchos::ArrayRCP<scalar_type>& myValues,
                  const Teuchos::RCP<const map_type>& rowMap,
                  Teuchos::RCP<const map_type>& colMap,
                  const Teuchos::RCP<const map_type>& domainMap,
                  const Teuchos::RCP<const map_type>& rangeMap,
                  const bool callFillComplete = true)
      {
        using Teuchos::ArrayView;
        using Teuchos::as;
        using Teuchos::null;
        using Teuchos::RCP;
        using Teuchos::rcp;
        typedef global_ordinal_type GO;

        // Construct the CrsMatrix.
        //
        // Create with DynamicProfile, so that the fillComplete() can
        // do first-touch reallocation.
        RCP<sparse_matrix_type> A; // the matrix to return.
        if (colMap.is_null ()) { // the user didn't provide a column Map
          A = rcp (new sparse_matrix_type (rowMap, myNumEntriesPerRow, DynamicProfile));
        } else { // the user provided a column Map
          A = rcp (new sparse_matrix_type (rowMap, colMap, myNumEntriesPerRow, DynamicProfile));
        }

        // List of the global indices of my rows.
        // They may or may not be contiguous.
        ArrayView<const GO> myRows = rowMap->getNodeElementList ();
        const size_type myNumRows = myRows.size ();

        // Add this process' matrix entries to the CrsMatrix.
        const GO indexBase = rowMap->getIndexBase ();
        for (size_type i = 0; i < myNumRows; ++i) {
          const size_type myCurPos = myRowPtr[i];
          const size_type curNumEntries = as<size_type> (myNumEntriesPerRow[i]);
          ArrayView<GO> curColInd = myColInd.view (myCurPos, curNumEntries);
          ArrayView<scalar_type> curValues = myValues.view (myCurPos, curNumEntries);

          // Modify the column indices in place to have the right index base.
          for (size_type k = 0; k < curNumEntries; ++k) {
            curColInd[k] += indexBase;
          }
          if (curNumEntries > 0) {
            A->insertGlobalValues (myRows[i], curColInd, curValues);
          }
        }
        // We've entered in all our matrix entries, so we can delete
        // the original data.  This will save memory when we call
        // fillComplete(), so that we never keep more than two copies
        // of the matrix's data in memory at once.
        myNumEntriesPerRow = null;
        myRowPtr = null;
        myColInd = null;
        myValues = null;

        if (callFillComplete) {
          A->fillComplete (domainMap, rangeMap);
          if (colMap.is_null ()) {
            colMap = A->getColMap ();
          }
        }
        return A;
      }

      template<class GO, class S>
      static void
      miniFE_get_row (size_t* rows, S* vals, GO* cols, size_t startrow,
                      size_t endrow, size_t& row,size_t o,size_t nx1,
                      size_t c1,size_t c2, size_t c3,size_t val,
                      size_t &miniFE_a,size_t &miniFE_b,size_t &miniFE_c)
      {
        // FIXME (mfh 25 Jun 2014) Seriously, "val27"???  Who writes
        // code like this???

        bool val27=false;
        if(c1*c2*c3==27) {
          val27=true;
        }

        if((row>=startrow)&&(row<endrow)) {
         size_t offset = rows[row-startrow];
         rows[row+1-startrow] = offset+c1*c2*c3;
         for(size_t i=0; i<c1; i++)
           for(size_t j=0;j<c2; j++)
             for(size_t k=0;k<c3;k++) {
               size_t m = i*c2*c3+j*c2+k;
               cols[offset+m] = o+i*nx1*nx1+j*nx1+k;
                 if(val27) {
                    bool doa = ((miniFE_a>0)&&(miniFE_a<nx1-3)) || ((miniFE_a==0)&&(m/9>=1)) || ((miniFE_a==nx1-3)&&(m/9<2));
                    bool dob = ((miniFE_b>0)&&(miniFE_b<nx1-3)) || ((miniFE_b==0)&&((m%9)/3>=1)) || ((miniFE_b==nx1-3)&&((m%9)/3<2));
                    bool doc = ((miniFE_c>0)&&(miniFE_c<nx1-3)) || ((miniFE_c==0)&&((m%3)>=1)) || ((miniFE_c==nx1-3)&&((m%3)<2));
                    if(doa&&dob&&doc) {
                       if(m==13)
                         vals[offset+m] = 8.0/3.0/(nx1-1);
                       else {
                         if(m%2==1)
                           vals[offset+m] = -5.0e-1/3.0/(nx1-1);
                         else {
                           if((m==4)||(m==22)|| ((m>9)&&(m<17)))
                             vals[offset+m] = -2.18960e-10/(nx1-1);
                           else
                             vals[offset+m] = -2.5e-1/3.0/(nx1-1);
                         }
                       }
                    } else vals[offset+m] = 0.0;
                 } else {
                   if(val==m)
                     vals[offset+m] = 1.0;
                   else
                     vals[offset+m] = 0.0;
                 }
            }
        }
         if(c1*c2*c3==27) {
           miniFE_c++;
           if(miniFE_c>nx1-3) {miniFE_c=0; miniFE_b++;}
           if(miniFE_b>nx1-3) {miniFE_b=0; miniFE_a++;}
         }

         row++;
       }

       template<class GO, class S>
       static void miniFE_get_block(size_t* rows, S* vals, GO* cols,size_t startrow, size_t endrow,size_t& row , size_t o, size_t nx1, size_t c1, size_t c2,size_t val1,size_t val2,size_t val3,size_t &miniFE_a,size_t &miniFE_b,size_t &miniFE_c) {
         miniFE_get_row(rows,vals,cols,startrow,endrow,row,o,nx1,c1,c2,2,val1,miniFE_a,miniFE_b,miniFE_c);
         for(size_t i=0;i<nx1-2;i++)
           miniFE_get_row(rows,vals,cols,startrow,endrow,row,o++,nx1,c1,c2,3,val2,miniFE_a,miniFE_b,miniFE_c);
         miniFE_get_row(rows,vals,cols,startrow,endrow,row,o++,nx1,c1,c2,2,val3,miniFE_a,miniFE_b,miniFE_c);
       }

       template<class GO, class S>
       static void miniFE_get_superblock(size_t* rows, S* vals, GO* cols,size_t startrow, size_t endrow,size_t& row , size_t o, size_t nx1, size_t c1,size_t val1,size_t val2,size_t val3,size_t &miniFE_a,size_t &miniFE_b,size_t &miniFE_c) {
         miniFE_get_block(rows,vals,cols,startrow,endrow,row,o,nx1,c1,2,val1+0,val1+val2+1,val1+1,miniFE_a,miniFE_b,miniFE_c);
         for(size_t i=0;i<nx1-2;i++){
           miniFE_get_block(rows,vals,cols,startrow,endrow,row,o,nx1,c1,3,val1+val2+3,val1+val2+val2+val3+4,val1+val2+4,miniFE_a,miniFE_b,miniFE_c);
           o+=nx1;
         }
         miniFE_get_block(rows,vals,cols,startrow,endrow,row,o,nx1,c1,2,val1+2,val1+val2+3,val1+3,miniFE_a,miniFE_b,miniFE_c);
       }

    public:


      static Teuchos::RCP<SparseMatrixType>
      generate_miniFE_matrix (int nx,
                              const Teuchos::RCP<const Comm<int> >& pComm,
                              const Teuchos::RCP<node_type>& pNode,
                              const bool callFillComplete=true,
                              const bool debug = false)
      {
        using Teuchos::ArrayRCP;
        using Teuchos::RCP;

        size_t miniFE_a = 0;
        size_t miniFE_b = 0;
        size_t miniFE_c = 0;

        const int myRank = pComm->getRank ();
        const int rootRank = 0;

        size_t nx1 = nx+1;

        int nrows_block = 1+nx-1+1;
        int nrows_superblock = (1+nx-1+1)*nrows_block;
        int nrows = (1+(nx-1)+1)*nrows_superblock;

        size_t nnz=0;
        nnz+=4*(8 + (nx-1)*12 + 8);
        nnz+=4*(nx-1)*(12 + (nx-1)*18 + 12);
        nnz+=(nx-1)*(nx-1)*(18 + (nx-1)*27 + 18);

        size_t dims[3];
        dims[0] = nrows;
        dims[1] = nrows;
        dims[2] = nnz;

        Teuchos::RCP<const map_type> pRangeMap = makeRangeMap (pComm, pNode, dims[0]);
        Teuchos::RCP<const map_type> pDomainMap = makeDomainMap (pRangeMap, dims[0], dims[1]);
        Teuchos::RCP<const map_type> pRowMap = makeRowMap (null, pComm, pNode, dims[0]);

        size_t startrow = pRowMap->getMinGlobalIndex();
        size_t endrow = pRowMap->getMaxGlobalIndex()+1;

        ArrayRCP<size_t> numEntriesPerRow(endrow-startrow);
        ArrayRCP<size_t> rowPtr(endrow-startrow+1);
        ArrayRCP<global_ordinal_type> colInd((endrow-startrow)*27);
        ArrayRCP<scalar_type> values((endrow-startrow)*27);

        size_t* rows = &rowPtr[0];
        scalar_type* vals = &values[0];
        global_ordinal_type* cols = &colInd[0];


        size_t row = 0;
        miniFE_get_superblock(rows,vals,cols,startrow,endrow,row,0,nx1,2,0,0,0,miniFE_a,miniFE_b,miniFE_c);
        for(size_t i=0;i<nx1-2;i++){
          miniFE_get_superblock(rows,vals,cols,startrow,endrow,row,i*nx1*nx1,nx1,3,4,2,1,miniFE_a,miniFE_b,miniFE_c);
        }
        miniFE_get_superblock(rows,vals,cols,startrow,endrow,row,(nx1-2)*nx1*nx1,nx1,2,4,2,1,miniFE_a,miniFE_b,miniFE_c);

        for(size_t i=0;i<endrow-startrow;i++)
          numEntriesPerRow[i]=rowPtr[i+1]-rowPtr[i];

        // Distribute the matrix data.  Each processor has to add the
        // rows that it owns.  If you try to make Proc 0 call
        // insertGlobalValues() for _all_ the rows, not just those it
        // owns, then fillComplete() will compute the number of
        // columns incorrectly.  That's why Proc 0 has to distribute
        // the matrix data and why we make all the processors (not
        // just Proc 0) call insertGlobalValues() on their own data.
        //
        // These arrays represent each processor's part of the matrix
        // data, in "CSR" format (sort of, since the row indices might
        // not be contiguous).
        /*ArrayRCP<size_t> myNumEntriesPerRow;
          ArrayRCP<size_t> myRowPtr;
          ArrayRCP<global_ordinal_type> myColInd;
          ArrayRCP<scalar_type> myValues;
          // Distribute the matrix data.  This is a collective operation.
          distribute (myNumEntriesPerRow, myRowPtr, myColInd, myValues, pRowMap,
          numEntriesPerRow, rowPtr, colInd, values, 0);
          RCP<sparse_matrix_type> pMatrix =
          makeMatrix (myNumEntriesPerRow, myRowPtr, myColInd, myValues,
          pRowMap, pRangeMap, pDomainMap, callFillComplete);*/
        RCP<sparse_matrix_type> pMatrix =
          makeMatrix (numEntriesPerRow, rowPtr, colInd, values,
                      pRowMap, pRangeMap, pDomainMap, callFillComplete);
        // Only use a reduce-all in debug mode to check if pMatrix is
        // null.  Otherwise, just throw an exception.  We never expect
        // a null pointer here, so we can save a communication.
        /*if (debug) {
          int localIsNull = pMatrix.is_null () ? 1 : 0;
          int globalIsNull = 0;
          reduceAll (*pComm, REDUCE_MAX, localIsNull, ptr (&globalIsNull));
          TEUCHOS_TEST_FOR_EXCEPTION(globalIsNull != 0, std::logic_error,
          "Reader::makeMatrix() returned a null pointer on at least one "
          "process.  Please report this bug to the Tpetra developers.");
          }
          else {*/
        TEUCHOS_TEST_FOR_EXCEPTION(pMatrix.is_null(), std::logic_error,
                                   "Reader::makeMatrix() returned a null pointer.  "
                                   "Please report this bug to the Tpetra developers.");
        //}

        // We can't get the dimensions of the matrix until after
        // fillComplete() is called.  Thus, we can't do the sanity
        // check (dimensions read from the Matrix Market data,
        // vs. dimensions reported by the CrsMatrix) unless the user
        // asked makeMatrix() to call fillComplete().
        //
        // Note that pMatrix->getGlobalNum{Rows,Cols}() does _not_ do
        // what one might think it does, so you have to ask the range
        // resp. domain map for the number of rows resp. columns.
        if (callFillComplete) {
          const int numProcs = pComm->getSize ();

          if (false && debug) {
            const size_t globalNumRows =
              pRangeMap->getGlobalNumElements();
            const size_t globalNumCols =
              pDomainMap->getGlobalNumElements();
            if (myRank == rootRank) {
              std::cerr << "-- Matrix is "
                        << globalNumRows << " x " << globalNumCols
                        << " with " << pMatrix->getGlobalNumEntries()
                        << " entries, and index base "
                        << pMatrix->getIndexBase() << "." << std::endl;
            }
            pComm->barrier ();
            for (int p = 0; p < numProcs; ++p) {
              if (myRank == p) {
                std::cerr << "-- Proc " << p << " owns "
                          << pMatrix->getNodeNumCols() << " columns, and "
                          << pMatrix->getNodeNumEntries() << " entries." << std::endl;
              }
              pComm->barrier ();
            }
          } // if (extraDebug && debug)
        } // if (callFillComplete)

        if (debug && myRank == rootRank) {
          std::cerr << "-- Done creating the CrsMatrix from the Matrix Market data"
                    << std::endl;
        }
        return pMatrix;
      }

    private:
       template<class S>
       static void miniFE_vector_generate_block(S* vec, int nx, S a, S b, int& count,int start,int end){
         if((count>=start) && (count<end))
         vec[count++ - start] = 0;
         for(int i=0; i<nx-2; i++)
           if((count>=start) && (count<end))
             vec[count++ - start] = a/nx/nx/nx;
         if((count>=start) && (count<end))
         vec[count++ - start] = a/nx/nx/nx + b/nx;
         if((count>=start) && (count<end))
         vec[count++ - start] = 1;
       }

       template<class S>
       static void miniFE_vector_generate_superblock(S* vec, int nx, S a,S b,S c, int& count,int start,int end){
         miniFE_vector_generate_block(vec,nx,0.0,0.0,count,start,end);
         miniFE_vector_generate_block(vec,nx,a,b,count,start,end);
         for(int i = 0;i<nx-3;i++)
           miniFE_vector_generate_block(vec,nx,a,c,count,start,end);
         miniFE_vector_generate_block(vec,nx,a,b,count,start,end);
         miniFE_vector_generate_block(vec,nx,0.0,0.0,count,start,end);
       }

    public:
       static Teuchos::RCP<Tpetra::Vector<scalar_type,
                                            local_ordinal_type,
                                            global_ordinal_type,
                                            node_type> >
       generate_miniFE_vector(
           int nx,
           const Teuchos::RCP<const Teuchos::Comm<int> >& pComm,
           const Teuchos::RCP<node_type>& pNode
           ) {
         typedef scalar_type ST;
         typedef local_ordinal_type LO;
         typedef global_ordinal_type GO;
         typedef node_type NT;
         //typedef Teuchos::ScalarTraits<ST> STS; // unused
         //typedef typename STS::magnitudeType MT; // unused
         //typedef Teuchos::ScalarTraits<MT> STM; // unused
         typedef Tpetra::Vector<ST, LO, GO, NT> MV;

         Tuple<GO, 3> dims;
         dims[0] = (nx+1)*(nx+1)*(nx+1);
         dims[1] = 1;

         const global_size_t numRows = static_cast<global_size_t> (dims[0]);
         // const size_t numCols = static_cast<size_t> (dims[1]);

         RCP<const map_type> map = createUniformContigMapWithNode<LO, GO, NT> (numRows, pComm, pNode);
         int start = map->getMinGlobalIndex();
         int end = map->getMaxGlobalIndex()+1;

         // Make a multivector X owned entirely by Proc 0.
         RCP<MV> X = createVector<ST, LO, GO, NT> (map);
         ArrayRCP<ST> X_view = X->get1dViewNonConst ();
         ST* vec = &X_view[0];
         int count = 0;
         miniFE_vector_generate_superblock(vec,nx,0.0,0.0,0.0,count,start,end);
         miniFE_vector_generate_superblock(vec,nx,1.0,5.0/12,8.0/12,count,start,end);
         for(int i = 0; i<nx-3; i++)
           miniFE_vector_generate_superblock(vec,nx,1.0,8.0/12,1.0,count,start,end);
         miniFE_vector_generate_superblock(vec,nx,1.0,5.0/12,8.0/12,count,start,end);
         miniFE_vector_generate_superblock(vec,nx,0.0,0.0,0.0,count,start,end);

         return X;
       }

    };
  }
}
#endif
