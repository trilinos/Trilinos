// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// Questions? Contact Jennifer A. Loe (jloe@sandia.gov)
//                or Heidi K. Thornquist (hkthorn@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
#ifndef BELOS_TEUCHOS_DENSE_MAT_TRAITS_HPP
#define BELOS_TEUCHOS_DENSE_MAT_TRAITS_HPP
//TODO This isn't really unique to epetra since we use arbitrary scalar type.
// Should not be deprecated with epetra. Needed for backward compat.
// Where should this file live? 

/*! \file BelosTeuchosDenseMatrixTraits.hpp
  \brief Full specialization of Belos::DenseMatTraits for Teuchos::SerialDenseMatrix
  with ordinal type int and arbitrary scalar type.
*/

//#include "BelosConfigDefs.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Belos {

  //! Full specialization of Belos::DenseMatTraits for Teuchos::SerialDenseMatrix<int,ST>.
  template<class ScalarType>
  class DenseMatTraits<ScalarType, Teuchos::SerialDenseMatrix<int,ScalarType>>{
  public:
    
    //@{ \name Creation methods

    /*! \brief Creates a new empty \c Teuchos::SerialDenseMatrix<int,ScalarType> containing \c numvecs columns.
     *         Will be initialized to zeros if last parameter is true.

    \return Reference-counted pointer to a new dense matrix of type \c Teuchos::SerialDenseMatrix<int,ScalarType>.
    */
    static Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType>> Create( const int numrows, const int numcols, bool initZero = true) { 
      return Teuchos::rcp(new Teuchos::SerialDenseMatrix<int,ScalarType>(numrows,numcols,initZero));
    }     

/* Kokkos Ex View-from-ptr constructor: 
View(const pointer_type &ptr, const IntType&... indices)

    Unmanaged data wrapping constructor.

        ptr: pointer to a user provided memory allocation. Must provide storage of size View::required_allocation_size(n0,...,nR)
        indices: Extents of the View.
        Requires: sizeof(IntType...)==rank_dynamic() or sizeof(IntType...)==rank(). In the latter case, the extents corresponding to compile-time dimensions must match the View type’s compile-time extents.
        Requires: array_layout::is_regular == true.
*/
    static Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType>> ViewFromPtr( ScalarType* &ptr, const int numrows, const int numcols)
    { return Teuchos::null; }     

    //! \brief Returns a raw pointer to the data on the host.
    static ScalarType* GetRawHostPtr(const Teuchos::SerialDenseMatrix<int,ScalarType> & dm )
    { return Teuchos::null; }     

    //! \brief Returns an RCP to a Teuchos::SerialDenseMatrix<int,ScalarType> which has a subview of the given Teuchos::SerialDenseMatrix<int,ScalarType>.
    //        Row and column indexing is zero-based.
    //        Row and column ranges include the first index and exclude the second index.
    //        So, for example, giving std::make_pair(5,7) for the rowRange will return rows 
    //        in range [5,7), so rows 5 and 6. 
    //        This follows the convention of Kokkos::subview.
    static Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType>> Subview( const Teuchos::SerialDenseMatrix<int,ScalarType> & dm, std::pair<int,int> rowRange, std::pair<int,int> colRange)
    { return Teuchos::null; }     
    //@}

    //@{ \name Attribute methods

    //! \brief Obtain the number of rows of \c dm.
    static int GetNumRows( const Teuchos::SerialDenseMatrix<int,ScalarType>& dm )
    { return dm.numRows(); }     

    //! \brief Obtain the number of columns of \c dm.
    static int GetNumCols( const Teuchos::SerialDenseMatrix<int,ScalarType>& dm )
    { return dm.numCols(); }     

    //@}

    //@{ \name Shaping methods

    /* \brief Reshaping method for changing the size of \c dm,
    *         keeping the entries. 
    \return Integer error code, set to 0 if successful.
    */
    static int Reshape( Teuchos::SerialDenseMatrix<int,ScalarType>& dm, const int numrows, const int numcols)
    { return 0; }     

    //@}

    //@{ \name Data access methods

    //! \brief Access a reference to the (i,j) entry of \c dm, \c e_i^T dm e_j.
    static ScalarType & Value( Teuchos::SerialDenseMatrix<int,ScalarType>& dm, const int i, const int j )
    { 
      return dm(i,j);
    }

    //! \brief Access a const reference to the (i,j) entry of \c dm, \c e_i^T dm e_j.
    static const ScalarType & Value( const Teuchos::SerialDenseMatrix<int,ScalarType>& dm, const int i, const int j )
    { 
      return dm(i,j);
    }

    static void SyncHostToDevice(Teuchos::SerialDenseMatrix<int,ScalarType> &)
    { }

    static void SyncDeviceToHost(Teuchos::SerialDenseMatrix<int,ScalarType> &)
    { }
    //@}
    //@{ \name Operator methods
    
    //!  \brief Adds sourceDM to thisDM and returns answer in thisDM.
    static void Add( Teuchos::SerialDenseMatrix<int,ScalarType>& thisDM, const Teuchos::SerialDenseMatrix<int,ScalarType>& sourceDM)
    { }

    //!  \brief Multiply all entries by a scalar. DM = value.*DM
    static void Scale( Teuchos::SerialDenseMatrix<int,ScalarType>& dm, ScalarType value)
    { }

    //!  \brief Fill the DM with random entries.
    //!   Entries are assumed to be the same on each MPI rank (each matrix copy). 
    static void Randomize( Teuchos::SerialDenseMatrix<int,ScalarType>& dm)
    { }

    //!  \brief Copies entries of source to dest (deep copy). 
    static void Assign( Teuchos::SerialDenseMatrix<int,ScalarType>& dest, const Teuchos::SerialDenseMatrix<int,ScalarType>& source)
    { }
    //@}
  };
  
} // namespace Belos

#endif // end file BELOS_TEUCHOS_DENSE_MAT_TRAITS_HPP
