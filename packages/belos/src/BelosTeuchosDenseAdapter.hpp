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

/*! \file BelosTeuchosDenseAdapter.hpp
  \brief Full specialization of Belos::DenseMatTraits for Teuchos::SerialDenseMatrix
  with ordinal type int and arbitrary scalar type.
*/

#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "BelosDenseMatTraits.hpp" 

namespace Belos {

  //! Full specialization of Belos::DenseMatTraits for Teuchos::SerialDenseMatrix<int,ST>.
  template<class ScalarType>
  class DenseMatTraits<ScalarType, Teuchos::SerialDenseMatrix<int,ScalarType>>{
  public:
    
    //@{ \name Creation methods

    /*! \brief Creates a new empty \c Teuchos::SerialDenseMatrix<int,ScalarType> with no dimension.

    \return Reference-counted pointer to a new dense matrix of type \c Teuchos::SerialDenseMatrix<int,ScalarType>.
    */
    static Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType>> Create() { 
      return Teuchos::rcp(new Teuchos::SerialDenseMatrix<int,ScalarType>());
    }     

    /*! \brief Creates a new empty \c Teuchos::SerialDenseMatrix<int,ScalarType> containing \c numvecs columns.
     *         Will be initialized to zeros if last parameter is true.

    \return Reference-counted pointer to a new dense matrix of type \c Teuchos::SerialDenseMatrix<int,ScalarType>.
    */
    static Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType>> Create( const int numrows, const int numcols, bool initZero = true) { 
      return Teuchos::rcp(new Teuchos::SerialDenseMatrix<int,ScalarType>(numrows,numcols,initZero));
    }     

    /*! \brief Create a new copy \c Teuchos::SerialDenseMatrix<int,ScalarType>, possibly transposed.
   
       \return Reference-counted pointer to a new dense matrix of type \c Teuchos::SerialDenseMatrix<int,ScalarType>.
    */
    static Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType>> CreateCopy(const Teuchos::SerialDenseMatrix<int,ScalarType> & dm, bool transpose=false) {
      if (transpose)    
        return Teuchos::rcp(new Teuchos::SerialDenseMatrix<int,ScalarType>(dm, Teuchos::CONJ_TRANS));
      else 
        return Teuchos::rcp(new Teuchos::SerialDenseMatrix<int,ScalarType>(dm, Teuchos::NO_TRANS));
    }

    //! \brief Returns a raw pointer to the (non-const) data on the host.
    static ScalarType* GetRawHostPtr(const Teuchos::SerialDenseMatrix<int,ScalarType> & dm )
    { return dm.values(); }     

    //! \brief Returns a raw pointer to const data on the host.
    static ScalarType const * GetConstRawHostPtr(const Teuchos::SerialDenseMatrix<int,ScalarType> & dm )  
    { return const_cast<ScalarType const *>(dm.values()); }     

    //! \brief Marks host data modified to avoid device sync errors. 
    /// \note Belos developers must call this function after EVERY
    ///   call to LAPACK!!!
    //static void RawPtrDataModified(Teuchos::SerialDenseMatrix<int,ScalarType> & dm ) { }

    //! \brief Returns an RCP to a DM which has a subview of the given DM.
    //        Row and column indexing is zero-based.
    //        Source  - Reference to another dense matrix from which values are to be copied.
    //        numRows - The number of rows in this matrix.
    //        numCols - The number of columns in this matrix.
    //        startRow  - The row of Source from which the submatrix copy should start.
    //        startCol  - The column of Source from which the submatrix copy should start.
    //  
    //        Should ints be const? Should they be ints or some other ordinal type?
    static Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType>> Subview(Teuchos::SerialDenseMatrix<int,ScalarType> & source, int numRows, int numCols, int startRow=0, int startCol=0)
    { return Teuchos::rcp(new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View, source, numRows, numCols, startRow, startCol)); }    

    static Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,ScalarType>> SubviewConst( 
                              const Teuchos::SerialDenseMatrix<int,ScalarType> & source, int numRows, int numCols, int startRow=0, int startCol=0)
    { return Teuchos::rcp(new const Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View, source, numRows, numCols, startRow, startCol)); }    

    //! \brief Returns a deep copy of the requested subview.
    static Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType>> SubviewCopy( const Teuchos::SerialDenseMatrix<int,ScalarType>& source, int numRows, int numCols, int startRow=0, int startCol=0)
    { return Teuchos::rcp(new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::Copy, source, numRows, numCols, startRow, startCol)); }    
    //@}

    //@{ \name Attribute methods

    //! \brief Obtain the number of rows of \c dm.
    static int GetNumRows( const Teuchos::SerialDenseMatrix<int,ScalarType>& dm )
    { return dm.numRows(); }     

    //! \brief Obtain the number of columns of \c dm.
    static int GetNumCols( const Teuchos::SerialDenseMatrix<int,ScalarType>& dm )
    { return dm.numCols(); }     

    //! \brief Obtain the stride between the columns of \c dm.
    static int GetStride( const Teuchos::SerialDenseMatrix<int,ScalarType>& dm )
    { return dm.stride(); }    

    //@}

    //@{ \name Shaping methods

    /* \brief Reshaping method for changing the size of \c dm to have \c numrows rows and \c numcols columns.
     *        All values will be initialized to zero if the final argument is true. 
     *        If the final argument is fale, the previous entries in 
     *        the matrix will be maintained. For new entries that did not exist in the previous matrix, values will
     *        contain noise from memory. 
    */
    static void Reshape( Teuchos::SerialDenseMatrix<int,ScalarType>& dm, const int numrows, const int numcols, bool initZero = false) { 
      if(initZero){
        int err =  dm.shape(numrows,numcols); 
        if(err != 0){throw std::runtime_error ("Error in DenseMatrixTraits::shape. Teuchos::SerialDenseMatrix.shape failed.");}
      }
      else{
        int err =  dm.reshape(numrows,numcols); 
        if(err != 0){throw std::runtime_error ("Error in DenseMatrixTraits::reshape. Teuchos::SerialDenseMatrix.reshape failed.");}
      }
    }     

    //@}

    //@{ \name Data access methods

    //! \brief Access a reference to the (i,j) entry of \c dm, \c e_i^T dm e_j.
    static ScalarType & Value( Teuchos::SerialDenseMatrix<int,ScalarType>& dm, const int i, const int j )
    { 
      return dm(i,j);
    }

    //! \brief Access a const reference to the (i,j) entry of \c dm, \c e_i^T dm e_j.
    static const ScalarType & ValueConst( const Teuchos::SerialDenseMatrix<int,ScalarType>& dm, const int i, const int j ) 
    { 
      return dm(i,j);
    }

    static void SyncDeviceToHost(Teuchos::SerialDenseMatrix<int,ScalarType> &){ }

    static void SyncHostToDevice(Teuchos::SerialDenseMatrix<int,ScalarType> &){ }
    //@}

    //@{ \name Operator methods
    
    //!  \brief Adds sourceDM to thisDM and returns answer in thisDM.
    static void Add( Teuchos::SerialDenseMatrix<int,ScalarType>& thisDM, const Teuchos::SerialDenseMatrix<int,ScalarType>& sourceDM){ 
      thisDM += sourceDM; 
    }

    //!  \brief Fill all entries with \c value. Value is zero if not specified.
    static void PutScalar( Teuchos::SerialDenseMatrix<int,ScalarType> & dm, ScalarType value = Teuchos::ScalarTraits<ScalarType>::zero()){ 
      dm.putScalar(value);
    }

    //!  \brief Multiply all entries by a scalar. DM = value.*DM
    static void Scale( Teuchos::SerialDenseMatrix<int,ScalarType>& dm, ScalarType value){
      dm.scale(value);
    }

    //!  \brief Fill the DM with random entries.
    //!   Entries are assumed to be the same on each MPI rank (each matrix copy). 
    //TODO What to do here? Kinda needs random synced version??
    static void Randomize( Teuchos::SerialDenseMatrix<int,ScalarType>& dm){ 
      dm.random();
    }

    //!  \brief Copies entries of source to dest (deep copy). 
    static void Assign( Teuchos::SerialDenseMatrix<int,ScalarType>& dest, const Teuchos::SerialDenseMatrix<int,ScalarType>& source){ 
      dest.assign(source);
    }

    //!  \brief Returns the Frobenius norm of the dense matrix.
    static typename Teuchos::ScalarTraits<ScalarType>::magnitudeType NormFrobenius(Teuchos::SerialDenseMatrix<int,ScalarType>& dm) { 
      return dm.normFrobenius(); 
    }

    //!  \brief Returns the one-norm of the dense matrix.
    static typename Teuchos::ScalarTraits<ScalarType>::magnitudeType NormOne( Teuchos::SerialDenseMatrix<int,ScalarType>& dm) { 
      return dm.normOne();
    }
    //@}
  };
  
} // namespace Belos

#endif // end file BELOS_TEUCHOS_DENSE_MAT_TRAITS_HPP
