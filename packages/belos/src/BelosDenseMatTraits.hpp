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
#ifndef BELOS_DENSE_MAT_TRAITS_HPP
#define BELOS_DENSE_MAT_TRAITS_HPP

/*! \file BelosDenseMatrixTraits.hpp
  \brief Virtual base class which defines basic traits for the dense matrix type
*/

//#include "BelosConfigDefs.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "BelosDenseSolver.hpp"

namespace Belos {

/*! \struct UndefinedDenseMatTraits
   \brief This is the default struct used by DenseMatrixTraits<OrdinalType, ScalarType> class to 
   produce a compile time error when the specialization does not exist for
   dense matrix type <tt>DM</tt>.
*/

  template< class ScalarType, class DM >
  struct UndefinedDenseMatTraits
  {
    //! This function should not compile if there is an attempt to instantiate!
    /*! \note Any attempt to compile this function results in a compile time error.  This means
      that the template specialization of Anasazi::DenseMatTraits class for type <tt>DM</tt> does 
      not exist, or is not complete.
    */
    static inline ScalarType notDefined() { return DM::this_type_is_missing_a_specialization(); };
  };
  
  /*! \class DenseMatTraits
    \brief Virtual base class which defines basic traits for the multi-vector type.
    
    An adapter for this traits class must exist for the <tt>DM</tt> type.
    If not, this class will produce a compile-time error.
  */

  template<class ScalarType, class DM>
  class DenseMatTraits 
  {
  public:
    
    //@{ \name Creation methods

    /*! \brief Creates a new empty \c DM with no dimension.

    \return Reference-counted pointer to a new dense matrix of type \c DM.
    */
    static Teuchos::RCP<DM> Create()
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return Teuchos::null; }     

    /*! \brief Creates a new empty \c DM containing \c numvecs columns.
     *         Will be initialized to zeros if last parameter is true.

    \return Reference-counted pointer to a new dense matrix of type \c DM.
    */
    static Teuchos::RCP<DM> Create( const int numrows, const int numcols, bool initZero = true)
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return Teuchos::null; }     

    /*! \brief Create a new copy \c DM, possibly transposed.
  
    \return Reference-counted pointer to a new dense matrix of type \c DM.
    */
    static Teuchos::RCP<DM> CreateCopy(const DM & dm, bool transpose=false)
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return Teuchos::null; }

/* Kokkos Ex View-from-ptr constructor: 
View(const pointer_type &ptr, const IntType&... indices)

    Unmanaged data wrapping constructor.

        ptr: pointer to a user provided memory allocation. Must provide storage of size View::required_allocation_size(n0,...,nR)
        indices: Extents of the View.
        Requires: sizeof(IntType...)==rank_dynamic() or sizeof(IntType...)==rank(). In the latter case, the extents corresponding to compile-time dimensions must match the View type’s compile-time extents.
        Requires: array_layout::is_regular == true.
*/

    //! \brief Returns a raw pointer to the (non-const) data on the host.
    /// \note We assume that return data in in a column-major format
    /// because this is what is expected by LAPACK.
    /// \note This raw pointer is intended only for passing data to LAPACK
    /// functions. Other operations on the raw data may result in undefined behavior!
    static ScalarType* GetRawHostPtr(DM & dm )
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return Teuchos::null; }     

    //! \brief Returns a raw pointer to const data on the host.
    static ScalarType const * GetConstRawHostPtr(const DM & dm )  
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return Teuchos::null; }     

    //! \brief Marks host data modified to avoid device sync errors. 
    /// \note Belos developers must call this function after EVERY
    ///   call to LAPACK that modifies dense matrix data accessed via raw pointer. 
    //static void RawPtrDataModified(DM & dm)
    //{ UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); }

    //! \brief Returns an RCP to a DM which has a subview of the given DM.
    //        Row and column indexing is zero-based.
    //        Source  - Reference to another dense matrix from which values are to be copied.
    //        numRows - The number of rows in this matrix.
    //        numCols - The number of columns in this matrix.
    //        startRow  - The row of Source from which the submatrix copy should start.
    //        startCol  - The column of Source from which the submatrix copy should start.
    //
    //        Should ints be const? Should they be ints or some other ordinal type?
    static Teuchos::RCP<DM> Subview( DM & source, int numRows, int numCols, int startRow=0, int startCol=0)
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return Teuchos::null; }     

    static Teuchos::RCP<const DM> SubviewConst( const DM & source, int numRows, int numCols, int startRow=0, int startCol=0)
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return Teuchos::null; }     

    //! \brief Returns a deep copy of the requested subview.
    static Teuchos::RCP<DM> SubviewCopy( const DM & source, int numRows, int numCols, int startRow=0, int startCol=0)
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return Teuchos::null; }     
    //@}

    //@{ \name Attribute methods

    //! \brief Obtain the number of rows of \c dm.
    static int GetNumRows( const DM& dm )
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return 0; }     

    //! \brief Obtain the number of columns of \c dm.
    static int GetNumCols( const DM& dm )
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return 0; }     

    //! \brief Obtain the stride between the columns of \c dm.
    static int GetStride( const DM& dm )
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return 0; }     

    //@}

    //@{ \name Shaping methods

    /* \brief Reshaping method for changing the size of \c dm,
    *         keeping the entries. 
    */
    
    /* \brief Reshaping method for changing the size of \c dm to have \c numrows rows and \c numcols columns.
     *        All values will be initialized to zero if the final argument is true. 
     *        If the final argument is fale, the previous entries in 
     *        the matrix will be maintained. For new entries that did not exist in the previous matrix, values will
     *        contain noise from memory. 
    */
    static void Reshape( DM& dm, const int numrows, const int numcols, bool initZero = true)
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); }     

    //@}

    //@{ \name Data access methods

    //! \brief Access a reference to the (i,j) entry of \c dm, \c e_i^T dm e_j.
    static ScalarType & Value( DM& dm, const int i, const int j )
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); }

    //! \brief Access a const reference to the (i,j) entry of \c dm, \c e_i^T dm e_j.
    static const ScalarType & ValueConst( const DM& dm, const int i, const int j ) 
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); }
    
    static void SyncDeviceToHost(DM & dm)
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); }

    static void SyncHostToDevice(DM & dm)
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); }

    //@}
    //@{ \name Operator methods
    
    //!  \brief Adds sourceDM to thisDM and returns answer in thisDM.
    static void Add( DM& thisDM, const DM& sourceDM)
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); }

    //!  \brief Fill all entries with \c value. Value is zero if not specified.
    static void PutScalar( DM& dm, ScalarType value = Teuchos::ScalarTraits<ScalarType>::zero()) 
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); }

    //!  \brief Multiply all entries by a scalar. DM = value.*DM
    static void Scale( DM& dm, ScalarType value)
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); }

    //!  \brief Fill the DM with random entries.
    //!   Entries are assumed to be the same on each MPI rank (each matrix copy). 
    static void Randomize( DM& dm)
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); }

    //!  \brief Copies entries of sourceDM to thisDM (deep copy). 
    static void Assign( DM& thisDM, const DM& sourceDM)
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); }

    //!  \brief Returns the Frobenius norm of the dense matrix.
    static typename Teuchos::ScalarTraits<ScalarType>::magnitudeType NormFrobenius( DM& dm)
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); }

    //!  \brief Returns the one-norm of the dense matrix.
    static typename Teuchos::ScalarTraits<ScalarType>::magnitudeType NormOne( DM& dm)
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); }

    //@}
    //@{ \name Solver methods

    //!  \brief Returns a dense solver object for the dense matrix.
    static Teuchos::RCP<DenseSolver<ScalarType, DM>> createDenseSolver()
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return Teuchos::null; }
    //@}
 
  };

} // namespace Belos

#endif // end file BELOS_DENSE_MAT_TRAITS_HPP
