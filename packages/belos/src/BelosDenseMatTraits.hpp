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

    /*! \brief Creates a new empty \c DM containing \c numvecs columns.
     *         Will be initialized to zeros if last parameter is true.

    \return Reference-counted pointer to a new dense matrix of type \c DM.
    */
    static Teuchos::RCP<DM> Create( const int numrows, const int numcols, bool initZero = true)
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return Teuchos::null; }     

/* Kokkos Ex View-from-ptr constructor: 
View(const pointer_type &ptr, const IntType&... indices)

    Unmanaged data wrapping constructor.

        ptr: pointer to a user provided memory allocation. Must provide storage of size View::required_allocation_size(n0,...,nR)
        indices: Extents of the View.
        Requires: sizeof(IntType...)==rank_dynamic() or sizeof(IntType...)==rank(). In the latter case, the extents corresponding to compile-time dimensions must match the View type’s compile-time extents.
        Requires: array_layout::is_regular == true.
*/
//TODO: What conditions should we have on this function? When does it work and not work?
// Do we need details about the stride and/or layout? 
// What corresponds to the Teuchos Copy or View parameter with Teuchos::SerialDenseMatrix? 
    static Teuchos::RCP<DM> ViewFromPtr( ScalarType* &ptr, const int numrows, const int numcols)
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return Teuchos::null; }     

    //! \brief Returns a raw pointer to the data on the host.
    // ....Expectations for data layout??
    //TODO: Will we need this permanently? Yes, for LAPACK....
    static ScalarType* GetRawHostPtr(const DM & dm )
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return Teuchos::null; }     

    //! \brief Returns an RCP to a DM which has a subview of the given DM.
    //        Row and column indexing is zero-based.
    //        Row and column ranges include the first index and exclude the second index.
    //        So, for example, giving std::make_pair(5,7) for the rowRange will return rows 
    //        in range [5,7), so rows 5 and 6. 
    //        This follows the convention of Kokkos::subview.
    static Teuchos::RCP<DM> Subview( const DM & dm, std::pair<int,int> rowRange, std::pair<int,int> colRange)
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return Teuchos::null; }     
    //@}

    //@{ \name Attribute methods

    //! \brief Obtain the number of rows of \c dm.
    static int GetNumRows( const DM& dm )
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return 0; }     

    //! \brief Obtain the number of columns of \c dm.
    static int GetNumCols( const DM& dm )
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return 0; }     

    //@}

    //@{ \name Shaping methods

    /* \brief Reshaping method for changing the size of \c dm,
    *         keeping the entries. 
    \return Integer error code, set to 0 if successful.
    */
    
    // TODO Specify what should happen if matrix shrinks? Should new entries init to zero?
    // No- new entries just get whatever is in memory. 
    // Are there Teuchos unit tests we can copy to verify this function?
    static int Reshape( DM& dm, const int numrows, const int numcols)
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return 0; }     

    //@}

    //@{ \name Data access methods

    //! \brief Access a reference to the (i,j) entry of \c dm, \c e_i^T dm e_j.
    static ScalarType & Value( DM& dm, const int i, const int j )
    { 
      UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); 
      ScalarType * ptrz = new ScalarType;
      return ptrz;
      //return Teuchos::ScalarTraits<ScalarType>::zero();
    }

    //! \brief Access a const reference to the (i,j) entry of \c dm, \c e_i^T dm e_j.
    static const ScalarType & Value( const DM& dm, const int i, const int j )
    { 
      UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); 
      ScalarType * ptrz = new ScalarType;
      return ptrz;
      //return Teuchos::ScalarTraits<ScalarType>::zero();
    }

    static void SyncHostToDevice(DM & dm)
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); }

    static void SyncDeviceToHost(DM & dm)
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); }
    //@}
    //@{ \name Operator methods
    
    //!  \brief Adds sourceDM to thisDM and returns answer in thisDM.
    static void Add( DM& thisDM, const DM& sourceDM)
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
    //@}
  };
  
} // namespace Belos

#endif // end file BELOS_DENSE_MAT_TRAITS_HPP
