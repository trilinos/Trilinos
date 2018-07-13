// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
#ifndef ANASAZI_DENSE_MAT_TRAITS_HPP
#define ANASAZI_DENSE_MAT_TRAITS_HPP

/*! \file AnasaziDenseMatrixTraits.hpp
  \brief Virtual base class which defines basic traits for the dense matrix type
*/

#include "AnasaziConfigDefs.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Anasazi {

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

    \return Reference-counted pointer to a new dense matrix of type \c DM.
    */
    static Teuchos::RCP<DM> Clone( const DM& dm, const int numrows, const int numcols )
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return Teuchos::null; }     

    /*! \brief Creates a new \c DM and copies the contents of \c dm into the new matrix (deep copy).
      
      \return Reference-counted pointer to the new matrix of type \c DM.
    */
    static Teuchos::RCP<DM> CloneCopy( const DM& dm )
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return Teuchos::null; }     

    /*! \brief Creates a new \c DM and copies the selected contents of \c dm into the new matrix (deep copy).
      
      \return Reference-counted pointer to the new matrix of type \c DM.
    */
    static Teuchos::RCP<DM> CloneCopy( const DM& dm, const int numrows, const int numcols, 
         const int firstrow, const int firstcol )
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return Teuchos::null; }     

    /*! \brief Creates a new \c DM that shares the selected contents of \c dm (shallow copy).

    \return Reference-counted pointer to the new matrix of type \c DM.
    */      
    static Teuchos::RCP<DM> CloneView( DM& dm, const int numrows, const int numcols, 
         const int firstrow, const int firstcol )
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); return Teuchos::null; }     

    /*! \brief Creates a new \c DM that shares the selected contents of \c dm (shallow copy).

    \return Reference-counted pointer to the new matrix of type \c DM.
    */      
    static Teuchos::RCP<const DM> CloneView( const DM& dm, const int numrows, const int numcols, 
         const int firstrow, const int firstcol )
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

    //@{ \name Data access methods

    //! \brief Access a reference to the (i,j) entry of \c dm, \c e_i^T dm e_j.
    static ScalarType & value( DM& dm, const int i, const int j )
    { 
      UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); 
      return Teuchos::ScalarTraits<ScalarType>::zero();
    }

    //! \brief Access a const reference to the (i,j) entry of \c dm, \c e_i^T dm e_j.
    static const ScalarType & value( const DM& dm, const int i, const int j )
    { 
      UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); 
      return Teuchos::ScalarTraits<ScalarType>::zero();
    }

    /*! \brief Access the pointers to the data in \c dm, information enabling
     * C-style access to the data.

        \return The return value is a pointer to the data, stored sequentially.
        \c *cor denotes whether the data is stored column-oriented (\f$*cor ==
        true$) or row-oriented (\f$cor == false$). \c *stride denotes the stride
        between columns/rows. 

    */
    static ScalarType * values( DM& dm, int *stride, bool *cor) {
      UndefinedDenseMatTraits<ScalarType, DM>::notDefined();
      return (ScalarType *)NULL;
    }

    //@}

    //@{ \name Update methods

    /*! \brief Update \c dm with \f$ \alpha op(A) op(B) + \beta dm \f$,
        where 
          op(A) := A^H if transA == true
          op(A) := A   if transA == false

          and 

          op(B) := B^H if transB == true
          op(B) := B   if transB == false
     */
    static void MvTimesMatAddMv( bool transA, bool transB,
         const ScalarType alpha, const DM& A, const DM& B, 
         const ScalarType beta, DM& dm )
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); }     

    /*! \brief Replace \c dm with \f$\alpha A + \beta B\f$.
     */
    static void MvAddMv( const ScalarType alpha, const DM& A, const ScalarType beta, const DM& B, DM& dm )
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); }     

    //@}

    //@{ \name Initialization methods

    /*! \brief Replace the entries of \c dm with random numbers.
     */
    static void DMRandom( DM& dm )
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); }     

    /*! \brief Replace each element of \c mv with \c alpha.
     */
    static void DMInit( DM& dm, const ScalarType alpha = Teuchos::ScalarTraits<ScalarType>::zero() )
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); }     

    //@}

    //@{ \name Print method

    /*! \brief Print the matrix \c dm to the output stream \c os.
     */
    static void DMPrint( const DM& dm, std::ostream& os )
    { UndefinedDenseMatTraits<ScalarType, DM>::notDefined(); }     

    //@}
  };
  
} // namespace Anasazi

#endif // end file ANASAZI_DENSE_MAT_TRAITS_HPP
