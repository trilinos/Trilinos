// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
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
// ***********************************************************************
// @HEADER

/*! \file AnasaziSpecializedEpetraAdapter.hpp
  \brief Declarations of specialized Anasazi multi-vector and operator classes using Epetra_MultiVector and Epetra_Operator
*/

#ifndef ANASAZI_SPECIALIZED_EPETRA_ADAPTER_HPP
#define ANASAZI_SPECIALIZED_EPETRA_ADAPTER_HPP

#include "AnasaziConfigDefs.hpp"
#include "Anasaziepetra_DLLExportMacro.h"
#include "AnasaziTypes.hpp"
#include "AnasaziMultiVec.hpp"
#include "AnasaziOperator.hpp"
#include "AnasaziEpetraAdapter.hpp"

#include "Teuchos_Assert.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"

#if defined(HAVE_ANASAZI_TPETRA) && defined(HAVE_ANASAZI_TSQR)
#  include <Tpetra_ConfigDefs.hpp> // HAVE_TPETRA_EPETRA 
#  if defined(HAVE_TPETRA_EPETRA)
#    include <Epetra_TsqrAdaptor.hpp>
#  endif // defined(HAVE_TPETRA_EPETRA)
#endif // defined(HAVE_ANASAZI_TPETRA) && defined(HAVE_ANASAZI_TSQR)

namespace Anasazi {

  //! @name Epetra Adapter Exceptions
  //@{

  /** \brief EpetraSpecializedMultiVecFailure is thrown when a return value from an Epetra
   * call on an Epetra_MultiVector is non-zero.
   */
  class EpetraSpecializedMultiVecFailure : public AnasaziError {public:
    EpetraSpecializedMultiVecFailure(const std::string& what_arg) : AnasaziError(what_arg)
    {}};

  //@}
  
  ///////////////////////////////////////////////////////////////
  //
  //--------template class AnasaziEpetraOpMultiVec-----------------
  //
  ///////////////////////////////////////////////////////////////
  
  /*! 
    \brief Specialized adapter class for Anasazi::MultiVec that uses Epetra_MultiVector and Epetra_Operator to define the inner-product

    \note The Epetra package performs double-precision arithmetic, so the use of Epetra with Anasazi will
    only provide a double-precision eigensolver.
  */
  class ANASAZIEPETRA_LIB_DLL_EXPORT EpetraOpMultiVec : public MultiVec<double>, public EpetraMultiVecAccessor {
  public:
    //! @name Constructors/Destructors
    //@{ 

    //! Basic EpetraOpMultiVec constructor.
    /*! @param Op [in] A reference-counted pointer to an existing fully constructed Epetra_Operator.
      @param Map [in] An Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
      @param numvecs [in] Number of vectors in multi-vector.

      \returns Pointer to an EpetraOpMultiVec
    */
    EpetraOpMultiVec(const Teuchos::RCP<Epetra_Operator> &Op, const Epetra_BlockMap& Map_in, const int numvecs);

    //! Create multi-vector with values from two dimensional array.
    /*! @param Op [in] A reference-counted pointer to an existing fully constructed Epetra_Operator.
      @param Map [in] An Epetra_LocalMap, Epetra_Map or Epetra_BlockMap
      @param array [in] Pointer to an array of double precision numbers.  The first vector starts at \c array, the
      second at \c array+stride, and so on.  This array is copied.
      @param numvecs [in] Number of vectors in the multi-vector.
      @param stride [in] The stride between vectors in memory of \c array.

      \returns Pointer to an EpetraOpMultiVec
    */
    EpetraOpMultiVec(const Teuchos::RCP<Epetra_Operator> &Op, const Epetra_BlockMap& Map_in, double * array, const int numvecs, const int stride=0);

    //! Create multi-vector from list of vectors in an existing EpetraOpMultiVec.
    /*! @param Op [in] A reference-counted pointer to an existing fully constructed Epetra_Operator.
      @param P_vec [in] A reference-counted pointer to an existing fully constructed Epetra_MultiVector.
      @param index [in] A integer vector containing the indices of the vectors to copy out of \c P_vec.

      \returns Pointer to an EpetraOpMultiVec
    */
    EpetraOpMultiVec(const Teuchos::RCP<Epetra_Operator> &Op, Epetra_DataAccess CV, const Epetra_MultiVector& P_vec, const std::vector<int>& index);

    //! Copy constructor. 
    EpetraOpMultiVec(const EpetraOpMultiVec& P_vec);

    //! Destructor
    virtual ~EpetraOpMultiVec() {};

    //@}

    //! @name Creation methods
    //@{ 

    /*! \brief Creates a new empty EpetraOpMultiVec containing \c numvecs columns.
      
    \returns Pointer to an EpetraOpMultiVec
    */
    MultiVec<double> * Clone ( const int numvecs ) const;

    /*! \brief Creates a new EpetraOpMultiVec and copies contents of \c *this into
      the new vector (deep copy).
      
      \returns Pointer to an EpetraOpMultiVec
    */
    MultiVec<double> * CloneCopy () const;

    /*! \brief Creates a new EpetraOpMultiVec and copies the selected contents of \c *this 
      into the new vector (deep copy).  
      
      The copied vectors from \c *this are indicated by the \c index.size() indices in \c index.
      
      \returns Pointer to an EpetraOpMultiVec
    */
    MultiVec<double> * CloneCopy ( const std::vector<int>& index ) const;
    
    /*! \brief Creates a new EpetraOpMultiVec that shares the selected contents of \c *this.
      
    The index of the \c numvecs vectors shallow copied from \c *this are indicated by the
    indices given in \c index.
    
    \returns Pointer to an EpetraOpMultiVec
    */
    MultiVec<double> * CloneViewNonConst ( const std::vector<int>& index );

    /*! \brief Creates a new EpetraOpMultiVec that shares the selected contents of \c *this.
      
    The index of the \c numvecs vectors shallow copied from \c *this are indicated by the
    indices given in \c index.
    
    \returns Pointer to an EpetraOpMultiVec
    */
    const MultiVec<double> * CloneView ( const std::vector<int>& index ) const;

    //@}

    //! @name Accessor methods
    Teuchos::RCP<Epetra_MultiVector> GetEpetraMultiVector() { return Epetra_MV; }

    //! @name Attribute methods
    //@{ 

    //! Obtain the number of vectors in *this.
    int GetVecLength () const { return Epetra_MV->GlobalLength(); }

    //! The number of rows in the multivector.
    //! \note This method supersedes GetVecLength, which will be deprecated.
    ptrdiff_t GetGlobalLength () const 
    {
       if ( Epetra_MV->Map().GlobalIndicesLongLong() )
          return static_cast<ptrdiff_t>( Epetra_MV->GlobalLength64() );       
       else
          return static_cast<ptrdiff_t>( Epetra_MV->GlobalLength() ); 
    }

    //! Obtain the vector length of *this.
    int GetNumberVecs () const { return Epetra_MV->NumVectors(); }

    //@}

    //! @name Update methods
    //@{ 
    /*! \brief Update \c *this with \f$\alpha AB + \beta (*this)\f$.
     */
    void MvTimesMatAddMv ( double alpha, const MultiVec<double>& A, 
                           const Teuchos::SerialDenseMatrix<int,double>& B, 
                           double beta );

    /*! \brief Replace \c *this with \f$\alpha A + \beta B\f$.
     */
    void MvAddMv ( double alpha, const MultiVec<double>& A, 
                   double beta, const MultiVec<double>& B);

    /*! \brief Compute a dense matrix \c B through the matrix-matrix multiply \f$\alpha A^T(*this)\f$.
    */
    void MvTransMv ( double alpha, const MultiVec<double>& A, Teuchos::SerialDenseMatrix<int,double>& B 
#ifdef HAVE_ANASAZI_EXPERIMENTAL
        , ConjType conj = Anasazi::CONJ
#endif
        ) const;
  
    /*! \brief Compute a vector \c b where the components are the individual dot-products, i.e. \f$ b[i] = A[i]^H(this[i])\f$ where \c A[i] is the i-th column of \c A.
    */
    void MvDot ( const MultiVec<double>& A, std::vector<double> &b
#ifdef HAVE_ANASAZI_EXPERIMENTAL
        , ConjType conj = Anasazi::CONJ
#endif
        ) const;

    /*! \brief Scale each element of the vectors in \c *this with \c alpha.
     */
    void MvScale ( double alpha ) { 
      TEUCHOS_TEST_FOR_EXCEPTION( Epetra_MV->Scale( alpha )!=0, EpetraSpecializedMultiVecFailure,
          "Anasazi::EpetraOpMultiVec::MvScale call to Epetra_MultiVector::Scale() returned a nonzero value.");
    }
    
    /*! \brief Scale each element of the \c i-th vector in \c *this with \c alpha[i].
     */
    void MvScale ( const std::vector<double>& alpha );

    //@}
    //! @name Norm method
    //@{ 
    
    /*! \brief Compute the 2-norm of each individual vector of \c *this.  
      Upon return, \c normvec[i] holds the 2-norm of the \c i-th vector of \c *this
    */
    void MvNorm ( std::vector<double> & normvec ) const;
    
    //@}
    
    //! @name Initialization methods
    //@{ 
    /*! \brief Copy the vectors in \c A to a set of vectors in \c *this.  
      
    The \c numvecs vectors in \c A are copied to a subset of vectors in \c *this
    indicated by the indices given in \c index.
    */
    void SetBlock ( const MultiVec<double>& A, const std::vector<int>& index );

    /*! \brief Fill the vectors in \c *this with random numbers.
     */
    void MvRandom() { 
      TEUCHOS_TEST_FOR_EXCEPTION( Epetra_MV->Random()!=0, EpetraSpecializedMultiVecFailure,
          "Anasazi::EpetraOpMultiVec::MvRandom call to Epetra_MultiVector::Random() returned a nonzero value.");
    }

    /*! \brief Replace each element of the vectors in \c *this with \c alpha.
     */
    void MvInit ( double alpha ) { 
      TEUCHOS_TEST_FOR_EXCEPTION( Epetra_MV->PutScalar( alpha )!=0, EpetraSpecializedMultiVecFailure,
          "Anasazi::EpetraOpMultiVec::MvInit call to Epetra_MultiVector::PutScalar() returned a nonzero value.");
    } 
   
    //! @name Accessor methods (inherited from EpetraMultiVecAccessor)
    //@{

    /*! \brief Return the pointer to the Epetra_MultiVector object. */
    Epetra_MultiVector* GetEpetraMultiVec() { return &*Epetra_MV; };

    /*! \brief Return the pointer to the Epetra_MultiVector object. */
    const Epetra_MultiVector* GetEpetraMultiVec() const { return &*Epetra_MV; };

    //@}
 
    //@}
    //! @name Print method
    //@{ 
    /*! \brief Print \c *this EpetraOpMultiVec.
     */
    void MvPrint( std::ostream& os ) const { Epetra_MV->Print( os ); }

    //@}

  private:
//use pragmas to disable some false-positive warnings for windows 
// sharedlibs export
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4251)
#endif
    Teuchos::RCP<Epetra_Operator> Epetra_OP;
    Teuchos::RCP<Epetra_MultiVector> Epetra_MV;
    Teuchos::RCP<Epetra_MultiVector> Epetra_MV_Temp;
#ifdef _MSC_VER
#pragma warning(pop)
#endif
  };

 
} // end of Anasazi namespace 

#endif // end of file ANASAZI_SPECIALIZED_EPETRA_ADAPTER_HPP
