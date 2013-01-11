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

/*! \file AnasaziEpetraAdapter.hpp
  \brief Declarations of Anasazi multi-vector and operator classes using Epetra_MultiVector and Epetra_Operator classes
*/

#ifndef ANASAZI_EPETRA_ADAPTER_HPP
#define ANASAZI_EPETRA_ADAPTER_HPP

#include "AnasaziConfigDefs.hpp"
#include "Anasaziepetra_DLLExportMacro.h"
#include "AnasaziTypes.hpp"
#include "AnasaziMultiVec.hpp"
#include "AnasaziOperator.hpp"

#include "Teuchos_Assert.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
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

  /** \brief EpetraMultiVecFailure is thrown when a return value from an Epetra
   * call on an Epetra_MultiVector is non-zero.
   */
  class EpetraMultiVecFailure : public AnasaziError {public:
    EpetraMultiVecFailure(const std::string& what_arg) : AnasaziError(what_arg)
    {}};

  /** \brief EpetraOpFailure is thrown when a return value from an Epetra
   * call on an Epetra_Operator is non-zero.
   */
  class EpetraOpFailure : public AnasaziError {public:
    EpetraOpFailure(const std::string& what_arg) : AnasaziError(what_arg)
    {}};

  //@}

  //! @name Epetra_MultiVector Accessor Interface
  //@{

  /** \brief EpetraMultiVecAccessor is an interfaceto allow any Anasazi::MultiVec implementation
   * that is based on Epetra_MultiVector to use the various Anasazi::Operator interfaces defined for Epetra_Operator.
   */
  class EpetraMultiVecAccessor {
 
  public:
    /*! \brief Return the pointer to the Epetra_MultiVector object. */
    virtual Epetra_MultiVector* GetEpetraMultiVec() { return 0; }

    /*! \brief Return the pointer to the Epetra_MultiVector object. */
    virtual const Epetra_MultiVector* GetEpetraMultiVec() const { return 0; }
  };

  //@}
  
  ///////////////////////////////////////////////////////////////
  //
  //--------template class AnasaziEpetraMultiVec-----------------
  //
  ///////////////////////////////////////////////////////////////
  
  /*! 
    \brief Basic adapter class for Anasazi::MultiVec that uses Epetra_MultiVector.

    \note The Epetra package performs double-precision arithmetic, so the use of Epetra with Anasazi will
    only provide a double-precision eigensolver.
  */
  class ANASAZIEPETRA_LIB_DLL_EXPORT EpetraMultiVec : public MultiVec<double>,  public Epetra_MultiVector, public EpetraMultiVecAccessor {
  public:
    //! @name Constructors/Destructors
    //@{ 

    //! Basic EpetraMultiVec constructor.
    /*! @param Map [in] An Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
      @param numvecs [in] Number of vectors in multi-vector.

      \returns Pointer to an EpetraMultiVec
    */
    EpetraMultiVec(const Epetra_BlockMap& Map_in, const int numvecs);

    //! Copy constructor.
    EpetraMultiVec(const Epetra_MultiVector & P_vec);
    
    //! Create multi-vector with values from two dimensional array.
    /*! @param Map [in] An Epetra_LocalMap, Epetra_Map or Epetra_BlockMap
      @param array [in] Pointer to an array of double precision numbers.  The first vector starts at \c array, the
      second at \c array+stride, and so on.  This array is copied.
      @param numvecs [in] Number of vectors in the multi-vector.
      @param stride [in] The stride between vectors in memory of \c array.

      \returns Pointer to an EpetraMultiVec
    */
    EpetraMultiVec(const Epetra_BlockMap& Map_in, double * array, const int numvecs, const int stride=0);

    //! Create multi-vector from list of vectors in an existing EpetraMultiVec.
    /*! @param CV [in] Enumerated type set to Copy or View.
      @param P_vec [in] An existing fully constructed Epetra_MultiVector.
      @param index [in] A integer vector containing the indices of the vectors to copy out of \c P_vec.

      \returns Pointer to an EpetraMultiVec
    */
    EpetraMultiVec(Epetra_DataAccess CV, const Epetra_MultiVector& P_vec, const std::vector<int>& index);

    //! Destructor
    virtual ~EpetraMultiVec() {};

    //@}

    //! @name Creation methods
    //@{ 

    /*! \brief Creates a new empty EpetraMultiVec containing \c numvecs columns.
      
    \returns Pointer to an EpetraMultiVec
    */
    MultiVec<double> * Clone ( const int numvecs ) const;

    /*! \brief Creates a new EpetraMultiVec and copies contents of \c *this into
      the new vector (deep copy).
      
      \returns Pointer to an EpetraMultiVec
    */
    MultiVec<double> * CloneCopy () const;

    /*! \brief Creates a new EpetraMultiVec and copies the selected contents of \c *this 
      into the new vector (deep copy).  
      
      The copied vectors from \c *this are indicated by the \c index.size() indices in \c index.
      
      \returns Pointer to an EpetraMultiVec
    */
    MultiVec<double> * CloneCopy ( const std::vector<int>& index ) const;
    
    /*! \brief Creates a new EpetraMultiVec that shares the selected contents of \c *this.
      
    The index of the \c numvecs vectors shallow copied from \c *this are indicated by the
    indices given in \c index.
    
    \returns Pointer to an EpetraMultiVec
    */
    MultiVec<double> * CloneViewNonConst ( const std::vector<int>& index );

    /*! \brief Creates a new EpetraMultiVec that shares the selected contents of \c *this.
      
    The index of the \c numvecs vectors shallow copied from \c *this are indicated by the
    indices given in \c index.
    
    \returns Pointer to an EpetraMultiVec
    */
    const MultiVec<double> * CloneView ( const std::vector<int>& index ) const;

    //@}

    //! Obtain the number of vectors in *this.
    int GetVecLength () const { return GlobalLength(); }

    //! The number of rows in the multivector.
    //! \note This method supersedes GetVecLength, which will be deprecated.
    ptrdiff_t GetGlobalLength () const 
    {
       if ( Map().GlobalIndicesLongLong() )
          return static_cast<ptrdiff_t>( GlobalLength64() );
       else 
          return static_cast<ptrdiff_t>( GlobalLength() );
    }

    //! @name Attribute methods
    //! Obtain the vector length of *this.
    int GetNumberVecs () const { return NumVectors(); }

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
      TEUCHOS_TEST_FOR_EXCEPTION( this->Scale( alpha )!=0, EpetraMultiVecFailure,
          "Anasazi::EpetraMultiVec::MvScale call to Epetra_MultiVector::Scale() returned a nonzero value.");
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
    void MvNorm ( std::vector<double> & normvec ) const {
      if (((int)normvec.size() >= GetNumberVecs()) ) {
        TEUCHOS_TEST_FOR_EXCEPTION( this->Norm2(&normvec[0])!=0, EpetraMultiVecFailure,
            "Anasazi::EpetraMultiVec::MvNorm call to Epetra_MultiVector::Norm2() returned a nonzero value.");
      }
    };
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
      TEUCHOS_TEST_FOR_EXCEPTION( this->Random()!=0, EpetraMultiVecFailure,
          "Anasazi::EpetraMultiVec::MvRandom call to Epetra_MultiVector::Random() returned a nonzero value.");
    };

    /*! \brief Replace each element of the vectors in \c *this with \c alpha.
     */
    void MvInit ( double alpha ) { 
      TEUCHOS_TEST_FOR_EXCEPTION( this->PutScalar( alpha )!=0, EpetraMultiVecFailure,
          "Anasazi::EpetraMultiVec::MvInit call to Epetra_MultiVector::PutScalar() returned a nonzero value.");
    };
   
    //! @name Accessor methods (inherited from EpetraMultiVecAccessor)
    //@{
   
    /*! \brief Return the pointer to the Epetra_MultiVector object. */
    Epetra_MultiVector* GetEpetraMultiVec() { return this; };
 
    /*! \brief Return the pointer to the Epetra_MultiVector object. */
    const Epetra_MultiVector* GetEpetraMultiVec() const { return this; };

    //@}
 
    //@}
    //! @name Print method
    //@{ 
    /*! \brief Print \c *this EpetraMultiVec.
     */
    void MvPrint( std::ostream& os ) const { os << *this << std::endl; };
    //@}

  private:
  };
  //-------------------------------------------------------------
  
  ///////////////////////////////////////////////////////////////
  //
  //--------template class AnasaziEpetraOp---------------------
  //
  ///////////////////////////////////////////////////////////////
  
  /*! 
    \brief Basic adapter class for Anasazi::Operator that uses Epetra_Operator.

    \note The Epetra package performs double-precision arithmetic, so the use of Epetra with Anasazi will
    only provide a double-precision eigensolver.
  */
  class ANASAZIEPETRA_LIB_DLL_EXPORT EpetraOp : public virtual Operator<double> {
  public:
    //! @name Constructor/Destructor
    //@{ 
    
    //! Basic constructor.  Accepts reference-counted pointer to an Epetra_Operator.
    EpetraOp(const Teuchos::RCP<Epetra_Operator> &Op );
    
    //! Destructor
    ~EpetraOp();
    //@}
    
    //! @name Operator application method
    //@{ 
    
    /*! \brief This method takes the Anasazi::MultiVec \c X and
      applies the operator to it resulting in the Anasazi::MultiVec \c Y.
    */
    void Apply ( const MultiVec<double>& X, MultiVec<double>& Y ) const;
    //@} 
    
  private:
//use pragmas to disable some false-positive warnings for windows 
// sharedlibs export
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4251)
#endif
    Teuchos::RCP<Epetra_Operator> Epetra_Op;
#ifdef _MSC_VER
#pragma warning(pop)
#endif
  };
  //-------------------------------------------------------------

  ///////////////////////////////////////////////////////////////
  //
  //--------template class AnasaziEpetraGenOp--------------------
  //
  ///////////////////////////////////////////////////////////////
  
  /*! 
    \brief Adapter class for creating an operators often used in solving generalized eigenproblems.

    This class will apply the operation \f$A^{-1}M\f$ [default] or \f$AM\f$, for the \c Apply method of the
    Epetra_Operator / Anasazi::Operator.  The Anasazi::EpetraGenOp operator is useful when spectral 
    transformations are used within eigensolvers.  For instance, \f$A^{-1}M\f$ is a shift and invert 
    spectral transformation commonly used with Anasazi::BlockKrylovSchur to compute the smallest-magnitude
    eigenvalues for the eigenproblem \f$Ax = \lambda Mx\f$.

    \note The Epetra package performs double-precision arithmetic, so the use of Epetra with Anasazi will
    only provide a double-precision eigensolver.
  */

  class ANASAZIEPETRA_LIB_DLL_EXPORT EpetraGenOp : public virtual Operator<double>, public virtual Epetra_Operator {
  public:
    //! Basic constructor for applying operator \f$A^{-1}M\f$ [default] or \f$AM\f$.
    /*! If \c isAInverse is true this operator will apply \f$A^{-1}M\f$, else
      it will apply \f$AM\f$.
    */
    EpetraGenOp(const Teuchos::RCP<Epetra_Operator> &AOp, 
                const Teuchos::RCP<Epetra_Operator> &MOp,
                bool isAInverse = true );

    //! Destructor
    ~EpetraGenOp();
    
    //! Apply method [inherited from Anasazi::Operator class]
    /*! This method will apply \f$A^{-1}M\f$ or \f$AM\f$ to \c X, returning \c Y.
     */
    void Apply ( const MultiVec<double>& X, MultiVec<double>& Y ) const; 

    //! Apply method [inherited from Epetra_Operator class]
    /*! This method will apply \f$A^{-1}M\f$ or \f$AM\f$ to \c X, returning \c Y.
     */
    int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    //! Apply inverse method [inherited from Epetra_Operator class]
    /*! This method will apply \f$(A^{-1}M)^{-1}\f$ or \f$(AM)^{-1}\f$ to \c X, returning \c Y.
     */
    int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    //! Returns a character string describing the operator.
    const char* Label() const { return "Epetra_Operator applying A^{-1}M"; };
    
    //! Returns the current UseTranspose setting [always false for this operator].
    bool UseTranspose() const { return (false); };

    //! If set true, the transpose of this operator will be applied [not functional for this operator].
    int SetUseTranspose(bool /*UseTranspose_in*/) { return 0; };
    
    //! Returns true if this object can provide an approximate inf-norm [always false for this operator].
    bool HasNormInf() const { return (false); };
    
    //! Returns the infinity norm of the global matrix [not functional for this operator].
    double NormInf() const  { return (-1.0); };
    
    //! Returns the Epetra_Comm communicator associated with this operator.
    const Epetra_Comm& Comm() const { return Epetra_AOp->Comm(); };

    //! Returns the Epetra_Map object associated with the domain of this operator.
    const Epetra_Map& OperatorDomainMap() const { return Epetra_AOp->OperatorDomainMap(); };

    //! Returns the Epetra_Map object associated with the range of this operator.
    const Epetra_Map& OperatorRangeMap() const { return Epetra_AOp->OperatorRangeMap(); };

  private:
    bool isAInverse;

//use pragmas to disable some false-positive warnings for windows 
// sharedlibs export
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4251)
#endif
    Teuchos::RCP<Epetra_Operator> Epetra_AOp;
    Teuchos::RCP<Epetra_Operator> Epetra_MOp;
#ifdef _MSC_VER
#pragma warning(pop)
#endif
  };
  
  ///////////////////////////////////////////////////////////////
  //
  //--------template class AnasaziEpetraSymOp--------------------
  //
  ///////////////////////////////////////////////////////////////

  /*! 
    \brief Adapter class for creating a symmetric operator from an Epetra_Operator.

    This class will apply the operation \f$A^TA\f$ [default] or \f$AA^T\f$, for the \c Apply method of the
    Epetra_Operator / Anasazi::Operator.  The Anasazi::EpetraSymOp operator is useful when trying to compute
    a few singular values of the operator \f$A\f$.  The singular values are the square-root of the eigenvalues
    of \f$A^TA\f$ and \f$AA^T\f$.

    \note The Epetra package performs double-precision arithmetic, so the use of Epetra with Anasazi will
    only provide a double-precision eigensolver.
  */

  class ANASAZIEPETRA_LIB_DLL_EXPORT EpetraSymOp : public virtual Operator<double>, public virtual Epetra_Operator {
  public:
    //! Basic constructor for applying operator \f$A^TA\f$ [default] or \f$AA^T\f$.
    /*! If \c isTrans is false this operator will apply \f$A^TA\f$, else it will apply \f$AA^T\f$.
    */
    EpetraSymOp(const Teuchos::RCP<Epetra_Operator> &Op, bool isTrans = false );

    //! Destructor
    ~EpetraSymOp();
    
    //! Apply method [inherited from Anasazi::Operator class]
    /*! This method will apply \f$A^TA\f$ or \f$AA^T\f$ to \c X, returning \c Y.
     */
    void Apply ( const MultiVec<double>& X, MultiVec<double>& Y ) const; 

    //! Apply method [inherited from Epetra_Operator class]
    /*! This method will apply \f$A^TA\f$ or \f$AA^T\f$ to \c X, returning \c Y.
     */
    int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    //! Apply inverse method [inherited from Epetra_Operator class]
    /*! This method will apply \f$(A^TA)^{-1}\f$ or \f$(AA^T)^{-1}\f$ to \c X, returning \c Y.
      \note This method is only defined if \f$A^{-1}\f$ is defined for the given Epetra_Operator.
     */
    int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    //! Returns a character string describing the operator.
    const char* Label() const { return "Epetra_Operator applying A^TA or AA^T"; };
    
    //! Returns the current UseTranspose setting [always false for this operator].
    bool UseTranspose() const { return (false); };

    //! If set true, the transpose of this operator will be applied [not functional for this operator].
    int SetUseTranspose(bool /*UseTranspose_in*/) { return 0; };
    
    //! Returns true if this object can provide an approximate inf-norm [always false for this operator].
    bool HasNormInf() const { return (false); };
    
    //! Returns the infinity norm of the global matrix [not functional for this operator].
    double NormInf() const  { return (-1.0); };
    
    //! Returns the Epetra_Comm communicator associated with this operator.
    const Epetra_Comm& Comm() const { return Epetra_Op->Comm(); };

    //! Returns the Epetra_Map object associated with the domain of this operator.
    const Epetra_Map& OperatorDomainMap() const { return Epetra_Op->OperatorDomainMap(); };

    //! Returns the Epetra_Map object associated with the range of this operator.
    const Epetra_Map& OperatorRangeMap() const { return Epetra_Op->OperatorRangeMap(); };

  private:

//use pragmas to disable false-positive warnings in generating windows sharedlib exports
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4251)
#endif
    Teuchos::RCP<Epetra_Operator> Epetra_Op;
#ifdef _MSC_VER
#pragma warning(pop)
#endif

    bool isTrans_;
  };


  //////////////////////////////////////////////////////////////////
  //
  //--------template class AnasaziEpetraSymMVOp---------------------
  //
  //////////////////////////////////////////////////////////////////

  /*! 
    \brief Adapter class for creating a symmetric operator from an Epetra_MultiVector.

    This class will apply the operation \f$A^TA\f$ [default] or \f$AA^T\f$, for the \c Apply method of the
    Epetra_Operator / Anasazi::Operator.  The Anasazi::EpetraSymMvOp operator is useful when trying to compute
    a few singular values of the Epetra_MultiVector \f$A\f$.  The singular values are the square-root of the 
    eigenvalues of \f$A^TA\f$ and \f$AA^T\f$.

    \note The Epetra package performs double-precision arithmetic, so the use of Epetra with Anasazi will
    only provide a double-precision eigensolver.
  */

  class ANASAZIEPETRA_LIB_DLL_EXPORT EpetraSymMVOp : public virtual Operator<double> {
  public:
    //! Basic constructor for applying operator \f$A^TA\f$ [default] or \f$AA^T\f$.
    /*! If \c isTrans is false this operator will apply \f$A^TA\f$, else it will apply \f$AA^T\f$.
    */
    EpetraSymMVOp(const Teuchos::RCP<const Epetra_MultiVector> &MV, 
                  bool isTrans = false );
    
    //! Destructor
    ~EpetraSymMVOp() {};
    
    //! Apply method 
    /*! This method will apply \f$A^TA\f$ or \f$AA^T\f$ to \c X, returning \c Y.
     */
    void Apply ( const MultiVec<double>& X, MultiVec<double>& Y ) const; 

  private:

//use pragmas to disable some false-positive warnings for windows 
// sharedlibs export
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4251)
#endif
    Teuchos::RCP<const Epetra_MultiVector> Epetra_MV;
    Teuchos::RCP<const Epetra_Map> MV_localmap;
    Teuchos::RCP<const Epetra_BlockMap> MV_blockmap;
#ifdef _MSC_VER
#pragma warning(pop)
#endif

    bool isTrans_;
  };

  //////////////////////////////////////////////////////////////////
  //
  //--------template class AnasaziEpetraWSymMVOp---------------------
  //
  //////////////////////////////////////////////////////////////////

  /*! 
    \brief Adapter class for creating a weighted operator from an Epetra_MultiVector and Epetra_Operator.

    This class will apply the operation \f$A^T*W*A\f$ for the \c Apply method of the
    Anasazi::Operator.  The Anasazi::EpetraWSymMvOp operator is useful when trying to compute
    a few singular values of the Epetra_MultiVector \f$A\f$ under the weighting matrix \f$W\f$.  
    The singular values are the square-root of the eigenvalues of \f$A^T*W*A\f$.

    \note The Epetra package performs double-precision arithmetic, so the use of Epetra with Anasazi will
    only provide a double-precision eigensolver.
  */

  class ANASAZIEPETRA_LIB_DLL_EXPORT EpetraWSymMVOp : public virtual Operator<double> {
  public:
    //! Basic constructor for applying operator \f$A^T*W*A\f$.
    EpetraWSymMVOp(const Teuchos::RCP<const Epetra_MultiVector> &MV, 
                   const Teuchos::RCP<Epetra_Operator> &OP );
    
    //! Destructor
    ~EpetraWSymMVOp() {};
    
    //! Apply method 
    /*! This method will apply \f$(WA)^T*WA\f$ to \c X, returning \c Y.
     */
    void Apply ( const MultiVec<double>& X, MultiVec<double>& Y ) const; 

  private:
//use pragmas to disable some false-positive warnings for windows 
// sharedlibs export
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4251)
#endif
    Teuchos::RCP<const Epetra_MultiVector> Epetra_MV;
    Teuchos::RCP<Epetra_Operator> Epetra_OP;
    Teuchos::RCP<Epetra_MultiVector> Epetra_WMV;
    Teuchos::RCP<const Epetra_Map> MV_localmap;
    Teuchos::RCP<const Epetra_BlockMap> MV_blockmap;
#ifdef _MSC_VER
#pragma warning(pop)
#endif
  };

  //////////////////////////////////////////////////////////////////
  //
  //--------template class AnasaziEpetraW2SymMVOp---------------------
  //
  //////////////////////////////////////////////////////////////////

  /*! 
    \brief Adapter class for creating a weighted symmetric operator from an Epetra_MultiVector and Epetra_Operator.

    This class will apply the operation \f$(WA)^T*WA\f$ for the \c Apply method of the
    Anasazi::Operator.  The Anasazi::EpetraW2SymMvOp operator is useful when trying to compute
    a few singular values of the Epetra_MultiVector \f$A\f$ under the weighting matrix \f$W\f$.  
    The singular values are the square-root of the eigenvalues of \f$(WA)^T*WA\f$.

    \note The Epetra package performs double-precision arithmetic, so the use of Epetra with Anasazi will
    only provide a double-precision eigensolver.
  */

  class ANASAZIEPETRA_LIB_DLL_EXPORT EpetraW2SymMVOp : public virtual Operator<double> {
  public:
    //! Basic constructor for applying operator \f$A^T*W*A\f$.
    EpetraW2SymMVOp(const Teuchos::RCP<const Epetra_MultiVector> &MV, 
                   const Teuchos::RCP<Epetra_Operator> &OP );
    
    //! Destructor
    ~EpetraW2SymMVOp() {};
    
    //! Apply method 
    /*! This method will apply \f$(WA)^T*WA\f$ to \c X, returning \c Y.
     */
    void Apply ( const MultiVec<double>& X, MultiVec<double>& Y ) const; 

  private:
//use pragmas to disable some false-positive warnings for windows 
// sharedlibs export
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4251)
#endif
    Teuchos::RCP<const Epetra_MultiVector> Epetra_MV;
    Teuchos::RCP<Epetra_Operator> Epetra_OP;
    Teuchos::RCP<Epetra_MultiVector> Epetra_WMV;
    Teuchos::RCP<const Epetra_Map> MV_localmap;
    Teuchos::RCP<const Epetra_BlockMap> MV_blockmap;
#ifdef _MSC_VER
#pragma warning(pop)
#endif
  };

  
  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Anasazi::MultiVecTraits for Epetra::MultiVector.
  //
  ////////////////////////////////////////////////////////////////////

  /*! 
    \brief Template specialization of Anasazi::MultiVecTraits class using the Epetra_MultiVector class.

    This interface will ensure that any Epetra_MultiVector will be accepted by the Anasazi
    templated solvers.  

    \note The Epetra package performs double-precision arithmetic, so the use of Epetra with Anasazi will
    only provide a double-precision eigensolver.
  */

  template<>
  class MultiVecTraits<double, Epetra_MultiVector>
  {
  public:

    //! @name Creation methods
    //@{ 

    /*! \brief Creates a new empty Epetra_MultiVector containing \c numVecs columns.
      
    \return Reference-counted pointer to the new Epetra_MultiVector.
    */
    static Teuchos::RCP<Epetra_MultiVector> 
    Clone (const Epetra_MultiVector& mv, const int outNumVecs)
    { 
      TEUCHOS_TEST_FOR_EXCEPTION(outNumVecs <= 0, std::invalid_argument,
			 "Belos::MultiVecTraits<double, Epetra_MultiVector>::"
			 "Clone(mv, outNumVecs = " << outNumVecs << "): "
			 "outNumVecs must be positive.");
      // FIXME (mfh 13 Jan 2011) Anasazi currently lets Epetra fill in
      // the entries of the returned multivector with zeros, but Belos
      // does not.  We retain this different behavior for now, but the
      // two versions will need to be reconciled.
      return Teuchos::rcp (new Epetra_MultiVector (mv.Map(), outNumVecs)); 
    }

    /*! \brief Creates a new Epetra_MultiVector and copies contents of \c mv into the new vector (deep copy).
      
      \return Reference-counted pointer to the new Epetra_MultiVector.
    */
    static Teuchos::RCP<Epetra_MultiVector> 
    CloneCopy (const Epetra_MultiVector& mv)
    { 
      return Teuchos::rcp (new Epetra_MultiVector (mv)); 
    }

    /*! \brief Creates a new Epetra_MultiVector and copies the selected contents of \c mv into the new vector (deep copy).  

      The copied vectors from \c mv are indicated by the \c indeX.size() indices in \c index.      
      \return Reference-counted pointer to the new Epetra_MultiVector.
    */
    static Teuchos::RCP<Epetra_MultiVector> 
    CloneCopy (const Epetra_MultiVector& mv, const std::vector<int>& index)
    { 
      const int inNumVecs = GetNumberVecs (mv);
      const int outNumVecs = index.size();

      // Simple, inexpensive tests of the index vector.
      TEUCHOS_TEST_FOR_EXCEPTION(outNumVecs == 0, std::invalid_argument,
			 "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::"
			 "CloneCopy(mv, index = {}): At least one vector must be"
			 " cloned from mv.");
      if (outNumVecs > inNumVecs)
	{
	  std::ostringstream os;
	  os << "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::"
	    "CloneCopy(mv, index = {";
	  for (int k = 0; k < outNumVecs - 1; ++k)
	    os << index[k] << ", ";
	  os << index[outNumVecs-1] << "}): There are " << outNumVecs 
	     << " indices to copy, but only " << inNumVecs << " columns of mv.";
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
	}
#ifdef TEUCHOS_DEBUG
      // In debug mode, we perform more expensive tests of the index
      // vector, to ensure all the elements are in range.
      // Dereferencing the iterator is valid because index has length
      // > 0.
      const int minIndex = *std::min_element (index.begin(), index.end());
      const int maxIndex = *std::max_element (index.begin(), index.end());

      if (minIndex < 0)
	{
	  std::ostringstream os;
	  os << "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::"
	    "CloneCopy(mv, index = {";
	  for (int k = 0; k < outNumVecs - 1; ++k)
	    os << index[k] << ", ";
	  os << index[outNumVecs-1] << "}): Indices must be nonnegative, but "
	    "the smallest index " << minIndex << " is negative.";
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
	}
      if (maxIndex >= inNumVecs)
	{
	  std::ostringstream os;
	  os << "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::"
	    "CloneCopy(mv, index = {";
	  for (int k = 0; k < outNumVecs - 1; ++k)
	    os << index[k] << ", ";
	  os << index[outNumVecs-1] << "}): Indices must be strictly less than "
	    "the number of vectors " << inNumVecs << " in mv; the largest index " 
	     << maxIndex << " is out of bounds.";
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
	}
#endif // TEUCHOS_DEBUG
      // Cast to nonconst, because Epetra_MultiVector's constructor
      // wants a nonconst int array argument.  It doesn't actually
      // change the entries of the array.
      std::vector<int>& tmpind = const_cast< std::vector<int>& > (index);
      return Teuchos::rcp (new Epetra_MultiVector (Copy, mv, &tmpind[0], index.size())); 
      // return Teuchos::rcp (new Epetra_MultiVector (::Copy, mv, &tmpind[0], index.size())); 
    }

    static Teuchos::RCP<Epetra_MultiVector> 
    CloneCopy (const Epetra_MultiVector& mv, const Teuchos::Range1D& index)
    { 
      const int inNumVecs = GetNumberVecs (mv);
      const int outNumVecs = index.size();
      const bool validRange = outNumVecs > 0 && index.lbound() >= 0 && 
	index.ubound() < inNumVecs;
      if (! validRange)
	{
	  std::ostringstream os;
	  os <<	"Anasazi::MultiVecTraits<double,Epetra_MultiVector>::Clone(mv,"
	    "index=[" << index.lbound() << ", " << index.ubound() << "]): ";
	  TEUCHOS_TEST_FOR_EXCEPTION(outNumVecs == 0, std::invalid_argument,
			     os.str() << "Column index range must be nonempty.");
	  TEUCHOS_TEST_FOR_EXCEPTION(index.lbound() < 0, std::invalid_argument,
			     os.str() << "Column index range must be nonnegative.");
	  TEUCHOS_TEST_FOR_EXCEPTION(index.ubound() >= inNumVecs, std::invalid_argument,
			     os.str() << "Column index range must not exceed "
			     "number of vectors " << inNumVecs << " in the "
			     "input multivector.");
	}
      return Teuchos::rcp (new Epetra_MultiVector (Copy, mv, index.lbound(), index.size()));
    }

    /*! \brief Creates a new Epetra_MultiVector that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new Epetra_MultiVector.
    */      
    static Teuchos::RCP<Epetra_MultiVector> 
    CloneViewNonConst (Epetra_MultiVector& mv, const std::vector<int>& index)
    { 
      const int inNumVecs = GetNumberVecs (mv);
      const int outNumVecs = index.size();

      // Simple, inexpensive tests of the index vector.
      TEUCHOS_TEST_FOR_EXCEPTION(outNumVecs == 0, std::invalid_argument,
			 "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::"
			 "CloneViewNonConst(mv, index = {}): The output view "
			 "must have at least one column.");
      if (outNumVecs > inNumVecs)
	{
	  std::ostringstream os;
	  os << "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::"
	    "CloneViewNonConst(mv, index = {";
	  for (int k = 0; k < outNumVecs - 1; ++k)
	    os << index[k] << ", ";
	  os << index[outNumVecs-1] << "}): There are " << outNumVecs 
	     << " indices to view, but only " << inNumVecs << " columns of mv.";
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
	}
#ifdef TEUCHOS_DEBUG
      // In debug mode, we perform more expensive tests of the index
      // vector, to ensure all the elements are in range.
      // Dereferencing the iterator is valid because index has length
      // > 0.
      const int minIndex = *std::min_element (index.begin(), index.end());
      const int maxIndex = *std::max_element (index.begin(), index.end());

      if (minIndex < 0)
	{
	  std::ostringstream os;
	  os << "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::"
	    "CloneViewNonConst(mv, index = {";
	  for (int k = 0; k < outNumVecs - 1; ++k)
	    os << index[k] << ", ";
	  os << index[outNumVecs-1] << "}): Indices must be nonnegative, but "
	    "the smallest index " << minIndex << " is negative.";
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
	}
      if (maxIndex >= inNumVecs)
	{
	  std::ostringstream os;
	  os << "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::"
	    "CloneViewNonConst(mv, index = {";
	  for (int k = 0; k < outNumVecs - 1; ++k)
	    os << index[k] << ", ";
	  os << index[outNumVecs-1] << "}): Indices must be strictly less than "
	    "the number of vectors " << inNumVecs << " in mv; the largest index " 
	     << maxIndex << " is out of bounds.";
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
	}
#endif // TEUCHOS_DEBUG
      // Cast to nonconst, because Epetra_MultiVector's constructor
      // wants a nonconst int array argument.  It doesn't actually
      // change the entries of the array.
      std::vector<int>& tmpind = const_cast< std::vector<int>& > (index);
      return Teuchos::rcp (new Epetra_MultiVector (View, mv, &tmpind[0], index.size()));
    }

    static Teuchos::RCP<Epetra_MultiVector> 
    CloneViewNonConst (Epetra_MultiVector& mv, const Teuchos::Range1D& index)
    { 
      const bool validRange = index.size() > 0 && 
	index.lbound() >= 0 && 
	index.ubound() < mv.NumVectors();
      if (! validRange)
	{
	  std::ostringstream os;
	  os <<	"Anasazi::MultiVecTraits<double,Epetra_MultiVector>::CloneView"
	    "NonConst(mv,index=[" << index.lbound() << ", " << index.ubound() 
	     << "]): ";
	  TEUCHOS_TEST_FOR_EXCEPTION(index.size() == 0, std::invalid_argument,
			     os.str() << "Column index range must be nonempty.");
	  TEUCHOS_TEST_FOR_EXCEPTION(index.lbound() < 0, std::invalid_argument,
			     os.str() << "Column index range must be nonnegative.");
	  TEUCHOS_TEST_FOR_EXCEPTION(index.ubound() >= mv.NumVectors(), 
			     std::invalid_argument,
			     os.str() << "Column index range must not exceed "
			     "number of vectors " << mv.NumVectors() << " in "
			     "the input multivector.");
	}
      return Teuchos::rcp (new Epetra_MultiVector (View, mv, index.lbound(), index.size()));
    }

    /*! \brief Creates a new const Epetra_MultiVector that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new const Epetra_MultiVector.
    */      
    static Teuchos::RCP<const Epetra_MultiVector> 
    CloneView (const Epetra_MultiVector& mv, const std::vector<int>& index)
    { 
      const int inNumVecs = GetNumberVecs (mv);
      const int outNumVecs = index.size();

      // Simple, inexpensive tests of the index vector.
      TEUCHOS_TEST_FOR_EXCEPTION(outNumVecs == 0, std::invalid_argument,
			 "Belos::MultiVecTraits<double,Epetra_MultiVector>::"
			 "CloneView(mv, index = {}): The output view "
			 "must have at least one column.");
      if (outNumVecs > inNumVecs)
	{
	  std::ostringstream os;
	  os << "Belos::MultiVecTraits<double,Epetra_MultiVector>::"
	    "CloneView(mv, index = {";
	  for (int k = 0; k < outNumVecs - 1; ++k)
	    os << index[k] << ", ";
	  os << index[outNumVecs-1] << "}): There are " << outNumVecs 
	     << " indices to view, but only " << inNumVecs << " columns of mv.";
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
	}
#ifdef TEUCHOS_DEBUG
      // In debug mode, we perform more expensive tests of the index
      // vector, to ensure all the elements are in range.
      // Dereferencing the iterator is valid because index has length
      // > 0.
      const int minIndex = *std::min_element (index.begin(), index.end());
      const int maxIndex = *std::max_element (index.begin(), index.end());

      if (minIndex < 0)
	{
	  std::ostringstream os;
	  os << "Belos::MultiVecTraits<double,Epetra_MultiVector>::"
	    "CloneView(mv, index = {";
	  for (int k = 0; k < outNumVecs - 1; ++k)
	    os << index[k] << ", ";
	  os << index[outNumVecs-1] << "}): Indices must be nonnegative, but "
	    "the smallest index " << minIndex << " is negative.";
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
	}
      if (maxIndex >= inNumVecs)
	{
	  std::ostringstream os;
	  os << "Belos::MultiVecTraits<double,Epetra_MultiVector>::"
	    "CloneView(mv, index = {";
	  for (int k = 0; k < outNumVecs - 1; ++k)
	    os << index[k] << ", ";
	  os << index[outNumVecs-1] << "}): Indices must be strictly less than "
	    "the number of vectors " << inNumVecs << " in mv; the largest index " 
	     << maxIndex << " is out of bounds.";
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
	}
#endif // TEUCHOS_DEBUG
      // Cast to nonconst, because Epetra_MultiVector's constructor
      // wants a nonconst int array argument.  It doesn't actually
      // change the entries of the array.
      std::vector<int>& tmpind = const_cast< std::vector<int>& > (index);
      return Teuchos::rcp (new Epetra_MultiVector (View, mv, &tmpind[0], index.size()));
    }

    static Teuchos::RCP<Epetra_MultiVector> 
    CloneView (const Epetra_MultiVector& mv, const Teuchos::Range1D& index)
    { 
      const bool validRange = index.size() > 0 && 
	index.lbound() >= 0 && 
	index.ubound() < mv.NumVectors();
      if (! validRange)
	{
	  std::ostringstream os;
	  os <<	"Anasazi::MultiVecTraits<double,Epetra_MultiVector>::CloneView"
	    "(mv,index=[" << index.lbound() << ", " << index.ubound() 
	     << "]): ";
	  TEUCHOS_TEST_FOR_EXCEPTION(index.size() == 0, std::invalid_argument,
			     os.str() << "Column index range must be nonempty.");
	  TEUCHOS_TEST_FOR_EXCEPTION(index.lbound() < 0, std::invalid_argument,
			     os.str() << "Column index range must be nonnegative.");
	  TEUCHOS_TEST_FOR_EXCEPTION(index.ubound() >= mv.NumVectors(), 
			     std::invalid_argument,
			     os.str() << "Column index range must not exceed "
			     "number of vectors " << mv.NumVectors() << " in "
			     "the input multivector.");
	}
      return Teuchos::rcp (new Epetra_MultiVector(View, mv, index.lbound(), index.size()));
    }

    //@}

    //! @name Attribute methods
    //@{ 

    //! Obtain the vector length of \c mv.
    static int GetVecLength( const Epetra_MultiVector& mv )
    { return mv.GlobalLength(); }

    //! Obtain the number of vectors in \c mv
    static int GetNumberVecs( const Epetra_MultiVector& mv )
    { return mv.NumVectors(); }

    static bool HasConstantStride( const Epetra_MultiVector& mv )
    { return mv.ConstantStride(); }
    //@}

    //! @name Update methods
    //@{ 

    /*! \brief Update \c mv with \f$ \alpha AB + \beta mv \f$.
     */
    static void MvTimesMatAddMv( double alpha, const Epetra_MultiVector& A, 
                                 const Teuchos::SerialDenseMatrix<int,double>& B, 
                                 double beta, Epetra_MultiVector& mv )
    { 
      Epetra_LocalMap LocalMap(B.numRows(), 0, mv.Map().Comm());
      Epetra_MultiVector B_Pvec(::View, LocalMap, B.values(), B.stride(), B.numCols());

      TEUCHOS_TEST_FOR_EXCEPTION( mv.Multiply( 'N', 'N', alpha, A, B_Pvec, beta )!=0, EpetraMultiVecFailure,
          "Anasazi::MultiVecTraits<double, Epetra_MultiVector>::MvNorm call to Epetra_MultiVector::Multiply() returned a nonzero value.");
    }

    /*! \brief Replace \c mv with \f$\alpha A + \beta B\f$.
     */
    static void MvAddMv( double alpha, const Epetra_MultiVector& A, double beta, const Epetra_MultiVector& B, Epetra_MultiVector& mv )
    { 
      // epetra mv.Update(alpha,A,beta,B,gamma) will check 
      //   alpha == 0.0 
      // and 
      //   beta == 0.0 
      // and based on this will call 
      //   mv.Update(beta,B,gamma) 
      // or 
      //   mv.Update(alpha,A,gamma)
      //
      // mv.Update(alpha,A,gamma) 
      // will then check for one of 
      //   gamma == 0
      // or 
      //   gamma == 1
      // or 
      //   alpha == 1 
      // in that order. however, it will not look for the combination
      //   alpha == 1 and gamma = 0
      // which is a common use case when we wish to assign 
      //   mv = A   (in which case alpha == 1, beta == gamma == 0)
      // or 
      //   mv = B   (in which case beta == 1, alpha == gamma == 0)
      // 
      // therefore, we will check for these use cases ourselves
      if (beta == 0.0) {
        if (alpha == 1.0) {
          // assign
          mv = A;
        }
        else {
          // single update
          TEUCHOS_TEST_FOR_EXCEPTION( mv.Update( alpha, A, 0.0 )!=0, EpetraMultiVecFailure,
              "Anasazi::MultiVecTraits<double, Epetra_MultiVector>::MvAddMv call to Epetra_MultiVector::Update(alpha,A,0.0) returned a nonzero value.");
        }
      }
      else if (alpha == 0.0) {
        if (beta == 1.0) {
          // assign 
          mv = B;
        }
        else {
          // single update
          TEUCHOS_TEST_FOR_EXCEPTION( mv.Update( beta, B, 0.0 )!=0, EpetraMultiVecFailure,
              "Anasazi::MultiVecTraits<double, Epetra_MultiVector>::MvAddMv call to Epetra_MultiVector::Update(beta,B,0.0) returned a nonzero value.");
        }
      }
      else {
        // double update
        TEUCHOS_TEST_FOR_EXCEPTION( mv.Update( alpha, A, beta, B, 0.0 )!=0, EpetraMultiVecFailure,
            "Anasazi::MultiVecTraits<double, Epetra_MultiVector>::MvAddMv call to Epetra_MultiVector::Update(alpha,A,beta,B,0.0) returned a nonzero value.");
      }
    }

    /*! \brief Compute a dense matrix \c B through the matrix-matrix multiply \f$ \alpha A^Tmv \f$.
    */
    static void MvTransMv( double alpha, const Epetra_MultiVector& A, const Epetra_MultiVector& mv, Teuchos::SerialDenseMatrix<int,double>& B
#ifdef HAVE_ANASAZI_EXPERIMENTAL
                          , ConjType conj = Anasazi::CONJ
#endif
                        )
    { 
      Epetra_LocalMap LocalMap(B.numRows(), 0, mv.Map().Comm());
      Epetra_MultiVector B_Pvec(::View, LocalMap, B.values(), B.stride(), B.numCols());
      
      TEUCHOS_TEST_FOR_EXCEPTION( B_Pvec.Multiply( 'T', 'N', alpha, A, mv, 0.0 )!=0, EpetraMultiVecFailure,
          "Anasazi::MultiVecTraits<double, Epetra_MultiVector>::MvTransMv call to Epetra_MultiVector::Multiply() returned a nonzero value.");
    }
    
    /*! \brief Compute a vector \c b where the components are the individual dot-products of the \c i-th columns of \c A and \c mv, i.e.\f$b[i] = A[i]^Tmv[i]\f$.
     */
    static void MvDot( const Epetra_MultiVector& A, const Epetra_MultiVector& B, std::vector<double> &b
#ifdef HAVE_ANASAZI_EXPERIMENTAL
                      , ConjType conj = Anasazi::CONJ
#endif
                      )
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(A.NumVectors() != B.NumVectors(),std::invalid_argument,
          "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::MvDot(A,B,b): A and B must have the same number of vectors.");
      TEUCHOS_TEST_FOR_EXCEPTION(b.size() != (unsigned int)A.NumVectors(),std::invalid_argument,
          "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::MvDot(A,B,b): b must have room for all dot products.");
#endif
      TEUCHOS_TEST_FOR_EXCEPTION( A.Dot( B, &b[0] )!=0, EpetraMultiVecFailure,
          "Anasazi::MultiVecTraits<double, Epetra_MultiVector>::MvDot(A,B,b) call to Epetra_MultiVector::Dot() returned a nonzero value.");     
    }

    //@}
    //! @name Norm method
    //@{ 

    /*! \brief Compute the 2-norm of each individual vector of \c mv.  
      Upon return, \c normvec[i] holds the value of \f$||mv_i||_2\f$, the \c i-th column of \c mv.
    */
    static void MvNorm( const Epetra_MultiVector& mv, std::vector<double> &normvec )
    { 
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION((unsigned int)mv.NumVectors() != normvec.size(),std::invalid_argument,
          "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::MvNorm(mv,normvec): normvec must be the same size of mv.");
#endif
      TEUCHOS_TEST_FOR_EXCEPTION( mv.Norm2(&normvec[0])!=0, EpetraMultiVecFailure,
          "Anasazi::MultiVecTraits<double, Epetra_MultiVector>::MvNorm call to Epetra_MultiVector::Norm2() returned a nonzero value."); 
    }

    //@}
    
    //! @name Initialization methods
    //@{ 
    /*! \brief Copy the vectors in \c A to a set of vectors in \c mv indicated by the indices given in \c index.
     */
    static void 
    SetBlock (const Epetra_MultiVector& A, 
	      const std::vector<int>& index, 
	      Epetra_MultiVector& mv)
    { 
      const int inNumVecs = GetNumberVecs (A);
      const int outNumVecs = index.size();

      // FIXME (mfh 13 Jan 2011) Belos allows A to have more columns
      // than index.size(), in which case we just take the first
      // index.size() columns of A.  Anasazi requires that A have the
      // same number of columns as index.size().  Changing Anasazi's
      // behavior should not break existing Anasazi solvers, but the
      // tests need to be done.
      if (inNumVecs != outNumVecs)
	{
	  std::ostringstream os;
	  os << "Belos::MultiVecTraits<double,Epetra_MultiVector>::"
	    "SetBlock(A, mv, index = {";
	  if (outNumVecs > 0)
	    {
	      for (int k = 0; k < outNumVecs - 1; ++k)
		os << index[k] << ", ";
	      os << index[outNumVecs-1];
	    }
	  os << "}): A has only " << inNumVecs << " columns, but there are "
	     << outNumVecs << " indices in the index vector.";
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str());
	}
      // Make a view of the columns of mv indicated by the index std::vector.
      Teuchos::RCP<Epetra_MultiVector> mv_view = CloneViewNonConst (mv, index);

      // View of columns [0, outNumVecs-1] of the source multivector A.
      // If A has fewer columns than mv_view, then create a view of
      // the first outNumVecs columns of A.
      Teuchos::RCP<const Epetra_MultiVector> A_view;
      if (outNumVecs == inNumVecs)
	A_view = Teuchos::rcpFromRef (A); // Const, non-owning RCP
      else
	A_view = CloneView (A, Teuchos::Range1D(0, outNumVecs - 1));

      // Assignment calls Epetra_MultiVector::Assign(), which deeply
      // copies the data directly, ignoring the underlying
      // Epetra_Map(s).  If A and mv don't have the same data
      // distribution (Epetra_Map), this may result in incorrect or
      // undefined behavior.  Epetra_MultiVector::Update() also
      // ignores the Epetra_Maps, so we might as well just use the
      // (perhaps slightly cheaper) Assign() method via operator=().
      *mv_view = *A_view;
    }

    static void 
    SetBlock (const Epetra_MultiVector& A, 
	      const Teuchos::Range1D& index, 
	      Epetra_MultiVector& mv)
    { 
      const int numColsA = A.NumVectors();
      const int numColsMv = mv.NumVectors();
      // 'index' indexes into mv; it's the index set of the target.
      const bool validIndex = index.lbound() >= 0 && index.ubound() < numColsMv;
      // We can't take more columns out of A than A has.
      const bool validSource = index.size() <= numColsA;

      if (! validIndex || ! validSource)
	{
	  std::ostringstream os;
	  os <<	"Anasazi::MultiVecTraits<double, Epetra_MultiVector>::SetBlock"
	    "(A, index=[" << index.lbound() << ", " << index.ubound() << "], "
	    "mv): ";
	  TEUCHOS_TEST_FOR_EXCEPTION(index.lbound() < 0, std::invalid_argument,
			     os.str() << "Range lower bound must be nonnegative.");
	  TEUCHOS_TEST_FOR_EXCEPTION(index.ubound() >= numColsMv, std::invalid_argument,
			     os.str() << "Range upper bound must be less than "
			     "the number of columns " << numColsA << " in the "
			     "'mv' output argument.");
	  TEUCHOS_TEST_FOR_EXCEPTION(index.size() > numColsA, std::invalid_argument,
			     os.str() << "Range must have no more elements than"
			     " the number of columns " << numColsA << " in the "
			     "'A' input argument.");
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
	}

      // View of columns [index.lbound(), index.ubound()] of the
      // target multivector mv.  We avoid view creation overhead by
      // only creating a view if the index range is different than [0,
      // (# columns in mv) - 1].
      Teuchos::RCP<Epetra_MultiVector> mv_view;
      if (index.lbound() == 0 && index.ubound()+1 == numColsMv)
	mv_view = Teuchos::rcpFromRef (mv); // Non-const, non-owning RCP
      else
	mv_view = CloneViewNonConst (mv, index);

      // View of columns [0, index.size()-1] of the source multivector
      // A.  If A has fewer columns than mv_view, then create a view
      // of the first index.size() columns of A.
      Teuchos::RCP<const Epetra_MultiVector> A_view;
      if (index.size() == numColsA)
	A_view = Teuchos::rcpFromRef (A); // Const, non-owning RCP
      else
	A_view = CloneView (A, Teuchos::Range1D(0, index.size()-1));

      // Assignment calls Epetra_MultiVector::Assign(), which deeply
      // copies the data directly, ignoring the underlying
      // Epetra_Map(s).  If A and mv don't have the same data
      // distribution (Epetra_Map), this may result in incorrect or
      // undefined behavior.  Epetra_MultiVector::Update() also
      // ignores the Epetra_Maps, so we might as well just use the
      // (perhaps slightly cheaper) Assign() method via operator=().
      *mv_view = *A_view; 
    }

    static void 
    Assign (const Epetra_MultiVector& A, 
	    Epetra_MultiVector& mv)
    {
      const int numColsA = GetNumberVecs (A);
      const int numColsMv = GetNumberVecs (mv);
      if (numColsA > numColsMv)
	{
	  std::ostringstream os;
	  os <<	"Anasazi::MultiVecTraits<double, Epetra_MultiVector>::Assign"
	    "(A, mv): ";
	  TEUCHOS_TEST_FOR_EXCEPTION(numColsA > numColsMv, std::invalid_argument,
			     os.str() << "Input multivector 'A' has " 
			     << numColsA << " columns, but output multivector "
			     "'mv' has only " << numColsMv << " columns.");
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
	}
      // View of the first [0, numColsA-1] columns of mv.
      Teuchos::RCP<Epetra_MultiVector> mv_view;
      if (numColsMv == numColsA)
	mv_view = Teuchos::rcpFromRef (mv); // Non-const, non-owning RCP
      else // numColsMv > numColsA
	mv_view = CloneView (mv, Teuchos::Range1D(0, numColsA - 1));
      
      // Assignment calls Epetra_MultiVector::Assign(), which deeply
      // copies the data directly, ignoring the underlying
      // Epetra_Map(s).  If A and mv don't have the same data
      // distribution (Epetra_Map), this may result in incorrect or
      // undefined behavior.  Epetra_MultiVector::Update() also
      // ignores the Epetra_Maps, so we might as well just use the
      // (perhaps slightly cheaper) Assign() method via operator=().
      *mv_view = A;
    }

    /*! \brief Scale each element of the vectors in \c mv with \c alpha.
     */
    static void MvScale ( Epetra_MultiVector& mv, double alpha ) 
    { 
      TEUCHOS_TEST_FOR_EXCEPTION( mv.Scale( alpha )!=0, EpetraMultiVecFailure,
          "Anasazi::MultiVecTraits<double, Epetra_MultiVector>::MvScale call to Epetra_MultiVector::Scale(mv,double alpha) returned a nonzero value."); 
    }

    /*! \brief Scale each element of the \c i-th vector in \c mv with \c alpha[i].
     */
    static void MvScale ( Epetra_MultiVector& mv, const std::vector<double>& alpha )
    { 
      // Check to make sure the vector is as long as the multivector has columns.
      int numvecs = mv.NumVectors();
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION( alpha.size() != (unsigned int)numvecs, std::invalid_argument,
                          "Anasazi::MultiVecTraits<double, Epetra_MultiVector>::MvScale(mv,vector alpha): size of alpha inconsistent with number of vectors in mv.")
#endif
      for (int i=0; i<numvecs; i++) {
        TEUCHOS_TEST_FOR_EXCEPTION( mv(i)->Scale(alpha[i])!=0, EpetraMultiVecFailure,
            "Anasazi::MultiVecTraits<double, Epetra_MultiVector>::MvScale call to Epetra_MultiVector::Scale() returned a nonzero value.");
      }
    }

    /*! \brief Replace the vectors in \c mv with random vectors.
     */
    static void MvRandom( Epetra_MultiVector& mv )
    { 
      TEUCHOS_TEST_FOR_EXCEPTION( mv.Random()!=0, EpetraMultiVecFailure,
          "Anasazi::MultiVecTraits<double, Epetra_MultiVector>::MvRandom call to Epetra_MultiVector::Random() returned a nonzero value.");
    }

    /*! \brief Replace each element of the vectors in \c mv with \c alpha.
     */
    static void MvInit( Epetra_MultiVector& mv, double alpha = Teuchos::ScalarTraits<double>::zero() )
    { 
      TEUCHOS_TEST_FOR_EXCEPTION( mv.PutScalar(alpha)!=0, EpetraMultiVecFailure,
          "Anasazi::MultiVecTraits<double, Epetra_MultiVector>::MvInit call to Epetra_MultiVector::PutScalar() returned a nonzero value.");
    }

    //@}

    //! @name Print method
    //@{ 

    /*! \brief Print the \c mv multi-vector to the \c os output stream.
     */
    static void MvPrint( const Epetra_MultiVector& mv, std::ostream& os )
    { os << mv << std::endl; }

    //@}

#if defined(HAVE_ANASAZI_TPETRA) && defined(HAVE_ANASAZI_TSQR)
#  if defined(HAVE_TPETRA_EPETRA)
    /// \typedef tsqr_adaptor_type
    /// \brief TsqrAdaptor specialization for Epetra_MultiVector
    ///
    /// \note This lives in Tpetra, for various hackish reasons.
    ///
    typedef Epetra::TsqrAdaptor tsqr_adaptor_type;
#  endif // defined(HAVE_TPETRA_EPETRA)
#endif // defined(HAVE_ANASAZI_TPETRA) && defined(HAVE_ANASAZI_TSQR)
  };        

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Anasazi::OperatorTraits for Epetra::Operator.
  //
  ////////////////////////////////////////////////////////////////////

  /*! 
    \brief Template specialization of Anasazi::OperatorTraits class using the Epetra_Operator virtual base class and 
    Epetra_MultiVector class.

    This interface will ensure that any Epetra_Operator and Epetra_MultiVector will be accepted by the Anasazi
    templated solvers.

    \note The Epetra package performs double-precision arithmetic, so the use of Epetra with Anasazi will
    only provide a double-precision eigensolver.
  */

  template <> 
  class OperatorTraits < double, Epetra_MultiVector, Epetra_Operator >
  {
  public:
    
    /*! \brief This method takes the Epetra_MultiVector \c x and
      applies the Epetra_Operator \c Op to it resulting in the Epetra_MultiVector \c y.
    */    
    static void Apply ( const Epetra_Operator& Op, 
                        const Epetra_MultiVector& x, 
                        Epetra_MultiVector& y )
    { 
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(x.NumVectors() != y.NumVectors(),std::invalid_argument,
          "Anasazi::OperatorTraits<double,Epetra_MultiVector,Epetra_Operator>::Apply(Op,x,y): x and y must have the same number of columns.");
#endif
      int ret = Op.Apply(x,y);
      TEUCHOS_TEST_FOR_EXCEPTION(ret != 0, OperatorError, 
          "Anasazi::OperatorTraits<double,Epetra_Multivector,Epetra_Operator>::Apply(): Error in Epetra_Operator::Apply(). Code " << ret);
    }
    
  };

  template<>
  class MultiVecTraitsExt<double, Epetra_MultiVector>
  {
  public:
    //! @name New attribute methods
    //@{

    //! Obtain the vector length of \c mv.
    //! \note This method supersedes GetVecLength, which will be deprecated.
    static ptrdiff_t GetGlobalLength( const Epetra_MultiVector& mv )
    { 
      if (mv.Map().GlobalIndicesLongLong())
        return static_cast<ptrdiff_t>( mv.GlobalLength64() );
      else
        return static_cast<ptrdiff_t>( mv.GlobalLength() );
    }

    //@}
  };
  
} // end of Anasazi namespace 

#endif 
// end of file ANASAZI_EPETRA_ADAPTER_HPP
