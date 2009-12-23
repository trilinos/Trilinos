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

#include "Teuchos_TestForException.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"

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
  class ANASAZIEPETRA_LIB_DLL_EXPORT EpetraMultiVec : public MultiVec<double>, public Epetra_MultiVector {
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
    MultiVec<double> * CloneView ( const std::vector<int>& index );

    //@}

    //! @name Attribute methods
    //@{ 

    //! Obtain the vector length of *this.
    int GetNumberVecs () const { return NumVectors(); }

    //! Obtain the number of vectors in *this.
    int GetVecLength () const { return GlobalLength(); }

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
      TEST_FOR_EXCEPTION( this->Scale( alpha )!=0, EpetraMultiVecFailure,
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
        TEST_FOR_EXCEPTION( this->Norm2(&normvec[0])!=0, EpetraMultiVecFailure,
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
      TEST_FOR_EXCEPTION( this->Random()!=0, EpetraMultiVecFailure,
          "Anasazi::EpetraMultiVec::MvRandom call to Epetra_MultiVector::Random() returned a nonzero value.");
    };

    /*! \brief Replace each element of the vectors in \c *this with \c alpha.
     */
    void MvInit ( double alpha ) { 
      TEST_FOR_EXCEPTION( this->PutScalar( alpha )!=0, EpetraMultiVecFailure,
          "Anasazi::EpetraMultiVec::MvInit call to Epetra_MultiVector::PutScalar() returned a nonzero value.");
    };
    
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
#pragma warning(push)
#pragma warning(disable:4251)
    Teuchos::RCP<Epetra_Operator> Epetra_Op;
#pragma warning(pop)
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
#pragma warning(push)
#pragma warning(disable:4251)
    Teuchos::RCP<Epetra_Operator> Epetra_AOp;
    Teuchos::RCP<Epetra_Operator> Epetra_MOp;
#pragma warning(pop)
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
#pragma warning(push)
#pragma warning(disable:4251)
    Teuchos::RCP<Epetra_Operator> Epetra_Op;
#pragma warning(pop)

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
#pragma warning(push)
#pragma warning(disable:4251)
    Teuchos::RCP<const Epetra_MultiVector> Epetra_MV;
    Teuchos::RCP<const Epetra_Map> MV_localmap;
    Teuchos::RCP<const Epetra_BlockMap> MV_blockmap;
#pragma warning(pop)

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
#pragma warning(push)
#pragma warning(disable:4251)
    Teuchos::RCP<const Epetra_MultiVector> Epetra_MV;
    Teuchos::RCP<Epetra_Operator> Epetra_OP;
    Teuchos::RCP<Epetra_MultiVector> Epetra_WMV;
    Teuchos::RCP<const Epetra_Map> MV_localmap;
    Teuchos::RCP<const Epetra_BlockMap> MV_blockmap;
#pragma warning(pop)
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
#pragma warning(push)
#pragma warning(disable:4251)
    Teuchos::RCP<const Epetra_MultiVector> Epetra_MV;
    Teuchos::RCP<Epetra_Operator> Epetra_OP;
    Teuchos::RCP<Epetra_MultiVector> Epetra_WMV;
    Teuchos::RCP<const Epetra_Map> MV_localmap;
    Teuchos::RCP<const Epetra_BlockMap> MV_blockmap;
#pragma warning(pop)
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

    /*! \brief Creates a new empty Epetra_MultiVector containing \c numvecs columns.
      
    \return Reference-counted pointer to the new Epetra_MultiVector.
    */
    static Teuchos::RCP<Epetra_MultiVector> Clone( const Epetra_MultiVector& mv, const int numvecs )
    { 
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(numvecs <= 0,std::invalid_argument,
          "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::Clone(mv,numvecs): numvecs must be greater than zero.");
#endif
      return Teuchos::rcp( new Epetra_MultiVector(mv.Map(), numvecs) ); 
    }

    /*! \brief Creates a new Epetra_MultiVector and copies contents of \c mv into the new vector (deep copy).
      
      \return Reference-counted pointer to the new Epetra_MultiVector.
    */
    static Teuchos::RCP<Epetra_MultiVector> CloneCopy( const Epetra_MultiVector& mv )
    { return Teuchos::rcp( new Epetra_MultiVector( mv ) ); }

    /*! \brief Creates a new Epetra_MultiVector and copies the selected contents of \c mv into the new vector (deep copy).  

      The copied vectors from \c mv are indicated by the \c indeX.size() indices in \c index.      
      \return Reference-counted pointer to the new Epetra_MultiVector.
    */
    static Teuchos::RCP<Epetra_MultiVector> CloneCopy( const Epetra_MultiVector& mv, const std::vector<int>& index )
    { 
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(index.size() == 0,std::invalid_argument,
          "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::Clone(mv,index): numvecs must be greater than zero.");
      TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::invalid_argument,
          "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::Clone(mv,index): indices must be >= zero.");
      TEST_FOR_EXCEPTION( *std::max_element(index.begin(),index.end()) >= mv.NumVectors(), std::invalid_argument,
          "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::Clone(mv,index): indices must be < mv.NumVectors().");
#endif
      std::vector<int>& tmp_index = const_cast<std::vector<int> &>( index );
      return Teuchos::rcp( new Epetra_MultiVector(::Copy, mv, &tmp_index[0], index.size()) ); 
    }

    /*! \brief Creates a new Epetra_MultiVector that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new Epetra_MultiVector.
    */      
    static Teuchos::RCP<Epetra_MultiVector> CloneView( Epetra_MultiVector& mv, const std::vector<int>& index )
    { 
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(index.size() == 0,std::invalid_argument,
          "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::CloneView(mv,index): numvecs must be greater than zero.");
      TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::invalid_argument,
          "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::CloneView(mv,index): indices must be >= zero.");
      TEST_FOR_EXCEPTION( *std::max_element(index.begin(),index.end()) >= mv.NumVectors(), std::invalid_argument,
          "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::CloneView(mv,index): indices must be < mv.NumVectors().");
#endif
      std::vector<int>& tmp_index = const_cast<std::vector<int> &>( index );
      return Teuchos::rcp( new Epetra_MultiVector(::View, mv, &tmp_index[0], index.size()) ); 
    }

    /*! \brief Creates a new const Epetra_MultiVector that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new const Epetra_MultiVector.
    */      
    static Teuchos::RCP<const Epetra_MultiVector> CloneView( const Epetra_MultiVector& mv, const std::vector<int>& index )
    { 
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(index.size() == 0,std::invalid_argument,
          "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::Clone(const mv,index): numvecs must be greater than zero.");
      TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::invalid_argument,
          "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::Clone(const mv,index): indices must be >= zero.");
      TEST_FOR_EXCEPTION( *std::max_element(index.begin(),index.end()) >= mv.NumVectors(), std::invalid_argument,
          "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::Clone(const mv,index): indices must be < mv.NumVectors().");
#endif
      std::vector<int>& tmp_index = const_cast<std::vector<int> &>( index );
      return Teuchos::rcp( new Epetra_MultiVector(::View, mv, &tmp_index[0], index.size()) ); 
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

      TEST_FOR_EXCEPTION( mv.Multiply( 'N', 'N', alpha, A, B_Pvec, beta )!=0, EpetraMultiVecFailure,
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
          TEST_FOR_EXCEPTION( mv.Update( alpha, A, 0.0 )!=0, EpetraMultiVecFailure,
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
          TEST_FOR_EXCEPTION( mv.Update( beta, B, 0.0 )!=0, EpetraMultiVecFailure,
              "Anasazi::MultiVecTraits<double, Epetra_MultiVector>::MvAddMv call to Epetra_MultiVector::Update(beta,B,0.0) returned a nonzero value.");
        }
      }
      else {
        // double update
        TEST_FOR_EXCEPTION( mv.Update( alpha, A, beta, B, 0.0 )!=0, EpetraMultiVecFailure,
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
      
      TEST_FOR_EXCEPTION( B_Pvec.Multiply( 'T', 'N', alpha, A, mv, 0.0 )!=0, EpetraMultiVecFailure,
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
      TEST_FOR_EXCEPTION(A.NumVectors() != B.NumVectors(),std::invalid_argument,
          "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::MvDot(A,B,b): A and B must have the same number of vectors.");
      TEST_FOR_EXCEPTION(b.size() != (unsigned int)A.NumVectors(),std::invalid_argument,
          "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::MvDot(A,B,b): b must have room for all dot products.");
#endif
      TEST_FOR_EXCEPTION( A.Dot( B, &b[0] )!=0, EpetraMultiVecFailure,
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
      TEST_FOR_EXCEPTION((unsigned int)mv.NumVectors() != normvec.size(),std::invalid_argument,
          "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::MvNorm(mv,normvec): normvec must be the same size of mv.");
#endif
      TEST_FOR_EXCEPTION( mv.Norm2(&normvec[0])!=0, EpetraMultiVecFailure,
          "Anasazi::MultiVecTraits<double, Epetra_MultiVector>::MvNorm call to Epetra_MultiVector::Norm2() returned a nonzero value."); 
    }

    //@}
    
    //! @name Initialization methods
    //@{ 
    /*! \brief Copy the vectors in \c A to a set of vectors in \c mv indicated by the indices given in \c index.
     */
    static void SetBlock( const Epetra_MultiVector& A, const std::vector<int>& index, Epetra_MultiVector& mv )
    { 
      TEST_FOR_EXCEPTION((unsigned int)A.NumVectors() != index.size(),std::invalid_argument,
          "Anasazi::MultiVecTraits<double,Epetra_MultiVector>::SetBlock(A,index,mv): index must be the same size as A");
      Teuchos::RCP<Epetra_MultiVector> mv_sub = CloneView(mv,index);
      *mv_sub = A;
    }

    /*! \brief Scale each element of the vectors in \c mv with \c alpha.
     */
    static void MvScale ( Epetra_MultiVector& mv, double alpha ) 
    { 
      TEST_FOR_EXCEPTION( mv.Scale( alpha )!=0, EpetraMultiVecFailure,
          "Anasazi::MultiVecTraits<double, Epetra_MultiVector>::MvScale call to Epetra_MultiVector::Scale(mv,double alpha) returned a nonzero value."); 
    }

    /*! \brief Scale each element of the \c i-th vector in \c mv with \c alpha[i].
     */
    static void MvScale ( Epetra_MultiVector& mv, const std::vector<double>& alpha )
    { 
      // Check to make sure the vector is as long as the multivector has columns.
      int numvecs = mv.NumVectors();
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION( alpha.size() != (unsigned int)numvecs, std::invalid_argument,
                          "Anasazi::MultiVecTraits<double, Epetra_MultiVector>::MvScale(mv,vector alpha): size of alpha inconsistent with number of vectors in mv.")
#endif
      for (int i=0; i<numvecs; i++) {
        TEST_FOR_EXCEPTION( mv(i)->Scale(alpha[i])!=0, EpetraMultiVecFailure,
            "Anasazi::MultiVecTraits<double, Epetra_MultiVector>::MvScale call to Epetra_MultiVector::Scale() returned a nonzero value.");
      }
    }

    /*! \brief Replace the vectors in \c mv with random vectors.
     */
    static void MvRandom( Epetra_MultiVector& mv )
    { 
      TEST_FOR_EXCEPTION( mv.Random()!=0, EpetraMultiVecFailure,
          "Anasazi::MultiVecTraits<double, Epetra_MultiVector>::MvRandom call to Epetra_MultiVector::Random() returned a nonzero value.");
    }

    /*! \brief Replace each element of the vectors in \c mv with \c alpha.
     */
    static void MvInit( Epetra_MultiVector& mv, double alpha = Teuchos::ScalarTraits<double>::zero() )
    { 
      TEST_FOR_EXCEPTION( mv.PutScalar(alpha)!=0, EpetraMultiVecFailure,
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
      TEST_FOR_EXCEPTION(x.NumVectors() != y.NumVectors(),std::invalid_argument,
          "Anasazi::OperatorTraits<double,Epetra_MultiVector,Epetra_Operator>::Apply(Op,x,y): x and y must have the same number of columns.");
#endif
      int ret = Op.Apply(x,y);
      TEST_FOR_EXCEPTION(ret != 0, OperatorError, 
          "Anasazi::OperatorTraits<double,Epetra_Multivector,Epetra_Operator>::Apply(): Error in Epetra_Operator::Apply(). Code " << ret);
    }
    
  };
  
} // end of Anasazi namespace 

#endif 
// end of file ANASAZI_EPETRA_ADAPTER_HPP
