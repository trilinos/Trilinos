// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ANASAZI_BASIC_EIGENPROBLEM_H
#define ANASAZI_BASIC_EIGENPROBLEM_H

/*! \file AnasaziBasicEigenproblem.hpp
  \brief Basic implementation of the Anasazi::Eigenproblem class
*/

#include "AnasaziEigenproblem.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"

/*! \class Anasazi::BasicEigenproblem
  \brief This provides a basic implementation for defining standard or 
  generalized eigenvalue problems.
*/

namespace Anasazi {
  
  template<class ScalarType, class MV, class OP>
  class BasicEigenproblem : public virtual Eigenproblem<ScalarType, MV, OP> {
    
  public:
    
    //! @name Constructors/Destructor
    //@{ 
    
    //! Empty constructor - allows Anasazi::BasicEigenproblem to be described at a later time through "Set Methods".
    BasicEigenproblem();
    
    //! Standard Eigenvalue Problem Constructor.
    BasicEigenproblem( const Teuchos::RCP<const OP>& Op, const Teuchos::RCP<MV>& InitVec );
    
    //! Generalized Eigenvalue Problem Constructor.
    BasicEigenproblem( const Teuchos::RCP<const OP>& Op, const Teuchos::RCP<const OP>& B, const Teuchos::RCP<MV>& InitVec );
    
    //! Copy Constructor.
    BasicEigenproblem( const BasicEigenproblem<ScalarType, MV, OP>& Problem );
    
    //! Destructor.
    virtual ~BasicEigenproblem() {};
    //@}
    
    //! @name Set Methods
    //@{ 
    
    /*! \brief Set the operator for which eigenvalues will be computed.  

    \note This may be different from the \c A if a spectral transformation is employed. 
    For example, this operator may apply the operation \f$(A-\sigma I)^{-1}\f$ if you are
    looking for eigenvalues of \c A around \f$\sigma\f$.  
    */
    void setOperator( const Teuchos::RCP<const OP>& Op ) { _Op = Op; _isSet=false; };
    
    /*! \brief Set the operator \c A of the eigenvalue problem \f$ Ax=Mx\lambda\f$.
    */
    void setA( const Teuchos::RCP<const OP>& A ) { _AOp = A; _isSet=false; };
    
    /*! \brief Set the operator \c M of the eigenvalue problem \f$ Ax = Mx\lambda\f$.
     */
    void setM( const Teuchos::RCP<const OP>& M ) { _MOp = M; _isSet=false; };
    
    /*! \brief Set the preconditioner for this eigenvalue problem \f$ Ax = Mx\lambda\f$.
     */
    void setPrec( const Teuchos::RCP<const OP>& Prec ) { _Prec = Prec; _isSet=false; };
    
    /*! \brief Set the initial guess.  

    This vector is required to create all the space needed 
    by Anasazi to solve the eigenvalue problem.  

    \note Even if an initial guess is not known by the user, an initial vector must be passed in.  
    */
    void setInitVec( const Teuchos::RCP<MV>& InitVec ) { _InitVec = InitVec; _isSet=false; };
    
    /*! \brief Set auxiliary vectors.

    \note This multivector can have any number of columns, and most likely will contain vectors that
    will be used by the eigensolver to orthogonalize against.
    */
    void setAuxVecs( const Teuchos::RCP<const MV>& AuxVecs ) { _AuxVecs = AuxVecs; _isSet=false; };

    //! Specify the number of eigenvalues (NEV) that are requested.
    void setNEV( int nev ){ _nev = nev; _isSet=false; };

    //! Specify the symmetry of this eigenproblem.
    /*! This knowledge may allow the solver to take advantage of the eigenproblems' symmetry.
      Some computational work can be avoided by setting this properly.
    */
    void setHermitian( bool isSym ){ _isSym = isSym; _isSet=false; };

    /*! \brief Specify that this eigenproblem is fully defined.
     *
     * This routine serves multiple purpose:
     *    - sanity check that the eigenproblem has been fully and consistently defined
     *    - opportunity for the eigenproblem to allocate internal storage for eigenvalues
     * and eigenvectors (to be used by eigensolvers and solver managers)
     * </ul>
     *
     * This method reallocates internal storage, so that any previously retrieved references to 
     * internal storage (eigenvectors or eigenvalues) are invalidated.
     *
     * \note The user MUST call this routine before they send the eigenproblem to any solver or solver manager.
     *
     * \returns \c true signifies success, \c false signifies error.
     */
    bool setProblem();

    /*! \brief Set the solution to the eigenproblem.
     *
     * This mechanism allows an Eigensolution struct to be associated with an Eigenproblem object.
     * setSolution() is usually called by a solver manager at the end of its SolverManager::solve() 
     * routine.
     */
    void setSolution(const Eigensolution<ScalarType,MV> &sol) {_sol = sol;}

    //@}
    
    //! @name Accessor Methods
    //@{ 
    
    //! Get a pointer to the operator for which eigenvalues will be computed.
    Teuchos::RCP<const OP> getOperator() const { return( _Op ); };
    
    //! Get a pointer to the operator \c A of the eigenproblem \f$ Ax=\lambda Mx\f$.
    Teuchos::RCP<const OP> getA() const { return( _AOp ); };
    
    //! Get a pointer to the operator \c M of the eigenproblem \f$ Ax=\lambda Mx\f$.
    Teuchos::RCP<const OP> getM() const { return( _MOp ); };
    
    //! Get a pointer to the preconditioner of the eigenproblem \f$ Ax=\lambda Mx\f$.
    Teuchos::RCP<const OP> getPrec() const { return( _Prec ); };
    
    //! Get a pointer to the initial vector
    Teuchos::RCP<const MV> getInitVec() const { return( _InitVec ); };
    
    //! Get a pointer to the auxiliary vector
    Teuchos::RCP<const MV> getAuxVecs() const { return( _AuxVecs ); };
    
    //! Get the number of eigenvalues (NEV) that are required by this eigenproblem.
    int getNEV() const { return( _nev ); }

    //! Get the symmetry information for this eigenproblem.
    bool isHermitian() const { return( _isSym ); }
    
    //! If the problem has been set, this method will return true.
    bool isProblemSet() const { return( _isSet ); }

    /*! \brief Get the solution to the eigenproblem.
     *
     * There is no computation associated with this method. It only provides a 
     * mechanism for associating an Eigensolution with a Eigenproblem.
     */
    const Eigensolution<ScalarType,MV> & getSolution() const { return(_sol); }

    //@}
    
    //! @name Apply / Compute Methods
    //@{

    /*! \brief Returns the residual vector corresponding to the computed solution
     *
     * If there is no computed solution, returns Teuchos::null.
     *
     * \note If A has not been set, this function will compute \f$ R = OpX-MX\Lambda\f$
     * If you are using a spectral transformation, this may not be what you want.
     */
    Teuchos::RCP<const MV> computeCurrResVec() const;

    //@}

  protected:
    
    //! Reference-counted pointer for \c A of the eigenproblem \f$ Ax=\lambda Mx\f$
    Teuchos::RCP<const OP> _AOp;

    //! Reference-counted pointer for \c M of the eigenproblem \f$ Ax=\lambda Mx\f$
    Teuchos::RCP<const OP> _MOp; 

    //! Reference-counted pointer for the operator of the eigenproblem \f$ Ax=\lambda Mx\f$
    Teuchos::RCP<const OP> _Op;

    //! Reference-counted pointer for the preconditioner of the eigenproblem \f$ Ax=\lambda Mx\f$
    Teuchos::RCP<const OP> _Prec;

    //! Reference-counted pointer for the initial vector of the eigenproblem \f$ Ax=\lambda Mx\f$
    Teuchos::RCP<MV> _InitVec;

    //! Reference-counted pointer for the auxiliary vector of the eigenproblem \f$ Ax=\lambda Mx\f$
    Teuchos::RCP<const MV> _AuxVecs;

    //! Number of eigenvalues requested
    int _nev;

    //! Symmetry of the eigenvalue problem
    /*! \note A generalized eigenvalue problem \f$ Ax= \lambda Mx\f$ is considered symmetric
      if the operator \c M is positive (semi) definite.
    */
    bool _isSym;

    //! Sanity Check Flag
    bool _isSet;

    //! Type-definition for the MultiVecTraits class corresponding to the \c MV type
    typedef MultiVecTraits<ScalarType,MV> MVT;
    //! Type-definition for the OperatorTraits class corresponding to the \c OP type
    typedef OperatorTraits<ScalarType,MV,OP> OPT;

    //! Solution to problem
    Eigensolution<ScalarType,MV> _sol;
  };


  //=============================================================================
  //     Implementations (Constructors / Destructors)
  //=============================================================================
  template <class ScalarType, class MV, class OP>
  BasicEigenproblem<ScalarType, MV, OP>::BasicEigenproblem() : 
    _nev(0), 
    _isSym(false),
    _isSet(false)
  {
  }


  //=============================================================================
  template <class ScalarType, class MV, class OP>
  BasicEigenproblem<ScalarType, MV, OP>::BasicEigenproblem( const Teuchos::RCP<const OP>& Op, const Teuchos::RCP<MV>& InitVec ) :    
    _Op(Op), 
    _InitVec(InitVec), 
    _nev(0), 
    _isSym(false),
    _isSet(false)
  {
  }


  //=============================================================================
  template <class ScalarType, class MV, class OP>
  BasicEigenproblem<ScalarType, MV, OP>::BasicEigenproblem( const Teuchos::RCP<const OP>& Op, const Teuchos::RCP<const OP>& M,
                                                            const Teuchos::RCP<MV>& InitVec ) :
    _MOp(M), 
    _Op(Op), 
    _InitVec(InitVec), 
    _nev(0), 
    _isSym(false),
    _isSet(false)
  {
  }


  //=============================================================================
  template <class ScalarType, class MV, class OP>
  BasicEigenproblem<ScalarType, MV, OP>::BasicEigenproblem( const BasicEigenproblem<ScalarType,MV,OP>& Problem ) :
    _AOp(Problem._AOp), 
    _MOp(Problem._MOp), 
    _Op(Problem._Op), 
    _Prec(Problem._Prec), 
    _InitVec(Problem._InitVec), 
    _nev(Problem._nev), 
    _isSym(Problem._isSym),
    _isSet(Problem._isSet),
    _sol(Problem._sol)
  {
  }


  //=============================================================================
  //    SetProblem (sanity check method)
  //=============================================================================
  template <class ScalarType, class MV, class OP>
  bool BasicEigenproblem<ScalarType, MV, OP>::setProblem() 
  {
    //----------------------------------------------------------------
    // Sanity Checks
    //----------------------------------------------------------------
    // If there is no operator, then we can't proceed.
    if ( !_AOp.get() && !_Op.get() ) { return false; }
    
    // If there is no initial vector, then we don't have anything to clone workspace from.
    if ( !_InitVec.get() ) { return false; }
    
    // If we don't need any eigenvalues, we don't need to continue.
    if (_nev == 0) { return false; }
    
    // If there is an A, but no operator, we can set them equal.
    if (_AOp.get() && !_Op.get()) { _Op = _AOp; }

    // Clear the storage from any previous call to setSolution()
    Eigensolution<ScalarType,MV> emptysol;
    _sol = emptysol;
    
    // mark the problem as set and return no-error
    _isSet=true;
    return true;
  }

  //=============================================================================
  //    Computes the residual vector
  //=============================================================================
  template <class ScalarType, class MV, class OP>
  Teuchos::RCP<const MV> BasicEigenproblem<ScalarType, MV, OP>::computeCurrResVec() const
  {
    using Teuchos::RCP;

    TEUCHOS_TEST_FOR_EXCEPTION(!isHermitian(), std::invalid_argument,
        "BasicEigenproblem::computeCurrResVec: This method is not currently "
        "implemented for non-Hermitian problems.  Sorry for any inconvenience.");

    const Eigensolution<ScalarType,MV> sol = getSolution();

    if(sol.numVecs <= 0)
      return Teuchos::null;

    // Copy the eigenvalues
    RCP<MV> X = sol.Evecs;
    std::vector<ScalarType> Lambda(sol.numVecs);
    for(int i = 0; i < sol.numVecs; i++)
      Lambda[i] = sol.Evals[i].realpart;

    // Compute AX
    RCP<MV> AX = MVT::Clone(*X,sol.numVecs);
    if(_AOp != Teuchos::null)
      OPT::Apply(*_AOp,*X,*AX);
    else
      OPT::Apply(*_Op,*X,*AX);

    // Compute MX if necessary
    RCP<MV> MX;
    if(_MOp != Teuchos::null)
    {
      MX = MVT::Clone(*X,sol.numVecs);
      OPT::Apply(*_MOp,*X,*MX);
    }
    else
    {
      MX = MVT::CloneCopy(*X);
    }

    // Compute R = AX - MX \Lambda
    RCP<MV> R = MVT::Clone(*X,sol.numVecs);
    MVT::MvScale(*MX,Lambda);
    MVT::MvAddMv(1.0,*AX,-1.0,*MX,*R);

    return R;
  }
  
} // end Anasazi namespace
#endif

// end AnasaziBasicEigenproblem.hpp
