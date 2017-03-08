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

#ifndef ANASAZI_EIGENPROBLEM_H
#define ANASAZI_EIGENPROBLEM_H

/*! \file AnasaziEigenproblem.hpp
  \brief Abstract base class which defines the interface required by an eigensolver and
  status test class to compute solutions to an eigenproblem
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_RCP.hpp"


/*! \class Anasazi::Eigenproblem
    \brief This class defines the interface required by an eigensolver and status
    test class to compute solutions to an eigenproblem.
*/

namespace Anasazi {
  
  template<class ScalarType, class MV, class OP>
  class Eigenproblem {

  public:

    //! @name Constructors/Destructor
    //@{ 
    
    //! Empty constructor 
    Eigenproblem() {};
    
    //! Destructor.
    virtual ~Eigenproblem() {};
    //@}
    
    //! @name Set Methods
    //@{ 
    
    /*! \brief Set the operator for which eigenvalues will be computed.  
     * 
     * \note This may be different from the \c A if a spectral transformation
     * is employed.  For example, this operator may apply the operation
     * \f$(A-\sigma I)^{-1}\f$ if you are looking for eigenvalues of \c A
     * around \f$\sigma\f$.  
     */
    virtual void setOperator( const Teuchos::RCP<const OP> &Op ) = 0;

    //! \brief Set the operator \c A of the eigenvalue problem \f$Ax=\lambda Mx\f$.
    virtual void setA( const Teuchos::RCP<const OP> &A ) = 0;

    //! \brief Set the operator \c M of the eigenvalue problem \f$Ax=\lambda Mx\f$.
    virtual void setM( const Teuchos::RCP<const OP> &M ) = 0;

    //! \brief Set the preconditioner for this eigenvalue problem \f$Ax=\lambda Mx\f$.
    virtual void setPrec( const Teuchos::RCP<const OP> &Prec ) = 0;

    /*! \brief Set the initial guess.  
     *
     * \note This multivector should have the same number of columns as the blocksize.
     */
    virtual void setInitVec( const Teuchos::RCP<MV> &InitVec ) = 0; 

    /*! \brief Set auxiliary vectors. 
     *
     * \note This multivector can have any number of columns, and most likely
     * will contain vectors that will be used by the eigensolver to
     * orthogonalize against.
     */
    virtual void setAuxVecs( const Teuchos::RCP<const MV> &AuxVecs ) = 0;

    //! The number of eigenvalues (NEV) that are requested.
    virtual void setNEV( int nev ) = 0;

    /*! \brief Specify the symmetry of the eigenproblem.
     *
     *  This knowledge may allow the solver to take advantage of the eigenproblems' symmetry.
     *  Some computational work may be avoided by setting this properly.
     */
    virtual void setHermitian( bool isSym ) = 0;

    /*! \brief Specify that this eigenproblem is fully defined.
     *
     * This routine serves multiple purpose:
     * <ul>
     * <li> sanity check that the eigenproblem has been fully and consistently defined
     * <li> opportunity for the eigenproblem to allocate internal storage for eigenvalues
     * and eigenvectors (to be used by eigensolvers and solver managers)
     * </ul>
     *
     * \note The user MUST call this routine before they send the eigenproblem to any solver or solver manager.
     *
     * \returns \c true signifies success, \c false signifies error.
     */
    virtual bool setProblem() = 0;

    /*! \brief Set the solution to the eigenproblem.
     *
     * This mechanism allows an Eigensolution struct to be associated with an Eigenproblem object.
     * setSolution() is usually called by a solver manager at the end of its SolverManager::solve() 
     * routine.
     */
    virtual void setSolution(const Eigensolution<ScalarType,MV> &sol) = 0;

    //@}
    
    //! @name Accessor Methods
    //@{ 

    //! Get a pointer to the operator for which eigenvalues will be computed.
    virtual Teuchos::RCP<const OP> getOperator() const = 0;

    //! Get a pointer to the operator \c A of the eigenproblem \f$AX=\lambda Mx\f$.
    virtual Teuchos::RCP<const OP> getA() const = 0;

    //! Get a pointer to the operator \c M of the eigenproblem \f$AX=\lambda Mx\f$.
    virtual Teuchos::RCP<const OP> getM() const = 0;

    //! Get a pointer to the preconditioner.
    virtual Teuchos::RCP<const OP> getPrec() const = 0;

    //! Get a pointer to the initial vector
    virtual Teuchos::RCP<const MV> getInitVec() const = 0;

    //! Get a pointer to the auxiliary vector
    virtual Teuchos::RCP<const MV> getAuxVecs() const = 0;

    //! Get the number of eigenvalues (NEV) that are required by this eigenproblem.
    virtual int getNEV() const = 0;

    //! Get the symmetry information for this eigenproblem.
    virtual bool isHermitian() const = 0;

    //! If the problem has been set, this method will return true.
    virtual bool isProblemSet() const = 0;

    /*! \brief Get the solution to the eigenproblem.
     *
     * There is no computation associated with this method. It only provides a 
     * mechanism for associating an Eigensolution with a Eigenproblem.
     */
    virtual const Eigensolution<ScalarType,MV> & getSolution() const = 0;

    //@}	
  };
   
} // end Anasazi namespace
#endif

// end AnasaziEigenproblem.hpp
