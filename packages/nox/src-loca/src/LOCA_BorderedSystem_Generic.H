// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_BORDEREDSYSTEM_GENERIC_H
#define LOCA_BORDEREDSYSTEM_GENERIC_H

#include "NOX_Abstract_Group.H"  // for ReturnType, Parameter::List, etc
#include "NOX_Abstract_MultiVector.H" // for DenseMatrix

namespace LOCA {

  //! %Strategies for solving bordered systems of equations
  namespace BorderedSystem {
    
    //! %Generic bordered solver interface.
    /*! 
     * %Generic interface for solving systems of equations of the form
     * \f[
     *     \begin{bmatrix}
     *          J & A    \\
     *        B^T & C
     *     \end{bmatrix}
     *     \begin{bmatrix}
     *        X \\
     *        Y 
     *     \end{bmatrix} = 
     *     \begin{bmatrix}
     *        F \\
     *        G
     *     \end{bmatrix}
     * \f]
     * where \f$J\f$ is an \f$n\times n\f$ matrix, \f$A\f$ and \f$B\f$ are
     * \f$n\times m\f$, \f$C\f$ is \f$m\times m\f$, \f$X\f$ and \f$F\f$ are
     * \f$n\times p\f$ and \f$Y\f$ and \f$G\f$ are \f$m\times p\f$.  The action
     * of \f$J\f$ and its inverse are represnted by a NOX::Abstract::Group
     * while \f$A\f$ and \f$B\f$ are NOX::Abstract::MultiVector's and
     * \f$C\f$ is a dense matrix.  All classes the implement a method for 
     * computing solutions to this system of equations should be derived from 
     * this class.
    */
    class Generic {

    public:

      //! Constructor. 
      /*! 
	Constructors of derived objects should look like reset.
      */
      Generic() {};

      //! Copy constructor
      Generic(const Generic& source) {}

      //! Destructor
      virtual ~Generic() {};

      //! Clone function
      virtual Generic* clone() const = 0;

      //! Assignment operator
      virtual Generic& operator = (const Generic& source) = 0;

      //! Reset parameters
      virtual NOX::Abstract::Group::ReturnType 
      reset(NOX::Parameter::List& params) = 0;

      //! Set flags indicating if any blocks are zero
      virtual void setIsZero(bool flagA, bool flagB, bool flagC, bool flagF,
			     bool flagG) = 0;

      //! Set flag indicating whether F and A are continguous
      virtual void setIsContiguous(bool flag) = 0;

      //! Set blocks
      /*!
       * Any of these block pointers (execpt the group pointer) may be NULL 
       * if the corresponding flag is set indicating the block is zero.
       */
      virtual void setMatrixBlocks(
		   const NOX::Abstract::Group* group,
		   const NOX::Abstract::MultiVector* blockA,
		   const NOX::Abstract::MultiVector* blockB,
		   const NOX::Abstract::MultiVector::DenseMatrix* blockC) = 0;
  
      /*! 
       * \brief Computed extended matrix-multivector product
       */
      /*!
       * Computes
       * \f[
       *     \begin{bmatrix}
       *        U \\
       *        V
       *     \end{bmatrix} = 
       *     \begin{bmatrix}
       *          J & A    \\
       *        B^T & C
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *        X \\
       *        Y 
       *     \end{bmatrix}
       * \f]
       * where \f$U\f$ is \f$n\times p\f$, \f$V\f$ is \f$m\times p\f$ and the
       * other blocks are as defined above.
       */
      virtual NOX::Abstract::Group::ReturnType 
      apply(const NOX::Abstract::MultiVector& X,
	    const NOX::Abstract::MultiVector::DenseMatrix& Y,
	    NOX::Abstract::MultiVector& U,
	    NOX::Abstract::MultiVector::DenseMatrix& V) const = 0;

      /*! 
       * \brief Computed extended matrix transpose-multivector product
       */
      /*!
       * Computes
       * \f[
       *     \begin{bmatrix}
       *        U \\
       *        V
       *     \end{bmatrix} = 
       *     \begin{bmatrix}
       *        J^T & B    \\
       *        A^T & C^T
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *        X \\
       *        Y 
       *     \end{bmatrix}
       * \f]
       * where \f$U\f$ is \f$n\times p\f$, \f$V\f$ is \f$m\times p\f$ and the
       * other blocks are as defined above.
       */
      virtual NOX::Abstract::Group::ReturnType 
      applyTranspose(const NOX::Abstract::MultiVector& X,
		     const NOX::Abstract::MultiVector::DenseMatrix& Y,
		     NOX::Abstract::MultiVector& U,
		     NOX::Abstract::MultiVector::DenseMatrix& V) const = 0;

      /*! 
       * \brief Solves the extended system as defined above
       */
      /*!
       * The \em params argument is the linear solver parameters. 
       */
      virtual NOX::Abstract::Group::ReturnType 
      applyInverse(NOX::Parameter::List& params,
		   const NOX::Abstract::MultiVector* F,
		   const NOX::Abstract::MultiVector::DenseMatrix* G,
		   NOX::Abstract::MultiVector& X,
		   NOX::Abstract::MultiVector::DenseMatrix& Y) const = 0;

    };
  } // namespace BorderedSystem
} // namespace LOCA

#endif
