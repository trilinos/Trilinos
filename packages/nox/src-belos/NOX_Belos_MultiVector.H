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

#ifndef NOX_BELOS_MULTIVEC_H
#define NOX_BELOS_MULTIVEC_H

#include "BelosMultiVec.hpp"                  // base class
#include "NOX_Abstract_MultiVector.H"         // data member
#include "NOX_Abstract_Group.H"               // for ReturnType

namespace NOX { 

  namespace Belos { 

    //! Implementation of Belos::MultiVec using NOX::Abstract::MultiVector.
    class MultiVector : public ::Belos::MultiVec<double> {

      public:

      //! Constructor.
      MultiVector(NOX::Abstract::MultiVector& noxMultiVec);

      //! Copy Constructor.
      MultiVector(const NOX::Belos::MultiVector& source);

      //! Destructor.
      virtual ~MultiVector();

      //! Assignment
      virtual NOX::Belos::MultiVector&
      operator=(const NOX::Belos::MultiVector& source);

      //! Creates a new empty multi-vec containing \c numvecs columns
      virtual ::Belos::MultiVec<double>* Clone(const int numvecs);

      /*! 
       * \brief Creates a new multi-vec and copies contents of \c *this into
       * the new vector (deep copy).
       */
      virtual ::Belos::MultiVec<double>* CloneCopy();
	
      /*! 
       * \brief Creates a new multi-vec and copies the selected contents of 
       * \c *this into the new vector (deep copy).  
       */
      virtual ::Belos::MultiVec<double>* CloneCopy(int index[], int numvecs);
	
      /*! 
       * \brief Creates a new multi-vec that shares the selected contents 
       * of \c *this (shallow copy).
       */
      virtual ::Belos::MultiVec<double>* CloneView(int index[], int numvecs);


      //! Obtain the vector length of *this multivector block.
      virtual int GetVecLength() const;

      //! Obtain the number of vectors in *this multivector block.
      virtual int GetNumberVecs() const;

      /*! 
       * \brief Update \c *this with 
       * \c alpha * \c A * \c B + \c beta * (\c *this).
       */
      virtual void MvTimesMatAddMv(double alpha, 
				   ::Belos::MultiVec<double>& A, 
				   Teuchos::SerialDenseMatrix<int,double>& B,
				   double beta);

      /*! 
       * \brief Replace \c *this with \c alpha * \c A + \c beta * \c B.
       */
      virtual void MvAddMv(double alpha, 
			   ::Belos::MultiVec<double>& A, 
			   double beta, 
			   ::Belos::MultiVec<double>& B);

      /*! 
       * \brief Compute a dense matrix \c B through the matrix-matrix multiply 
       \c alpha * \c A^T * (\c *this).
      */
      virtual void MvTransMv(double alpha, 
			     ::Belos::MultiVec<double>& A, 
			     Teuchos::SerialDenseMatrix<int,double>& B);
      
      //! Compute the norm of each individual vector of \c *this.  
      virtual ::Belos::ReturnType 
      MvNorm(double *normvec, ::Belos::NormType norm_type = ::Belos::TwoNorm);

      //! Copy the vectors in \c A to a set of vectors in \c *this.  
      virtual void SetBlock(::Belos::MultiVec<double>& A, int index[], 
			    int numvecs);
	
      //! Replace the vectors in \c *this with random vectors.
      virtual void MvRandom();

      //! Replace each element of the vectors in \c *this with \c alpha.
      virtual void 
      MvInit(double alpha = Teuchos::ScalarTraits<double>::zero());

      //! Print the \c *this multivector.  
      virtual void MvPrint(ostream& os);

      //! Returns reference to underlying NOX::Abstract::MultiVector
      NOX::Abstract::MultiVector& getNoxMultiVector();

      //! Returns const reference to underlying NOX::Abstract::MultiVector
      const NOX::Abstract::MultiVector& getNoxMultiVector() const;

    protected:

      //! Constructor that sets ownsVec
      MultiVector(NOX::Abstract::MultiVector& noxMultiVec, bool ownsVecFlag);

      /*! 
       * \brief Converts NOX::Abstract::Group::ReturnType to the corresponding
       * Belos type.
       */
      ::Belos::ReturnType 
      noxReturnTypeToBelos(NOX::Abstract::Group::ReturnType noxStatus) const;

    protected:

      //! Underlying NOX multivector
      NOX::Abstract::MultiVector* vecPtr;

      //! Flag indicating whether we own the NOX::Abstract::MultiVector
      bool ownsVec;
  
    };
  } // namespace Belos
} // namespace LOCA

#endif
