// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
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

#ifndef ANASAZI_LOCA_MULTIVEC_H
#define ANASAZI_LOCA_MULTIVEC_H

// Anasazi includes
#include "AnasaziMultiVec.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziReturnType.hpp"

// NOX/LOCA includes
#include "NOX_Abstract_Vector.H"

namespace Anasazi {

  enum DataAccess {Copy, View};

  namespace LOCA {
    
    /*!
     * \brief Implementation of the Anasazi::MultiVec for 
     * NOX::Abstract::Vectors.
     */
    /*!
     * This class provides an implementation of the Anasazi::MultiVec class
     * for NOX::Abstract::Vector's, to be used in conjunction with the
     * Anasazi::LOCA::Matrix class.  In its current implementation, this class
     * holds an array of NOX::Abstract::Vector pointers, an implements the
     * multi-vector operations defined by the Anasazi::MultiVec interface
     * serially.
     */
    class MultiVec : public Anasazi::MultiVec<double> {
    public:
      
      //! Constructor
      MultiVec(const NOX::Abstract::Vector& N_vec, int NumVecs );

      //! Constructor
      MultiVec(const vector<NOX::Abstract::Vector*> N_vecPtrs, 
	       Anasazi::DataAccess type = Anasazi::Copy );

      //! Constructor
      MultiVec(const MultiVec& source, 
	       Anasazi::DataAccess type = Anasazi::Copy );

      //! Constructor
      MultiVec(Anasazi::DataAccess type, const MultiVec& source, 
	      const std::vector<int>& index); 
    
      //! Destructor
      ~MultiVec();
      
      //! Create a new MultiVec with \em n columns
      virtual Anasazi::MultiVec<double>* Clone(const int n) const;
      
      //! Create a new MultiVec as a copy of \em this
      virtual Anasazi::MultiVec<double>* CloneCopy() const;
      
      //! Create a new MultiVec as a copy of selected columns
      virtual Anasazi::MultiVec<double>* CloneCopy( const std::vector<int>& index ) const;
      
      //! Create a new MultiVec as a view of selected columns
      virtual Anasazi::MultiVec<double>* CloneView( const std::vector<int>& index );
      
      //! Returns number of vectors
      virtual int GetNumberVecs () const;

      //! Returns vector length
      virtual int GetVecLength () const;

      /*! 
       * \brief Set a block of this multivec with the multivecs specified by
       *the index.
       */
      virtual void SetBlock(const Anasazi::MultiVec<double>& A, const std::vector<int>& index);

      //! *this <- alpha * A * B + beta * (*this)
      virtual void MvTimesMatAddMv(const double alpha, const Anasazi::MultiVec<double>& A, 
				   const Teuchos::SerialDenseMatrix<int,double>& B,
				   const double beta);

      //!*this <- alpha * A + beta * B
      virtual void MvAddMv(const double alpha, const Anasazi::MultiVec<double>& A, 
			   const double beta, const Anasazi::MultiVec<double>& B);
      
      //! B <- alpha * A^T * (*this)
      virtual void MvTransMv(const double alpha, const Anasazi::MultiVec<double>& A, 
			     Teuchos::SerialDenseMatrix<int,double>& B) const;
      
      //! b[i] = \c A[i]^T*\c this[i] where \c A[i] is the i-th column of A.
      virtual void MvDot ( const Anasazi::MultiVec<double>& A, std::vector<double>* b ) const;

      //! alpha[i] = norm of i-th column of (*this)
      virtual void MvNorm(std::vector<double>* normvec) const;
      
      //! random vectors in i-th column of (*this)
      virtual void MvRandom();

      //! initializes each element of (*this) with alpha
      virtual void MvInit (const double alpha);

      //! print (*this)
      virtual void MvPrint() const;

      //! Return a reference to NOX::Abstract::Vector
      virtual NOX::Abstract::Vector& GetNOXVector(int index);

      //! Return a const reference to NOX::Abstract::Vector
      virtual const NOX::Abstract::Vector& GetNOXVector(int index) const;
  
    protected:

      //! Array of vectors
      vector< NOX::Abstract::Vector* > mvPtrs;

      //! View or Copy
      Anasazi::DataAccess CV;

    }; // class MultiVec
  } // namespace LOCA
} // namespace Anasazi

#endif 
