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

/*! \file NOX_LAPACK_Vector.H - Header file for NOX::LAPACK::Vector. */

#ifndef NOX_LAPACK_VECTOR_H
#define NOX_LAPACK_VECTOR_H

#include "NOX_Abstract_Vector.H" // base class
#include "NOX_Common.H" // for #include<vector>

namespace NOX {

  //! NOX BLAS/LAPACK support
  namespace LAPACK {

  //! 1.0
  const double d_one = 1.0;
  //! -1.0
  const double d_mone = -1.0;
  //! 0.0
  const double d_zero = 0.0;
  //! 1
  const int i_one = 1;
  //! 0
  const int i_zero = 0;
  
    //! Implementation of NOX::Abstract::Vector for STL vector<double> (using LAPACK for some computations)
    class Vector : public NOX::Abstract::Vector {

    public:	

      //! Construct an empty vector
      Vector();

      //! Construct a zero vector of length n
      Vector(int n);

      //! Construct a vector of length n from given array
      Vector(int n, double *v);
    
      //! Copy constructor
      Vector(const NOX::LAPACK::Vector& source, CopyType type = DeepCopy);
      
      //! Destruct Vector.
      ~Vector();

      //@{ \name Initialization methods.

      NOX::Abstract::Vector& init(double value);

      //! Initialize every element of this vector with random values
      virtual NOX::Abstract::Vector& random(bool useSeed = false, int seed =  1);

      //! Replace this vector with STL double vector y
      NOX::Abstract::Vector& operator=(const vector<double>& y);

      // derived
      NOX::Abstract::Vector& operator=(const NOX::LAPACK::Vector& y);
      NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector& y);
  
      // derived
      NOX::Abstract::Vector& abs(const NOX::LAPACK::Vector& y);
      NOX::Abstract::Vector& abs(const NOX::Abstract::Vector& y);
  
      // derived
      NOX::Abstract::Vector& reciprocal(const NOX::LAPACK::Vector& y);
      NOX::Abstract::Vector& reciprocal(const NOX::Abstract::Vector& y);
  
      //@}
  
      //@{ \name Update methods.
  
      // derived
      NOX::Abstract::Vector& scale(double gamma);
  
      // derived
      NOX::Abstract::Vector& scale(const NOX::LAPACK::Vector& a);
      NOX::Abstract::Vector& scale(const NOX::Abstract::Vector& a);
  
      // derived
      NOX::Abstract::Vector& update(double alpha, const NOX::LAPACK::Vector& a, double gamma = 0.0);
      NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a, double gamma = 0.0);
  
      // derived
      NOX::Abstract::Vector& update(double alpha, const NOX::LAPACK::Vector& a, 
			       double beta, const NOX::LAPACK::Vector& b,
			       double gamma = 0.0);
      NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a, 
			       double beta, const NOX::Abstract::Vector& b,
			       double gamma = 0.0);
  
      //@}
  
      //@{ \name Creating new Vectors. 
  
      // derived
      NOX::Abstract::Vector* clone(NOX::CopyType type = NOX::DeepCopy) const;
  
      //@}
  
      //@{ \name Norms.
  
      // derived
      double norm(NOX::Abstract::Vector::NormType type = NOX::Abstract::Vector::TwoNorm) const;
  
      // derived
      double norm(const NOX::LAPACK::Vector& weights) const;
      double norm(const NOX::Abstract::Vector& weights) const;
  
      //@}
  
      //@{ \name Dot products
  
      // derived
      double dot(const NOX::LAPACK::Vector& y) const;
      double dot(const NOX::Abstract::Vector& y) const;
  
      //@}
  
      // derived
      int length() const;
  
      //! Return the i-th element
      double& operator() (int i);

      //! Return the i-th element (const version)
      const double& operator() (int i) const;

      //! Prints out the vector to the specified stream. 
      /*! 
	For example, a vector would appear as
	\f[ \left[ \; 0.1 \; 2.34 \; 5 \; \right] \f] 
	It will be all on one line, with a single space between each entry, bracketed on either side.
      */
      ostream& leftshift(ostream& stream) const;

      // derived
      void print() const;

    private:

      //! Return the i-th element
      double& operator[] (int i);

      //! Return the i-th element
      const double& operator[] (int i) const;

      //! The used length of vector
      int n;

      //! The vector owned by this object
      vector<double> x;
    };

  } // namespace LAPACK
} // namespace NOX

//! Function for printing
ostream& operator<<(ostream& stream, const NOX::LAPACK::Vector& v);



#endif
