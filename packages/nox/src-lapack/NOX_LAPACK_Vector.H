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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef NOX_LAPACK_VECTOR_H
#define NOX_LAPACK_VECTOR_H

#include "NOX_Abstract_Vector.H" // base class
#include "NOX_Common.H" // for #include<vector>
#include "NOX_Random.H" // for Random class

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
  
    //! Implementation of NOX::Abstract::Vector for vector<double> using LAPACK
    class Vector : public Abstract::Vector {

    public:			
    
      //! Construct a zero vector of length n
      Vector(int n);
    
      //! Copy constructor
      Vector(const Vector& source, CopyType type = DeepCopy);
      
      //! Destruct Vector.
      ~Vector();

      //@{ \name Initialization methods.

      Abstract::Vector& init(double value);

      //! Initialize every element of this vector with random values
      virtual NOX::Abstract::Vector& random(bool useSeed = false, 
					    double seed = 0.0);

      //! Copies y into "this".
      Abstract::Vector& operator=(const vector<double>& y);

      Abstract::Vector& operator=(const Vector& y);

      //! See above.
      Abstract::Vector& operator=(const Abstract::Vector& y);
  
      Abstract::Vector& abs(const Vector& y);

      //! See above.
      Abstract::Vector& abs(const Abstract::Vector& y);
  
      Abstract::Vector& reciprocal(const Vector& y);

      //! See above.
      Abstract::Vector& reciprocal(const Abstract::Vector& y);
  
      //@}
  
      //@{ \name Update methods.
  
      Abstract::Vector& scale(double gamma);
  
      Abstract::Vector& scale(const Vector& a);
  
      //! See above
      Abstract::Vector& scale(const Abstract::Vector& a);
  
      Abstract::Vector& update(double alpha, const Vector& a, double gamma = 0.0);
      //! See above.
      Abstract::Vector& update(double alpha, const Abstract::Vector& a, double gamma = 0.0);
  
      Abstract::Vector& update(double alpha, const Vector& a, 
			       double beta, const Vector& b,
			       double gamma = 0.0);
      //! See above.
      Abstract::Vector& update(double alpha, const Abstract::Vector& a, 
			       double beta, const Abstract::Vector& b,
			       double gamma = 0.0);
  
      //@}
  
      //@{ \name Creating new Vectors. 
  
      Abstract::Vector* clone(CopyType type = DeepCopy) const;
  
      //@}
  
      //@{ \name Norms.
  
      double norm(NormType type = TwoNorm) const;
  
      double norm(const Vector& weights) const;

      //! See above.
      double norm(const Abstract::Vector& weights) const;
  
      //@}
  
      //@{ \name Dot products
  
      double dot(const Vector& y) const;

      //! See above.
      double dot(const Abstract::Vector& y) const;
  
      //@}
  
      int length() const;
  
      //! Return the i-th element
      double& operator() (int i);

      //! Return the i-th element (const version)
      const double& operator() (int i) const;

      //! Prints out the vector to the specified stream. 
      /*! For example, a vector would appear as
	\f[ \left[ \; 0.1 \; 2.34 \; 5 \; \right] \f] 
	It will be all on one line, with a single space between each entry, bracketed on either side.
      */
      ostream& leftshift(ostream& stream) const;

      //! Prints out the vector to the cout stream. 
      void print() const;

    private:

      //! Return the i-th element
      double& operator[] (int i);

      //! Return the i-th element
      const double& operator[] (int i) const;

      //! Pointer to vector owned by this object
      vector<double> x;

      //! Random number generator
      NOX::Random urand;

    };

  } // namespace LAPACK
} // namespace NOX

ostream& operator<<(ostream& stream, const NOX::LAPACK::Vector& v);



#endif
