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

/*! \file NOX_TSF_Vector.H - Header file for NOX::TSF::Vector. */

#ifndef NOX_TSF_VECTOR_H
#define NOX_TSF_VECTOR_H

#include "NOX_Abstract_Vector.H" // base class
#include "NOX_Common.H" // for #include<vector> 
#include "TSFVector.hpp" // definitions for TSFExtended Vectors
#include "TSFVectorSpace.hpp"

namespace NOX 
{

  //! NOX BLAS/TSF support
  namespace TSF 
  {

  
    /** Implementation of the NOX::Abstract::Vector interface for 
     * TSFExtended::Vector<double> (using TSFExtended handle to TSFCore).
     *
     * @author Jill Reese and Kevin Long
     */
    class Vector : public NOX::Abstract::Vector 
    {

    public:	

      //! Copy constructor 
      Vector(const NOX::TSF::Vector& source, NOX::CopyType type = DeepCopy);

      //! Copy constructor
      Vector(const TSFExtended::Vector<double>& source, NOX::CopyType type = DeepCopy);

      /* There's no need to write an explicit dtor for this class, because
       * the TSFExtended vector cleans up its own memory. I've commented out 
       * this method for now, but it should probably just be deleted. 
       * - KL 12 July 2004 */ 
      ////! Destruct Vector.
      //~Vector();

      /** Access to a writable reference to the underlying TSF vector*/
      TSFExtended::Vector<double>& NOX::TSF::Vector::getTSFVector();

      /** Access to a read-only reference to the underlying TSF vector*/
      const TSFExtended::Vector<double>& NOX::TSF::Vector::getTSFVector() const;
      //@{ \name Initialization methods.

      /** Initialize all elements to a constant */
      NOX::Abstract::Vector& init(double value);

      /** Assignment operator */
      NOX::Abstract::Vector& operator=(const NOX::TSF::Vector& y);
      /** Assignment operator */
      NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector& y);

      /* Element access should not be implemented in this class; it should
       * be done through the LoadableVector and AccessibleVector 
       * interfaces of TSFExtended.
       * -KL 12 July 2004  */
      // Store value in the i-th element of a vector 
      //void setElement(int i, const double& value);

      // Return the i-th element (const version).
      //const double& getElement(int i) const;
  
      /** Element-wise absolute value */
      NOX::Abstract::Vector& abs(const NOX::TSF::Vector& y);
      /** Element-wise absolute value */
      NOX::Abstract::Vector& abs(const NOX::Abstract::Vector& y);
  
      /** Element-wise reciprocal */
      NOX::Abstract::Vector& reciprocal(const NOX::TSF::Vector& y);
      /** Element-wise reciprocal */
      NOX::Abstract::Vector& reciprocal(const NOX::Abstract::Vector& y);
  
      //@}
  
      //@{ \name Update methods.
  
      // derived
      NOX::Abstract::Vector& scale(double gamma);
  
      // derived
      NOX::Abstract::Vector& scale(const NOX::TSF::Vector& a);
      NOX::Abstract::Vector& scale(const NOX::Abstract::Vector& a);
  
      // derived
      NOX::Abstract::Vector& update(double alpha, const NOX::TSF::Vector& a, double gamma = 0.0);
      NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a, double gamma = 0.0);
  
      // derived
      NOX::Abstract::Vector& update(double alpha, const NOX::TSF::Vector& a, 
                                    double beta, const NOX::TSF::Vector& b,
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
      double norm(const NOX::TSF::Vector& weights) const;
      double norm(const NOX::Abstract::Vector& weights) const;
  
      //@}
  
      //@{ \name Dot products
  
      // derived
      double dot(const NOX::TSF::Vector& y) const;
      double dot(const NOX::Abstract::Vector& y) const;
  
      //@}
  
      // derived
      //! The length of a TSFExtended vector is known only by the vector space it belongs to
      //! We query the vector space and ask its dimension  
      int length() const;

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
      /** The TSF vector owned by this object. */
      TSFExtended::Vector<double> x;
    };

  } // namespace TSF
} // namespace NOX

//! Function for printing
ostream& operator<<(ostream& stream, const NOX::TSF::Vector& v);

#endif
