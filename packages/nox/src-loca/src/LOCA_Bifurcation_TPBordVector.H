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

#ifndef LOCA_BIFURCATION_TPBORDVECTOR_H
#define LOCA_BIFURCATION_TPBORDVECTOR_H

#include "LOCA_MultiVector.H"  // Base class

namespace LOCA {

  namespace Bifurcation {

    /*! 
      \brief %Vector class to hold solution vectors, %Newton vectors, etc.
      for turning point equations.
    */
    /*! 
       This class uses the LOCA::MultiVector implementation to store the
       solution, null, and parameter components of the turning point vector
       and merely provides an interface for naming which components of the
       multivector these quantities correspond to.
    */

    class TPBordVector : public LOCA::MultiVector {

    public:
  
      //! %Constructor.
      TPBordVector(const NOX::Abstract::Vector& xVec,
		   const NOX::Abstract::Vector& nullVec,
		   double bifParam);

      //! Copy constructor.
      TPBordVector(const TPBordVector& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~TPBordVector();

      //! Assignment operator
      virtual NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector& y);

      //! Assignment operator
      virtual LOCA::MultiVector& operator=(const LOCA::MultiVector& y);
  
      //! Assignment operator
      virtual TPBordVector& operator=(const TPBordVector& y);

      //! Cloning function
      virtual NOX::Abstract::Vector* clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Sets the %TPBordVector by setting its three components 
      virtual void setVec(const NOX::Abstract::Vector& xVec,
			  const NOX::Abstract::Vector& nullVec,
			  double bifPar);

      //! Returns the solution vector component of extended vector
      virtual const NOX::Abstract::Vector& getXVec() const;

      //! Returns the null vector component of extended vector
      virtual const NOX::Abstract::Vector& getNullVec() const;

      //! Get %Bifurcation parameter
      virtual double getBifParam() const;

      //! Returns the solution vector component of extended vector
      virtual NOX::Abstract::Vector& getXVec();

      //! Returns the null vector component of extended vector
      virtual NOX::Abstract::Vector& getNullVec();

      //! Get %Bifurcation parameter
      virtual double& getBifParam();

    }; // class Vector
  } // namespace Bifurcation
} // namespace LOCA

#endif
