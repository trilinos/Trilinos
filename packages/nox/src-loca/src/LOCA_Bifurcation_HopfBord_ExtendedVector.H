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

#ifndef LOCA_BIFURCATION_HOPFBORD_EXTENDEDVECTOR_H
#define LOCA_BIFURCATION_HOPFBORD_EXTENDEDVECTOR_H

#include "LOCA_Extended_Vector.H"  // Base class

namespace LOCA {

  namespace Bifurcation {

    namespace HopfBord {

      /*! 
	\brief %Vector class to hold solution vectors, %Newton vectors, etc.
	for Hopf equations.
      */
      /*! 
	This class uses the LOCA::Extended::Vector implementation to store the
	solution, real and imaginary eigenvector, frequency and parameter 
	components of the Hopf vector and merely provides an interface for 
	naming which components of the multivector these quantities 
	correspond to.
      */

      class ExtendedVector : public LOCA::Extended::Vector {

      public:
  
	//! %Constructor.
	ExtendedVector(const NOX::Abstract::Vector& xVec,
		       const NOX::Abstract::Vector& realEigenVec,
		       const NOX::Abstract::Vector& imagEigenVec,
		       double frequency, double bifParam);

	//! Copy constructor.
	ExtendedVector(const ExtendedVector& source, 
		       NOX::CopyType type = NOX::DeepCopy);

	//! Destructor.
	virtual ~ExtendedVector();

	//! Assignment operator
	virtual NOX::Abstract::Vector& 
	operator=(const NOX::Abstract::Vector& y);

	//! Assignment operator
	virtual LOCA::Extended::Vector& 
	operator=(const LOCA::Extended::Vector& y);
  
	//! Assignment operator
	virtual ExtendedVector& 
	operator=(const ExtendedVector& y);

	//! Cloning function
	virtual NOX::Abstract::Vector* 
	clone(NOX::CopyType type = NOX::DeepCopy) const;

	//! Sets the %HopfBordVector by setting its five components 
	virtual void setVec(const NOX::Abstract::Vector& xVec,
			    const NOX::Abstract::Vector& realEigenVec,
			    const NOX::Abstract::Vector& imagEigenVec,
			    double frequency,
			    double bifPar);

	//! Returns the solution vector component of extended vector
	virtual const NOX::Abstract::Vector& getXVec() const;

	//! Returns the real part of the eigenvector component of extended vector
	virtual const NOX::Abstract::Vector& getRealEigenVec() const;

	//! Returns the imaginary part of the eigenvector component of extended vector
	virtual const NOX::Abstract::Vector& getImagEigenVec() const;

	//! Returns the frequency component of the extended vector
	virtual double getFrequency() const;

	//! Get %Bifurcation parameter
	virtual double getBifParam() const;

	//! Returns the solution vector component of extended vector
	virtual NOX::Abstract::Vector& getXVec();

	//! Returns the real part of the eigenvector component of extended vector
	virtual NOX::Abstract::Vector& getRealEigenVec();

	//! Returns the imaginary part of the eigenvector component of extended vector
	virtual NOX::Abstract::Vector& getImagEigenVec();

	//! Returns the frequency component of the extended vector
	virtual double& getFrequency();

	//! Get %Bifurcation parameter
	virtual double& getBifParam();

      }; // class Vector
    } // namespace HopfBord
  } // namespace Bifurcation
} // namespace LOCA

#endif
