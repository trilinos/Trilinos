// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef LOCA_HOPF_COMPLEXVECTOR_H
#define LOCA_HOPF_COMPLEXVECTOR_H

#include "LOCA_Extended_Vector.H"  // Base class

// Forward declarations
namespace LOCA {
  namespace Hopf {
    class ComplexMultiVector;
  }
}

namespace LOCA {

  namespace Hopf {

    /*! 
     * \brief %Vector class to hold two vectors to represent a complex vector.
     */
    /*!
     * This is not a true complex vector.  Operations like innerProduct()
     * are not correct for complex vectors.  This class
     * exists to make some aspects of the real-equivalent formulation of
     * complex linear algebra simpler to implement.
     */
    class ComplexVector : public LOCA::Extended::Vector {

      /*! 
       * \brief Declare LOCA::Hopf::ComplexMultiVector 
       * as a friend class so it can call protected methods.
       */
      friend class ComplexMultiVector;

    public:
  
      //! Constructor.
      ComplexVector(const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		    const NOX::Abstract::Vector& realVec,
		    const NOX::Abstract::Vector& imagVec);

      //! Copy constructor.
      ComplexVector(const ComplexVector& source, 
		    NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~ComplexVector();

      //! Assignment operator
      virtual NOX::Abstract::Vector& 
      operator=(const NOX::Abstract::Vector& y);

      //! Assignment operator
      virtual LOCA::Extended::Vector& 
      operator=(const LOCA::Extended::Vector& y);
  
      //! Assignment operator
      virtual ComplexVector& 
      operator=(const ComplexVector& y);

      //! Cloning function
      virtual Teuchos::RefCountPtr<NOX::Abstract::Vector> 
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Sets the vector by setting its components 
      virtual void setVec(const NOX::Abstract::Vector& realVec,
			  const NOX::Abstract::Vector& imagVec);

      //! Returns the real component of extended vector
      virtual Teuchos::RefCountPtr<const NOX::Abstract::Vector> 
      getRealVec() const;

      //! Returns the imaginary component of extended vector
      virtual Teuchos::RefCountPtr<const NOX::Abstract::Vector> 
      getImagVec() const;

      //! Returns the real component of extended vector
      virtual Teuchos::RefCountPtr<NOX::Abstract::Vector> 
      getRealVec();

      //! Returns the imaginary component of extended vector
      virtual Teuchos::RefCountPtr<NOX::Abstract::Vector> 
      getImagVec();

    protected:

      //! Default Constructor
      ComplexVector(
		   const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data);

      //! generate an extended multi-vector
      /*!
       * Returns an empty multivector of type 
       * LOCA::Hopf::ComplexMultiVector
       */
      virtual Teuchos::RefCountPtr<LOCA::Extended::MultiVector>
      generateMultiVector(int nColumns, int nVectorRows, 
			  int nScalarRows) const;

    }; // class ComplexVector
  } // namespace Hopf
} // namespace LOCA

#endif
