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

#ifndef LOCA_BORDEREDSOLVER_COMPLEXOPERATOR_H
#define LOCA_BORDEREDSOLVER_COMPLEXOPERATOR_H

#include "Teuchos_RefCountPtr.hpp"

#include "LOCA_BorderedSolver_AbstractOperator.H" // base class

// forward declarations
namespace LOCA {
  namespace Hopf {
    namespace MooreSpence {
      class AbstractGroup;
    }
  }
}

namespace LOCA {

  namespace BorderedSolver {
    
    /*! 
     * \brief Bordered solver operator representing the \f$J + i\omega M\f$ as 
     * implemented in the LOCA::Hopf::MooreSpence::AbstractGroup.
     */
    class ComplexOperator : public LOCA::BorderedSolver::AbstractOperator {

    public:

      //! Constructor
      ComplexOperator(const Teuchos::RefCountPtr<const LOCA::Hopf::MooreSpence::AbstractGroup>& grp,
		      double omega);

      //! Destructor
      virtual ~ComplexOperator();

      //! Get group pointer
      virtual Teuchos::RefCountPtr<const NOX::Abstract::Group>
      getGroup() const;

      //! Get frequency
      virtual double getFrequency() const;
  
      //! Apply the operator
      virtual NOX::Abstract::Group::ReturnType 
      apply(const NOX::Abstract::MultiVector& X,
	    NOX::Abstract::MultiVector& Y) const;

      //! Apply transpose of the operator
      /*!
       * Group must be of type LOCA::Hopf::MinimallyAugmented::AbstractGroup
       * for this method to be implemented.
       */
      virtual NOX::Abstract::Group::ReturnType 
      applyTranspose(const NOX::Abstract::MultiVector& X,
		     NOX::Abstract::MultiVector& Y) const;

      //! Apply inverse of the operator
      virtual NOX::Abstract::Group::ReturnType 
      applyInverse(Teuchos::ParameterList& params,
		   const NOX::Abstract::MultiVector& B,
		   NOX::Abstract::MultiVector& X) const;

      //! Apply inverse transpose of the operator
      /*!
       * Group must be of type LOCA::Hopf::MinimallyAugmented::AbstractGroup
       * for this method to be implemented.
       */
      virtual NOX::Abstract::Group::ReturnType 
      applyInverseTranspose(Teuchos::ParameterList& params,
			    const NOX::Abstract::MultiVector& B,
			    NOX::Abstract::MultiVector& X) const;

    private:

      //! Private to prohibit copying
      ComplexOperator(const ComplexOperator&);

      //! Private to prohibit copying
      ComplexOperator& operator = (const ComplexOperator&);

    protected:

      //! Pointer to group
      Teuchos::RefCountPtr<const LOCA::Hopf::MooreSpence::AbstractGroup> grpPtr;

      //! Complex frequency
      double omega;

    }; // Class ComplexOperator

  } // Namespace BorderedSolver

} // Namespace LOCA

#endif // LOCA_BORDEREDSOLVER_COMPLEXOPERATOR_H
