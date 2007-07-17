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

#ifndef LOCA_THYRA_GROUP_H
#define LOCA_THYRA_GROUP_H

#include "LOCA_Abstract_Group.H"               // base class
#include "NOX_Thyra_Group.H"                   // base class
#include "LOCA_Parameter_Vector.H"             // class data element
#include "Thyra_VectorBase.hpp"                // class data element

// Forward declares
namespace Teuchos {
  class ParameterList;
}

namespace LOCA {

  //! %LOCA %Thyra support
  namespace Thyra {

    //! Extension of the NOX::Thyra::Group to %LOCA.  
    class Group : public virtual NOX::Thyra::Group, 
		  public virtual LOCA::Abstract::Group {

    public:
      
      //! Constructor
      Group(const Teuchos::RCP<LOCA::GlobalData>& global_data,
	    const NOX::Thyra::Vector& initial_guess,
	    const Teuchos::RCP< ::Thyra::ModelEvaluator<double> >& model,
	    const LOCA::ParameterVector& p,
	    int p_index);

      /*! 
       * \brief Copy constructor. If type is DeepCopy, takes ownership of
       * valid shared Jacobian. 
       */
      Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~Group();

      //! Assignment operator
      virtual Group& operator=(const Group& source);

      /*! 
       * @name Overloaded NOX::Thyra::Group  methods.
       */
      //@{

      //! Assignment operator.
      virtual NOX::Abstract::Group& 
      operator=(const NOX::Abstract::Group& source);

      //! Assignment operator.
      virtual NOX::Abstract::Group& 
      operator=(const NOX::Thyra::Group& source);

      //! Cloning function
      virtual Teuchos::RCP<NOX::Abstract::Group>
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Overloaded computeF()
      virtual NOX::Abstract::Group::ReturnType 
      computeF();

      //! Overloaded computeJacobian()
      virtual NOX::Abstract::Group::ReturnType 
      computeJacobian();

      //@}

      /*! 
       * @name Implementation of LOCA::MultiContinuation::AbstractGroup virtual methods.
       */
      //@{

      //! Copy
      virtual void copy(const NOX::Abstract::Group& source);

      //! Set the parameters
      virtual void setParams(const ParameterVector& p);

      //! Set parameter indexed by paramID
      virtual void setParam(int paramID, double val);
  
      //! Set parameter indexed by paramID
      virtual void setParam(string paramID, double val);

      //! Return a const reference to the ParameterVector owned by the group. 
      const LOCA::ParameterVector& getParams() const;

      //! Return copy of parameter indexed by paramID
      virtual double getParam(int paramID) const;

      //! Return copy of parameter indexed by paramID
      virtual double getParam(string paramID) const;

      //! Call the user interface print() routine, solution vector
      virtual void printSolution(const double conParam) const;

      //! Call the user interface print() routine, any vector
      virtual void printSolution(const NOX::Abstract::Vector& x,
                                 const double conParam) const;

      //@}

      /*! 
       * @name Implementation of LOCA::TimeDependent::AbstractGroup virtual methods.
       */
      //@{

      //! Compute the shifted matrix
      virtual NOX::Abstract::Group::ReturnType
      computeShiftedMatrix(double alpha, double beta);

      //! Multiply the shifted matrix by a vector.
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrix(const NOX::Abstract::Vector& input,
			 NOX::Abstract::Vector& result) const;

      //! Multiply the shifted matrix by a multi-vector.
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrixMultiVector(
				const NOX::Abstract::MultiVector& input,
				NOX::Abstract::MultiVector& result) const;
      
      /*! 
       * \brief Apply the inverse of the shifted matrix by a multi-vector, as 
       * needed by the shift-and-invert and generalized Cayley transformations.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrixInverseMultiVector(
			        Teuchos::ParameterList& params, 
				const NOX::Abstract::MultiVector& input,
				NOX::Abstract::MultiVector& result) const;
      //@}

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Parameter vector
      LOCA::ParameterVector params;

      //! Parameter index in Thyra::ModelEvaluator
      int param_index;

      //! Pointer to Thyra vector storing parameters
      Teuchos::RCP< ::Thyra::VectorBase<double> > param_thyra_vec;

      //! Pointer to x_dot vector of zeros
      Teuchos::RCP< ::Thyra::VectorBase<double> > x_dot_vec;

    };

  } // namespace Thyra
} // namespace LOCA


#endif
