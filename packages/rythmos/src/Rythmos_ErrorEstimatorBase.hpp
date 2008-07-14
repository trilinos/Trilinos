//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_ERROR_ESTIMATOR_BASE_H
#define Rythmos_ERROR_ESTIMATOR_BASE_H

#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Thyra_ModelEvaluator.hpp"

#include "Rythmos_ErrorEstimateBase.hpp"

namespace Rythmos {

enum ERROR_QUANTITY_OF_INTEREST { INVALID_ERROR_QTY, AVERAGE_ERROR_QTY, ENDPOINT_ERROR_QTY };

template<class Scalar> 
class ErrorEstimatorBase
  : virtual public Teuchos::Describable
  , virtual public Teuchos::ParameterListAcceptor
  , virtual public Teuchos::VerboseObject<ErrorEstimatorBase<Scalar> >
{
public:

  // setModel
  virtual void setModel( const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model ) = 0;
  
  // setQuantityOfInterest
  virtual void setQuantityOfInterest( ERROR_QUANTITY_OF_INTEREST qtyOfInterest ) = 0;

  // getErrorEstimate
  virtual Teuchos::RCP<const ErrorEstimateBase<Scalar> > getErrorEstimate() = 0;

};


} // namespace Rythmos


#endif // Rythmos_ERROR_ESTIMATOR_BASE_H
