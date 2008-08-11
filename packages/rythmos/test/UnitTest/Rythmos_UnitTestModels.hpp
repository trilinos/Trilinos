//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
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

#ifndef Rythmos_UNITTEST_MODELS_H
#define Rythmos_UNITTEST_MODELS_H

#include "Teuchos_DefaultComm.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "EpetraExt_DiagonalTransientModel.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Epetra_SerialComm.h"

namespace Rythmos {

const std::string Stratimikos_name = "Stratimikos";
const std::string DiagonalTransientModel_name = "DiagonalTransientModel";

template<class Scalar>
Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > getWFactory(Teuchos::RCP<Teuchos::ParameterList> paramList) {
  Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
  linearSolverBuilder.setParameterList(sublist(paramList,Stratimikos_name));
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
    W_factory = createLinearSolveStrategy(linearSolverBuilder);
  return(W_factory);
}
  
template<class Scalar>
Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > 
  getDiagonalModel( Teuchos::RCP<Teuchos::ParameterList> paramList ) {

  Teuchos::RCP<Epetra_Comm> epetra_comm = Teuchos::rcp(new Epetra_SerialComm);
  Teuchos::RCP<EpetraExt::DiagonalTransientModel>
    epetraDiagonalModel = EpetraExt::diagonalTransientModel(
        epetra_comm, 
        sublist(paramList,DiagonalTransientModel_name)
        );

  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
    W_factory = getWFactory<Scalar>(paramList);

  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >
    thyraDiagonalModel = Thyra::epetraModelEvaluator(epetraDiagonalModel,W_factory);

  return(thyraDiagonalModel);
}


} // namespace Rythmos
#endif // Rythmos_UNITTEST_MODELS_H

