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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_UNITTEST_MODELS_H
#define Rythmos_UNITTEST_MODELS_H

#include "Rythmos_Types.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_ObjectBuilder.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "EpetraExt_DiagonalTransientModel.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_ModelEvaluator.hpp"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

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
Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > getWFactory() {
  RCP<Teuchos::ParameterList> paramList  = Teuchos::parameterList();
  RCP<Teuchos::ParameterList> stratPl = sublist(paramList,Stratimikos_name);
  RCP<Teuchos::ParameterList> modelPl = sublist(paramList,DiagonalTransientModel_name);
  stratPl->set("Linear Solver Type","AztecOO");
  stratPl->set("Preconditioner Type","None");
  modelPl->set("NumElements",2);
  return(getWFactory<Scalar>(paramList));
}

  
template<class Scalar>
Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > 
  getDiagonalModel( Teuchos::RCP<Teuchos::ParameterList> paramList ) {

#ifdef HAVE_MPI
  Teuchos::RCP<Epetra_Comm> epetra_comm = 
    Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD) );
#else
  Teuchos::RCP<Epetra_Comm> epetra_comm = Teuchos::rcp(new Epetra_SerialComm);
#endif
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

template<class Scalar>
Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > 
  getDiagonalModel() {

  RCP<Teuchos::ParameterList> paramList  = Teuchos::parameterList();
  RCP<Teuchos::ParameterList> stratPl = sublist(paramList,Stratimikos_name);
  RCP<Teuchos::ParameterList> modelPl = sublist(paramList,DiagonalTransientModel_name);
  stratPl->set("Linear Solver Type","AztecOO");
  stratPl->set("Preconditioner Type","None");
  modelPl->set("NumElements",2);

  return(getDiagonalModel<Scalar>(paramList));
}

//// Class for Unit Testing which builds Thyra::ModelEvaluators
//class TestModelBuilder : virtual public Teuchos::ParameterListAcceptor
//{
//  public:
//
//    TestModelBuilder();
//
//    virtual ~TestModelBuilder();
//
//    void setModelFactory(
//      const RCP<const Teuchos::AbstractFactory<Thyra::ModelEvaluator<double> > > &modelFactory,
//      const std::string &modelFactoryName
//      );
//
//    std::string getModelName() const;
//
//    RCP<Thyra::ModelEvaluator<double> > create(
//      const std::string &modelName = ""
//      ) const;
//
//    /** \name Overridden from Teuchos::ParameterListAcceptor */
//    //@{
//
//    /** \brief . */
//    void setParameterList(const RCP<Teuchos::ParameterList> & paramList);
//    
//    /** \brief . */
//    RCP<Teuchos::ParameterList> getNonconstParameterList();
//    
//    /** \brief . */
//    RCP<Teuchos::ParameterList> unsetParameterList();
//    
//    /** \brief. */
//    RCP<const ParameterList> getParameterList() const;
//
//    /** \brief. */
//    RCP<const Teuchos::ParameterList> getValidParameters() const;
//  
//    //@}
//
//  private:
//    Teuchos::ObjectBuilder<Thyra::ModelEvaluator<double> > builder_;
//    void initializeDefaults_();
//};
//
//// Nonmember constructor
//RCP<TestModelBuilder> testModelBuilder();
//// Nonmember helper
//RCP<Thyra::ModelEvaluator<double> > createTestModel(
//    const std::string& modelName,
//    const RCP<ParameterList>& modelPL = Teuchos::null
//    );


} // namespace Rythmos
#endif // Rythmos_UNITTEST_MODELS_H

