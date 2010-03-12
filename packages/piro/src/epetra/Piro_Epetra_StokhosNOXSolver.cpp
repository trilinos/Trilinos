// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include "Piro_Epetra_StokhosNOXSolver.hpp"
#include "Piro_Epetra_NOXSolver.hpp"
#include "Piro_ValidPiroParameters.hpp"
#include "Stokhos.hpp"
#include "Stokhos_SGModelEvaluator.hpp"
#include "Stokhos_SGQuadModelEvaluator.hpp"
#include "Sacado_PCE_OrthogPoly.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Ifpack.h"

// The preconditioner we will use for PCE
class IfpackPreconditionerFactory : public Stokhos::PreconditionerFactory {
public:
  IfpackPreconditionerFactory(const Teuchos::RCP<Teuchos::ParameterList>& p) :
    precParams(p) {}
  virtual ~IfpackPreconditionerFactory() {}
  virtual Teuchos::RCP<Epetra_Operator> 
  compute(const Teuchos::RCP<Epetra_Operator>& op) {
    Teuchos::RCP<Epetra_RowMatrix> mat = 
      Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(op, true);
    Ifpack Factory;
    std::string prec = precParams->get("Ifpack Preconditioner", "ILU");
    int overlap = precParams->get("Overlap", 0);
    ifpackPrec = Teuchos::rcp(Factory.Create(prec, mat.get(), overlap));
    ifpackPrec->SetParameters(*precParams);
    int err = ifpackPrec->Initialize();   
    err = ifpackPrec->Compute();
    return ifpackPrec;
  }
protected:
  Teuchos::RCP<Teuchos::ParameterList> precParams;
  Teuchos::RCP<Ifpack_Preconditioner> ifpackPrec;
};

Piro::Epetra::StokhosNOXSolver::
StokhosNOXSolver(const Teuchos::RCP<Teuchos::ParameterList>& appParams,
	    const Teuchos::RCP<EpetraExt::ModelEvaluator>& model,
//	    const Teuchos::RCP<FEApp::Application>& app,
	    const Teuchos::RCP<const Epetra_Comm>& comm,
            Teuchos::RCP<Piro::Epetra::NOXObserver> noxObserver)
{
  //appParams->validateParameters(*Piro::getValidPiroParameters(),0);

  Teuchos::ParameterList& problemParams = appParams->sublist("Problem");
  Teuchos::ParameterList& sgParams =
    problemParams.sublist("Stochastic Galerkin");
  sgParams.validateParameters(*getValidSGParameters(),0);

  // Get SG expansion type
  std::string sg_type = sgParams.get("SG Method", "AD");
  SG_METHOD sg_method;
  if (sg_type == "AD")
    sg_method = SG_AD;
  else if (sg_type == "Global")
    sg_method = SG_GLOBAL;
  else if (sg_type == "Non-intrusive")
    sg_method = SG_NI;
  else
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
		       std::endl << "Error!  Piro::Epetra::StokhosNOXSolver():  " <<
		       "Invalid SG Method  " << sg_type << std::endl);

  // Create SG basis
  Teuchos::ParameterList& sg_parameterParams =
    sgParams.sublist("SG Parameters");
  int numParameters = sg_parameterParams.get("Number", 0);
  Teuchos::ParameterList& basisParams = sgParams.sublist("Basis");
  Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(numParameters);
  for (int i=0; i<numParameters; i++) {
    std::ostringstream ss;
    ss << "Basis " << i;
    Teuchos::ParameterList& bp = basisParams.sublist(ss.str());
    std::string type = bp.get("Type","Legendre");
    int order = bp.get("Order", 3);
    if (type == "Legendre")
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<int,double>(order));
    else if (type == "Hermite")
      bases[i] = Teuchos::rcp(new Stokhos::HermiteBasis<int,double>(order));
    else
      TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
			 std::endl << "Error!  Piro::Epetra_StokhosNOXSolver():  " <<
			 "Invalid basis type  " << type << std::endl);

    
  }
  basis = 
    Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
  std::cout << "Basis size = " << basis->size() << std::endl;

  // Set up stochastic parameters
  Teuchos::Array< Teuchos::Array< Teuchos::RCP<Epetra_Vector> > > sg_p(1);
  unsigned int sz = basis->size();
  sg_p[0].resize(sz);
  Epetra_LocalMap p_sg_map(numParameters, 0, *comm);
  for (unsigned int i=0; i<sz; i++)
    sg_p[0][i] = Teuchos::rcp(new Epetra_Vector(p_sg_map));
  for (int i=0; i<numParameters; i++) {
    std::ostringstream ss;
    ss << "Basis " << i;
    Teuchos::ParameterList& bp = basisParams.sublist(ss.str());
    Teuchos::Array<double> initial_p_vals = 
      Teuchos::getArrayFromStringParameter<double>(bp, 
						   std::string("Initial Expansion Coefficients"), 
						   -1, false);
    if (initial_p_vals.size() == 0)
      (*(sg_p[0][i+1]))[i] = 1.0;  // Set order 1 coeff to 1 for this RV
    else
      for (unsigned int j = 0; j<initial_p_vals.size(); j++)
	(*(sg_p[0][j]))[i] = initial_p_vals[j];
  }

  // SG Quadrature
  std::string quad_type = sgParams.get("Quadrature Type", "Tensor Product");
  if (quad_type == "Tensor Product")
    quad = 
      Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
  else if (quad_type == "Sparse Grid") {
#ifdef HAVE_DAKOTA
    if (sgParams.isType<int>("Sparse Grid Level")) {
      int level = sgParams.get<int>("Sparse Grid Level");
      quad = 
	Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(basis,
								   level));
    }
    else
      quad = 
	Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(basis));
#else
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
		       std::endl << "Error!  Piro::Epetra_StokhosNOXSolver():  " <<
		       "SparseGrid quadrature requires Dakota" << std::endl);
#endif
  }
  else
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
		       std::endl << "Error!  Piro::Epetra_StokhosNOXSolver():  " <<
		       "Invalid quadrature type  " << quad_type << std::endl);

  // SG AD Expansion
  std::string exp_type = sgParams.get("AD Expansion Type", "Quadrature");
  Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion;
  if (exp_type == "Quadrature")
    expansion = 
      Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<int,double>(basis, 
								    quad));
#ifdef HAVE_STOKHOS_FORUQTK
  else if (exp_type == "For UQTK") {
    if (sgParams.isType<double>("Taylor Expansion Tolerance")) {
      double rtol = sgParams.get<double>("Taylor Expansion Tolerance");
      expansion = 
	Teuchos::rcp(new Stokhos::ForUQTKOrthogPolyExpansion<int,double>(basis,
									 Stokhos::ForUQTKOrthogPolyExpansion<int,double>::TAYLOR,
									 rtol));
    }
    else
      expansion = 
	Teuchos::rcp(new Stokhos::ForUQTKOrthogPolyExpansion<int,double>(basis));
  }
#endif
  else if (exp_type == "Derivative")
    expansion = 
      Teuchos::rcp(new Stokhos::DerivOrthogPolyExpansion<int,double>(basis));
  else
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
		       std::endl << "Error!  Piro::Epetra_StokhosNOXSolver():  " <<
		       "Invalid expansion type  " << exp_type << std::endl);

  // Set up stochastic Galerkin model
  Teuchos::Array<int> sg_p_index(1);
  Teuchos::Array<int> sg_g_index(1);
  sg_p_index[0] = 1; // 2nd parameter vector are the stochastic parameters
  sg_g_index[0] = 0; // 1st response vector are the stochastic responses
  Teuchos::RCP<EpetraExt::ModelEvaluator> sg_model;
  if (sg_method == SG_AD) {
    sg_model = model;
cout << " AGS!! COMMENTED OUT CRITICAL CALL TO ALLOW StokhosNOX SOLVER TO COMPILE!!!" << endl;
 // app->init_sg(expansion);
  }
  else {
    Teuchos::RCP<EpetraExt::ModelEvaluator> underlying_model;
    if (sg_method == SG_GLOBAL)
      underlying_model = model;
    else 
      underlying_model =
	Teuchos::rcp(new Piro::Epetra::NOXSolver(appParams, model));
    sg_model =
      Teuchos::rcp(new Stokhos::SGQuadModelEvaluator(underlying_model, 
						     quad, sg_p_index,
						     sg_g_index));
  }

  // Set up SG nonlinear model
  std::string solver_method;
  Teuchos::RCP<Teuchos::ParameterList> sgSolverParams = 
    Teuchos::rcp(&(sgParams.sublist("SG Solver Parameters")),false);
  if (sg_method != SG_NI) {
    solver_method = sgSolverParams->get("Jacobian Method", "Matrix Free");
    if (solver_method == "Matrix Free" || solver_method == "Matrix Free Orig"){
      Teuchos::RCP<Teuchos::ParameterList> precParams = 
	Teuchos::rcp(&(sgSolverParams->sublist("SG Preconditioner")),false);
      Teuchos::RCP<Stokhos::PreconditionerFactory> sg_prec = 
	Teuchos::rcp(new IfpackPreconditionerFactory(precParams));
      sgSolverParams->set("Preconditioner Factory", sg_prec);
    }
  }
  Teuchos::RCP<Stokhos::SGModelEvaluator> sg_nonlin_model =
    Teuchos::rcp(new Stokhos::SGModelEvaluator(sg_model, basis, sg_p_index,
					       sg_g_index, sg_p, sgSolverParams,
					       comm));

  // Set up Observer to call noxObserver for each vector block
  Teuchos::RCP<Piro::Epetra::NOXObserver> sgnoxObserver;
  if (noxObserver != Teuchos::null)
    sgnoxObserver = Teuchos::rcp(new Piro::Epetra::StokhosNOXObserver(noxObserver,
                                                          *(model->get_x_map()),
                                                          sz));

  // Create SG NOX solver
  if (sg_method != SG_NI) {
    // Will find preconditioner for Matrix-Free method
    sg_solver = Teuchos::rcp(new Piro::Epetra::NOXSolver(appParams, sg_nonlin_model, 
						 sgnoxObserver));
  }
  else 
    sg_solver = sg_nonlin_model;
}

Piro::Epetra::StokhosNOXSolver::~StokhosNOXSolver()
{
}

Teuchos::RCP<const Epetra_Map> 
Piro::Epetra::StokhosNOXSolver::get_x_map() const
{
  return sg_solver->get_x_map();
}

Teuchos::RCP<const Epetra_Map> 
Piro::Epetra::StokhosNOXSolver::get_f_map() const
{
  return sg_solver->get_f_map();
}

Teuchos::RCP<const Epetra_Map> 
Piro::Epetra::StokhosNOXSolver::get_p_map(int l) const
{
  return sg_solver->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map> 
Piro::Epetra::StokhosNOXSolver::get_g_map(int j) const
{
  return sg_solver->get_g_map(j);
}

Teuchos::RCP<const Epetra_Vector> 
Piro::Epetra::StokhosNOXSolver::get_x_init() const
{
  return sg_solver->get_x_init();
}

Teuchos::RCP<const Epetra_Vector> 
Piro::Epetra::StokhosNOXSolver::get_p_init(int l) const
{
  return sg_solver->get_p_init(l);
}

EpetraExt::ModelEvaluator::InArgs 
Piro::Epetra::StokhosNOXSolver::createInArgs() const
{
  return sg_solver->createInArgs();
}

EpetraExt::ModelEvaluator::OutArgs 
Piro::Epetra::StokhosNOXSolver::createOutArgs() const
{
  return sg_solver->createOutArgs();
}

void 
Piro::Epetra::StokhosNOXSolver::evalModel(const InArgs& inArgs,
			     const OutArgs& outArgs ) const
{
  sg_solver->evalModel(inArgs, outArgs);
}

Teuchos::RCP<const Teuchos::ParameterList>
Piro::Epetra::StokhosNOXSolver::getValidSGParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
     rcp(new Teuchos::ParameterList("ValidSGParams"));;
  validPL->sublist("SG Parameters", false, "");
  validPL->sublist("SG Solver Parameters", false, "");
  validPL->sublist("Basis", false, "");
  validPL->set<std::string>("SG Method", "","");
  validPL->set<std::string>("SG Method", "","");
  validPL->set<std::string>("AD Expansion Type", "","");
  validPL->set<std::string>("Quadrature Type", "","");

  return validPL;
}

