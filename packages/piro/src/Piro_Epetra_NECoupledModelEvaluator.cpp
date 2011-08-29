#include "Piro_Epetra_NECoupledModelEvaluator.hpp"
#include "Piro_Epetra_Factory.hpp"
#include "Piro_Epetra_StokhosSolver.hpp"

#include "Epetra_LocalMap.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_TestForException.hpp"

#include "Teuchos_TimeMonitor.hpp"


#include "Stokhos.hpp"
#include "Stokhos_StieltjesGramSchmidtBuilder.hpp"

Piro::Epetra::NECoupledModelEvaluator::
NECoupledModelEvaluator(
  const Teuchos::RCP<EpetraExt::ModelEvaluator>& modelA_, 
  const Teuchos::RCP<EpetraExt::ModelEvaluator>& modelB_,
  const Teuchos::RCP<Teuchos::ParameterList>& piroParamsA_,
  const Teuchos::RCP<Teuchos::ParameterList>& piroParamsB_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_,
  const Teuchos::RCP<const Epetra_Comm>& comm_):
  modelA(modelA_),
  modelB(modelB_),
  piroParamsA(piroParamsA_),
  piroParamsB(piroParamsB_),
  params(params_),
  comm(comm_),

  basis(),
  DnumRandVar(dnumRandVar),
  HnumRandVar(hnumRandVar),
  reduce_dimension(do_dim_reduction),
  orthogonalize_bases(do_basis_orthog),
  eval_W_with_f(eval_W_with_f_)
{
  // Create solvers for models A and B
  bool stochastic = params->get("Stochastic", false);
  if (stochastic) {
    sgSolverA = Teuchos::rcp(new Piro::Epetra::StokhosSolver(piroParamsA, 
							     comm));
    sgSolverB = Teuchos::rcp(new Piro::Epetra::StokhosSolver(piroParamsB, 
							     comm));
    sgSolverA->setup(modelA);
    sgSolverB->setup(modelB);
    solverA = sgSolverA;
    solverB = sgSolverB;
  }
  else {
    solverA = Piro::Epetra::Factory(piroParamsA, modelA);
    solverB = Piro::Epetra::Factory(piroParamsB, modelB);
  }

  // Get connectivity information
  pIndexA = params->get("Model A Parameter Index", 0);
  pIndexB = params->get("Model B Parameter Index", 0);
  gIndexA = params->get("Model A Response Index", 0);
  gIndexB = params->get("Model B Response Index", 0);

  // Get number of parameter and response vectors
  EpetraExt::ModelEvaluator::InArgs modelA_inargs = modelA->createInArgs();
  EpetraExt::ModelEvaluator::InArgs modelB_inargs = modelB->createInArgs();
  EpetraExt::ModelEvaluator::OutArgs modelA_outargs = modelA->createOutArgs();
  EpetraExt::ModelEvaluator::OutArgs modelB_outargs = modelB->createOutArgs();
  num_params_A = modelA_inargs.Np();
  num_params_B = modelB_inargs.Np();
  num_responses_A = modelA_outargs.Ng();
  num_responses_B = modelB_outargs.Ng();

  // Building indexing maps between coupled system parameters/responses and
  // individual components
  param_map.resize(num_params_A + num_params_B - 2);
  for (int i=0; i<pIndexA; i++)
    param_map[i] = i;
  for (int i=pIndexA+1; i<num_params_A; i++)
    param_map[i-1] = i;
  for (int i=0; i<pIndexB; i++)
    param_map[num_params_A-1+i] = i;
  for (int i=pIndexB+1; i<num_params_B; i++)
    param_map[num_params_A-2+i] = i;
  response_map.resize(num_responses_A + num_responses_B - 2);
  for (int i=0; i<pIndexA; i++)
    response_map[i] = i;
  for (int i=pIndexA+1; i<num_responses_A; i++)
    response_map[i-1] = i;
  for (int i=0; i<pIndexB; i++)
    response_map[num_responses_A-1+i] = i;
  for (int i=pIndexB+1; i<num_responses_B; i++)
    response_map[num_responses_A-2+i] = i;

  //
  // The network equations look like:
  //    p_1 - g_2(x_2,p_2) = 0 s.t. f_1(x_1,p_1) = 0
  //    p_2 - g_1(x_1,p_1) = 0 s.t. f_2(x_2,p_2) = 0
  //
  // We define x = [p_1; p_2] and f = [ p_1 - g_2(x_2,p_2); p_2 - g_1(x_1,p_1)]
  //

  // Build x map, which is the product of the pIndexA and pIndexB parameter maps
  // For the time being, we will assume local maps, in the future we need to
  // build proper product maps
  Teuchos::RCP<const Epetra_Map> p_map_A = modelA->get_p_map(pIndexA);
  Teuchos::RCP<const Epetra_Map> p_map_B = modelB->get_p_map(pIndexB);
  int n_p_A = p_map_A->NumGlobalElements();
  int n_p_B = p_map_B->NumGlobalElements();
  int nx = n_p_A + n_p_B;
  x_map = Teuchos::rcp(new Epetra_Map(nx, 0, *comm));
  x_overlap_map = Teuchos::rcp(new Epetra_LocalMap(nx, 0, *comm));
  x_importer = Teuchos::rcp(new Epetra_Import(*x_overlap_map, *x_map));

  // Build f map, which is the product of the gIndexA and gIndexB response maps
  // For the time being, we will assume local maps, in the future we need to
  // build proper product maps
  Teuchos::RCP<const Epetra_Map> g_map_A = modelA->get_g_map(gIndexA);
  Teuchos::RCP<const Epetra_Map> g_map_B = modelB->get_g_map(gIndexB);
  int n_g_A = g_map_A->NumGlobalElements();
  int n_g_B = g_map_B->NumGlobalElements();
  int nf = n_g_A + n_g_B;
  f_map = Teuchos::rcp(new Epetra_Map(nf, 0, *comm));
  f_overlap_map = Teuchos::rcp(new Epetra_LocalMap(nf, 0, *comm));
  f_exporter = Teuchos::rcp(new Epetra_Export(*f_overlap_map, *f_map));

  // Do some consistency checking
  TEST_FOR_EXCEPTION(
    !p_map_A->IsSameAs(*g_map_B), std::logic_error,
    "Model A parameter map for index " << pIndexA << " must be the same " <<
    "map as model B response map for index " << gIndexB << "!");
  TEST_FOR_EXCEPTION(
    !p_map_B->IsSameAs(*g_map_A), std::logic_error,
    "Model B parameter map for index " << pIndexB << " must be the same " <<
    "map as model A response map for index " << gIndexA << "!");

  // Build the Jacobian graph, which looks like
  //   [     I      dg_2/dp_2 ]
  //   [ dg_1/dp_1     I      ]
  W_graph = Teuchos::rcp(new Epetra_CrsGraph(Copy, *f_map, nx));
  for (int i=0; i<f_map->NumMyElements(); i++) {
    int row = f_map->GID(i);
    // Diagonal part
    int index = row;
    W_graph->InsertGlobalIndices(row, 1, &index);
    if (row < n_p_A) {
      // dg_2/dp_2 part
      for (int j=0; j<n_p_B; j++) {
	index = n_p_A + j;
	W_graph->InsertGlobalIndices(row, 1, &index);
      }
    }
    else {
      // dg_1/dp_1 part
      for (int j=0; j<n_p_A; j++) {
	index = j;
	W_graph->InsertGlobalIndices(row, 1, &index);
      }
    }
  }
  W_graph->FillComplete();

  // Build initial guess
  Epetra_Vector x_init_overlap(*x_overlap_map);
  Teuchos::RCP<const Epetra_Vector> p_init_A = modelA->get_p_init(pIndexA);
  Teuchos::RCP<const Epetra_Vector> p_init_B = modelB->get_p_init(pIndexB);
  for (int i=0; i<n_p_A; i++)
    x_init_overlap[i] = (*p_init_A)[i];
  for (int i=0; i<n_p_B; i++)
    x_init_overlap[n_p_A + i] = (*p_init_B)[i];
  x_init = Teuchos::rcp(new Epetra_Vector(*x_map));
  x_init->Export(*x_init_overlap, *x_importer, Insert);

  // Create storage for parameters, responses, and derivatives
  p_A = Teuchos::rcp(new Epetra_Vector(*p_map_A));
  p_B = Teuchos::rcp(new Epetra_Vector(*p_map_B));
  g_A = Teuchos::rcp(new Epetra_Vector(*g_map_A));
  g_B = Teuchos::rcp(new Epetra_Vector(*g_map_B));
  Teuchos::RCP<Epetra_MultiVector> dgdp_vec_A = 
    Teuchos::rcp(new Epetra_MultiVector(*g_map_A, n_p_A));
  Teuchos::RCP<Epetra_MultiVector> dgdp_vec_B = 
    Teuchos::rcp(new Epetra_MultiVector(*g_map_B, n_p_B));
  dgdp_A = EpetraExt::ModelEvaluator::Derivative(dgdp_vec_A, DERIV_MV_BY_COL);
  dgdp_B = EpetraExt::ModelEvaluator::Derivative(dgdp_vec_B, DERIV_MV_BY_COL);

  // Get the process ID and the total number of processors
  MyPID = Comm->MyPID();
  NumProc = Comm->NumProc();

}

/*
void
Piro::Epetra::NECoupledModelEvaluator::init_sg(
  const Teuchos::RCP<const Stokhos::ProductBasis<int,double> >& Basis, 
  const Teuchos::RCP<const Stokhos::Quadrature<int,double> >& Quad,
  const Teuchos::RCP<const EpetraExt::MultiComm>& multiComm)
{
  basis = Basis;
  quad = Quad;

  if (basis != Teuchos::null) {
    sg_overlap_map =
      Teuchos::rcp(new Epetra_LocalMap(basis->size(), 0, 
				       multiComm->TimeDomainComm()));
    p1_sg = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		     basis, sg_overlap_map, map_p, multiComm));
    p2_sg = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		     basis, sg_overlap_map, map_p, multiComm));
    rVars_sg = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		     basis, sg_overlap_map, map_rvar, multiComm));
    g1_sg = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		     basis, sg_overlap_map, map_g, multiComm));
    g2_sg = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		     basis, sg_overlap_map, map_g, multiComm));
    dgdp1_sg = 
      Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
		     basis, sg_overlap_map, map_g, multiComm, 
		     map_p->NumGlobalElements()));
    dgdp2_sg = 
      Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
		     basis, sg_overlap_map, map_g, multiComm,
		     map_p->NumGlobalElements()));
  }
}
*/

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RefCountPtr<const Epetra_Map>
Piro::Epetra::NECoupledModelEvaluator::
get_x_map() const
{
  return x_map;
}

Teuchos::RefCountPtr<const Epetra_Map>
Piro::Epetra::NECoupledModelEvaluator::
get_f_map() const
{
  return f_map;
}

Teuchos::RefCountPtr<const Epetra_Vector>
Piro::Epetra::NECoupledModelEvaluator::
get_x_init() const
{
  return x_init;
}

Teuchos::RefCountPtr<const Epetra_Map>
Piro::Epetra::NECoupledModelEvaluator::
get_p_map(int j) const
{
  TEST_FOR_EXCEPTION(
    j >= num_params_total || j < 0, Teuchos::Exceptions::InvalidParameter,
    std::endl <<
    "Error in Piro::Epetra::NECoupledModelEvaluator::get_p_map():  " <<
    "Invalid parameter index j = " << j << std::endl);

  if (j < num_params_A-1)
    return modelA->get_p_map(param_map[j]);
  return modelB->get_p_map(param_map[j]);
}

Teuchos::RefCountPtr<const Epetra_Map>
Piro::Epetra::NECoupledModelEvaluator::
get_g_map(int j) const
{
  TEST_FOR_EXCEPTION(
    j >= num_responses_total || j < 0, Teuchos::Exceptions::InvalidParameter,
    std::endl <<
    "Error in Piro::Epetra::NECoupledModelEvaluator::get_g_map():  " <<
    "Invalid response index j = " << j << std::endl);

  if (j < num_responses_A-1)
    return modelA->get_g_map(response_map[j]);
  return modelB->get_g_map(response_map[j]);
}

Teuchos::RCP<const Teuchos::Array<std::string> > 
Piro::Epetra::NECoupledModelEvaluator::
get_p_names(int j) const
{
  TEST_FOR_EXCEPTION(
    j >= num_params_total || j < 0, Teuchos::Exceptions::InvalidParameter,
    std::endl <<
    "Error in Piro::Epetra::NECoupledModelEvaluator::get_p_names():  " <<
    "Invalid parameter index j = " << j << std::endl);

  if (j < num_params_A-1)
    return modelA->get_p_names(param_map[j]);
  return modelB->get_p_names(param_map[j]);
}

Teuchos::RefCountPtr<const Epetra_Vector>
Piro::Epetra::NECoupledModelEvaluator::
get_p_init(int j) const
{
  TEST_FOR_EXCEPTION(
    j >= num_params_total || j < 0, Teuchos::Exceptions::InvalidParameter,
    std::endl <<
    "Error in Piro::Epetra::NECoupledModelEvaluator::get_p_init():  " <<
    "Invalid parameter index j = " << j << std::endl);

  if (j < num_params_A-1)
    return modelA->get_p_init(param_map[j]);
  return modelB->get_p_init(param_map[j]);
}

Teuchos::RefCountPtr<Epetra_Operator>
Piro::Epetra::NECoupledModelEvaluator::
create_W() const
{
  Teuchos::RCP<Epetra_CrsMatrix> mat =
    Teuchos::rcp(new Epetra_CrsMatrix(Copy, *W_graph));
  mat->FillComplete();
  return mat;
}

EpetraExt::ModelEvaluator::InArgs
Piro::Epetra::NECoupledModelEvaluator::
createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());

  EpetraExt::ModelEvaluator::InArgs modelA_inargs = modelA->createInArgs();
  EpetraExt::ModelEvaluator::InArgs modelB_inargs = modelB->createInArgs();

  // Deterministic InArgs
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.set_Np(num_params_total);
  
  // Stochastic InArgs
  inArgs.setSupports(IN_ARG_x_sg,true);
  for (int i=0; i<num_params_A-1; i++)
    inArgs.setSupports(IN_ARG_p_sg, i, 
		       modelA_inargs.supports(IN_ARG_p_sg, param_map[i])); 
  for (int i=0; i<num_params_B-1; i++)
    inArgs.setSupports(IN_ARG_p_sg, i, 
		       modelB_inargs.supports(IN_ARG_p_sg, 
					      param_map[num_params_A-1+i])); 
  inArgs.setSupports(IN_ARG_sg_basis,true);
  inArgs.setSupports(IN_ARG_sg_quadrature,true);
  inArgs.setSupports(IN_ARG_sg_expansion,true);

  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
Piro::Epetra::NECoupledModelEvaluator::
createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());

  EpetraExt::ModelEvaluator::OutArgs modelA_outargs = modelA->createOutArgs();
  EpetraExt::ModelEvaluator::OutArgs modelB_outargs = modelB->createOutArgs();

  // Deterministic OutArgs
  outArgs.setSupports(OUT_ARG_f, true);
  outArgs.setSupports(OUT_ARG_W, true);
  outArgs.set_W_properties(
    DerivativeProperties(DERIV_LINEARITY_NONCONST, DERIV_RANK_FULL, true));
  outArgs.set_Np_Ng(num_params_total, num_responses_total);

  // Stochastic OutArgs
  outArgs.setSupports(OUT_ARG_f_sg, true);
  outArgs.setSupports(OUT_ARG_W_sg, true);
  for (int i=0; i<num_responses_A-1; i++)
    outArgs.setSupports(OUT_ARG_g_sg, i, 
		       modelA_outargs.supports(OUT_ARG_g_sg, response_map[i])); 
  for (int i=0; i<num_responses_B-1; i++)
    outArgs.setSupports(OUT_ARG_g_sg, i, 
		       modelB_outargs.supports(OUT_ARG_g_sg, 
					       response_map[num_responses_A-1+i]));
  
  return outArgs;
}

void Piro::Epetra::NECoupledModelEvaluator::
evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const
{
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;

  //
  // Get the input arguments
  //
  const Epetra_Vector* x = inArgs.get_x().get();
  const Epetra_Vector* p = inArgs.get_p(0).get();
  InArgs::sg_const_vector_t x_sg = inArgs.get_x_sg();
  InArgs::sg_const_vector_t p_sg = inArgs.get_p_sg(0);

  //
  // Get the output arguments
  //
  Epetra_Vector       *f_out = outArgs.get_f().get();
  Epetra_Operator     *W_out = outArgs.get_W().get();
  OutArgs::sg_vector_t f_sg = outArgs.get_f_sg();
  OutArgs::sg_operator_t W_sg = outArgs.get_W_sg();

  // See if there is anything to do
  if (f_out == NULL && W_out == NULL &&
      f_sg == Teuchos::null && W_sg == Teuchos::null)
    return;

  if (eval_W_with_f && f_out == NULL && W_out != NULL) {
    Epetra_CrsMatrix* W_crs = dynamic_cast<Epetra_CrsMatrix*>(W_out);
    *W_crs = *A;
    return;
  }
  if (eval_W_with_f && f_out != NULL && W_out == NULL)
    W_out = A.get();

  DiffProbLocal = DiffProb;
  HeatProbLocal = HeatProb;

  if (x_sg != Teuchos::null || p_sg != Teuchos::null) {
    basis = 
      Teuchos::rcp_dynamic_cast<const Stokhos::ProductBasis<int,double> >(
	inArgs.get_sg_basis(), true);
    quad = inArgs.get_sg_quadrature();
    expansion = inArgs.get_sg_expansion();

    if (p1_sg == Teuchos::null) {
      Teuchos::RCP<const EpetraExt::MultiComm> multiComm;
      if (x_sg != Teuchos::null)
	multiComm = x_sg->productComm();
      else
	multiComm = p_sg->productComm();
      sg_overlap_map =
	Teuchos::rcp(new Epetra_LocalMap(basis->size(), 0, 
					 multiComm->TimeDomainComm()));
      p1_sg = 
	Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		       basis, sg_overlap_map, map_p, multiComm));
      p2_sg = 
	Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		       basis, sg_overlap_map, map_p, multiComm));
      rVars_sg = 
	Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		       basis, sg_overlap_map, map_rvar, multiComm));
      g1_sg = 
	Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		       basis, sg_overlap_map, map_g, multiComm));
      g2_sg = 
	Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		       basis, sg_overlap_map, map_g, multiComm));
      dgdp1_sg = 
	Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
		       basis, sg_overlap_map, map_g, multiComm, 
		       map_p->NumGlobalElements()));
      dgdp2_sg = 
	Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
		       basis, sg_overlap_map, map_g, multiComm,
		       map_p->NumGlobalElements()));
    }
  }

  if (p!=NULL)
    *rVars = *p;
  if (p_sg != Teuchos::null)
    for (int i=0; i<p_sg->size(); i++)
      (*rVars_sg)[i] = (*p_sg)[i];


  double tempx0;
  double tempx1;
  double xval0;
  double xval1;

  //Load the parameter vectors from all of the processors
  if (x) {
    if (MyPID==0){
      xval0=(*x)[0];
      xval1=(*x)[1];
    }
    else {
      xval0=0.0;
      xval1=0.0;
    }
    (*Comm).SumAll(&xval0,&tempx0,1);
    (*Comm).SumAll(&xval1,&tempx1,1);
    
    (*p1)[0]=tempx1;
    (*p2)[0]=tempx0;
  }

  if (x_sg != Teuchos::null) {
    for (int i=0; i<x_sg->size(); i++) {
      if (MyPID==0){
	xval0=(*x_sg)[i][0];
	xval1=(*x_sg)[i][1];
      }
      else {
	xval0=0.0;
	xval1=0.0;
      }
      (*Comm).SumAll(&xval0,&tempx0,1);
      (*Comm).SumAll(&xval1,&tempx1,1);
    
      (*p1_sg)[i][0]=tempx1;
      (*p2_sg)[i][0]=tempx0;
    }
  }

  //Prepare to split the vector of random values into the number required 
  //for the two different problems
  Epetra_LocalMap DrVarMap(DnumRandVar,0,*Comm);
  Epetra_LocalMap HrVarMap(HnumRandVar,0,*Comm);

  Teuchos::RCP<Epetra_Vector> DrVars=Teuchos::rcp(new Epetra_Vector(DrVarMap));
  Teuchos::RCP<Epetra_Vector> HrVars=Teuchos::rcp(new Epetra_Vector(HrVarMap));

  if (x) {
    for (int i=0;i<DnumRandVar;++i) 
      (*DrVars)[i]=(*rVars)[i];
    
    for (int i=0;i<HnumRandVar;++i) 
      (*HrVars)[i]=(*rVars)[(*rVars).MyLength()-DnumRandVar+i];
  }

  // Orthogonalize SG bases for each problem
  OutArgs::sg_vector_t p1_gs;
  OutArgs::sg_vector_t p2_gs;
  OutArgs::sg_vector_t DrVars_sg;
  OutArgs::sg_vector_t HrVars_sg;
  OutArgs::sg_vector_t g1_gs;
  OutArgs::sg_vector_t g2_gs;
  Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > dgdp1_gs;
  Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > dgdp2_gs;
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > diff_gs_basis;
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > heat_gs_basis;
  Teuchos::RCP<const Stokhos::Quadrature<int,double> > diff_gs_quad;
  Teuchos::RCP<const Stokhos::Quadrature<int,double> > heat_gs_quad;
  Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > diff_gs_expansion;
  Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > heat_gs_expansion;
  Teuchos::RCP<const Teuchos::Array< Teuchos::Array<double> > > diff_basis_vals;
  Teuchos::RCP<const Teuchos::Array< Teuchos::Array<double> > > heat_basis_vals;
  if (x_sg != Teuchos::null && reduce_dimension) {
    TEUCHOS_FUNC_TIME_MONITOR("NECoupledModelEvaluator -- dimension reduction");
    int order = basis->order();
    int new_order = order;
    int sz = p1_sg->size();

    // Copy Epetra PCEs into Stokhos PCE objects
    Teuchos::Array< Stokhos::OrthogPolyApprox<int,double> > diff_opa(DnumRandVar + 1);
    Teuchos::Array< Stokhos::OrthogPolyApprox<int,double> > heat_opa(HnumRandVar + 1);
    for (int j=0; j<DnumRandVar+1; j++)
      diff_opa[j].reset(basis);
    for (int j=0; j<HnumRandVar+1; j++)
      heat_opa[j].reset(basis);
    for (int i=0; i<sz; i++) {
      diff_opa[0][i] = (*p1_sg)[i][0];
      for (int j=0; j<DnumRandVar; j++)
	diff_opa[j+1][i] = (*rVars_sg)[i][j];

      heat_opa[0][i] = (*p2_sg)[i][0];
      for (int j=0; j<HnumRandVar; j++)
	heat_opa[j+1][i] = (*rVars_sg)[i][j+DnumRandVar];
    }

    // Build Stieltjes basis, quadrature, and new PCEs
    Teuchos::Array<Stokhos::OrthogPolyApprox<int,double> > diff_gs_pces;
    if (orthogonalize_bases) {
      Stokhos::StieltjesGramSchmidtBuilder<int,double> diff_gs_builder(
      quad, diff_opa, new_order, true, false);
      diff_gs_basis = diff_gs_builder.getReducedBasis();
      diff_gs_quad = diff_gs_builder.getReducedQuadrature();
      diff_basis_vals = Teuchos::rcp(&(diff_gs_quad->getBasisAtQuadPoints()),
				     false);
      diff_gs_builder.computeReducedPCEs(diff_opa, diff_gs_pces);
    }
    else {
      Teuchos::RCP<const Stokhos::ProductBasis<int,double> > product_basis = 
	Teuchos::rcp_dynamic_cast<const Stokhos::ProductBasis<int,double> >(
	  basis, true);
      Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double > > >
	coordinate_bases = product_basis->getCoordinateBases();
      Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double > > >
	new_coordinate_bases(diff_opa.size());
       Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk = 
	 DiffParams->sublist("Problem").sublist("Stochastic Galerkin").get< Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > >(
	   "Triple Product Tensor");
      if (st_quad == Teuchos::null) {
	st_quad = quad;
	// st_quad =
	//   Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(
	// 		 basis, new_order+1));
      }
      new_coordinate_bases[0] = Teuchos::rcp(
	new Stokhos::StieltjesPCEBasis<int,double>(
	  new_order, Teuchos::rcp(&(diff_opa[0]),false), st_quad, 
	  false, false, true, Cijk));
      for (int j=0; j<DnumRandVar; j++)
	new_coordinate_bases[j+1] = coordinate_bases[j];
      Teuchos::RCP<const Stokhos::ProductBasis<int,double> > tensor_basis = 
	Teuchos::rcp(
	  new Stokhos::CompletePolynomialBasis<int,double>(new_coordinate_bases)
	  );
      diff_gs_basis = tensor_basis;
      if (diff_gs_basis->dimension() <= 0)
      	diff_gs_quad = 
      	  Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(
      			 tensor_basis));
      else
#ifdef HAVE_STOKHOS_DAKOTA
	diff_gs_quad = 
	  Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(
			 tensor_basis, new_order));
#else
        diff_gs_quad = 
      	  Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(
      			 tensor_basis));
#endif
      const Teuchos::Array< Teuchos::Array<double> >& points = 
	  quad->getQuadPoints();
      const Teuchos::Array< Teuchos::Array<double> >& basis_vals = 
	  quad->getBasisAtQuadPoints();
      int nqp = points.size();
      Teuchos::Array<double> diff_opa_val(diff_opa.size());
      Teuchos::RCP< Teuchos::Array< Teuchos::Array<double> > > ncdiff_basis_vals
	= Teuchos::rcp(new Teuchos::Array< Teuchos::Array<double> >(nqp));
      for (int i=0; i<nqp; i++) {
	for (int j=0; j<diff_opa_val.size(); j++)
	  diff_opa_val[j] = diff_opa[j].evaluate(points[i], basis_vals[i]);
	(*ncdiff_basis_vals)[i].resize(diff_gs_basis->size());
	diff_gs_basis->evaluateBases(diff_opa_val, (*ncdiff_basis_vals)[i]);
      }
      diff_basis_vals = ncdiff_basis_vals;
      diff_gs_pces.resize(diff_opa.size());
      for (int k=0; k<diff_opa.size(); k++) {
	diff_gs_pces[k].reset(diff_gs_basis);
	diff_gs_pces[k].term(k, 0) = diff_opa[k].mean();
	diff_gs_pces[k].term(k, 1) = 1.0; 
      }
    }
    
    Teuchos::Array<Stokhos::OrthogPolyApprox<int,double> > heat_gs_pces;
    if (orthogonalize_bases) {
      Stokhos::StieltjesGramSchmidtBuilder<int,double> heat_gs_builder(
      quad, heat_opa, new_order, true, false);
      heat_gs_basis = heat_gs_builder.getReducedBasis();
      heat_gs_quad = heat_gs_builder.getReducedQuadrature();
      heat_basis_vals = Teuchos::rcp(&(heat_gs_quad->getBasisAtQuadPoints()),
				     false);
      heat_gs_builder.computeReducedPCEs(heat_opa, heat_gs_pces);
    }
    else {
      Teuchos::RCP<const Stokhos::ProductBasis<int,double> > product_basis = 
	Teuchos::rcp_dynamic_cast<const Stokhos::ProductBasis<int,double> >(
	  basis, true);
      Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double > > >
	coordinate_bases = product_basis->getCoordinateBases();
      Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double > > >
	new_coordinate_bases(heat_opa.size());
      Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk = 
	 HeatParams->sublist("Problem").sublist("Stochastic Galerkin").get< Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > >(
	   "Triple Product Tensor");
      new_coordinate_bases[0] = Teuchos::rcp(
	new Stokhos::StieltjesPCEBasis<int,double>(
	  new_order, Teuchos::rcp(&(heat_opa[0]),false), st_quad, 
	  false, false, true, Cijk));
      for (int j=0; j<HnumRandVar; j++)
	new_coordinate_bases[j+1] = coordinate_bases[j];
      Teuchos::RCP<const Stokhos::ProductBasis<int,double> > tensor_basis = 
	Teuchos::rcp(
	  new Stokhos::CompletePolynomialBasis<int,double>(new_coordinate_bases)
	  );
      heat_gs_basis = tensor_basis;
      if (heat_gs_basis->dimension() <= 0)
      	heat_gs_quad = 
      	  Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(
      			 tensor_basis));
      else
#ifdef HAVE_STOKHOS_DAKOTA
	heat_gs_quad = 
	  Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(
			 tensor_basis, new_order));
#else
        heat_gs_quad = 
      	  Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(
      			 tensor_basis));
#endif
      const Teuchos::Array< Teuchos::Array<double> >& points = 
	  quad->getQuadPoints();
      const Teuchos::Array< Teuchos::Array<double> >& basis_vals = 
	  quad->getBasisAtQuadPoints();
      int nqp = points.size();
      Teuchos::Array<double> heat_opa_val(heat_opa.size());
      Teuchos::RCP< Teuchos::Array< Teuchos::Array<double> > > ncheat_basis_vals
	= Teuchos::rcp(new Teuchos::Array< Teuchos::Array<double> >(nqp));
      for (int i=0; i<nqp; i++) {
	for (int j=0; j<heat_opa_val.size(); j++)
	  heat_opa_val[j] = heat_opa[j].evaluate(points[i], basis_vals[i]);
	(*ncheat_basis_vals)[i].resize(heat_gs_basis->size());
	heat_gs_basis->evaluateBases(heat_opa_val, (*ncheat_basis_vals)[i]);
      }
      heat_basis_vals = ncheat_basis_vals;
      heat_gs_pces.resize(heat_opa.size());
      for (int k=0; k<heat_opa.size(); k++) {
	heat_gs_pces[k].reset(heat_gs_basis);
	heat_gs_pces[k].term(k, 0) = heat_opa[k].mean();
	heat_gs_pces[k].term(k, 1) = 1.0; 
      }
    }

    Teuchos::RCP<const EpetraExt::MultiComm> multiComm = x_sg->productComm();
    
    // Copy into Epetra objects
    int diff_sz = diff_gs_basis->size();
    Teuchos::RCP<const Epetra_BlockMap> diff_gs_overlap_map =
      Teuchos::rcp(new Epetra_LocalMap(diff_sz, 0, 
				       multiComm->TimeDomainComm()));
    p1_gs = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		     diff_gs_basis, diff_gs_overlap_map, map_p, multiComm));
    DrVars_sg = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		     diff_gs_basis, diff_gs_overlap_map, 
		     Teuchos::rcp(&DrVarMap, false), multiComm));
    for (int k=0; k<diff_sz; k++) {
      (*p1_gs)[k][0] = diff_gs_pces[0][k];
      for (int i=0; i<DnumRandVar; i++)
	(*DrVars_sg)[k][i] = diff_gs_pces[i+1][k];
    }

    int heat_sz = heat_gs_basis->size();
    Teuchos::RCP<const Epetra_BlockMap> heat_gs_overlap_map =
      Teuchos::rcp(new Epetra_LocalMap(heat_sz, 0, 
				       multiComm->TimeDomainComm()));
    p2_gs = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		     heat_gs_basis, heat_gs_overlap_map, map_p, multiComm));
    HrVars_sg = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		     heat_gs_basis, heat_gs_overlap_map, 
		     Teuchos::rcp(&HrVarMap, false), multiComm));
    for (int k=0; k<heat_sz; k++) {
      (*p2_gs)[k][0] = heat_gs_pces[0][k];
      for (int i=0;i<HnumRandVar;++i)
	(*HrVars_sg)[k][i] = heat_gs_pces[i+1][k];
    }
    
    g1_gs = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		     diff_gs_basis, diff_gs_overlap_map, map_g, multiComm));
    g2_gs = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		     heat_gs_basis, heat_gs_overlap_map, map_g, multiComm));
    dgdp1_gs = 
      Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
		     diff_gs_basis, diff_gs_overlap_map, map_g, multiComm,
		     map_p->NumGlobalElements()));
    dgdp2_gs = 
      Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
		     heat_gs_basis, heat_gs_overlap_map, map_g, multiComm,
		     map_p->NumGlobalElements()));

    // Setup new model evaluators
    Teuchos::RCP<Teuchos::ParameterList> DiffParamsGS = 
      Teuchos::rcp(new Teuchos::ParameterList(*DiffParams));
    Teuchos::ParameterList& diff_sg_params = 
      DiffParamsGS->sublist("Problem").sublist("Stochastic Galerkin");
    diff_sg_params.sublist("Basis").set("Stochastic Galerkin Basis",
					diff_gs_basis);
    diff_sg_params.sublist("Quadrature").set("Stochastic Galerkin Quadrature",
					     diff_gs_quad);
    std::string diff_sg_method = diff_sg_params.get("SG Method", "AD");
    if (diff_sg_params.sublist("Expansion").isParameter("Stochastic Galerkin Expansion"))
      diff_sg_params.sublist("Expansion").remove("Stochastic Galerkin Expansion");
    if (diff_sg_params.isParameter("Triple Product Tensor"))
      diff_sg_params.remove("Triple Product Tensor");
    if (diff_sg_method != "NI" && diff_sg_method != "Non-intrusive") {
      Stokhos::ExpansionFactory<int,double>::create(diff_sg_params);
      diff_gs_expansion = diff_sg_params.sublist("Expansion").get< Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > >(
	"Stochastic Galerkin Expansion");
    }
    Teuchos::RCP<Stokhos::SGModelEvaluator> diff_sg_model;
    DiffProbLocal = 
      SGModelEvaluatorFactory::createModelEvaluator(DiffModel, DiffParamsGS, 
						    Comm, diff_sg_model);

    Teuchos::RCP<Teuchos::ParameterList> HeatParamsGS = 
      Teuchos::rcp(new Teuchos::ParameterList(*HeatParams));
    Teuchos::ParameterList& heat_sg_params = 
      HeatParamsGS->sublist("Problem").sublist("Stochastic Galerkin");
    heat_sg_params.sublist("Basis").set("Stochastic Galerkin Basis",
					heat_gs_basis);
    heat_sg_params.sublist("Quadrature").set("Stochastic Galerkin Quadrature",
					     heat_gs_quad);
    std::string heat_sg_method = heat_sg_params.get("SG Method", "AD");
    if (heat_sg_params.sublist("Expansion").isParameter("Stochastic Galerkin Expansion"))
      heat_sg_params.sublist("Expansion").remove("Stochastic Galerkin Expansion");
    if (heat_sg_params.isParameter("Triple Product Tensor"))
      heat_sg_params.remove("Triple Product Tensor");
    if (heat_sg_method != "NI" && heat_sg_method != "Non-intrusive") {
      Stokhos::ExpansionFactory<int,double>::create(heat_sg_params);
      heat_gs_expansion = diff_sg_params.sublist("Expansion").get< Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > >(
	"Stochastic Galerkin Expansion");
    }
    Teuchos::RCP<Stokhos::SGModelEvaluator> heat_sg_model;
    HeatProbLocal = 
      SGModelEvaluatorFactory::createModelEvaluator(HeatModel, HeatParamsGS, 
						    Comm, heat_sg_model);
  }
  else if (x_sg != Teuchos::null) {
    diff_gs_basis = basis;
    heat_gs_basis = basis;
    diff_gs_quad = quad;
    heat_gs_quad = quad;
    diff_gs_expansion = expansion;
    heat_gs_expansion = expansion;
    p1_gs = p1_sg;
    p2_gs = p2_sg;
    Teuchos::RCP<const EpetraExt::MultiComm> multiComm = x_sg->productComm();
    DrVars_sg = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		     basis, sg_overlap_map, Teuchos::rcp(&DrVarMap,false),
		     multiComm));
    HrVars_sg = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		     basis, sg_overlap_map, Teuchos::rcp(&HrVarMap,false),
		     multiComm));
    for (int k=0; k<x_sg->size(); k++) {
      for (int i=0;i<DnumRandVar;++i) 
	(*DrVars_sg)[k][i]=(*rVars_sg)[k][i];
      
      for (int i=0;i<HnumRandVar;++i) 
	(*HrVars_sg)[k][i]=(*rVars_sg)[k][(*rVars).MyLength()-DnumRandVar+i];
    }
    g1_gs = g1_sg;
    g2_gs = g2_sg;
    dgdp1_gs = dgdp1_sg;
    dgdp2_gs = dgdp2_sg;
  }

  //Setup the Diffusion Problem
  EpetraExt::ModelEvaluator::InArgs inArgs1 = DiffProb->createInArgs();
  EpetraExt::ModelEvaluator::OutArgs outArgs1 = DiffProb->createOutArgs();
  if (x) {
    inArgs1.set_p(0,p1);
    inArgs1.set_p(1,DrVars);
  }
  if (x_sg != Teuchos::null) {
    if (inArgs1.supports(IN_ARG_sg_basis))
      inArgs1.set_sg_basis(diff_gs_basis);
    if (inArgs1.supports(IN_ARG_sg_quadrature))
      inArgs1.set_sg_quadrature(diff_gs_quad);
    if (inArgs1.supports(IN_ARG_sg_expansion))
      inArgs1.set_sg_expansion(diff_gs_expansion);
    inArgs1.set_p_sg(0, p1_gs);
    inArgs1.set_p_sg(1, DrVars_sg);
  }

  //Setup the Heat Problem
  EpetraExt::ModelEvaluator::InArgs inArgs2 = HeatProb->createInArgs();
  EpetraExt::ModelEvaluator::OutArgs outArgs2 = HeatProb->createOutArgs();
  if (x) {
    inArgs2.set_p(0,p2);
    inArgs2.set_p(1,HrVars);
  }
  if (x_sg != Teuchos::null) {
    if (inArgs2.supports(IN_ARG_sg_basis))
      inArgs2.set_sg_basis(heat_gs_basis);
    if (inArgs2.supports(IN_ARG_sg_quadrature))
      inArgs2.set_sg_quadrature(heat_gs_quad);
    if (inArgs2.supports(IN_ARG_sg_expansion))
      inArgs2.set_sg_expansion(heat_gs_expansion);
    inArgs2.set_p_sg(0, p2_gs);
    inArgs2.set_p_sg(1, HrVars_sg);
  }

  
  if(f_out) {
    outArgs1.set_g(0,g1);
    outArgs2.set_g(0,g2);
  }
  if(W_out) {
    outArgs1.set_DgDp(0,0,*dgdp1);
    outArgs2.set_DgDp(0,0,*dgdp2);
  }
  if(f_sg != Teuchos::null) {
    outArgs1.set_g_sg(0,g1_gs);
    outArgs2.set_g_sg(0,g2_gs);
  }
  if(W_sg != Teuchos::null) {
    outArgs1.set_DgDp_sg(0,0,dgdp1_gs);
    outArgs2.set_DgDp_sg(0,0,dgdp2_gs);
  }


  {
  TEUCHOS_FUNC_TIME_MONITOR("NECoupledModelEvaluator -- Diffusion nonlinear elimination");
  std::cout << "eliminating diffusion states..." << std::endl;
  DiffProbLocal->evalModel(inArgs1,outArgs1);
  }

  {
  TEUCHOS_FUNC_TIME_MONITOR("NECoupledModelEvaluator -- Heat nonlinear elimination");
  std::cout << "eliminating heat states..." << std::endl;
  HeatProbLocal->evalModel(inArgs2,outArgs2);
  }

  if(f_out){  

    //Residual on the first processor
    f_out->PutScalar(0.0);
    if (MyPID==0){
      (*f_out)[0]=(*x)[0]-(*g1)[0];
      (*f_out)[1]=(*x)[1]-(*g2)[0];
    }
  }


  int* index=new int[1];
  double* valu=new double[1];

  //Jacobian on the first processor
  if(W_out){
    (*dirvec1).Scale(-1);
    (*dirvec2).Scale(-1);
    Epetra_CrsMatrix &DfDx = dyn_cast<Epetra_CrsMatrix>(*W_out);
    if (MyPID==0){
      (index)[0]=0;
      valu[0]=1.0;
      DfDx.InsertGlobalValues(0,1,valu,index);
      DfDx.InsertGlobalValues(1,1,(*dirvec2)[0],index);
      DfDx.ReplaceGlobalValues(0,1,valu,index);
      DfDx.ReplaceGlobalValues(1,1,(*dirvec2)[0],index);
      (index)[0]=1;
      DfDx.InsertGlobalValues(1,1,valu,index);
      DfDx.InsertGlobalValues(0,1,(*dirvec1)[0],index);
      DfDx.ReplaceGlobalValues(1,1,valu,index);
      DfDx.ReplaceGlobalValues(0,1,(*dirvec1)[0],index);
    }
    DfDx.FillComplete();
  }

  if(f_sg != Teuchos::null){  

    // Project back to original basis
    if (reduce_dimension) {
      TEUCHOS_FUNC_TIME_MONITOR("NECoupledModelEvaluator -- dimension reduction g integration");
      const Teuchos::Array<double>& weights = quad->getQuadWeights();
      const Teuchos::Array< Teuchos::Array<double> >& basis_vals = 
	quad->getBasisAtQuadPoints();
      int nqp = weights.size();
      const Teuchos::Array<double>& norms = basis->norm_squared();
      Epetra_Vector g1_val(*map_g);
      Epetra_Vector g2_val(*map_g);
      g1_sg->init(0.0);
      g2_sg->init(0.0);
      for (int i=0; i<nqp; i++) {
	g1_gs->evaluate((*diff_basis_vals)[i], g1_val);
	g2_gs->evaluate((*heat_basis_vals)[i], g2_val);
	g1_sg->sumIntoAllTerms(weights[i], basis_vals[i], norms, g1_val);
	g2_sg->sumIntoAllTerms(weights[i], basis_vals[i], norms, g2_val);
      }
    }

    //Residual on the first processor
    f_sg->init(0.0);
    if (MyPID==0){
      for (int i=0; i<f_sg->size(); i++) {
	(*f_sg)[i][0]=(*x_sg)[i][0]-(*g1_sg)[i][0];
	(*f_sg)[i][1]=(*x_sg)[i][1]-(*g2_sg)[i][0];
      }
    }
    //std::cout << "f_sg = " << std::endl << *f_sg << std::endl;
  }

  //Jacobian on the first processor
  if(W_sg != Teuchos::null){
    
    // Project back to original basis
    if (reduce_dimension) {
      TEUCHOS_FUNC_TIME_MONITOR("NECoupledModelEvaluator -- dimension reduction dg/dp integration");
      const Teuchos::Array<double>& weights = quad->getQuadWeights();
      const Teuchos::Array< Teuchos::Array<double> >& basis_vals = 
	quad->getBasisAtQuadPoints();
      int nqp = weights.size();
      const Teuchos::Array<double>& norms = basis->norm_squared();
      Epetra_MultiVector dgdp1_val(*map_g, 1);
      Epetra_MultiVector dgdp2_val(*map_g, 1);
      dgdp1_sg->init(0.0);
      dgdp2_sg->init(0.0);
      for (int i=0; i<nqp; i++) {
	dgdp1_gs->evaluate((*diff_basis_vals)[i], dgdp1_val);
	dgdp2_gs->evaluate((*heat_basis_vals)[i], dgdp2_val);
	dgdp1_sg->sumIntoAllTerms(weights[i], basis_vals[i], norms, dgdp1_val);
	dgdp2_sg->sumIntoAllTerms(weights[i], basis_vals[i], norms, dgdp2_val);
      }
    }

    for (int i=0; i<W_sg->size(); i++) {
      Epetra_MultiVector& dgdp1_mv = (*dgdp1_sg)[i];
      Epetra_MultiVector& dgdp2_mv = (*dgdp2_sg)[i];
      dgdp1_mv.Scale(-1);
      dgdp2_mv.Scale(-1);
      Epetra_CrsMatrix &DfDx = dyn_cast<Epetra_CrsMatrix>((*W_sg)[i]);
      if (MyPID==0){
	(index)[0]=0;
	if (i == 0)
	  valu[0]=1.0;
	else
	  valu[0]=0.0;
	DfDx.InsertGlobalValues(0,1,valu,index);
	DfDx.InsertGlobalValues(1,1,dgdp2_mv[0],index);
	DfDx.ReplaceGlobalValues(0,1,valu,index);
	DfDx.ReplaceGlobalValues(1,1,dgdp2_mv[0],index);
	(index)[0]=1;
	DfDx.InsertGlobalValues(1,1,valu,index);
	DfDx.InsertGlobalValues(0,1,dgdp1_mv[0],index);
	DfDx.ReplaceGlobalValues(1,1,valu,index);
	DfDx.ReplaceGlobalValues(0,1,dgdp1_mv[0],index);
	
      }
      DfDx.FillComplete();
    }
  }

  delete [] index;
  delete [] valu;
}
