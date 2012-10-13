// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Stokhos_SGModelEvaluator_Interlaced.hpp"

#include <algorithm>
#include "Teuchos_Assert.hpp"
#include "EpetraExt_BlockUtility.h"
#include "EpetraExt_BlockMultiVector.h"
#include "Stokhos_InterlacedOperator.hpp"
#include "Stokhos_SGPreconditionerFactory.hpp"
#include "Stokhos_MatrixFreeOperator.hpp"
#include "Stokhos_EpetraMultiVectorOperator.hpp"
#include "Epetra_LocalMap.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"

Stokhos::SGModelEvaluator_Interlaced::SGModelEvaluator_Interlaced(
  const Teuchos::RCP<EpetraExt::ModelEvaluator>& me_,
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis_,
  const Teuchos::RCP<const Stokhos::Quadrature<int,double> >& sg_quad_,
  const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> >& sg_exp_,
  const Teuchos::RCP<const Stokhos::ParallelData>& sg_parallel_data_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_,
  bool scaleOP_) 
  : me(me_),
    sg_basis(sg_basis_),
    sg_quad(sg_quad_),
    sg_exp(sg_exp_),
    params(params_),
    num_sg_blocks(sg_basis->size()),
    num_W_blocks(sg_basis->size()),
    num_p_blocks(sg_basis->size()),
    supports_x(false),
    x_map(me->get_x_map()),
    f_map(me->get_f_map()),
    sg_parallel_data(sg_parallel_data_),
    sg_comm(sg_parallel_data->getMultiComm()),
    epetraCijk(sg_parallel_data->getEpetraCijk()),
    serialCijk(),
    num_p(0),
    num_p_sg(0),
    sg_p_index_map(),
    sg_p_map(),
    sg_p_names(),
    num_g(0),
    num_g_sg(0),
    sg_g_index_map(),
    sg_g_map(),
    x_dot_sg_blocks(),
    x_sg_blocks(),
    f_sg_blocks(),
    W_sg_blocks(),
    dgdx_dot_sg_blocks(),
    dgdx_sg_blocks(),
    sg_x_init(),
    sg_p_init(),
    eval_W_with_f(false),
    scaleOP(scaleOP_)
{
  if (x_map != Teuchos::null)
    supports_x = true;

  overlapped_stoch_row_map = 
    Teuchos::rcp(new Epetra_LocalMap(
		   static_cast<int>(num_sg_blocks), 0, *(sg_parallel_data->getStochasticComm())));
  stoch_row_map = overlapped_stoch_row_map;
  if (epetraCijk != Teuchos::null)
    stoch_row_map = epetraCijk->getStochasticRowMap();

  if (supports_x) {

    // Create interlace SG x and f maps
    interlace_x_map = buildInterlaceMap(*x_map,*stoch_row_map);
    interlace_overlapped_x_map = buildInterlaceMap(*x_map,*overlapped_stoch_row_map);

    interlace_f_map = buildInterlaceMap(*f_map,*stoch_row_map);
    interlace_overlapped_f_map = buildInterlaceMap(*f_map,*overlapped_stoch_row_map);

    // Create importer/exporter from/to overlapped distribution
    interlace_overlapped_x_importer = 
      Teuchos::rcp(new Epetra_Import(*interlace_overlapped_x_map, *get_x_map()));
    interlace_overlapped_f_exporter = 
      Teuchos::rcp(new Epetra_Export(*interlace_overlapped_f_map, *interlace_f_map));
    
    // now we create the underlying Epetra block vectors
    // that will be used by the model evaluator to construct
    // the solution of the deterministic problem.
    ////////////////////////////////////////////////////////////

    // Create vector blocks
    x_dot_sg_blocks = create_x_sg_overlap();
    x_sg_blocks = create_x_sg_overlap();
    f_sg_blocks = create_f_sg_overlap();

    // Create default sg_x_init
    sg_x_init = create_x_sg();
    if (sg_x_init->myGID(0))
      (*sg_x_init)[0] = *(me->get_x_init());
    
    // Preconditioner needs an x: This is interlaced
    my_x = Teuchos::rcp(new Epetra_Vector(*get_x_map()));


    // setup storage for W, these are blocked in Stokhos
    // format
    ///////////////////////////////////////////////////////
 
    // Determine W expansion type
    std::string W_expansion_type = 
      params->get("Jacobian Expansion Type", "Full");
    if (W_expansion_type == "Linear")
      num_W_blocks = sg_basis->dimension() + 1;
    else
      num_W_blocks = num_sg_blocks;

    Teuchos::RCP<Epetra_BlockMap> W_overlap_map =
      Teuchos::rcp(new Epetra_LocalMap(
		     static_cast<int>(num_W_blocks), 0, 
		     *(sg_parallel_data->getStochasticComm())));
    W_sg_blocks = 
      Teuchos::rcp(new Stokhos::EpetraOperatorOrthogPoly(
		     sg_basis, W_overlap_map, x_map, f_map, interlace_f_map,
		     sg_comm));
    for (unsigned int i=0; i<num_W_blocks; i++)
      W_sg_blocks->setCoeffPtr(i, me->create_W()); // allocate a bunch of matrices

    eval_W_with_f = params->get("Evaluate W with F", false);
  }
    
  // Parameters -- The idea here is to add new parameter vectors
  // for the stochastic Galerkin components of the parameters

  InArgs me_inargs = me->createInArgs();
  OutArgs me_outargs = me->createOutArgs();
  num_p = me_inargs.Np();
  
  // Get the p_sg's supported and build index map
  for (int i=0; i<num_p; i++) {
    if (me_inargs.supports(IN_ARG_p_sg, i))
      sg_p_index_map.push_back(i);
  }
  num_p_sg = sg_p_index_map.size();

  sg_p_map.resize(num_p_sg);
  sg_p_names.resize(num_p_sg);
  sg_p_init.resize(num_p_sg);

  // Determine parameter expansion type
  std::string p_expansion_type = 
    params->get("Parameter Expansion Type", "Full");
  if (p_expansion_type == "Linear")
    num_p_blocks = sg_basis->dimension() + 1;
  else
    num_p_blocks = num_sg_blocks;
  
  // Create parameter maps, names, and initial values
  overlapped_stoch_p_map =
    Teuchos::rcp(new Epetra_LocalMap(
		   static_cast<int>(num_p_blocks), 0, 
		   *(sg_parallel_data->getStochasticComm())));
  for (int i=0; i<num_p_sg; i++) {
    Teuchos::RCP<const Epetra_Map> p_map = me->get_p_map(sg_p_index_map[i]);
    sg_p_map[i] = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
		     *p_map, *overlapped_stoch_p_map, *sg_comm));
    
    Teuchos::RCP<const Teuchos::Array<std::string> > p_names = 
      me->get_p_names(sg_p_index_map[i]);
    if (p_names != Teuchos::null) {
      sg_p_names[i] = 
	Teuchos::rcp(new Teuchos::Array<std::string>(num_sg_blocks*(p_names->size())));
      for (int j=0; j<p_names->size(); j++) {
	std::stringstream ss;
	ss << (*p_names)[j] << " -- SG Coefficient " << i;
	(*sg_p_names[i])[j] = ss.str();
      }
    }

    // Create default sg_p_init
    sg_p_init[i] = create_p_sg(sg_p_index_map[i]);
    sg_p_init[i]->init(0.0);
  }
  
  // Responses -- The idea here is to add new parameter vectors
  // for the stochastic Galerkin components of the respones

  // Get number of SG parameters model supports derivatives w.r.t.
  num_g = me_outargs.Ng();

  // Get the g_sg's supported and build index map
  for (int i=0; i<num_g; i++) {
    if (me_outargs.supports(OUT_ARG_g_sg, i))
      sg_g_index_map.push_back(i);
  }
  num_g_sg = sg_g_index_map.size();

  sg_g_map.resize(num_g_sg);
  dgdx_dot_sg_blocks.resize(num_g_sg);
  dgdx_sg_blocks.resize(num_g_sg);

  // Create response maps
  for (int i=0; i<num_g_sg; i++) {
    Teuchos::RCP<const Epetra_Map> g_map = me->get_g_map(sg_g_index_map[i]);
    sg_g_map[i] = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
		     *g_map, *overlapped_stoch_row_map, *sg_comm));
    
    // Create dg/dxdot, dg/dx SG blocks
    if (supports_x) {
      dgdx_dot_sg_blocks[i] = 
	Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
		       sg_basis,  overlapped_stoch_row_map));
      dgdx_sg_blocks[i] = 
	Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
		       sg_basis, overlapped_stoch_row_map));
    }
  }
  
  // We don't support parallel for dgdx yet, so build a new EpetraCijk
  if (supports_x) {
    serialCijk =
      Teuchos::rcp(new Stokhos::EpetraSparse3Tensor(sg_basis, 
						    epetraCijk->getCijk(), 
						    sg_comm,
						    overlapped_stoch_row_map));
  }

}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
Stokhos::SGModelEvaluator_Interlaced::get_x_map() const
{
  return interlace_x_map;
}

Teuchos::RCP<const Epetra_Map>
Stokhos::SGModelEvaluator_Interlaced::get_f_map() const
{
  return interlace_f_map;
}

Teuchos::RCP<const Epetra_Map>
Stokhos::SGModelEvaluator_Interlaced::get_p_map(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l < 0 || l >= num_p + num_p_sg, std::logic_error,
		     "Error!  Invalid p map index " << l);
  if (l < num_p)
    return me->get_p_map(l);
  else 
    return sg_p_map[l-num_p];

  return Teuchos::null;
}

Teuchos::RCP<const Epetra_Map>
Stokhos::SGModelEvaluator_Interlaced::get_g_map(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l < 0 || l >= num_g_sg, std::logic_error,
		     "Error!  Invalid g map index " << l);
  return sg_g_map[l];
}

Teuchos::RCP<const Teuchos::Array<std::string> >
Stokhos::SGModelEvaluator_Interlaced::get_p_names(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l < 0 || l >= num_p + num_p_sg, std::logic_error,
		     "Error!  Invalid p map index " << l);
  if (l < num_p)
    return me->get_p_names(l);
  else 
    return sg_p_names[l-num_p];
  
  return Teuchos::null;
}

Teuchos::RCP<const Epetra_Vector>
Stokhos::SGModelEvaluator_Interlaced::get_x_init() const
{
  return sg_x_init->getBlockVector();
}

Teuchos::RCP<const Epetra_Vector>
Stokhos::SGModelEvaluator_Interlaced::get_p_init(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l < 0 || l >= num_p + num_p_sg, std::logic_error,
		     "Error!  Invalid p map index " << l);
  if (l < num_p)
    return me->get_p_init(l);
  else 
    return sg_p_init[l-num_p]->getBlockVector();

  return Teuchos::null;
}

Teuchos::RCP<Epetra_Operator>
Stokhos::SGModelEvaluator_Interlaced::create_W() const
{
  if (supports_x) {
    Teuchos::RCP<Teuchos::ParameterList> sgOpParams =
      Teuchos::rcp(&(params->sublist("SG Operator")), false);
    Teuchos::RCP<Epetra_CrsMatrix> W_crs;
    W_crs = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(me->create_W(), true);
    Teuchos::RCP<const Epetra_CrsGraph> W_graph =
      Teuchos::rcp(&(W_crs->Graph()), false);

    my_W = Teuchos::rcp(new Stokhos::InterlacedOperator(sg_comm, sg_basis, epetraCijk, W_graph,sgOpParams));
    my_W->setupOperator(W_sg_blocks);

    return my_W;
  }
  
  return Teuchos::null;
}

EpetraExt::ModelEvaluator::InArgs
Stokhos::SGModelEvaluator_Interlaced::createInArgs() const
{
  InArgsSetup inArgs;
  InArgs me_inargs = me->createInArgs();

  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p + num_p_sg); 
  inArgs.setSupports(IN_ARG_x_dot, me_inargs.supports(IN_ARG_x_dot_sg));
  inArgs.setSupports(IN_ARG_x, me_inargs.supports(IN_ARG_x_sg));
  inArgs.setSupports(IN_ARG_t, me_inargs.supports(IN_ARG_t));
  inArgs.setSupports(IN_ARG_alpha, me_inargs.supports(IN_ARG_alpha));
  inArgs.setSupports(IN_ARG_beta, me_inargs.supports(IN_ARG_beta));
  inArgs.setSupports(IN_ARG_sg_basis, me_inargs.supports(IN_ARG_sg_basis));
  inArgs.setSupports(IN_ARG_sg_quadrature, 
		     me_inargs.supports(IN_ARG_sg_quadrature));
  inArgs.setSupports(IN_ARG_sg_expansion, 
		     me_inargs.supports(IN_ARG_sg_expansion));
  
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
Stokhos::SGModelEvaluator_Interlaced::createOutArgs() const
{
  OutArgsSetup outArgs;
  OutArgs me_outargs = me->createOutArgs();

  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(num_p+num_p_sg, num_g_sg);
  outArgs.setSupports(OUT_ARG_f, me_outargs.supports(OUT_ARG_f_sg));
  outArgs.setSupports(OUT_ARG_W, me_outargs.supports(OUT_ARG_W_sg));
  outArgs.setSupports(OUT_ARG_WPrec, false);
  for (int j=0; j<num_p; j++)
    outArgs.setSupports(OUT_ARG_DfDp, j, 
			me_outargs.supports(OUT_ARG_DfDp_sg, j));
  for (int i=0; i<num_g_sg; i++) {
    int ii = sg_g_index_map[i];
//    if (!me_outargs.supports(OUT_ARG_DgDx_dot_sg, ii).none())
//      outArgs.setSupports(OUT_ARG_DgDx_dot, i,  DERIV_LINEAR_OP);
//    if (!me_outargs.supports(OUT_ARG_DgDx_sg, i).none())
//      outArgs.setSupports(OUT_ARG_DgDx, i,  DERIV_LINEAR_OP);
    for (int j=0; j<num_p; j++)
      outArgs.setSupports(OUT_ARG_DgDp, i, j, 
			  me_outargs.supports(OUT_ARG_DgDp_sg, ii, j));
  }

  // We do not support derivatives w.r.t. the new SG parameters, so their
  // support defaults to none.
  
  return outArgs;
}

void 
Stokhos::SGModelEvaluator_Interlaced::evalModel(const InArgs& inArgs,
				     const OutArgs& outArgs) const
{
  // Get the input arguments
  Teuchos::RCP<const Epetra_Vector> x;
  if (inArgs.supports(IN_ARG_x)) {
    x = inArgs.get_x();
    if (x != Teuchos::null)
      *my_x = *x;
  }
  Teuchos::RCP<const Epetra_Vector> x_dot;
  if (inArgs.supports(IN_ARG_x_dot))
    x_dot = inArgs.get_x_dot();

  // Get the output arguments
  EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> f_out;
  if (outArgs.supports(OUT_ARG_f))
    f_out = outArgs.get_f();
  Teuchos::RCP<Epetra_Operator> W_out;
  if (outArgs.supports(OUT_ARG_W))
    W_out = outArgs.get_W();

  // Create underlying inargs
  InArgs me_inargs = me->createInArgs();
  if (x != Teuchos::null) {
    Teuchos::RCP<Epetra_Vector> overlapped_x 
          = Teuchos::rcp(new Epetra_Vector(*interlace_overlapped_x_map));
    overlapped_x->Import(*x,*interlace_overlapped_x_importer,Insert);

    // x_sg_blocks->getBlockVector()->Import(*x, *interlace_overlapped_x_importer, 
    // 					  Insert);

    copyToPolyOrthogVector(*overlapped_x,*x_sg_blocks);
    me_inargs.set_x_sg(x_sg_blocks);
  }
  if (x_dot != Teuchos::null) {
    Teuchos::RCP<Epetra_Vector> overlapped_x_dot 
          = Teuchos::rcp(new Epetra_Vector(*interlace_overlapped_x_map));
    overlapped_x_dot->Import(*x_dot,*interlace_overlapped_x_importer,Insert);

    // x_dot_sg_blocks->getBlockVector()->Import(*x_dot, *interlace_overlapped_x_importer,
    // 					      Insert);

    copyToPolyOrthogVector(*overlapped_x_dot,*x_dot_sg_blocks);
    me_inargs.set_x_dot_sg(x_dot_sg_blocks);
  }
  if (me_inargs.supports(IN_ARG_alpha))
    me_inargs.set_alpha(inArgs.get_alpha());
  if (me_inargs.supports(IN_ARG_beta))
    me_inargs.set_beta(inArgs.get_beta());
  if (me_inargs.supports(IN_ARG_t))
    me_inargs.set_t(inArgs.get_t());
  if (me_inargs.supports(IN_ARG_sg_basis)) {
    if (inArgs.get_sg_basis() != Teuchos::null)
      me_inargs.set_sg_basis(inArgs.get_sg_basis());
    else
      me_inargs.set_sg_basis(sg_basis);
  }
  if (me_inargs.supports(IN_ARG_sg_quadrature)) {
    if (inArgs.get_sg_quadrature() != Teuchos::null)
      me_inargs.set_sg_quadrature(inArgs.get_sg_quadrature());
    else
      me_inargs.set_sg_quadrature(sg_quad);
  }
  if (me_inargs.supports(IN_ARG_sg_expansion)) {
    if (inArgs.get_sg_expansion() != Teuchos::null)
      me_inargs.set_sg_expansion(inArgs.get_sg_expansion());
    else
      me_inargs.set_sg_expansion(sg_exp);
  }

  // Pass parameters
  for (int i=0; i<num_p; i++)
    me_inargs.set_p(i, inArgs.get_p(i));
  for (int i=0; i<num_p_sg; i++) {
    Teuchos::RCP<const Epetra_Vector> p = inArgs.get_p(i+num_p);

    // We always need to pass in the SG parameters, so just use
    // the initial parameters if it is NULL
    if (p == Teuchos::null)
      p = sg_p_init[i]->getBlockVector();

    // Convert block p to SG polynomial
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> p_sg = 
      create_p_sg(sg_p_index_map[i], View, p.get());
    me_inargs.set_p_sg(sg_p_index_map[i], p_sg);
  }

  // Create underlying outargs
  OutArgs me_outargs = me->createOutArgs();

  // f
  if (f_out != Teuchos::null) {
    me_outargs.set_f_sg(f_sg_blocks);
    if (eval_W_with_f)
      me_outargs.set_W_sg(W_sg_blocks);
  }

  // W
  if (W_out != Teuchos::null && !eval_W_with_f )
     me_outargs.set_W_sg(W_sg_blocks);

  // df/dp -- deterministic p
  for (int i=0; i<outArgs.Np(); i++) {
    if (!outArgs.supports(OUT_ARG_DfDp, i).none()) {
      Derivative dfdp = outArgs.get_DfDp(i);
      if (dfdp.getMultiVector() != Teuchos::null) {
	Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> dfdp_sg;
	if (dfdp.getMultiVectorOrientation() == DERIV_MV_BY_COL)
	  dfdp_sg = 
	    Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			   sg_basis, overlapped_stoch_row_map, 
			   me->get_f_map(), interlace_overlapped_f_map, sg_comm, 
			   me->get_p_map(i)->NumMyElements()));
	else if (dfdp.getMultiVectorOrientation() == DERIV_TRANS_MV_BY_ROW)
	  dfdp_sg = 
	    Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			   sg_basis, overlapped_stoch_row_map, 
			   me->get_p_map(i), sg_comm, 
			   me->get_f_map()->NumMyElements()));
	me_outargs.set_DfDp_sg(i, 
			       SGDerivative(dfdp_sg,
					    dfdp.getMultiVectorOrientation()));
      }
      TEUCHOS_TEST_FOR_EXCEPTION(dfdp.getLinearOp() != Teuchos::null, std::logic_error,
			 "Error!  Stokhos::SGModelEvaluator_Interlaced::evalModel " << 
			 "cannot handle operator form of df/dp!");
    }
  }

  // Responses (g, dg/dx, dg/dp, ...)
  for (int i=0; i<num_g_sg; i++) {
    int ii = sg_g_index_map[i];

    // g
    Teuchos::RCP<Epetra_Vector> g = outArgs.get_g(i);
    if (g != Teuchos::null) {
      Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> g_sg =
	create_g_sg(sg_g_index_map[i], View, g.get());
      me_outargs.set_g_sg(i, g_sg);
    }

    // dg/dxdot
    if (outArgs.supports(OUT_ARG_DgDx_dot, i).supports(DERIV_LINEAR_OP)) {
      Derivative dgdx_dot = outArgs.get_DgDx_dot(i);
      if (dgdx_dot.getLinearOp() != Teuchos::null) {
	Teuchos::RCP<Stokhos::SGOperator> op =
	  Teuchos::rcp_dynamic_cast<Stokhos::SGOperator>(
	    dgdx_dot.getLinearOp(), true);
	Teuchos::RCP< Stokhos::EpetraOperatorOrthogPoly > sg_blocks =
	  op->getSGPolynomial();
	if (me_outargs.supports(OUT_ARG_DgDx, ii).supports(DERIV_LINEAR_OP))
	  me_outargs.set_DgDx_dot_sg(ii, sg_blocks);
	else {
	  for (unsigned int k=0; k<num_sg_blocks; k++) {
	    Teuchos::RCP<Epetra_MultiVector> mv = 
	      Teuchos::rcp_dynamic_cast<Stokhos::EpetraMultiVectorOperator>(
		sg_blocks->getCoeffPtr(k), true)->getMultiVector();
	    dgdx_dot_sg_blocks[i]->setCoeffPtr(k, mv);
	  }
	  if (me_outargs.supports(OUT_ARG_DgDx_dot_sg, ii).supports(DERIV_MV_BY_COL))
	    me_outargs.set_DgDx_dot_sg(ii, SGDerivative(dgdx_dot_sg_blocks[i],
					 	        DERIV_MV_BY_COL));
	  else
	    me_outargs.set_DgDx_dot_sg(ii, SGDerivative(dgdx_dot_sg_blocks[i],
						        DERIV_TRANS_MV_BY_ROW));
	}
      }
      TEUCHOS_TEST_FOR_EXCEPTION(dgdx_dot.getLinearOp() == Teuchos::null &&
			 dgdx_dot.isEmpty() == false, 
			 std::logic_error,
			 "Error!  Stokhos::SGModelEvaluator_Interlaced::evalModel: " << 
			 "Operator form of dg/dxdot is required!");
    }

    // dg/dx
    if (outArgs.supports(OUT_ARG_DgDx, i).supports(DERIV_LINEAR_OP)) {
      Derivative dgdx = outArgs.get_DgDx(i);
      if (dgdx.getLinearOp() != Teuchos::null) {
	Teuchos::RCP<Stokhos::SGOperator> op =
	  Teuchos::rcp_dynamic_cast<Stokhos::SGOperator>(
	    dgdx.getLinearOp(), true);
	Teuchos::RCP< Stokhos::EpetraOperatorOrthogPoly > sg_blocks =
	  op->getSGPolynomial();
	if (me_outargs.supports(OUT_ARG_DgDx, ii).supports(DERIV_LINEAR_OP))
	  me_outargs.set_DgDx_sg(i, sg_blocks);
	else {
	  for (unsigned int k=0; k<num_sg_blocks; k++) {
	    Teuchos::RCP<Epetra_MultiVector> mv = 
	      Teuchos::rcp_dynamic_cast<Stokhos::EpetraMultiVectorOperator>(
		sg_blocks->getCoeffPtr(k), true)->getMultiVector();
	    dgdx_sg_blocks[i]->setCoeffPtr(k, mv);
	  }
	  if (me_outargs.supports(OUT_ARG_DgDx_sg, ii).supports(DERIV_MV_BY_COL))
	    me_outargs.set_DgDx_sg(ii, SGDerivative(dgdx_sg_blocks[i],
						    DERIV_MV_BY_COL));
	  else
	    me_outargs.set_DgDx_sg(ii, SGDerivative(dgdx_sg_blocks[i],
						    DERIV_TRANS_MV_BY_ROW));
	}
      }
      TEUCHOS_TEST_FOR_EXCEPTION(dgdx.getLinearOp() == Teuchos::null &&
			 dgdx.isEmpty() == false, 
			 std::logic_error,
			 "Error!  Stokhos::SGModelEvaluator_Interlaced::evalModel: " << 
			 "Operator form of dg/dxdot is required!");
    }

    // dg/dp -- deterministic p
    // Rembember, no derivatives w.r.t. sg parameters
    for (int j=0; j<num_p; j++) {
      if (!outArgs.supports(OUT_ARG_DgDp, i, j).none()) {
	Derivative dgdp = outArgs.get_DgDp(i,j);
	if (dgdp.getMultiVector() != Teuchos::null) {
	  Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> dgdp_sg; 
	  if (dgdp.getMultiVectorOrientation() == DERIV_MV_BY_COL)
	    dgdp_sg = 
	      Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			     sg_basis, overlapped_stoch_row_map, 
			     me->get_g_map(ii), sg_g_map[i], sg_comm, 
			     View, *(dgdp.getMultiVector())));
	  else if (dgdp.getMultiVectorOrientation() == DERIV_TRANS_MV_BY_ROW) {
	    Teuchos::RCP<const Epetra_BlockMap> product_map =
	      Teuchos::rcp(&(dgdp.getMultiVector()->Map()),false);
	    dgdp_sg = 
	      Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			     sg_basis, overlapped_stoch_row_map, 
			     me->get_p_map(j), product_map, sg_comm, 
			     View, *(dgdp.getMultiVector())));
	  }
	  me_outargs.set_DgDp_sg(ii, j, 
				 SGDerivative(dgdp_sg,
					      dgdp.getMultiVectorOrientation()));
	}
	TEUCHOS_TEST_FOR_EXCEPTION(dgdp.getLinearOp() != Teuchos::null, 
			   std::logic_error,
			   "Error!  Stokhos::SGModelEvaluator_Interlaced::evalModel " << 
			   "cannot handle operator form of dg/dp!");
      }
    }

  }

  // Compute the functions
  me->evalModel(me_inargs, me_outargs);

  // Copy block SG components for W
  if ((W_out != Teuchos::null || (eval_W_with_f && f_out != Teuchos::null)) ) {
      
    Teuchos::RCP<Epetra_Operator> W;
    if (W_out != Teuchos::null)
      W = W_out;
    else
      W = my_W;
    Teuchos::RCP<Stokhos::SGOperator> W_sg = 
      Teuchos::rcp_dynamic_cast<Stokhos::SGOperator>(W, true);
    W_sg->setupOperator(W_sg_blocks);
  }

  // f
  if (f_out!=Teuchos::null){
    if (!scaleOP)
      for (int i=0; i<sg_basis->size(); i++)
	(*f_sg_blocks)[i].Scale(sg_basis->norm_squared(i));

    Teuchos::RCP<Epetra_Vector> overlapped_f 
          = Teuchos::rcp(new Epetra_Vector(*interlace_overlapped_f_map));
    copyToInterlacedVector(*f_sg_blocks,*overlapped_f);
    f_out->Export(*overlapped_f,*interlace_overlapped_f_exporter,Insert);

    // f_out->Export(*(f_sg_blocks->getBlockVector()), *interlace_overlapped_f_exporter, 
    // 		  Insert);
  }

  // df/dp -- determinsistic p
  for (int i=0; i<num_p; i++) {
    if (!outArgs.supports(OUT_ARG_DfDp, i).none()) {
      Derivative dfdp = outArgs.get_DfDp(i);
      SGDerivative dfdp_sg = me_outargs.get_DfDp_sg(i);
      if (dfdp.getMultiVector() != Teuchos::null) {
	dfdp.getMultiVector()->Export(
	  *(dfdp_sg.getMultiVector()->getBlockMultiVector()), 
	  *interlace_overlapped_f_exporter, Insert);
      }
    }
  }
}

void 
Stokhos::SGModelEvaluator_Interlaced::set_x_sg_init(
  const Stokhos::EpetraVectorOrthogPoly& x_sg_in)
{
  *sg_x_init = x_sg_in;
}

Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly>
Stokhos::SGModelEvaluator_Interlaced::get_x_sg_init() const
{
  return sg_x_init;
}

void 
Stokhos::SGModelEvaluator_Interlaced::set_p_sg_init(
  int i, const Stokhos::EpetraVectorOrthogPoly& p_sg_in)
{
  *sg_p_init[i] = p_sg_in;
}

Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly>
Stokhos::SGModelEvaluator_Interlaced::get_p_sg_init(int l) const
{
  return sg_p_init[l];
}

Teuchos::Array<int> 
Stokhos::SGModelEvaluator_Interlaced::get_p_sg_map_indices() const
{
  return sg_p_index_map;
}

Teuchos::Array<int> 
Stokhos::SGModelEvaluator_Interlaced::get_g_sg_map_indices() const
{
  return sg_g_index_map;
}

Teuchos::Array< Teuchos::RCP<const Epetra_Map> > 
Stokhos::SGModelEvaluator_Interlaced::get_g_sg_base_maps() const
{
  Teuchos::Array< Teuchos::RCP<const Epetra_Map> > base_maps(num_g);
  for (int i=0; i<num_g; i++)
    base_maps[i] = me->get_g_map(i);
  return base_maps;
 }

Teuchos::RCP<const Epetra_BlockMap> 
Stokhos::SGModelEvaluator_Interlaced::get_overlap_stochastic_map() const
{
  return overlapped_stoch_row_map;
}

Teuchos::RCP<const Epetra_BlockMap> 
Stokhos::SGModelEvaluator_Interlaced::get_x_sg_overlap_map() const
{
  return interlace_overlapped_x_map;
}

Teuchos::RCP<const Epetra_Import> 
Stokhos::SGModelEvaluator_Interlaced::get_x_sg_importer() const
{
  return interlace_overlapped_x_importer;
}

Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly>
Stokhos::SGModelEvaluator_Interlaced::create_x_sg(Epetra_DataAccess CV, 
				       const Epetra_Vector* v) const
{
  Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_x;
  if (v == NULL)
    sg_x = Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
			  sg_basis, stoch_row_map, x_map, get_x_map(), sg_comm));
  else
    sg_x = Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
			  sg_basis, stoch_row_map, x_map, get_x_map(), sg_comm,
			  CV, *v));
  return sg_x;
}

Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly>
Stokhos::SGModelEvaluator_Interlaced::create_x_sg_overlap(Epetra_DataAccess CV, 
					       const Epetra_Vector* v) const
{
  Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_x;
  if (v == NULL)
    sg_x = Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
			  sg_basis, overlapped_stoch_row_map, x_map, 
			  get_x_sg_overlap_map(), sg_comm));
  else
    sg_x = Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
			  sg_basis, overlapped_stoch_row_map, x_map, 
			  get_x_sg_overlap_map(), sg_comm, CV, *v));
  return sg_x;
}

Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly>
Stokhos::SGModelEvaluator_Interlaced::create_x_mv_sg(int num_vecs, Epetra_DataAccess CV, 
					  const Epetra_MultiVector* v) const
{
  Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> sg_x;
  if (v == NULL)
    sg_x = Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			  sg_basis, stoch_row_map, x_map, get_x_map(), sg_comm,
			  num_vecs));
  else
    sg_x = Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			  sg_basis, stoch_row_map, x_map, get_x_map(), sg_comm,
			  CV, *v));
  return sg_x;
}

Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly>
Stokhos::SGModelEvaluator_Interlaced::create_x_mv_sg_overlap(
  int num_vecs, 
  Epetra_DataAccess CV, 
  const Epetra_MultiVector* v) const
{
  Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> sg_x;
  if (v == NULL)
    sg_x = Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			  sg_basis, overlapped_stoch_row_map, x_map, 
			  get_x_sg_overlap_map(), sg_comm, num_vecs));
  else
    sg_x = Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			  sg_basis, overlapped_stoch_row_map, x_map, 
			  get_x_sg_overlap_map(), sg_comm, CV, *v));
  return sg_x;
}

Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly>
Stokhos::SGModelEvaluator_Interlaced::create_p_sg(int l, Epetra_DataAccess CV, 
				       const Epetra_Vector* v) const
{
  Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_p;
  Teuchos::Array<int>::const_iterator it = std::find(sg_p_index_map.begin(),
						     sg_p_index_map.end(), 
						     l);
  TEUCHOS_TEST_FOR_EXCEPTION(it == sg_p_index_map.end(), std::logic_error,
		     "Error!  Invalid p map index " << l);
  int ll = it - sg_p_index_map.begin();
  if (v == NULL)
    sg_p = Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
			  sg_basis, overlapped_stoch_p_map, me->get_p_map(l),
			  sg_p_map[ll], sg_comm));
  else
    sg_p = Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
			  sg_basis, overlapped_stoch_p_map, me->get_p_map(l),
			  sg_p_map[ll], sg_comm, CV, *v));
  return sg_p;
}

Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly>
Stokhos::SGModelEvaluator_Interlaced::create_f_sg(Epetra_DataAccess CV, 
				       const Epetra_Vector* v) const
{
  Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_f;
  if (v == NULL)
    sg_f = Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
			  sg_basis, stoch_row_map, f_map, interlace_f_map, sg_comm));
  else
    sg_f = Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
			  sg_basis, stoch_row_map, f_map, interlace_f_map, sg_comm,
			  CV, *v));
  return sg_f;
}

Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly>
Stokhos::SGModelEvaluator_Interlaced::create_f_sg_overlap(Epetra_DataAccess CV, 
					       const Epetra_Vector* v) const
{
  Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_f;
  if (v == NULL)
    sg_f = Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
			  sg_basis, overlapped_stoch_row_map, f_map, 
			  interlace_overlapped_f_map, sg_comm));
  else
    sg_f = Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
			  sg_basis, overlapped_stoch_row_map, f_map, 
			  interlace_overlapped_f_map, sg_comm, CV, *v));
  return sg_f;
}

Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly>
Stokhos::SGModelEvaluator_Interlaced::create_f_mv_sg(
  int num_vecs, 
  Epetra_DataAccess CV, 
  const Epetra_MultiVector* v) const
{
  Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> sg_f;
  if (v == NULL)
    sg_f = Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			  sg_basis, stoch_row_map, f_map, interlace_f_map, sg_comm,
			  num_vecs));
  else
    sg_f = Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			  sg_basis, stoch_row_map, f_map, interlace_f_map, sg_comm,
			  CV, *v));
  return sg_f;
}

Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly>
Stokhos::SGModelEvaluator_Interlaced::create_f_mv_sg_overlap(
  int num_vecs, 
  Epetra_DataAccess CV, 
  const Epetra_MultiVector* v) const
{
  Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> sg_f;
  if (v == NULL)
    sg_f = Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			  sg_basis, overlapped_stoch_row_map, f_map, 
			  interlace_overlapped_f_map, sg_comm, num_vecs));
  else
    sg_f = Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			  sg_basis, overlapped_stoch_row_map, f_map, 
			  interlace_overlapped_f_map, sg_comm, CV, *v));
  return sg_f;
}

Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly>
Stokhos::SGModelEvaluator_Interlaced::create_g_sg(int l, Epetra_DataAccess CV, 
				       const Epetra_Vector* v) const
{
  Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_g;
  Teuchos::Array<int>::const_iterator it = std::find(sg_g_index_map.begin(),
						     sg_g_index_map.end(), 
						     l);
  TEUCHOS_TEST_FOR_EXCEPTION(it == sg_g_index_map.end(), std::logic_error,
		     "Error!  Invalid g map index " << l);
  int ll = it - sg_g_index_map.begin();
  if (v == NULL)
    sg_g = Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
			  sg_basis, overlapped_stoch_row_map, 
			  me->get_g_map(l), 
			  sg_g_map[ll], sg_comm));
  else
    sg_g = Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
			  sg_basis, overlapped_stoch_row_map, 
			  me->get_g_map(l), 
			  sg_g_map[ll], sg_comm, CV, *v));
  return sg_g;
}

Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly>
Stokhos::SGModelEvaluator_Interlaced::create_g_mv_sg(int l, int num_vecs,
					  Epetra_DataAccess CV, 
					  const Epetra_MultiVector* v) const
{
  Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> sg_g;
  Teuchos::Array<int>::const_iterator it = std::find(sg_g_index_map.begin(),
						     sg_g_index_map.end(), 
						     l);
  TEUCHOS_TEST_FOR_EXCEPTION(it == sg_g_index_map.end(), std::logic_error,
		     "Error!  Invalid g map index " << l);
  int ll = it - sg_g_index_map.begin();
  if (v == NULL)
    sg_g = Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			  sg_basis, overlapped_stoch_row_map, 
			  me->get_g_map(l), 
			  sg_g_map[ll], sg_comm, num_vecs));
  else
    sg_g = Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			  sg_basis, overlapped_stoch_row_map, 
			  me->get_g_map(l), 
			  sg_g_map[ll], sg_comm, CV, *v));
  return sg_g;
}

Teuchos::RCP<Epetra_Map>
Stokhos::SGModelEvaluator_Interlaced::buildInterlaceMap(const Epetra_BlockMap & determ_map,const Epetra_BlockMap & stocha_map)
{
   int stochaUnks = stocha_map.NumMyElements();
   int determUnks = determ_map.NumMyElements();

   // these must be equal for the "adapt" model evaluator
   TEUCHOS_ASSERT(stocha_map.NumGlobalElements()==stochaUnks);

   // build interlaced unknowns
   std::vector<int> interlacedUnks(stochaUnks*determUnks);
   std::size_t i=0;
   for(int d=0;d<determUnks;d++)
      for(int s=0;s<stochaUnks;s++,i++)
         interlacedUnks[i] = determ_map.GID(d)*stochaUnks+s;

   return Teuchos::rcp(new Epetra_Map(-1,interlacedUnks.size(),&interlacedUnks[0],0,determ_map.Comm()));
}

void Stokhos::SGModelEvaluator_Interlaced::copyToInterlacedVector(const Stokhos::EpetraVectorOrthogPoly & x_sg,
                                                             Epetra_Vector & x)
{
   std::size_t numBlocks = x_sg.size();
   Teuchos::RCP<const EpetraExt::BlockVector> bv_x = x_sg.getBlockVector();

   // loop over all blocks
   for(std::size_t blk=0;blk<numBlocks;blk++) { 
      const Epetra_Vector & v = *bv_x->GetBlock(blk);
      
      for(int dof=0;dof<v.MyLength();dof++)
         x[dof*numBlocks+blk] = v[dof]; 
   }
}

void Stokhos::SGModelEvaluator_Interlaced::copyToPolyOrthogVector(const Epetra_Vector & x,
                                                                  Stokhos::EpetraVectorOrthogPoly & x_sg)
{
   std::size_t numBlocks = x_sg.size();
   Teuchos::RCP<EpetraExt::BlockVector> bv_x = x_sg.getBlockVector();

   // loop over all blocks
   for(std::size_t blk=0;blk<numBlocks;blk++) { 
      Epetra_Vector & v = *bv_x->GetBlock(blk);
      
      for(int dof=0;dof<v.MyLength();dof++)
         v[dof] = x[dof*numBlocks+blk];
   }
}
