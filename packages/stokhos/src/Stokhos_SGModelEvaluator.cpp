// $Id$ 
// $Source$ 
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

#include "Stokhos_SGModelEvaluator.hpp"

#include <algorithm>
#include "Teuchos_TestForException.hpp"
#include "EpetraExt_BlockUtility.h"
#include "EpetraExt_BlockMultiVector.h"
#include "Stokhos_SGOperatorFactory.hpp"
#include "Stokhos_SGPreconditionerFactory.hpp"
#include "Stokhos_MatrixFreeOperator.hpp"
#include "Stokhos_EpetraMultiVectorOperator.hpp"
#include "Epetra_LocalMap.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"

Stokhos::SGModelEvaluator::SGModelEvaluator(
  const Teuchos::RCP<EpetraExt::ModelEvaluator>& me_,
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis_,
  const Teuchos::RCP<const Stokhos::Quadrature<int,double> >& sg_quad_,
  const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> >& sg_exp_,
  const Teuchos::RCP<const Stokhos::ParallelData>& sg_parallel_data_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_,
  const Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly>& initial_x_sg,
  const Teuchos::Array< Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> >& initial_p_sg,
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
    sg_x_map(),
    sg_overlapped_x_map(),
    sg_f_map(),
    sg_overlapped_f_map(),
    sg_overlapped_x_importer(),
    sg_overlapped_f_exporter(),
    sg_overlapped_x(),
    sg_overlapped_x_dot(),
    sg_overlapped_f(),
    num_p(0),
    num_p_sg(0),
    sg_p_map(),
    sg_p_names(),
    num_g(0),
    num_g_sg(0),
    sg_g_map(),
    x_dot_sg_blocks(),
    x_sg_blocks(),
    p_sg_blocks(),
    f_sg_blocks(),
    W_sg_blocks(),
    g_sg_blocks(),
    dfdp_sg_blocks(),
    dgdx_dot_sg_blocks(),
    dgdx_sg_blocks(),
    dgdp_sg_blocks(),
    sg_p_init(),
    eval_W_with_f(false),
    scaleOP(scaleOP_)
{
  if (x_map != Teuchos::null)
    supports_x = true;

  Teuchos::RCP<const Epetra_BlockMap> stoch_row_map =
    epetraCijk->getStochasticRowMap();
  Teuchos::RCP<const Epetra_BlockMap> overlapped_stoch_row_map = 
    Teuchos::rcp(new Epetra_LocalMap(
		   num_sg_blocks, 0, *(sg_parallel_data->getStochasticComm())));

  if (supports_x) {

    // Create block SG x and f maps
    sg_x_map = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
		     *x_map, *stoch_row_map, *sg_comm));
    sg_overlapped_x_map = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
		     *x_map, *overlapped_stoch_row_map, *sg_comm));

    sg_f_map = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
		     *f_map, *stoch_row_map, *sg_comm));
    sg_overlapped_f_map = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
		     *f_map, *overlapped_stoch_row_map, *sg_comm));

    // Create importer/exporter from/to overlapped distribution
    sg_overlapped_x_importer = 
      Teuchos::rcp(new Epetra_Import(*sg_overlapped_x_map, *sg_x_map));
    sg_overlapped_f_exporter = 
      Teuchos::rcp(new Epetra_Export(*sg_overlapped_f_map, *sg_f_map));

    // Create overlapped vectors
    sg_overlapped_x = Teuchos::rcp(new Epetra_Vector(*sg_overlapped_x_map));
    sg_overlapped_x_dot = Teuchos::rcp(new Epetra_Vector(*sg_overlapped_x_map));
    sg_overlapped_f = Teuchos::rcp(new Epetra_Vector(*sg_overlapped_f_map));
    
    // Create vector blocks
    x_dot_sg_blocks = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(sg_basis));
    x_sg_blocks = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(sg_basis));
    f_sg_blocks = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(sg_basis));

    // Preconditioner needs an x
    my_x = Teuchos::rcp(new Epetra_Vector(*sg_x_map));

    // Determine W expansion type
    std::string W_expansion_type = 
      params->get("Jacobian Expansion Type", "Full");
    if (W_expansion_type == "Linear")
      num_W_blocks = sg_basis->dimension() + 1;
    else
      num_W_blocks = num_sg_blocks;
    W_sg_blocks = 
      Teuchos::rcp(new Stokhos::VectorOrthogPoly<Epetra_Operator>(sg_basis,
								  num_W_blocks));
    for (unsigned int i=0; i<num_W_blocks; i++)
      W_sg_blocks->setCoeffPtr(i, me->create_W());
  }
    
  // Parameters -- The idea here is to add new parameter vectors
  // for the stochastic Galerkin components of the parameters

  InArgs me_inargs = me->createInArgs();
  num_p = me_inargs.Np();
  num_p_sg = me_inargs.Np_sg();
  sg_p_map.resize(num_p_sg);
  sg_p_names.resize(num_p_sg);
  sg_p_init.resize(num_p_sg);
  p_sg_blocks.resize(num_p_sg);

  // Determine parameter expansion type
  std::string p_expansion_type = 
    params->get("Parameter Expansion Type", "Full");
  if (p_expansion_type == "Linear")
    num_p_blocks = sg_basis->dimension() + 1;
  else
    num_p_blocks = num_sg_blocks;
  
  // Create parameter maps, names, and initial values
  for (int i=0; i<num_p_sg; i++) {
    Teuchos::RCP<const Epetra_Map> p_map = me->get_p_sg_map(i);
    sg_p_map[i] = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
		     *p_map, *overlapped_stoch_row_map, *sg_comm));
    
    Teuchos::RCP<const Teuchos::Array<std::string> > p_names = 
      me->get_p_sg_names(i);
    if (p_names != Teuchos::null) {
      sg_p_names[i] = 
	Teuchos::rcp(new Teuchos::Array<std::string>(num_sg_blocks*(p_names->size())));
      for (int j=0; j<p_names->size(); j++) {
	std::stringstream ss;
	ss << (*p_names)[j] << " -- SG Coefficient " << i;
	(*sg_p_names[i])[j] = ss.str();
      }
    }
    
    // Create initial p
    sg_p_init[i] = 
      Teuchos::rcp(new EpetraExt::BlockVector(*p_map, *sg_p_map[i]));
    if (initial_p_sg.size() == num_p_sg && initial_p_sg[i] != Teuchos::null) {
      // Use provided initial parameters
      initial_p_sg[i]->assignToBlockVector(*sg_p_init[i]);
    }
    else {
      Teuchos::RCP<const EpetraVectorOrthogPoly> init_p_sg = 
	me->get_p_sg_init(i);
      if (init_p_sg != Teuchos::null) {
	// Use initial parameters provided by underlying model evaluator
	init_p_sg->assignToBlockVector(*sg_p_init[i]);
      }
      else {
	TEST_FOR_EXCEPTION(true, std::logic_error,
			   "Error!  Must provide initial parameters through " <<
			   "constructor or undelying model evaluator's " <<
			   "get_p_sg_init() method.");
      }
    }
    
    // Create p SG blocks
    p_sg_blocks[i] = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(sg_basis,
						       num_p_blocks));
    
  }

  // Responses -- The idea here is to add new parameter vectors
  // for the stochastic Galerkin components of the respones

  // Get number of SG parameters model supports derivatives w.r.t.
  OutArgs me_outargs = me->createOutArgs();
  num_g = me_outargs.Ng();
  num_g_sg = me_outargs.Ng_sg();
  sg_g_map.resize(num_g_sg);
  g_sg_blocks.resize(num_g_sg);
  dgdx_dot_sg_blocks.resize(num_g_sg);
  dgdx_sg_blocks.resize(num_g_sg);
  dgdp_sg_blocks.resize(num_g_sg);

  // Create response maps
  for (int i=0; i<num_g_sg; i++) {
    Teuchos::RCP<const Epetra_Map> g_map = me->get_g_sg_map(i);
    sg_g_map[i] = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
		     *g_map, *overlapped_stoch_row_map, *sg_comm));
    
    // Create g SG blocks
    g_sg_blocks[i] = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(sg_basis));
    
    // Create dg/dxdot, dg/dx SG blocks
    dgdx_dot_sg_blocks[i] = 
      Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(sg_basis));
    dgdx_sg_blocks[i] = 
      Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(sg_basis));
    
    // Create dg/dp SG blocks
    dgdp_sg_blocks[i].resize(num_p_sg);
    for (int j=0; j<num_p_sg; j++)
      dgdp_sg_blocks[i][j] = 
	Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(sg_basis));
  }

  
  if (supports_x) {

    // Create df/dp SG blocks
    dfdp_sg_blocks.resize(num_p_sg);
    sg_overlapped_dfdp.resize(num_p_sg);
    for (int i=0; i<num_p_sg; i++)
      dfdp_sg_blocks[i] = 
	Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(sg_basis));

    // Create initial x
    sg_x_init = Teuchos::rcp(new EpetraExt::BlockVector(*x_map, *sg_x_map));
    if (initial_x_sg != Teuchos::null) {
      // Use supplied initial guess
      EpetraExt::BlockVector sg_overlapped_x_init(*x_map, 
						  *sg_overlapped_x_map);
      initial_x_sg->assignToBlockVector(sg_overlapped_x_init);
      sg_x_init->Export(sg_overlapped_x_init, *sg_overlapped_x_importer, 
			Insert);
    }
    else {
      Teuchos::RCP<const EpetraVectorOrthogPoly> init_x_sg = 
	me->get_x_sg_init();
      if (init_x_sg != Teuchos::null) {
	// Use initial guess provided by underlying model evaluator
	EpetraExt::BlockVector sg_overlapped_x_init(*x_map, 
						    *sg_overlapped_x_map);
	init_x_sg->assignToBlockVector(sg_overlapped_x_init);
	sg_x_init->Export(sg_overlapped_x_init, *sg_overlapped_x_importer, 
			  Insert);
      }
      else {
	// Use default choice of deterministic initial guess for mean
	if (stoch_row_map->MyGID(0))
	  sg_x_init->LoadBlockValues(*(me->get_x_init()), 0);
      }
    }
   
  }
  
  eval_W_with_f = params->get("Evaluate W with F", false);

  // We don't support parallel for dgdx yet, so build a new EpetraCijk
  serialCijk =
    Teuchos::rcp(new Stokhos::EpetraSparse3Tensor(sg_basis, 
						  epetraCijk->getCijk(), 
						  sg_comm,
						  overlapped_stoch_row_map));
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
Stokhos::SGModelEvaluator::get_x_map() const
{
  return sg_x_map;
}

Teuchos::RCP<const Epetra_Map>
Stokhos::SGModelEvaluator::get_f_map() const
{
  return sg_f_map;
}

Teuchos::RCP<const Epetra_Map>
Stokhos::SGModelEvaluator::get_p_map(int l) const
{
  if (l < num_p)
    return me->get_p_map(l);
  else if (l < num_p + num_p_sg)
    return sg_p_map[l-num_p];
  else
    TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  Invalid p map index " << l);
  return Teuchos::null;
}

Teuchos::RCP<const Epetra_Map>
Stokhos::SGModelEvaluator::get_g_map(int l) const
{
  if (l < num_g)
    return me->get_g_map(l);
  else if (l < num_g + num_g_sg)
    return sg_g_map[l-num_g];
  else
    TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  Invalid g map index " << l);
  return Teuchos::null;
}

Teuchos::RCP<const Teuchos::Array<std::string> >
Stokhos::SGModelEvaluator::get_p_names(int l) const
{
  if (l < num_p)
    return me->get_p_names(l);
  else if (l < num_p + num_p_sg)
    return sg_p_names[l-num_p];
  else
    TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  Invalid p names index " << l);
  return Teuchos::null;
}

Teuchos::RCP<const Epetra_Vector>
Stokhos::SGModelEvaluator::get_x_init() const
{
  return sg_x_init;
}

Teuchos::RCP<const Epetra_Vector>
Stokhos::SGModelEvaluator::get_p_init(int l) const
{
  if (l < num_p)
    return me->get_p_init(l);
  else if (l < num_p + num_p_sg)
    return sg_p_init[l-num_p];
  else
    TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  Invalid p names index " << l);
  return Teuchos::null;
}

Teuchos::RCP<Epetra_Operator>
Stokhos::SGModelEvaluator::create_W() const
{
  if (supports_x) {
    Teuchos::RCP<Teuchos::ParameterList> sgOpParams =
      Teuchos::rcp(&(params->sublist("SG Operator")), false);
    Teuchos::RCP<Epetra_CrsMatrix> W_crs;
    if (sgOpParams->get("Operator Method", "Matrix Free") == 
	"Fully Assembled") {
      W_crs = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(me->create_W(), true);
      Teuchos::RCP<const Epetra_CrsGraph> W_graph =
	Teuchos::rcp(&(W_crs->Graph()), false);
      sgOpParams->set< Teuchos::RCP<const Epetra_CrsGraph> >("Base Graph", 
	W_graph);
    }
    Stokhos::SGOperatorFactory sg_op_factory(sgOpParams);
    my_W = sg_op_factory.build(sg_comm, sg_basis, epetraCijk, x_map, f_map, 
			       sg_x_map, sg_f_map);
    my_W->setupOperator(W_sg_blocks);

    return my_W;
  }
  
  return Teuchos::null;
}

Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner>
Stokhos::SGModelEvaluator::create_WPrec() const
{
  if (supports_x) {
    Teuchos::RCP<Teuchos::ParameterList> sgPrecParams =
      Teuchos::rcp(&(params->sublist("SG Preconditioner")), false);
    Stokhos::SGPreconditionerFactory sg_prec_factory(sgPrecParams);
    Teuchos::RCP<Epetra_Operator> precOp = 
      sg_prec_factory.build(sg_comm, sg_basis, epetraCijk, x_map, sg_x_map);
    return Teuchos::rcp(new EpetraExt::ModelEvaluator::Preconditioner(precOp,
								      true));
  }
  
  return Teuchos::null;
}

Teuchos::RCP<Epetra_Operator>
Stokhos::SGModelEvaluator::create_DgDx_op(int j) const
{
  if (j < num_g)
    return me->create_DgDx_op(j);
  else if (j < num_g + num_g_sg && supports_x) {
    int jj = j-num_g;
    Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_Operator> > sg_blocks =
      Teuchos::rcp(new Stokhos::VectorOrthogPoly<Epetra_Operator>(sg_basis));
    OutArgs me_outargs = me->createOutArgs();
    if (me_outargs.supports(OUT_ARG_DgDx_sg, jj).supports(DERIV_LINEAR_OP))
      for (unsigned int i=0; i<num_sg_blocks; i++)
	sg_blocks->setCoeffPtr(i, me->create_DgDx_sg_op(jj));
    else if (me_outargs.supports(OUT_ARG_DgDx_sg, jj).supports(DERIV_MV_BY_COL))
      for (unsigned int i=0; i<num_sg_blocks; i++) {
	Teuchos::RCP<Epetra_MultiVector> mv = 
	  Teuchos::rcp(new Epetra_MultiVector(*(me->get_g_sg_map(jj)), 
					      me->get_x_map()->NumMyElements()));
	Teuchos::RCP<Epetra_Operator> block = 
	  Teuchos::rcp(new Stokhos::EpetraMultiVectorOperator(mv));
	sg_blocks->setCoeffPtr(i, block);
      }
    else if (me_outargs.supports(OUT_ARG_DgDx_sg, jj).supports(DERIV_TRANS_MV_BY_ROW))
      for (unsigned int i=0; i<num_sg_blocks; i++) {
	Teuchos::RCP<Epetra_MultiVector> mv = 
	  Teuchos::rcp(new Epetra_MultiVector(*(me->get_x_map()), 
					      me->get_g_sg_map(jj)->NumMyElements()));
	Teuchos::RCP<Epetra_Operator> block = 
	  Teuchos::rcp(new Stokhos::EpetraMultiVectorOperator(mv));
	sg_blocks->setCoeffPtr(i, block);
      }
    else
      TEST_FOR_EXCEPTION(true, std::logic_error,
			 "Error!  me_outargs.supports(OUT_ARG_DgDx_sg, " << j
			 << ").none() is true!");

    Teuchos::RCP<Teuchos::ParameterList> pl = 
      Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::RCP<Stokhos::SGOperator> dgdx_sg = 
      Teuchos::rcp(new Stokhos::MatrixFreeOperator(
		     sg_comm, sg_basis, serialCijk, x_map, me->get_g_sg_map(jj),
		     sg_x_map, sg_g_map[jj], pl));
    dgdx_sg->setupOperator(sg_blocks);
    return dgdx_sg;
  }
  else 
    TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  DgDx_op is not supported for index " << j
		       << "!");
  
  return Teuchos::null;
}

Teuchos::RCP<Epetra_Operator>
Stokhos::SGModelEvaluator::create_DgDx_dot_op(int j) const
{
  if (j < num_g)
    return me->create_DgDx_dot_op(j);
  else if (j < num_g + num_g_sg && supports_x) {
    int jj = j-num_g;
    Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_Operator> > sg_blocks =
      Teuchos::rcp(new Stokhos::VectorOrthogPoly<Epetra_Operator>(sg_basis));
    OutArgs me_outargs = me->createOutArgs();
    if (me_outargs.supports(OUT_ARG_DgDx_dot_sg, jj).supports(DERIV_LINEAR_OP))
      for (unsigned int i=0; i<num_sg_blocks; i++)
	sg_blocks->setCoeffPtr(i, me->create_DgDx_dot_sg_op(jj));
    else if (me_outargs.supports(OUT_ARG_DgDx_dot_sg, jj).supports(DERIV_MV_BY_COL))
      for (unsigned int i=0; i<num_sg_blocks; i++) {
	Teuchos::RCP<Epetra_MultiVector> mv = 
	  Teuchos::rcp(new Epetra_MultiVector(*(me->get_g_sg_map(jj)), 
					      me->get_x_map()->NumMyElements()));
	Teuchos::RCP<Epetra_Operator> block = 
	  Teuchos::rcp(new Stokhos::EpetraMultiVectorOperator(mv));
	sg_blocks->setCoeffPtr(i, block);
      }
    else if (me_outargs.supports(OUT_ARG_DgDx_dot_sg, jj).supports(DERIV_TRANS_MV_BY_ROW))
      for (unsigned int i=0; i<num_sg_blocks; i++) {
	Teuchos::RCP<Epetra_MultiVector> mv = 
	  Teuchos::rcp(new Epetra_MultiVector(*(me->get_x_map()), 
					      me->get_g_sg_map(jj)->NumMyElements()));
	
	Teuchos::RCP<Epetra_Operator> block = 
	  Teuchos::rcp(new Stokhos::EpetraMultiVectorOperator(mv));
	sg_blocks->setCoeffPtr(i, block);
      }
    else
      TEST_FOR_EXCEPTION(true, std::logic_error,
			 "Error!  me_outargs.supports(OUT_ARG_DgDx_dot_sg, " 
			 << j << ").none() is true!");

    Teuchos::RCP<Teuchos::ParameterList> pl = 
      Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::RCP<Stokhos::SGOperator> dgdx_dot_sg = 
      Teuchos::rcp(new Stokhos::MatrixFreeOperator(
		     sg_comm, sg_basis, serialCijk, x_map, me->get_g_sg_map(jj),
		     sg_x_map, sg_g_map[jj], pl));
    dgdx_dot_sg->setupOperator(sg_blocks);
  }
  else 
    TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  DgDx_dot_op is not supported for index " << j
		       << "!");

  return Teuchos::null;
}

EpetraExt::ModelEvaluator::InArgs
Stokhos::SGModelEvaluator::createInArgs() const
{
  InArgsSetup inArgs;
  InArgs me_inargs = me->createInArgs();

  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p + num_p_sg); 
  inArgs.set_Np_sg(0);
  inArgs.setSupports(IN_ARG_x_dot, me_inargs.supports(IN_ARG_x_dot));
  inArgs.setSupports(IN_ARG_x, me_inargs.supports(IN_ARG_x));
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
Stokhos::SGModelEvaluator::createOutArgs() const
{
  OutArgsSetup outArgs;
  OutArgs me_outargs = me->createOutArgs();

  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(me_outargs.Np() + num_p_sg, num_g + num_g_sg);
  outArgs.set_Np_Ng_sg(0, 0);
  outArgs.setSupports(OUT_ARG_f, me_outargs.supports(OUT_ARG_f));
  outArgs.setSupports(OUT_ARG_W, me_outargs.supports(OUT_ARG_W));
  outArgs.setSupports(OUT_ARG_WPrec, me_outargs.supports(OUT_ARG_W));
  for (int j=0; j<me_outargs.Np(); j++)
    outArgs.setSupports(OUT_ARG_DfDp, j, 
			me_outargs.supports(OUT_ARG_DfDp, j));
  for (int i=0; i<num_g; i++) {
    outArgs.setSupports(OUT_ARG_DgDx_dot, i,
			me_outargs.supports(OUT_ARG_DgDx_dot, i));
    outArgs.setSupports(OUT_ARG_DgDx, i,
			me_outargs.supports(OUT_ARG_DgDx, i));
    for (int j=0; j<me_outargs.Np(); j++)
      outArgs.setSupports(OUT_ARG_DgDp, i, j, 
			  me_outargs.supports(OUT_ARG_DgDp, i, j));
  }
  for (int i=0; i<num_g_sg; i++) {
    if (!me_outargs.supports(OUT_ARG_DgDx_dot_sg, i).none())
      outArgs.setSupports(OUT_ARG_DgDx_dot, i+num_g,  DERIV_LINEAR_OP);
    if (!me_outargs.supports(OUT_ARG_DgDx_sg, i).none())
      outArgs.setSupports(OUT_ARG_DgDx, i+num_g,  DERIV_LINEAR_OP);
    for (int j=0; j<me_outargs.Np(); j++)
      outArgs.setSupports(OUT_ARG_DgDp, i+num_g, j, 
			  me_outargs.supports(OUT_ARG_DgDp_sg, i, j));
  }

  // We do not support derivatives w.r.t. the new SG parameters, so their
  // support defaults to none.
  
  return outArgs;
}

void 
Stokhos::SGModelEvaluator::evalModel(const InArgs& inArgs,
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
  Teuchos::RCP<Epetra_Operator> WPrec_out;
  if (outArgs.supports(OUT_ARG_WPrec))
    WPrec_out = outArgs.get_WPrec();

  // Check if we are using the "matrix-free" method for W and we are 
  // computing a preconditioner.  
  bool eval_prec = (W_out == Teuchos::null && WPrec_out != Teuchos::null);

  // Here we are assuming a full W fill occurred previously which we can use
  // for the preconditioner.  Given the expense of computing the SG W blocks
  // this saves significant computational cost
  if (eval_prec) {
    Teuchos::RCP<Stokhos::SGPreconditioner> W_prec = 
      Teuchos::rcp_dynamic_cast<Stokhos::SGPreconditioner>(WPrec_out, true);
    W_prec->setupPreconditioner(my_W, *my_x);
    
    // We can now quit unless a fill of f, g, or dg/dp was also requested
    bool done = (f_out == Teuchos::null);
    for (int i=0; i<outArgs.Ng(); i++) {
      done = done && (outArgs.get_g(i) == Teuchos::null);
      for (int j=0; j<outArgs.Np(); j++)
	if (!outArgs.supports(OUT_ARG_DgDp, i, j).none())
	done = done && (outArgs.get_DgDp(i,j).isEmpty());
    }
    if (done)
      return;
  }

  // Create underlying inargs
  InArgs me_inargs = me->createInArgs();
  if (x != Teuchos::null) {
    sg_overlapped_x->Import(*x, *sg_overlapped_x_importer, Insert);
    x_sg_blocks->resetCoefficients(View, *x_map, *sg_overlapped_x);
    me_inargs.set_x_sg(x_sg_blocks);
  }
  if (x_dot != Teuchos::null) {
    sg_overlapped_x_dot->Import(*x_dot, *sg_overlapped_x_importer, Insert);
    x_dot_sg_blocks->resetCoefficients(View, *x_map, *sg_overlapped_x_dot);
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
      p = sg_p_init[i];

    // Convert block p to SG polynomial
    p_sg_blocks[i]->resetCoefficients(View, sg_p_init[i]->GetBaseMap(), *p);
    me_inargs.set_p_sg(i, p_sg_blocks[i]);
  }

  // Create underlying outargs
  OutArgs me_outargs = me->createOutArgs();

  // f
  if (f_out != Teuchos::null) {
    f_sg_blocks->resetCoefficients(View, *f_map, *sg_overlapped_f);
    me_outargs.set_f_sg(f_sg_blocks);
    if (eval_W_with_f)
      me_outargs.set_W_sg(W_sg_blocks);
  }

  // W
  if (W_out != Teuchos::null && !eval_W_with_f && !eval_prec)
     me_outargs.set_W_sg(W_sg_blocks);

  // df/dp
  for (int i=0; i<outArgs.Np(); i++) {
    if (!outArgs.supports(OUT_ARG_DfDp, i).none()) {
      Derivative dfdp = outArgs.get_DfDp(i);
      if (dfdp.getMultiVector() != Teuchos::null) {
	if (sg_overlapped_dfdp[i] == Teuchos::null)
	  sg_overlapped_dfdp[i] = 
	    Teuchos::rcp(new Epetra_MultiVector(
			   *sg_overlapped_f_map, 
			   dfdp.getMultiVector()->NumVectors()));
	
	dfdp_sg_blocks[i]->resetCoefficients(View, *(me->get_f_map()), 
					     *sg_overlapped_dfdp[i]);
	me_outargs.set_DfDp_sg(i, 
			       SGDerivative(dfdp_sg_blocks[i],
					    dfdp.getMultiVectorOrientation()));
      }
      TEST_FOR_EXCEPTION(dfdp.getLinearOp() != Teuchos::null, std::logic_error,
			 "Error!  Stokhos::SGModelEvaluator::evalModel " << 
			 "cannot handle operator form of df/dp!");
    }
  }

  // Responses (g, dg/dx, dg/dp, ...)
  for (int i=0; i<num_g; i++) {
    // g
    me_outargs.set_g(i, outArgs.get_g(i));

    // dg/dxdot
    if (!outArgs.supports(OUT_ARG_DgDx_dot, i).none())
      me_outargs.set_DgDx_dot(i, outArgs.get_DgDx_dot(i));

    // dg/dx
    if (!outArgs.supports(OUT_ARG_DgDx, i).none())
      me_outargs.set_DgDx(i, outArgs.get_DgDx(i));

    // dg/dp
    for (int j=0; j<outArgs.Np(); j++)
      if (!outArgs.supports(OUT_ARG_DgDp, i, j).none())
	me_outargs.set_DgDp(i, j, outArgs.get_DgDp(i,j));
  }
  for (int i=0; i<num_g_sg; i++) {
    // g
    Teuchos::RCP<Epetra_Vector> g = outArgs.get_g(i+num_g);
    if (g != Teuchos::null) {
      g_sg_blocks[i]->resetCoefficients(View, *(me->get_g_sg_map(i)), *g);
      me_outargs.set_g_sg(i, g_sg_blocks[i]);
    }

    // dg/dxdot
    if (outArgs.supports(OUT_ARG_DgDx_dot, i+num_g).supports(DERIV_LINEAR_OP)) {
      Derivative dgdx_dot = outArgs.get_DgDx_dot(i+num_g);
      if (dgdx_dot.getLinearOp() != Teuchos::null) {
	Teuchos::RCP<Stokhos::SGOperator> op =
	  Teuchos::rcp_dynamic_cast<Stokhos::SGOperator>(
	    dgdx_dot.getLinearOp(), true);
	Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_Operator> > sg_blocks =
	  op->getSGPolynomial();
	if (me_outargs.supports(OUT_ARG_DgDx, i).supports(DERIV_LINEAR_OP))
	  me_outargs.set_DgDx_dot_sg(i, sg_blocks);
	else {
	  for (unsigned int k=0; k<num_sg_blocks; k++) {
	    Teuchos::RCP<Epetra_MultiVector> mv = 
	      Teuchos::rcp_dynamic_cast<Stokhos::EpetraMultiVectorOperator>(
		sg_blocks->getCoeffPtr(k), true)->getMultiVector();
	    dgdx_dot_sg_blocks[i]->setCoeffPtr(k, mv);
	  }
	  if (me_outargs.supports(OUT_ARG_DgDx_dot_sg, i).supports(DERIV_MV_BY_COL))
	    me_outargs.set_DgDx_dot_sg(i, SGDerivative(dgdx_dot_sg_blocks[i],
						       DERIV_MV_BY_COL));
	  else
	    me_outargs.set_DgDx_dot_sg(i, SGDerivative(dgdx_dot_sg_blocks[i],
						       DERIV_TRANS_MV_BY_ROW));
	}
      }
      TEST_FOR_EXCEPTION(dgdx_dot.getLinearOp() == Teuchos::null &&
			 dgdx_dot.isEmpty() == false, 
			 std::logic_error,
			 "Error!  Stokhos::SGModelEvaluator::evalModel: " << 
			 "Operator form of dg/dxdot is required!");
    }

    // dg/dx
    if (outArgs.supports(OUT_ARG_DgDx, i+num_g).supports(DERIV_LINEAR_OP)) {
      Derivative dgdx = outArgs.get_DgDx(i+num_g);
      if (dgdx.getLinearOp() != Teuchos::null) {
	Teuchos::RCP<Stokhos::SGOperator> op =
	  Teuchos::rcp_dynamic_cast<Stokhos::SGOperator>(
	    dgdx.getLinearOp(), true);
	Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_Operator> > sg_blocks =
	  op->getSGPolynomial();
	if (me_outargs.supports(OUT_ARG_DgDx, i).supports(DERIV_LINEAR_OP))
	  me_outargs.set_DgDx_sg(i, sg_blocks);
	else {
	  for (unsigned int k=0; k<num_sg_blocks; k++) {
	    Teuchos::RCP<Epetra_MultiVector> mv = 
	      Teuchos::rcp_dynamic_cast<Stokhos::EpetraMultiVectorOperator>(
		sg_blocks->getCoeffPtr(k), true)->getMultiVector();
	    dgdx_sg_blocks[i]->setCoeffPtr(k, mv);
	  }
	  if (me_outargs.supports(OUT_ARG_DgDx_sg, i).supports(DERIV_MV_BY_COL))
	    me_outargs.set_DgDx_sg(i, SGDerivative(dgdx_sg_blocks[i],
						   DERIV_MV_BY_COL));
	  else
	    me_outargs.set_DgDx_sg(i, SGDerivative(dgdx_sg_blocks[i],
						   DERIV_TRANS_MV_BY_ROW));
	}
      }
      TEST_FOR_EXCEPTION(dgdx.getLinearOp() == Teuchos::null &&
			 dgdx.isEmpty() == false, 
			 std::logic_error,
			 "Error!  Stokhos::SGModelEvaluator::evalModel: " << 
			 "Operator form of dg/dxdot is required!");
    }

    // dg/dp
    // Rembember, no derivatives w.r.t. sg parameters
    for (int j=0; j<outArgs.Np(); j++) {
      if (!outArgs.supports(OUT_ARG_DgDp, i+num_g, j).none()) {
	Derivative dgdp = outArgs.get_DgDp(i+num_g,j);
	if (dgdp.getMultiVector() != Teuchos::null) {
	  dgdp_sg_blocks[i][j]->resetCoefficients(
	    View, *(me->get_g_sg_map(i)), *(dgdp.getMultiVector()));
	  me_outargs.set_DgDp_sg(i, j, 
				 SGDerivative(dgdp_sg_blocks[i][j],
					      dgdp.getMultiVectorOrientation()));
	}
	TEST_FOR_EXCEPTION(dgdp.getLinearOp() != Teuchos::null, 
			   std::logic_error,
			   "Error!  Stokhos::SGModelEvaluator::evalModel " << 
			   "cannot handle operator form of dg/dp!");
      }
    }

  }

  // Compute the functions
  me->evalModel(me_inargs, me_outargs);

  // Copy block SG components for W
  if ((W_out != Teuchos::null || (eval_W_with_f && f_out != Teuchos::null)) 
      && !eval_prec) {
    Teuchos::RCP<Epetra_Operator> W;
    if (W_out != Teuchos::null)
      W = W_out;
    else
      W = my_W;
    Teuchos::RCP<Stokhos::SGOperator> W_sg = 
      Teuchos::rcp_dynamic_cast<Stokhos::SGOperator>(W, true);
    W_sg->setupOperator(W_sg_blocks);
      
    if (WPrec_out != Teuchos::null) {
      Teuchos::RCP<Stokhos::SGPreconditioner> W_prec = 
	Teuchos::rcp_dynamic_cast<Stokhos::SGPreconditioner>(WPrec_out, true);
      W_prec->setupPreconditioner(W_sg, *my_x);
    }
  }

  // f
  if (f_out!=Teuchos::null){
    if (!scaleOP)
      for (int i=0; i<sg_basis->size(); i++)
	(*f_sg_blocks)[i].Scale(sg_basis->norm_squared(i));
    f_out->Export(*sg_overlapped_f, *sg_overlapped_f_exporter, Insert);
  }

  // df/dp
  for (int i=0; i<outArgs.Np(); i++) {
    if (!outArgs.supports(OUT_ARG_DfDp, i).none()) {
      Derivative dfdp = outArgs.get_DfDp(i);
      if (dfdp.getMultiVector() != Teuchos::null) {
	dfdp.getMultiVector()->Export(*sg_overlapped_dfdp[i], 
				      *sg_overlapped_f_exporter, Insert);
      }
    }
  }
}

void 
Stokhos::SGModelEvaluator::set_x_init(const Epetra_Vector& x_in)
{
  sg_x_init->Scale(1.0, x_in);
}

Teuchos::Array<int> 
Stokhos::SGModelEvaluator::get_p_sg_indices() const
{
  Teuchos::Array<int> sg_p_index(num_p_sg);
  for (int i=0; i<num_p_sg; i++)
    sg_p_index[i] = num_p + i;
  return sg_p_index;
}

Teuchos::Array<int> 
Stokhos::SGModelEvaluator::get_non_p_sg_indices() const
{
  Teuchos::Array<int> non_sg_p_index(num_p);
  for (int i=0; i<num_p; i++)
    non_sg_p_index[i] = i;
  return non_sg_p_index;
}

Teuchos::Array<int> 
Stokhos::SGModelEvaluator::get_g_sg_indices() const
{
  Teuchos::Array<int> sg_g_index(num_g_sg);
  for (int i=0; i<num_g_sg; i++)
    sg_g_index[i] = num_g + i;
  return sg_g_index;
}

Teuchos::Array<int> 
Stokhos::SGModelEvaluator::get_non_g_sg_indices() const
{
  Teuchos::Array<int> non_sg_g_index(num_g);
  for (int i=0; i<num_g; i++)
    non_sg_g_index[i] = i;
  return non_sg_g_index;
}

Teuchos::Array< Teuchos::RCP<const Epetra_Map> > 
Stokhos::SGModelEvaluator::get_p_sg_base_maps() const
{
  Teuchos::Array< Teuchos::RCP<const Epetra_Map> > base_maps(num_p_sg);
  for (int i=0; i<num_p_sg; i++)
    base_maps[i] = me->get_p_sg_map(i);
  return base_maps;
}

Teuchos::Array< Teuchos::RCP<const Epetra_Map> > 
Stokhos::SGModelEvaluator::get_g_sg_base_maps() const
{
  Teuchos::Array< Teuchos::RCP<const Epetra_Map> > base_maps(num_g_sg);
  for (int i=0; i<num_g_sg; i++)
    base_maps[i] = me->get_g_sg_map(i);
  return base_maps;
}

Teuchos::RCP<EpetraExt::BlockVector> 
Stokhos::SGModelEvaluator::import_solution(const Epetra_Vector& x) const
{
  Teuchos::RCP<EpetraExt::BlockVector> x_overlapped = 
    Teuchos::rcp(new EpetraExt::BlockVector(*x_map, *sg_overlapped_x_map));
  x_overlapped->Import(x, *sg_overlapped_x_importer, Insert);
  return x_overlapped;
}

Teuchos::RCP<EpetraExt::BlockVector> 
Stokhos::SGModelEvaluator::export_solution(const Epetra_Vector& x_overlapped) const
{
  Teuchos::RCP<EpetraExt::BlockVector> x = 
    Teuchos::rcp(new EpetraExt::BlockVector(*x_map, *sg_x_map));
  x->Export(x_overlapped, *sg_overlapped_x_importer, Insert);
  return x;
}

Teuchos::RCP<EpetraExt::BlockVector> 
Stokhos::SGModelEvaluator::import_residual(const Epetra_Vector& f) const
{
  Teuchos::RCP<EpetraExt::BlockVector> f_overlapped = 
    Teuchos::rcp(new EpetraExt::BlockVector(*f_map, *sg_overlapped_f_map));
  f_overlapped->Import(f, *sg_overlapped_f_exporter, Insert);
  return f_overlapped;
}

Teuchos::RCP<EpetraExt::BlockVector> 
Stokhos::SGModelEvaluator::export_residual(const Epetra_Vector& f_overlapped) const
{
  Teuchos::RCP<EpetraExt::BlockVector> f = 
    Teuchos::rcp(new EpetraExt::BlockVector(*f_map, *sg_f_map));
  f->Export(f_overlapped, *sg_overlapped_f_exporter, Insert);
  return f;
}
