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
#ifdef HAVE_MPI
#include "EpetraExt_MultiMpiComm.h"
#else
#include "EpetraExt_MultiSerialComm.h"
#endif
#include "EpetraExt_BlockUtility.h"
#include "EpetraExt_BlockCrsMatrix.h"
#include "EpetraExt_BlockMultiVector.h"
#include "Stokhos_MatrixFreeEpetraOp.hpp"
#include "Stokhos_KLMatrixFreeEpetraOp.hpp"
#include "Stokhos_KLReducedMatrixFreeEpetraOp.hpp"
#include "Stokhos_MeanEpetraOp.hpp"
#include "Stokhos_IfpackPreconditionerFactory.hpp"
#include "Stokhos_MLPreconditionerFactory.hpp"
#include "Stokhos_EpetraMultiVectorOperator.hpp"

Stokhos::SGModelEvaluator::SGModelEvaluator(
  const Teuchos::RCP<EpetraExt::ModelEvaluator>& me_,
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis_,
  const Teuchos::RCP<const Stokhos::Quadrature<int,double> >& sg_quad_,
  const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> >& sg_exp_,
  const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_,
  const Teuchos::RCP<const Epetra_Comm>& comm,
  const Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly>& initial_x_sg,
  const Teuchos::Array< Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> >& initial_p_sg) 
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
    sg_comm(),
    sg_x_map(),
    sg_f_map(),
    num_p(0),
    num_p_sg(0),
    sg_p_map(),
    sg_p_names(),
    num_g(0),
    num_g_sg(0),
    sg_g_map(),
    Cijk(Cijk_),
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
    jacobianMethod(MATRIX_FREE),
    sg_p_init(),
    eval_W_with_f(false)
{
  if (x_map != Teuchos::null)
    supports_x = true;

  Teuchos::RCP<EpetraExt::MultiComm> multiComm;
#ifdef HAVE_MPI
  // No parallelism over blocks, so spatial partition is unchanged 
  // as comm->NumProc()
  multiComm =
    Teuchos::rcp(new EpetraExt::MultiMpiComm(MPI_COMM_WORLD, 
					     comm->NumProc(), 
					     num_sg_blocks));
#else
  multiComm = Teuchos::rcp(new EpetraExt::MultiSerialComm(num_sg_blocks));
#endif
  sg_comm = multiComm;
  int numBlockRows =  multiComm->NumTimeSteps();
  int myBlockRows  =  multiComm->NumTimeStepsOnDomain();
  int myFirstBlockRow = multiComm->FirstTimeStepOnDomain();

  // DENSE STENCIL for Stochastic Galerkin
  // For 3 blocks on 2 procs, this should be:
  // Proc  nBR  mBR  mFBR     Stencil      Index
  //  0     3    2    0       0  1  2        0
  //                         -1  0  1        1
  //  1     3    1    2      -2 -1  0        2
  //
  rowStencil.resize(myBlockRows);
  rowIndex.resize(myBlockRows);
  for (int i=0; i < myBlockRows; i++) {
    for (int j=0; j < numBlockRows; j++) 
      rowStencil[i].push_back(-myFirstBlockRow - i + j);
    rowIndex[i] = (i + myFirstBlockRow);
  }

  if (supports_x) {

    // Create block SG x and f maps
    sg_x_map = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(*x_map,
							     rowIndex,
							     *sg_comm));
    sg_f_map = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(*f_map,
							     rowIndex,
							     *sg_comm));
    
    // Create vector blocks
    x_dot_sg_blocks = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(sg_basis));
    x_sg_blocks = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(sg_basis));
    f_sg_blocks = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(sg_basis));

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

    // Get block Jacobian method
    std::string method = params->get("Jacobian Method", "Matrix Free");
    if (method == "Matrix Free")
      jacobianMethod = MATRIX_FREE;
    else if (method == "KL Matrix Free")
      jacobianMethod = KL_MATRIX_FREE;
    else if (method == "KL Reduced Matrix Free")
      jacobianMethod = KL_REDUCED_MATRIX_FREE;
    else if (method == "Fully Assembled")
      jacobianMethod = FULLY_ASSEMBLED;
    else
      TEST_FOR_EXCEPTION(true, 
			 std::logic_error,
			 std::endl << 
			 "Error!  Stokhos::SGModelEvaluator():  " <<
			 "Unknown Jacobian Method " << method << "!" << 
			 std::endl);
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
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(*p_map,
							     rowIndex,
							     *sg_comm));
    
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
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(*g_map,
							     rowIndex,
							     *sg_comm));
    
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
    for (int i=0; i<num_p_sg; i++)
      dfdp_sg_blocks[i] = 
	Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(sg_basis));

    // Create initial x
    sg_x_init = Teuchos::rcp(new EpetraExt::BlockVector(*x_map, *sg_x_map));
    if (initial_x_sg != Teuchos::null) {
      // Use supplied initial guess
      initial_x_sg->assignToBlockVector(*sg_x_init);
    }
    else {
      Teuchos::RCP<const EpetraVectorOrthogPoly> init_x_sg = 
	me->get_x_sg_init();
      if (init_x_sg != Teuchos::null) {
	// Use initial guess provided by underlying model evaluator
	init_x_sg->assignToBlockVector(*sg_x_init);
      }
      else {
	// Use default choice of deterministic initial guess for mean
	sg_x_init->LoadBlockValues(*(me->get_x_init()), 0);
      }
    }
    
    // Get preconditioner factory for matrix-free
    if (jacobianMethod == MATRIX_FREE || 
	jacobianMethod == KL_MATRIX_FREE || 
	jacobianMethod == KL_REDUCED_MATRIX_FREE) {
      std::string prec_name = 
	params->get("Mean Preconditioner Type", "Ifpack");
      Teuchos::RCP<Teuchos::ParameterList> precParams = 
	Teuchos::rcp(&params->sublist("Preconditioner Parameters"),false);
      if (prec_name == "Ifpack")
	precFactory = 
	  Teuchos::rcp(new Stokhos::IfpackPreconditionerFactory(precParams));
      else if (prec_name == "ML")
	precFactory = 
	  Teuchos::rcp(new Stokhos::MLPreconditionerFactory(precParams));
      else
	TEST_FOR_EXCEPTION(true, std::logic_error,
			   "Error!  Unknown preconditioner type " << prec_name
			   << ".  Valid choices are \"Ifpack\" and \"ML\".");
    }
  }
  
  eval_W_with_f = params->get("Evaluate W with F", false);
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
    if (jacobianMethod == MATRIX_FREE)
      my_W = Teuchos::rcp(new Stokhos::MatrixFreeEpetraOp(x_map, f_map,
							  sg_x_map, sg_f_map, 
							  sg_basis, Cijk, 
							  W_sg_blocks));
    else if (jacobianMethod == KL_MATRIX_FREE)
      my_W = Teuchos::rcp(new Stokhos::KLMatrixFreeEpetraOp(
			    x_map, sg_x_map, sg_basis, Cijk, 
			    W_sg_blocks->getCoefficients()));
    else if (jacobianMethod == KL_REDUCED_MATRIX_FREE) {
      int num_KL = params->get<int>("Number of KL Terms");
      my_W = Teuchos::rcp(new Stokhos::KLReducedMatrixFreeEpetraOp(
			    x_map, sg_x_map, sg_basis, Cijk, W_sg_blocks,
			    num_KL));
    }
    else if (jacobianMethod == FULLY_ASSEMBLED) {
      Teuchos::RCP<Epetra_Operator> W = me->create_W();
      Teuchos::RCP<Epetra_RowMatrix> W_row = 
	Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(W, true);
      my_W =
	Teuchos::rcp(new EpetraExt::BlockCrsMatrix(*W_row, rowStencil, 
						   rowIndex, *sg_comm));
    }

    return my_W;
  }
  
  return Teuchos::null;
}

Teuchos::RCP<Epetra_Operator>
Stokhos::SGModelEvaluator::create_M() const
{
  if (supports_x) {
    if (jacobianMethod == MATRIX_FREE || 
	jacobianMethod == KL_MATRIX_FREE || 
	jacobianMethod == KL_REDUCED_MATRIX_FREE)
      return Teuchos::rcp(new Stokhos::MeanEpetraOp(x_map, sg_x_map, 
						    num_sg_blocks, 
						    Teuchos::null));
    else if (jacobianMethod == FULLY_ASSEMBLED) {
      Teuchos::RCP<Epetra_Operator> W = me->create_W();
      Teuchos::RCP<Epetra_RowMatrix> W_row = 
	Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(W, true);
      return
	Teuchos::rcp(new EpetraExt::BlockCrsMatrix(*W_row, rowStencil, 
						   rowIndex, *sg_comm));
    }
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

    return Teuchos::rcp(new Stokhos::MatrixFreeEpetraOp(
			  x_map, me->get_g_sg_map(jj), sg_x_map, 
			  sg_g_map[jj], sg_basis, Cijk, sg_blocks));
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

    return Teuchos::rcp(new Stokhos::MatrixFreeEpetraOp(
			  x_map, me->get_g_sg_map(jj), sg_x_map, 
			  sg_g_map[jj], sg_basis, Cijk, sg_blocks));
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
  if (jacobianMethod == MATRIX_FREE ||
      jacobianMethod == KL_MATRIX_FREE || 
      jacobianMethod == KL_REDUCED_MATRIX_FREE)
    outArgs.setSupports(OUT_ARG_M, me_outargs.supports(OUT_ARG_W));
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
  // Get the output arguments
  EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> f_out;
  if (outArgs.supports(OUT_ARG_f))
    f_out = outArgs.get_f();
  Teuchos::RCP<Epetra_Operator> W_out;
  if (outArgs.supports(OUT_ARG_W))
    W_out = outArgs.get_W();
  Teuchos::RCP<Epetra_Operator> M_out;
  if (outArgs.supports(OUT_ARG_M))
    M_out = outArgs.get_M();

  // Check if we are using the "matrix-free" method for W and we are 
  // computing a preconditioner.  
  bool eval_mean = (W_out == Teuchos::null && M_out != Teuchos::null);

  // Here we are assuming a full W fill occurred previously which we can use
  // for the preconditioner.  Given the expense of computing the SG W blocks
  // this saves significant computational cost
  if (eval_mean) {
    Teuchos::RCP<Stokhos::MeanEpetraOp> W_mean = 
      Teuchos::rcp_dynamic_cast<Stokhos::MeanEpetraOp>(M_out, true);
    Teuchos::RCP<Epetra_Operator> prec = 
      precFactory->compute(W_sg_blocks->getCoeffPtr(0));    
    W_mean->setMeanOperator(prec);

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

  // Get the input arguments
  Teuchos::RCP<const Epetra_Vector> x;
  if (inArgs.supports(IN_ARG_x))
    x = inArgs.get_x();
  Teuchos::RCP<const Epetra_Vector> x_dot;
  if (inArgs.supports(IN_ARG_x_dot))
    x_dot = inArgs.get_x_dot();

  // Create underlying inargs
  InArgs me_inargs = me->createInArgs();
  if (x != Teuchos::null) {
    x_sg_blocks->resetCoefficients(View, *x_map, *x);
    me_inargs.set_x_sg(x_sg_blocks);
  }
  if (x_dot != Teuchos::null) {
    me_inargs.set_x_dot_sg(x_dot_sg_blocks);
    x_dot_sg_blocks->resetCoefficients(View, *x_map, *x_dot);
  }
  if (me_inargs.supports(IN_ARG_alpha))
    me_inargs.set_alpha(inArgs.get_alpha());
  if (me_inargs.supports(IN_ARG_beta))
    me_inargs.set_beta(inArgs.get_beta());
  if (me_inargs.supports(IN_ARG_t))
    me_inargs.set_t(inArgs.get_t());
  if (me_inargs.supports(IN_ARG_sg_basis))
    me_inargs.set_sg_basis(sg_basis);
  if (me_inargs.supports(IN_ARG_sg_quadrature))
    me_inargs.set_sg_quadrature(sg_quad);
  if (me_inargs.supports(IN_ARG_sg_expansion))
    me_inargs.set_sg_expansion(sg_exp);

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
    f_sg_blocks->resetCoefficients(View, *f_map, *f_out);
    me_outargs.set_f_sg(f_sg_blocks);
    if (eval_W_with_f)
      me_outargs.set_W_sg(W_sg_blocks);
  }

  // W
  if (W_out != Teuchos::null && !eval_W_with_f && !eval_mean)
     me_outargs.set_W_sg(W_sg_blocks);

  // df/dp
  for (int i=0; i<outArgs.Np(); i++) {
    if (!outArgs.supports(OUT_ARG_DfDp, i).none()) {
      Derivative dfdp = outArgs.get_DfDp(i);
      if (dfdp.getMultiVector() != Teuchos::null) {
	dfdp_sg_blocks[i]->resetCoefficients(View, *(me->get_f_map()), 
					     *(dfdp.getMultiVector()));
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
	Teuchos::RCP<Stokhos::MatrixFreeEpetraOp> op =
	  Teuchos::rcp_dynamic_cast<Stokhos::MatrixFreeEpetraOp>(
	    dgdx_dot.getLinearOp(), true);
	Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_Operator> > sg_blocks =
	  op->getOperatorBlocks();
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
	Teuchos::RCP<Stokhos::MatrixFreeEpetraOp> op =
	  Teuchos::rcp_dynamic_cast<Stokhos::MatrixFreeEpetraOp>(
	    dgdx.getLinearOp(), true);
	Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_Operator> > sg_blocks =
	  op->getOperatorBlocks();
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

  // Copy block SG components for f and W into f and W
  if ((W_out != Teuchos::null || (eval_W_with_f && f_out != Teuchos::null)) 
      && !eval_mean) {
    Teuchos::RCP<Epetra_Operator> W;
    if (W_out != Teuchos::null)
      W = W_out;
    else
      W = my_W;
    if (jacobianMethod == MATRIX_FREE ||
	jacobianMethod == KL_MATRIX_FREE || 
	jacobianMethod == KL_REDUCED_MATRIX_FREE) {
      if (jacobianMethod == MATRIX_FREE) {
	Teuchos::RCP<Stokhos::MatrixFreeEpetraOp> W_mf = 
	  Teuchos::rcp_dynamic_cast<Stokhos::MatrixFreeEpetraOp>(W, true);
	W_mf->reset(W_sg_blocks);
      }
      else if (jacobianMethod == KL_MATRIX_FREE) {
	Teuchos::RCP<Stokhos::KLMatrixFreeEpetraOp> W_mf = 
	  Teuchos::rcp_dynamic_cast<Stokhos::KLMatrixFreeEpetraOp>(W, true);
	W_mf->reset(W_sg_blocks->getCoefficients());
      }
      else if (jacobianMethod == KL_REDUCED_MATRIX_FREE) {
	Teuchos::RCP<Stokhos::KLReducedMatrixFreeEpetraOp> W_mf = 
	  Teuchos::rcp_dynamic_cast<Stokhos::KLReducedMatrixFreeEpetraOp>(W, 
									  true);
	W_mf->reset(W_sg_blocks);
      }

      if (M_out != Teuchos::null) {
	Teuchos::RCP<Stokhos::MeanEpetraOp> W_mean = 
	  Teuchos::rcp_dynamic_cast<Stokhos::MeanEpetraOp>(M_out, true);
	Teuchos::RCP<Epetra_Operator> prec = 
	  precFactory->compute(W_sg_blocks->getCoeffPtr(0));    
	W_mean->setMeanOperator(prec);
      }
    }
    else if (jacobianMethod == FULLY_ASSEMBLED) {
      Teuchos::RCP<EpetraExt::BlockCrsMatrix> W_sg = 
	Teuchos::rcp_dynamic_cast<EpetraExt::BlockCrsMatrix>(W, true);
      Teuchos::RCP<Epetra_RowMatrix> W_row;
      int i,j;
      double cijk;
      W_sg->PutScalar(0.0);
      for (unsigned int k=0; k<num_W_blocks; k++) {
	W_row = 
	  Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(
	    W_sg_blocks->getCoeffPtr(k), true);
	int nl = Cijk->num_values(k);
	for (int l=0; l<nl; l++) {
	  Cijk->value(k,l,i,j,cijk);
	  W_sg->SumIntoBlock(cijk/sg_basis->norm_squared(i), *W_row, i, j);
	}
      }
    }
  }
}

void 
Stokhos::SGModelEvaluator::set_x_init(const Epetra_Vector& x_in)
{
  sg_x_init->Scale(1.0, x_in);
}
