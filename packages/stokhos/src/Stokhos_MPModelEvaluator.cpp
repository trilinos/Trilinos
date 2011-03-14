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

#include "Stokhos_MPModelEvaluator.hpp"

#include <algorithm>
#include "Teuchos_TestForException.hpp"
#include "EpetraExt_BlockUtility.h"
#include "EpetraExt_BlockMultiVector.h"
#include "Stokhos_MPPreconditionerFactory.hpp"
#include "Stokhos_EpetraMultiVectorOperator.hpp"


Stokhos::MPModelEvaluator::MPModelEvaluator(
  const Teuchos::RCP<EpetraExt::ModelEvaluator>& me_,
  const Teuchos::RCP<const EpetraExt::MultiComm>& mp_comm_,
  const Teuchos::RCP<const Epetra_Map>& mp_block_map_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_) 
  : me(me_),
    num_mp_blocks(mp_block_map_->NumMyElements()),
    mp_comm(mp_comm_),
    mp_block_map(mp_block_map_),
    params(params_),
    supports_x(false),
    x_map(me->get_x_map()),
    f_map(me->get_f_map()),
    mp_x_map(),
    mp_f_map(),
    num_p(0),
    num_p_mp(0),
    mp_p_map(),
    mp_p_names(),
    num_g(0),
    num_g_mp(0),
    mp_g_map(),
    W_mp_blocks(),
    dgdx_dot_mp_blocks(),
    dgdx_mp_blocks(),
    mp_p_init(),
    my_W(),
    my_x()
{
  if (x_map != Teuchos::null)
    supports_x = true;
  
  if (supports_x) {

    // Create block MP x and f maps
    mp_x_map = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
		     *x_map, *mp_block_map, *mp_comm));

    mp_f_map = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
		     *f_map, *mp_block_map, *mp_comm));

    // Create default mp_x_init
    mp_x_init = Teuchos::rcp(new ProductEpetraVector(
			       mp_block_map, x_map, mp_x_map, mp_comm));
    for (unsigned int i=0; i<num_mp_blocks; i++)
      (*mp_x_init)[i] = *(me->get_x_init());

    // Preconditioner needs an x
    my_x = Teuchos::rcp(new Epetra_Vector(*mp_x_map));

    W_mp_blocks = 
      Teuchos::rcp(new Stokhos::ProductContainer<Epetra_Operator>(
		     mp_block_map));
    for (unsigned int i=0; i<num_mp_blocks; i++)
      W_mp_blocks->setCoeffPtr(i, me->create_W());
  }
    
  // Parameters -- The idea here is to add new parameter vectors
  // for the stochastic Galerkin components of the parameters

  InArgs me_inargs = me->createInArgs();
  num_p = me_inargs.Np();
  num_p_mp = me_inargs.Np_mp();
  mp_p_map.resize(num_p_mp);
  mp_p_names.resize(num_p_mp);
  mp_p_init.resize(num_p_mp);
  
  // Create parameter maps, names, and initial values
  for (int i=0; i<num_p_mp; i++) {
    Teuchos::RCP<const Epetra_Map> p_map = me->get_p_mp_map(i);
    mp_p_map[i] = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
		     *p_map, *mp_block_map, *mp_comm));
    
    Teuchos::RCP<const Teuchos::Array<std::string> > p_names = 
      me->get_p_mp_names(i);
    if (p_names != Teuchos::null) {
      mp_p_names[i] = 
	Teuchos::rcp(new Teuchos::Array<std::string>(num_mp_blocks*(p_names->size())));
      for (int j=0; j<p_names->size(); j++) {
	std::stringstream ss;
	ss << (*p_names)[j] << " -- MP Coefficient " << i;
	(*mp_p_names[i])[j] = ss.str();
      }
    }
    
    // Create default mp_p_init
    mp_p_init[i] = 
      Teuchos::rcp(new ProductEpetraVector(
		     mp_block_map, p_map, mp_p_map[i], mp_comm));
    mp_p_init[i]->init(0.0);
    // for (unsigned int j=0; j<num_mp_blocks; j++)
    //   (*mp_p_init[i])[j] = *(me->get_p_init(i));
    
  }

  // Responses -- The idea here is to add new parameter vectors
  // for the stochastic Galerkin components of the respones

  // Get number of MP parameters model supports derivatives w.r.t.
  OutArgs me_outargs = me->createOutArgs();
  num_g = me_outargs.Ng();
  num_g_mp = me_outargs.Ng_mp();
  mp_g_map.resize(num_g_mp);
  dgdx_dot_mp_blocks.resize(num_g_mp);
  dgdx_mp_blocks.resize(num_g_mp);

  // Create response maps
  for (int i=0; i<num_g_mp; i++) {
    Teuchos::RCP<const Epetra_Map> g_map = me->get_g_mp_map(i);
    mp_g_map[i] = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
		     *g_map, *mp_block_map, *mp_comm));
    
    // Create dg/dxdot, dg/dx MP blocks
    dgdx_dot_mp_blocks[i] = 
      Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(mp_block_map));
    dgdx_mp_blocks[i] = 
      Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(mp_block_map));
  }
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
Stokhos::MPModelEvaluator::get_x_map() const
{
  return mp_x_map;
}

Teuchos::RCP<const Epetra_Map>
Stokhos::MPModelEvaluator::get_f_map() const
{
  return mp_f_map;
}

Teuchos::RCP<const Epetra_Map>
Stokhos::MPModelEvaluator::get_p_map(int l) const
{
  if (l < num_p)
    return me->get_p_map(l);
  else if (l < num_p + num_p_mp)
    return mp_p_map[l-num_p];
  else
    TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  Invalid p map index " << l);
  return Teuchos::null;
}

Teuchos::RCP<const Epetra_Map>
Stokhos::MPModelEvaluator::get_g_map(int l) const
{
  if (l < num_g)
    return me->get_g_map(l);
  else if (l < num_g + num_g_mp)
    return mp_g_map[l-num_g];
  else
    TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  Invalid g map index " << l);
  return Teuchos::null;
}

Teuchos::RCP<const Teuchos::Array<std::string> >
Stokhos::MPModelEvaluator::get_p_names(int l) const
{
  if (l < num_p)
    return me->get_p_names(l);
  else if (l < num_p + num_p_mp)
    return mp_p_names[l-num_p];
  else
    TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  Invalid p names index " << l);
  return Teuchos::null;
}

Teuchos::RCP<const Epetra_Vector>
Stokhos::MPModelEvaluator::get_x_init() const
{
  return mp_x_init->getBlockVector();
}

Teuchos::RCP<const Epetra_Vector>
Stokhos::MPModelEvaluator::get_p_init(int l) const
{
  if (l < num_p)
    return me->get_p_init(l);
  else if (l < num_p + num_p_mp)
    return mp_p_init[l-num_p]->getBlockVector();
  else
    TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  Invalid p names index " << l);
  return Teuchos::null;
}

Teuchos::RCP<Epetra_Operator>
Stokhos::MPModelEvaluator::create_W() const
{
  if (supports_x) {
    my_W = 
      Teuchos::rcp(new Stokhos::BlockDiagonalOperator(mp_comm, num_mp_blocks, 
						      x_map, f_map, 
						      mp_x_map, mp_f_map));
    my_W->setupOperator(W_mp_blocks);
    
    return my_W;
  }
  
  return Teuchos::null;
}

Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner>
Stokhos::MPModelEvaluator::create_WPrec() const
{
  if (supports_x) {
    Teuchos::RCP<Teuchos::ParameterList> mpPrecParams =
      Teuchos::rcp(&(params->sublist("MP Preconditioner")), false);
    Stokhos::MPPreconditionerFactory mp_prec_factory(mpPrecParams);
    Teuchos::RCP<Epetra_Operator> precOp = 
      mp_prec_factory.build(mp_comm, num_mp_blocks, x_map, mp_x_map);
    return Teuchos::rcp(new EpetraExt::ModelEvaluator::Preconditioner(precOp,
								      true));
  }
  
  return Teuchos::null;
}

Teuchos::RCP<Epetra_Operator>
Stokhos::MPModelEvaluator::create_DgDx_op(int j) const
{
  if (j < num_g)
    return me->create_DgDx_op(j);
  else if (j < num_g + num_g_mp && supports_x) {
    int jj = j-num_g;
    Teuchos::RCP< Stokhos::ProductContainer<Epetra_Operator> > mp_blocks =
      Teuchos::rcp(new Stokhos::ProductContainer<Epetra_Operator>(
		     mp_block_map));
    OutArgs me_outargs = me->createOutArgs();
    if (me_outargs.supports(OUT_ARG_DgDx_mp, jj).supports(DERIV_LINEAR_OP))
      for (unsigned int i=0; i<num_mp_blocks; i++)
	mp_blocks->setCoeffPtr(i, me->create_DgDx_mp_op(jj));
    else if (me_outargs.supports(OUT_ARG_DgDx_mp, jj).supports(DERIV_MV_BY_COL))
      for (unsigned int i=0; i<num_mp_blocks; i++) {
	Teuchos::RCP<Epetra_MultiVector> mv = 
	  Teuchos::rcp(new Epetra_MultiVector(*(me->get_g_mp_map(jj)), 
					      me->get_x_map()->NumMyElements()));
	Teuchos::RCP<Epetra_Operator> block = 
	  Teuchos::rcp(new Stokhos::EpetraMultiVectorOperator(mv));
	mp_blocks->setCoeffPtr(i, block);
      }
    else if (me_outargs.supports(OUT_ARG_DgDx_mp, jj).supports(DERIV_TRANS_MV_BY_ROW))
      for (unsigned int i=0; i<num_mp_blocks; i++) {
	Teuchos::RCP<Epetra_MultiVector> mv = 
	  Teuchos::rcp(new Epetra_MultiVector(*(me->get_x_map()), 
					      me->get_g_mp_map(jj)->NumMyElements()));
	Teuchos::RCP<Epetra_Operator> block = 
	  Teuchos::rcp(new Stokhos::EpetraMultiVectorOperator(mv));
	mp_blocks->setCoeffPtr(i, block);
      }
    else
      TEST_FOR_EXCEPTION(true, std::logic_error,
			 "Error!  me_outargs.supports(OUT_ARG_DgDx_mp, " << j
			 << ").none() is true!");

    Teuchos::RCP<Stokhos::BlockDiagonalOperator> dgdx_mp = 
      Teuchos::rcp(new Stokhos::BlockDiagonalOperator(
		     mp_comm, num_mp_blocks, x_map, me->get_g_mp_map(jj),
		     mp_x_map, mp_g_map[jj]));
    dgdx_mp->setupOperator(mp_blocks);
    return dgdx_mp;
  }
  else 
    TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  DgDx_op is not supported for index " << j
		       << "!");
  
  return Teuchos::null;
}

Teuchos::RCP<Epetra_Operator>
Stokhos::MPModelEvaluator::create_DgDx_dot_op(int j) const
{
  if (j < num_g)
    return me->create_DgDx_dot_op(j);
  else if (j < num_g + num_g_mp && supports_x) {
    int jj = j-num_g;
    Teuchos::RCP< Stokhos::ProductContainer<Epetra_Operator> > mp_blocks =
      Teuchos::rcp(new Stokhos::ProductContainer<Epetra_Operator>(
		     mp_block_map));
    OutArgs me_outargs = me->createOutArgs();
    if (me_outargs.supports(OUT_ARG_DgDx_dot_mp, jj).supports(DERIV_LINEAR_OP))
      for (unsigned int i=0; i<num_mp_blocks; i++)
	mp_blocks->setCoeffPtr(i, me->create_DgDx_dot_mp_op(jj));
    else if (me_outargs.supports(OUT_ARG_DgDx_dot_mp, jj).supports(DERIV_MV_BY_COL))
      for (unsigned int i=0; i<num_mp_blocks; i++) {
	Teuchos::RCP<Epetra_MultiVector> mv = 
	  Teuchos::rcp(new Epetra_MultiVector(*(me->get_g_mp_map(jj)), 
					      me->get_x_map()->NumMyElements()));
	Teuchos::RCP<Epetra_Operator> block = 
	  Teuchos::rcp(new Stokhos::EpetraMultiVectorOperator(mv));
	mp_blocks->setCoeffPtr(i, block);
      }
    else if (me_outargs.supports(OUT_ARG_DgDx_dot_mp, jj).supports(DERIV_TRANS_MV_BY_ROW))
      for (unsigned int i=0; i<num_mp_blocks; i++) {
	Teuchos::RCP<Epetra_MultiVector> mv = 
	  Teuchos::rcp(new Epetra_MultiVector(*(me->get_x_map()), 
					      me->get_g_mp_map(jj)->NumMyElements()));
	
	Teuchos::RCP<Epetra_Operator> block = 
	  Teuchos::rcp(new Stokhos::EpetraMultiVectorOperator(mv));
	mp_blocks->setCoeffPtr(i, block);
      }
    else
      TEST_FOR_EXCEPTION(true, std::logic_error,
			 "Error!  me_outargs.supports(OUT_ARG_DgDx_dot_mp, " 
			 << j << ").none() is true!");

    Teuchos::RCP<Teuchos::ParameterList> pl = 
      Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::RCP<Stokhos::BlockDiagonalOperator> dgdx_dot_mp = 
      Teuchos::rcp(new Stokhos::BlockDiagonalOperator(
		     mp_comm, num_mp_blocks, x_map, me->get_g_mp_map(jj),
		     mp_x_map, mp_g_map[jj]));
    dgdx_dot_mp->setupOperator(mp_blocks);
  }
  else 
    TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  DgDx_dot_op is not supported for index " << j
		       << "!");

  return Teuchos::null;
}

EpetraExt::ModelEvaluator::InArgs
Stokhos::MPModelEvaluator::createInArgs() const
{
  InArgsSetup inArgs;
  InArgs me_inargs = me->createInArgs();

  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p + num_p_mp); 
  inArgs.set_Np_mp(0);
  inArgs.setSupports(IN_ARG_x_dot, me_inargs.supports(IN_ARG_x_dot));
  inArgs.setSupports(IN_ARG_x, me_inargs.supports(IN_ARG_x));
  inArgs.setSupports(IN_ARG_t, me_inargs.supports(IN_ARG_t));
  inArgs.setSupports(IN_ARG_alpha, me_inargs.supports(IN_ARG_alpha));
  inArgs.setSupports(IN_ARG_beta, me_inargs.supports(IN_ARG_beta));
  
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
Stokhos::MPModelEvaluator::createOutArgs() const
{
  OutArgsSetup outArgs;
  OutArgs me_outargs = me->createOutArgs();

  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(me_outargs.Np() + num_p_mp, num_g + num_g_mp);
  outArgs.set_Np_Ng_mp(0, 0);
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
  for (int i=0; i<num_g_mp; i++) {
    if (!me_outargs.supports(OUT_ARG_DgDx_dot_mp, i).none())
      outArgs.setSupports(OUT_ARG_DgDx_dot, i+num_g,  DERIV_LINEAR_OP);
    if (!me_outargs.supports(OUT_ARG_DgDx_mp, i).none())
      outArgs.setSupports(OUT_ARG_DgDx, i+num_g,  DERIV_LINEAR_OP);
    for (int j=0; j<me_outargs.Np(); j++)
      outArgs.setSupports(OUT_ARG_DgDp, i+num_g, j, 
			  me_outargs.supports(OUT_ARG_DgDp_mp, i, j));
  }

  // We do not support derivatives w.r.t. the new MP parameters, so their
  // support defaults to none.
  
  return outArgs;
}

void 
Stokhos::MPModelEvaluator::evalModel(const InArgs& inArgs,
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
  // for the preconditioner.  Given the expense of computing the MP W blocks
  // this saves significant computational cost
  if (eval_prec) {
    Teuchos::RCP<Stokhos::MPPreconditioner> W_prec = 
      Teuchos::rcp_dynamic_cast<Stokhos::MPPreconditioner>(WPrec_out, true);
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
   Teuchos::RCP<Stokhos::ProductEpetraVector> x_mp = 
      create_x_mp(View, x.get());
    me_inargs.set_x_mp(x_mp);
  }
  if (x_dot != Teuchos::null) {
    Teuchos::RCP<Stokhos::ProductEpetraVector> x_dot_mp = 
      create_x_mp(View, x_dot.get());
    me_inargs.set_x_dot_mp(x_dot_mp);
  }
  if (me_inargs.supports(IN_ARG_alpha))
    me_inargs.set_alpha(inArgs.get_alpha());
  if (me_inargs.supports(IN_ARG_beta))
    me_inargs.set_beta(inArgs.get_beta());
  if (me_inargs.supports(IN_ARG_t))
    me_inargs.set_t(inArgs.get_t());

  // Pass parameters
  for (int i=0; i<num_p; i++)
    me_inargs.set_p(i, inArgs.get_p(i));
  for (int i=0; i<num_p_mp; i++) {
    Teuchos::RCP<const Epetra_Vector> p = inArgs.get_p(i+num_p);

    // We always need to pass in the MP parameters, so just use
    // the initial parameters if it is NULL
    if (p == Teuchos::null)
      p = mp_p_init[i]->getBlockVector();

    // Convert block p to MP polynomial
    Teuchos::RCP<Stokhos::ProductEpetraVector> p_mp = 
      create_p_mp(i, View, p.get());
    me_inargs.set_p_mp(i, p_mp);
  }

  // Create underlying outargs
  OutArgs me_outargs = me->createOutArgs();

  // f
  if (f_out != Teuchos::null) {
    Teuchos::RCP<Stokhos::ProductEpetraVector> f_mp = 
      create_f_mp(View, f_out.get());
    me_outargs.set_f_mp(f_mp);
  }

  // W
  if (W_out != Teuchos::null && !eval_prec)
     me_outargs.set_W_mp(W_mp_blocks);

  // df/dp
  for (int i=0; i<outArgs.Np(); i++) {
    if (!outArgs.supports(OUT_ARG_DfDp, i).none()) {
      Derivative dfdp = outArgs.get_DfDp(i);
      if (dfdp.getMultiVector() != Teuchos::null) {
	Teuchos::RCP<Stokhos::ProductEpetraMultiVector> dfdp_mp; 
	if (dfdp.getMultiVectorOrientation() == DERIV_MV_BY_COL)
	  dfdp_mp = 
	    Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
			   mp_block_map, me->get_f_map(), mp_f_map, mp_comm, 
			   View, *(dfdp.getMultiVector())));
	else if (dfdp.getMultiVectorOrientation() == DERIV_TRANS_MV_BY_ROW)
	  dfdp_mp = 
	    Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
			   mp_block_map, me->get_p_map(i), mp_p_map[i], mp_comm,
			   View, *(dfdp.getMultiVector())));
	me_outargs.set_DfDp_mp(i, 
			       MPDerivative(dfdp_mp,
					    dfdp.getMultiVectorOrientation()));
      }
      TEST_FOR_EXCEPTION(dfdp.getLinearOp() != Teuchos::null, std::logic_error,
			 "Error!  Stokhos::MPModelEvaluator::evalModel " << 
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
  for (int i=0; i<num_g_mp; i++) {
    // g
    Teuchos::RCP<Epetra_Vector> g = outArgs.get_g(i+num_g);
    if (g != Teuchos::null) {
      Teuchos::RCP<Stokhos::ProductEpetraVector> g_mp =
	create_g_mp(i, View, g.get());
      me_outargs.set_g_mp(i, g_mp);
    }

    // dg/dxdot
    if (outArgs.supports(OUT_ARG_DgDx_dot, i+num_g).supports(DERIV_LINEAR_OP)) {
      Derivative dgdx_dot = outArgs.get_DgDx_dot(i+num_g);
      if (dgdx_dot.getLinearOp() != Teuchos::null) {
	Teuchos::RCP<Stokhos::BlockDiagonalOperator> op =
	  Teuchos::rcp_dynamic_cast<Stokhos::BlockDiagonalOperator>(
	    dgdx_dot.getLinearOp(), true);
	Teuchos::RCP< Stokhos::ProductContainer<Epetra_Operator> > mp_blocks =
	  op->getMPOps();
	if (me_outargs.supports(OUT_ARG_DgDx, i).supports(DERIV_LINEAR_OP))
	  me_outargs.set_DgDx_dot_mp(i, mp_blocks);
	else {
	  for (unsigned int k=0; k<num_mp_blocks; k++) {
	    Teuchos::RCP<Epetra_MultiVector> mv = 
	      Teuchos::rcp_dynamic_cast<Stokhos::EpetraMultiVectorOperator>(
		mp_blocks->getCoeffPtr(k), true)->getMultiVector();
	    dgdx_dot_mp_blocks[i]->setCoeffPtr(k, mv);
	  }
	  if (me_outargs.supports(OUT_ARG_DgDx_dot_mp, i).supports(DERIV_MV_BY_COL))
	    me_outargs.set_DgDx_dot_mp(i, MPDerivative(dgdx_dot_mp_blocks[i],
						       DERIV_MV_BY_COL));
	  else
	    me_outargs.set_DgDx_dot_mp(i, MPDerivative(dgdx_dot_mp_blocks[i],
						       DERIV_TRANS_MV_BY_ROW));
	}
      }
      TEST_FOR_EXCEPTION(dgdx_dot.getLinearOp() == Teuchos::null &&
			 dgdx_dot.isEmpty() == false, 
			 std::logic_error,
			 "Error!  Stokhos::MPModelEvaluator::evalModel: " << 
			 "Operator form of dg/dxdot is required!");
    }

    // dg/dx
    if (outArgs.supports(OUT_ARG_DgDx, i+num_g).supports(DERIV_LINEAR_OP)) {
      Derivative dgdx = outArgs.get_DgDx(i+num_g);
      if (dgdx.getLinearOp() != Teuchos::null) {
	Teuchos::RCP<Stokhos::BlockDiagonalOperator> op =
	  Teuchos::rcp_dynamic_cast<Stokhos::BlockDiagonalOperator>(
	    dgdx.getLinearOp(), true);
	Teuchos::RCP< Stokhos::ProductContainer<Epetra_Operator> > mp_blocks =
	  op->getMPOps();
	if (me_outargs.supports(OUT_ARG_DgDx, i).supports(DERIV_LINEAR_OP))
	  me_outargs.set_DgDx_mp(i, mp_blocks);
	else {
	  for (unsigned int k=0; k<num_mp_blocks; k++) {
	    Teuchos::RCP<Epetra_MultiVector> mv = 
	      Teuchos::rcp_dynamic_cast<Stokhos::EpetraMultiVectorOperator>(
		mp_blocks->getCoeffPtr(k), true)->getMultiVector();
	    dgdx_mp_blocks[i]->setCoeffPtr(k, mv);
	  }
	  if (me_outargs.supports(OUT_ARG_DgDx_mp, i).supports(DERIV_MV_BY_COL))
	    me_outargs.set_DgDx_mp(i, MPDerivative(dgdx_mp_blocks[i],
						   DERIV_MV_BY_COL));
	  else
	    me_outargs.set_DgDx_mp(i, MPDerivative(dgdx_mp_blocks[i],
						   DERIV_TRANS_MV_BY_ROW));
	}
      }
      TEST_FOR_EXCEPTION(dgdx.getLinearOp() == Teuchos::null &&
			 dgdx.isEmpty() == false, 
			 std::logic_error,
			 "Error!  Stokhos::MPModelEvaluator::evalModel: " << 
			 "Operator form of dg/dxdot is required!");
    }

    // dg/dp
    // Rembember, no derivatives w.r.t. mp parameters
    for (int j=0; j<outArgs.Np(); j++) {
      if (!outArgs.supports(OUT_ARG_DgDp, i+num_g, j).none()) {
	Derivative dgdp = outArgs.get_DgDp(i+num_g,j);
	if (dgdp.getMultiVector() != Teuchos::null) {
	  Teuchos::RCP<Stokhos::ProductEpetraMultiVector> dgdp_mp; 
	  if (dgdp.getMultiVectorOrientation() == DERIV_MV_BY_COL)
	    dgdp_mp = 
	      Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
			     mp_block_map, me->get_g_mp_map(i), mp_g_map[i], 
			     mp_comm, View, *(dgdp.getMultiVector())));
	  else if (dgdp.getMultiVectorOrientation() == DERIV_TRANS_MV_BY_ROW)
	    dgdp_mp = 
	      Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
			     mp_block_map, me->get_p_map(j), mp_p_map[j], 
			     mp_comm, View, *(dgdp.getMultiVector())));
	  me_outargs.set_DgDp_mp(i, j, 
				 MPDerivative(dgdp_mp,
					      dgdp.getMultiVectorOrientation()));
	}
	TEST_FOR_EXCEPTION(dgdp.getLinearOp() != Teuchos::null, 
			   std::logic_error,
			   "Error!  Stokhos::MPModelEvaluator::evalModel " << 
			   "cannot handle operator form of dg/dp!");
      }
    }

  }

  // Compute the functions
  me->evalModel(me_inargs, me_outargs);

  // Copy block MP components for W
  if (W_out != Teuchos::null && !eval_prec) {
    Teuchos::RCP<Epetra_Operator> W;
    if (W_out != Teuchos::null)
      W = W_out;
    else
      W = my_W;
    Teuchos::RCP<Stokhos::BlockDiagonalOperator> W_mp = 
      Teuchos::rcp_dynamic_cast<Stokhos::BlockDiagonalOperator>(W, true);
    W_mp->setupOperator(W_mp_blocks);
      
    if (WPrec_out != Teuchos::null) {
      Teuchos::RCP<Stokhos::MPPreconditioner> W_prec = 
	Teuchos::rcp_dynamic_cast<Stokhos::MPPreconditioner>(WPrec_out, true);
      W_prec->setupPreconditioner(W_mp, *my_x);
    }
  }
}

void 
Stokhos::MPModelEvaluator::set_x_mp_init(
  const Stokhos::ProductEpetraVector& x_mp_in)
{
  *mp_x_init = x_mp_in;
}

void 
Stokhos::MPModelEvaluator::set_p_mp_init(
  int i, const Stokhos::ProductEpetraVector& p_mp_in)
{
  *mp_p_init[i] = p_mp_in;
}

Teuchos::Array<int> 
Stokhos::MPModelEvaluator::get_p_mp_indices() const
{
  Teuchos::Array<int> mp_p_index(num_p_mp);
  for (int i=0; i<num_p_mp; i++)
    mp_p_index[i] = num_p + i;
  return mp_p_index;
}

Teuchos::Array<int> 
Stokhos::MPModelEvaluator::get_non_p_mp_indices() const
{
  Teuchos::Array<int> non_mp_p_index(num_p);
  for (int i=0; i<num_p; i++)
    non_mp_p_index[i] = i;
  return non_mp_p_index;
}

Teuchos::Array<int> 
Stokhos::MPModelEvaluator::get_g_mp_indices() const
{
  Teuchos::Array<int> mp_g_index(num_g_mp);
  for (int i=0; i<num_g_mp; i++)
    mp_g_index[i] = num_g + i;
  return mp_g_index;
}

Teuchos::Array<int> 
Stokhos::MPModelEvaluator::get_non_g_mp_indices() const
{
  Teuchos::Array<int> non_mp_g_index(num_g);
  for (int i=0; i<num_g; i++)
    non_mp_g_index[i] = i;
  return non_mp_g_index;
}

Teuchos::Array< Teuchos::RCP<const Epetra_Map> > 
Stokhos::MPModelEvaluator::get_p_mp_base_maps() const
{
  Teuchos::Array< Teuchos::RCP<const Epetra_Map> > base_maps(num_p_mp);
  for (int i=0; i<num_p_mp; i++)
    base_maps[i] = me->get_p_mp_map(i);
  return base_maps;
}

Teuchos::Array< Teuchos::RCP<const Epetra_Map> > 
Stokhos::MPModelEvaluator::get_g_mp_base_maps() const
{
  Teuchos::Array< Teuchos::RCP<const Epetra_Map> > base_maps(num_g_mp);
  for (int i=0; i<num_g_mp; i++)
    base_maps[i] = me->get_g_mp_map(i);
  return base_maps;
}

Teuchos::RCP<Stokhos::ProductEpetraVector>
Stokhos::MPModelEvaluator::create_x_mp(Epetra_DataAccess CV, 
				       const Epetra_Vector* v) const
{
  Teuchos::RCP<Stokhos::ProductEpetraVector> mp_x;
  if (v == NULL)
    mp_x = Teuchos::rcp(new Stokhos::ProductEpetraVector(
			  mp_block_map, x_map, mp_x_map, mp_comm));
  else
    mp_x = Teuchos::rcp(new Stokhos::ProductEpetraVector(
			  mp_block_map, x_map, mp_x_map, mp_comm,
			  CV, *v));
  return mp_x;
}

Teuchos::RCP<Stokhos::ProductEpetraMultiVector>
Stokhos::MPModelEvaluator::create_x_mv_mp(int num_vecs, Epetra_DataAccess CV, 
					  const Epetra_MultiVector* v) const
{
  Teuchos::RCP<Stokhos::ProductEpetraMultiVector> mp_x;
  if (v == NULL)
    mp_x = Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
			  mp_block_map, x_map, mp_x_map, mp_comm,
			  num_vecs));
  else
    mp_x = Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
			  mp_block_map, x_map, mp_x_map, mp_comm,
			  CV, *v));
  return mp_x;
}

Teuchos::RCP<Stokhos::ProductEpetraVector>
Stokhos::MPModelEvaluator::create_p_mp(int l, Epetra_DataAccess CV, 
				       const Epetra_Vector* v) const
{
  Teuchos::RCP<Stokhos::ProductEpetraVector> mp_p;
  if (v == NULL)
    mp_p = Teuchos::rcp(new Stokhos::ProductEpetraVector(
			  mp_block_map, me->get_p_mp_map(l),
			  mp_p_map[l], mp_comm));
  else
    mp_p = Teuchos::rcp(new Stokhos::ProductEpetraVector(
			  mp_block_map, me->get_p_mp_map(l),
			  mp_p_map[l], mp_comm, CV, *v));
  return mp_p;
}

Teuchos::RCP<Stokhos::ProductEpetraVector>
Stokhos::MPModelEvaluator::create_f_mp(Epetra_DataAccess CV, 
				       const Epetra_Vector* v) const
{
  Teuchos::RCP<Stokhos::ProductEpetraVector> mp_f;
  if (v == NULL)
    mp_f = Teuchos::rcp(new Stokhos::ProductEpetraVector(
			  mp_block_map, f_map, mp_f_map, mp_comm));
  else
    mp_f = Teuchos::rcp(new Stokhos::ProductEpetraVector(
			  mp_block_map, f_map, mp_f_map, mp_comm,
			  CV, *v));
  return mp_f;
}

Teuchos::RCP<Stokhos::ProductEpetraMultiVector>
Stokhos::MPModelEvaluator::create_f_mv_mp(
  int num_vecs, 
  Epetra_DataAccess CV, 
  const Epetra_MultiVector* v) const
{
  Teuchos::RCP<Stokhos::ProductEpetraMultiVector> mp_f;
  if (v == NULL)
    mp_f = Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
			  mp_block_map, f_map, mp_f_map, mp_comm,
			  num_vecs));
  else
    mp_f = Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
			  mp_block_map, f_map, mp_f_map, mp_comm,
			  CV, *v));
  return mp_f;
}

Teuchos::RCP<Stokhos::ProductEpetraVector>
Stokhos::MPModelEvaluator::create_g_mp(int l, Epetra_DataAccess CV, 
				       const Epetra_Vector* v) const
{
  Teuchos::RCP<Stokhos::ProductEpetraVector> mp_g;
  if (v == NULL)
    mp_g = Teuchos::rcp(new Stokhos::ProductEpetraVector(
			  mp_block_map, 
			  me->get_g_mp_map(l), 
			  mp_g_map[l], mp_comm));
  else
    mp_g = Teuchos::rcp(new Stokhos::ProductEpetraVector(
			  mp_block_map, 
			  me->get_g_mp_map(l), 
			  mp_g_map[l], mp_comm, CV, *v));
  return mp_g;
}

Teuchos::RCP<Stokhos::ProductEpetraMultiVector>
Stokhos::MPModelEvaluator::create_g_mv_mp(int l, int num_vecs,
					  Epetra_DataAccess CV, 
					  const Epetra_MultiVector* v) const
{
  Teuchos::RCP<Stokhos::ProductEpetraMultiVector> mp_g;
  if (v == NULL)
    mp_g = Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
			  mp_block_map, 
			  me->get_g_mp_map(l), 
			  mp_g_map[l], mp_comm, num_vecs));
  else
    mp_g = Teuchos::rcp(new Stokhos::ProductEpetraMultiVector(
			  mp_block_map, 
			  me->get_g_mp_map(l), 
			  mp_g_map[l], mp_comm, CV, *v));
  return mp_g;
}
