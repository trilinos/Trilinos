// @HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER


#include "EpetraExt_DiagonalTransientModel.hpp"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LocalMap.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"


namespace {


const std::string Implicit_name = "Implicit";
const bool Implicit_default = true;

const std::string Lambda_min_name = "Lambda_min";
const double Lambda_min_default = -0.9;

const std::string Lambda_max_name = "Lambda_max";
const double Lambda_max_default = -0.01;

const std::string Coeff_s_name = "Coeff_s";
const string Coeff_s_default = "{ 0.0 }";

const std::string Lambda_fit_name = "Lambda_fit";
const std::string Lambda_fit_default = "Linear"; // Will be validated at runtime!

const std::string NumElements_name = "NumElements";
const int NumElements_default = 1;

const std::string x0_name = "x0";
const double x0_default = 10.0;


inline
double evalR( const double& t, const double& lambda, const double& s )
{
  return (exp(lambda*t)*sin(s*t));
}


inline
double d_evalR_d_s( const double& t, const double& lambda, const double& s )
{
  return (exp(lambda*t)*cos(s*t));
}


inline
double f_func( const double& x, const double& t, const double& lambda, const double& s )
{
  return ( lambda*x + evalR(t,lambda,s) );
}


} // namespace


namespace EpetraExt {


// Constructors


DiagonalTransientModel::DiagonalTransientModel(
  Teuchos::RefCountPtr<Epetra_Comm> const& epetra_comm
  )
  : epetra_comm_(epetra_comm.assert_not_null()),
    implicit_(Implicit_default),
    numElements_(NumElements_default),
    lambda_min_(Lambda_min_default),
    lambda_max_(Lambda_max_default),
    coeff_s_(Teuchos::fromStringToArray<double>(Coeff_s_default)),
    lambda_fit_(LAMBDA_FIT_LINEAR), // Must be the same as Lambda_fit_default!
    x0_(x0_default),
    isIntialized_(false)
{
  initialize();
}


Teuchos::RefCountPtr<const Epetra_Vector>
DiagonalTransientModel::getExactSolution(const double t) const
{
  Teuchos::RefCountPtr<Epetra_Vector>
    x_star_ptr = Teuchos::rcp(new Epetra_Vector(*epetra_map_,false));
  Epetra_Vector& x_star = *x_star_ptr;
  Epetra_Vector& lambda = *lambda_;
  int myN = x_star.MyLength();
  for ( int i=0 ; i<myN ; ++i ) {
    const double coeff_s_i = coeff_s(i);
    const double t_lambda_i = t*lambda[i];
    if (coeff_s_i == 0.0)
      x_star[i] = x0_*exp(t_lambda_i);
    else
      x_star[i] = 
        (x0_+(1.0/coeff_s_i))*exp(t_lambda_i)
        - exp(t_lambda_i)*cos(t*coeff_s_i)/coeff_s_i;
  }
  return(x_star_ptr);
}


// Overridden from ParameterListAcceptor


void DiagonalTransientModel::setParameterList(
  Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList
  )
{
  using Teuchos::get; using Teuchos::getIntegralValue;
  using Teuchos::getArrayFromStringParameter;
  TEST_FOR_EXCEPT( is_null(paramList) );
  paramList->validateParametersAndSetDefaults(*this->getValidParameters());
  isIntialized_ = false;
  paramList_ = paramList;
  implicit_ = get<bool>(*paramList_,Implicit_name);
  numElements_ = get<int>(*paramList_,NumElements_name);
  lambda_min_ = get<double>(*paramList_,Lambda_min_name);
  lambda_max_ = get<double>(*paramList_,Lambda_max_name);
  coeff_s_ = getArrayFromStringParameter<double>(*paramList_,Coeff_s_name);
  lambda_fit_ = getIntegralValue<ELambdaFit>(*paramList_,Lambda_fit_name);
  x0_ = get<double>(*paramList_,x0_name);
  initialize();
}


Teuchos::RefCountPtr<Teuchos::ParameterList> 
DiagonalTransientModel::getParameterList()
{
  return paramList_;
}


Teuchos::RefCountPtr<Teuchos::ParameterList>
DiagonalTransientModel::unsetParameterList()
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}


Teuchos::RefCountPtr<const Teuchos::ParameterList>
DiagonalTransientModel::getParameterList() const
{
  return paramList_;
}


Teuchos::RefCountPtr<const Teuchos::ParameterList>
DiagonalTransientModel::getValidParameters() const
{
  using Teuchos::RefCountPtr; using Teuchos::ParameterList;
  using Teuchos::tuple;
  using Teuchos::setIntParameter; using Teuchos::setDoubleParameter;
  using Teuchos::setStringToIntegralParameter;
  static RefCountPtr<const ParameterList> validPL;
  if (is_null(validPL)) {
    RefCountPtr<ParameterList> pl = Teuchos::parameterList();
    pl->set(Implicit_name,true);
    setDoubleParameter(
      Lambda_min_name, Lambda_min_default, "",
      &*pl
      );
    setDoubleParameter(
      Lambda_max_name, Lambda_max_default, "",
      &*pl
      );
    setStringToIntegralParameter(
      Lambda_fit_name, Lambda_fit_default, "",
      tuple<std::string>("Linear","Random"),
      tuple<ELambdaFit>(LAMBDA_FIT_LINEAR,LAMBDA_FIT_RANDOM),
      &*pl
      );
    setIntParameter(
      NumElements_name, NumElements_default, "",
      &*pl
      );
    setDoubleParameter(
      x0_name, x0_default, "",
      &*pl
      );
    pl->set( Coeff_s_name, Coeff_s_default );
    validPL = pl;
  }
  return validPL;
}


// Overridden from EpetraExt::ModelEvaluator


Teuchos::RefCountPtr<const Epetra_Map>
DiagonalTransientModel::get_x_map() const
{
  return epetra_map_;
}


Teuchos::RefCountPtr<const Epetra_Map>
DiagonalTransientModel::get_f_map() const
{
  return epetra_map_;
}


Teuchos::RefCountPtr<const Epetra_Map>
DiagonalTransientModel::get_p_map(int l) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT( l == 0 );
#endif
  return map_p_[l];
}


Teuchos::RefCountPtr<const Epetra_Vector>
DiagonalTransientModel::get_x_init() const
{
  return x_init_;
}


Teuchos::RefCountPtr<const Epetra_Vector>
DiagonalTransientModel::get_x_dot_init() const
{
  return x_dot_init_;
}


Teuchos::RefCountPtr<const Epetra_Vector>
DiagonalTransientModel::get_p_init(int l) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT( l == 0 );
#endif
  return p_init_[l];
}


Teuchos::RefCountPtr<Epetra_Operator>
DiagonalTransientModel::create_W() const
{
  if(implicit_)
    return Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*W_graph_));
  return Teuchos::null;
}


EpetraExt::ModelEvaluator::InArgs
DiagonalTransientModel::createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.set_Np(Np_);
  inArgs.setSupports(IN_ARG_t,true);
  inArgs.setSupports(IN_ARG_x,true);
  if(implicit_) {
    inArgs.setSupports(IN_ARG_x_dot,true);
    inArgs.setSupports(IN_ARG_alpha,true);
    inArgs.setSupports(IN_ARG_beta,true);
  }
  return inArgs;
}


EpetraExt::ModelEvaluator::OutArgs
DiagonalTransientModel::createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.set_Np_Ng(Np_,0);
  outArgs.setSupports(OUT_ARG_f,true);
  if(implicit_) {
    outArgs.setSupports(OUT_ARG_W,true);
    outArgs.set_W_properties(
      DerivativeProperties(
        DERIV_LINEARITY_NONCONST
        ,DERIV_RANK_FULL
        ,true // supportsAdjoint
        )
      );
  }
  outArgs.setSupports(OUT_ARG_DfDp,0,DERIV_MV_BY_COL);
  outArgs.set_DfDp_properties(
    0,DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_DEFICIENT
      ,true // supportsAdjoint
      )
    );
  return outArgs;
}


void DiagonalTransientModel::evalModel(
  const InArgs& inArgs, const OutArgs& outArgs
  ) const
{

  using Teuchos::RefCountPtr;
  using Teuchos::null;
  using Teuchos::dyn_cast;

  const Epetra_Vector &x = *(inArgs.get_x());
  const double t = inArgs.get_t();
  const Epetra_Vector *p = 0;
  if (Np_) p = inArgs.get_p(0).get();

  Epetra_Operator *W_out = ( implicit_ ? outArgs.get_W().get() : 0 );
  Epetra_Vector *f_out = outArgs.get_f().get();
  Epetra_MultiVector *DfDp_out = 0;
  if (Np_) DfDp_out = get_DfDp_mv(0,outArgs).get();

  const Epetra_Vector &lambda = *lambda_;

  int localNumElements = x.MyLength();

  if (f_out) {
    Epetra_Vector &f = *f_out;
    if (implicit_) {
      const Epetra_Vector *x_dot_in = inArgs.get_x_dot().get();
      if (x_dot_in) {
        const Epetra_Vector &x_dot = *x_dot_in;
        for ( int i=0 ; i<localNumElements ; ++i )
          f[i] = x_dot[i] - f_func(x[i],t,lambda[i],coeff_s(i));
      }
      else {
        for ( int i=0 ; i<localNumElements ; ++i )
          f[i] = - f_func(x[i],t,lambda[i],coeff_s(i));
      }
    }
    else {
      for ( int i=0 ; i<localNumElements ; ++i )
        f[i] = f_func(x[i],t,lambda[i],coeff_s(i));
    }
  }

  if ( W_out ) {
    // If we get here then we are in implicit mode!
    const double alpha = inArgs.get_alpha();
    const double beta = inArgs.get_beta();
    Epetra_CrsMatrix &crsW = dyn_cast<Epetra_CrsMatrix>(*W_out);
    double values[1];
    int indices[1];
    const int offset
      = epetra_comm_->MyPID()*localNumElements + epetra_map_->IndexBase(); 
    for( int i = 0; i < localNumElements; ++i ) {
      values[0] = alpha - beta*lambda[i];
      indices[0] = i + offset;  // global column
      crsW.ReplaceGlobalValues(
        i + offset // GlobalRow
        ,1 // NumEntries
        ,values // Values
        ,indices // Indices
        );
    }
  }

  if (DfDp_out) {
    Epetra_MultiVector &DfDp = *DfDp_out;
    DfDp.PutScalar(0.0);
    if (implicit_) {
      for( int i = 0; i < localNumElements; ++i ) {
        DfDp[coeff_s_idx(i)][i]
          = - d_evalR_d_s(t,lambda[i],coeff_s(i));
      }
    }
    else {
      for( int i = 0; i < localNumElements; ++i ) {
        DfDp[coeff_s_idx(i)][i]
          = + d_evalR_d_s(t,lambda[i],coeff_s(i));
      }
    }
  }
  
}


// private


void DiagonalTransientModel::initialize()
{

  using Teuchos::rcp;
  using Teuchos::as;

  //
  // Setup map
  //

  const int numProcs = epetra_comm_->NumProc();
  const int procRank = epetra_comm_->MyPID();
  epetra_map_ = rcp( new Epetra_Map( numElements_ * numProcs, 0, *epetra_comm_ ) );

  //
  // Create lambda
  //

  lambda_ = Teuchos::rcp(new Epetra_Vector(*epetra_map_));
  Epetra_Vector &lambda = *lambda_;
  switch(lambda_fit_) {
    case LAMBDA_FIT_LINEAR: {
      const int N = lambda.GlobalLength();
      const double slope = (lambda_max_ - lambda_min_)/(N-1);
      const int MyLength = lambda.MyLength();
      for ( int i=0 ; i<MyLength ; ++i )
      {
        lambda[i] = slope*(procRank*MyLength+i)+lambda_min_;
      }
      break;
    }
    case LAMBDA_FIT_RANDOM: {
      unsigned int seed = Teuchos::ScalarTraits<int>::random()+10*procRank; 
      seed *= seed;
      lambda.SetSeed(seed);
      lambda.Random(); // fill with random numbers in (-1,1)
      // Scale random numbers to (lambda_min_,lambda_max_)
      const double slope = (lambda_min_ - lambda_max_)/2.0;
      const double tmp = (lambda_max_ + lambda_min_)/2.0;
      int MyLength = lambda.MyLength();
      for (int i=0 ; i<MyLength ; ++i)
      {
        lambda[i] = slope*lambda[i] + tmp;
      }
      break;
    }
    default:
      TEST_FOR_EXCEPT("Should never get here!");
  }

  //
  // Setup for parameters
  //

  Np_ = 1;
  np_ = coeff_s_.size();
  map_p_.clear();
  map_p_.push_back(
    rcp( new Epetra_LocalMap(np_,0,*epetra_comm_) )
    );

  coeff_s_idx_.resize(0);
  const int num_func_per_coeff_s = numElements_ / np_;
  const int num_func_per_coeff_s_rem = numElements_ % np_;
  for ( int coeff_s_idx_i = 0; coeff_s_idx_i < np_; ++coeff_s_idx_i ) {
    const int num_func_per_coeff_s_idx_i
      = num_func_per_coeff_s
      + ( coeff_s_idx_i < num_func_per_coeff_s_rem ? 1 : 0 );
    for (
      int coeff_s_idx_i_j = 0;
      coeff_s_idx_i_j < num_func_per_coeff_s_idx_i; 
      ++coeff_s_idx_i_j
      )
    {
      coeff_s_idx_.push_back(coeff_s_idx_i);
    }
  }
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(
    ( as<int>(coeff_s_idx_.size()) != numElements_ ) && "Internal programming error!" );
#endif

  //
  // Setup graph for W
  //

  if(implicit_) {

    W_graph_ = rcp(
      new Epetra_CrsGraph(
        ::Copy,*epetra_map_,
        1 // elements per row
        )
      );

    int indices[1];
    const int offset = procRank*numElements_ + epetra_map_->IndexBase(); 

    for( int i = 0; i < numElements_; ++i ) {
      indices[0] = i + offset;  // global column
      W_graph_->InsertGlobalIndices(
        i + offset // GlobalRow
        ,1 // NumEntries
        ,indices // Indices
        );
    }

    W_graph_->FillComplete();

  }

  //
  // Setup initial guess
  //

  // Set x_init
  x_init_ = Teuchos::rcp(new Epetra_Vector(*epetra_map_,false));
  x_init_->PutScalar(x0_);

  // Set x_dot_init to provide for a consistent inital guess for implicit mode
  // such that f(x_dot,x) = 0!
  if (implicit_) {
    x_dot_init_ = Teuchos::rcp(new Epetra_Vector(*epetra_map_,false));
    InArgs inArgs = this->createInArgs();
    inArgs.set_x(x_init_);
    inArgs.set_t(0.0);
    OutArgs outArgs = this->createOutArgs();
    outArgs.set_f(x_dot_init_);
    this->evalModel(inArgs,outArgs);
    x_dot_init_->Scale(-1.0);
  }

  // Set p_init
  p_init_.push_back(
    rcp( new Epetra_Vector( ::Copy, *map_p_[0], &coeff_s_[0] ) )
    );

}


} // namespace EpetraExt


// Nonmembers


Teuchos::RefCountPtr<EpetraExt::DiagonalTransientModel>
EpetraExt::diagonalTransientModel(
  Teuchos::RefCountPtr<Epetra_Comm> const& epetra_comm,
  Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList
  )
{
  Teuchos::RefCountPtr<DiagonalTransientModel>
    diagonalTransientModel
    = Teuchos::rcp(new DiagonalTransientModel(epetra_comm));
  if (!is_null(paramList))
    diagonalTransientModel->setParameterList(paramList);
  return diagonalTransientModel;
}

