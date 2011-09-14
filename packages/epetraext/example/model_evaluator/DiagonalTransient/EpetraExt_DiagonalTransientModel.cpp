//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER


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


using Teuchos::RCP;


const std::string Implicit_name = "Implicit";
const bool Implicit_default = true;

const std::string Gamma_min_name = "Gamma_min";
const double Gamma_min_default = -0.9;

const std::string Gamma_max_name = "Gamma_max";
const double Gamma_max_default = -0.01;

const std::string Coeff_s_name = "Coeff_s";
const string Coeff_s_default = "{ 0.0 }";

const std::string Gamma_fit_name = "Gamma_fit";
const std::string Gamma_fit_default = "Linear"; // Will be validated at runtime!

const std::string NumElements_name = "NumElements";
const int NumElements_default = 1;

const std::string x0_name = "x0";
const double x0_default = 10.0;

const std::string ExactSolutionAsResponse_name = "Exact Solution as Response";
const bool ExactSolutionAsResponse_default = false;


inline
double evalR( const double& t, const double& gamma, const double& s )
{
  return (exp(gamma*t)*sin(s*t));
}


inline
double d_evalR_d_s( const double& t, const double& gamma, const double& s )
{
  return (exp(gamma*t)*cos(s*t)*t);
}


inline
double f_func(
  const double& x, const double& t, const double& gamma, const double& s
  )
{
  return ( gamma*x + evalR(t,gamma,s) );
}


inline
double x_exact(
  const double& x0, const double& t, const double& gamma, const double& s
  )
{
  if ( s == 0.0 )
    return ( x0 * exp(gamma*t) );
  return ( exp(gamma*t) * (x0 + (1.0/s) * ( 1.0 - cos(s*t) ) ) );
  // Note that the limit of (1.0/s) * ( 1.0 - cos(s*t) ) as s goes to zero is
  // zero.  This limit is neeed to avoid the 0/0 that would occur if floating
  // point was used to evaluate this term.  This means that cos(t*s) goes to
  // one at the same rate as s goes to zero giving 1-1=0..
}


inline
double dxds_exact(
  const double& t, const double& gamma, const double& s
  )
{
  if ( s == 0.0 )
    return 0.0;
  return ( -exp(gamma*t)/(s*s) * ( 1.0 - sin(s*t)*(s*t) - cos(s*t) ) );
}


class UnsetParameterVector {
public:
  ~UnsetParameterVector()
    {
      if (!is_null(vec_))
        *vec_ = Teuchos::null;
    }
  UnsetParameterVector(
    const RCP<RCP<const Epetra_Vector> > &vec
    )
    {
      setVector(vec);
    }
  void setVector( const RCP<RCP<const Epetra_Vector> > &vec )
    {
      vec_ = vec;
    }
private:
  RCP<RCP<const Epetra_Vector> > vec_;
};


} // namespace


namespace EpetraExt {


// Constructors


DiagonalTransientModel::DiagonalTransientModel(
  Teuchos::RCP<Epetra_Comm> const& epetra_comm
  )
  : epetra_comm_(epetra_comm.assert_not_null()),
    implicit_(Implicit_default),
    numElements_(NumElements_default),
    gamma_min_(Gamma_min_default),
    gamma_max_(Gamma_max_default),
    coeff_s_(Teuchos::fromStringToArray<double>(Coeff_s_default)),
    gamma_fit_(GAMMA_FIT_LINEAR), // Must be the same as Gamma_fit_default!
    x0_(x0_default),
    exactSolutionAsResponse_(ExactSolutionAsResponse_default),
    isIntialized_(false)
{
  initialize();
}


Teuchos::RCP<const Epetra_Vector>
DiagonalTransientModel::get_gamma() const
{
  return gamma_;
}


Teuchos::RCP<const Epetra_Vector>
DiagonalTransientModel::getExactSolution(
  const double t, const Epetra_Vector *coeff_s_p
  ) const
{
  set_coeff_s_p(Teuchos::rcp(coeff_s_p,false));
  UnsetParameterVector unsetParameterVector(Teuchos::rcp(&coeff_s_p_,false));
  Teuchos::RCP<Epetra_Vector>
    x_star_ptr = Teuchos::rcp(new Epetra_Vector(*epetra_map_,false));
  Epetra_Vector& x_star = *x_star_ptr;
  Epetra_Vector& gamma = *gamma_;
  int myN = x_star.MyLength();
  for ( int i=0 ; i<myN ; ++i ) {
    x_star[i] = x_exact( x0_, t, gamma[i], coeff_s(i) );
  }
  return(x_star_ptr);
}


Teuchos::RCP<const Epetra_MultiVector>
DiagonalTransientModel::getExactSensSolution(
  const double t, const Epetra_Vector *coeff_s_p
  ) const
{
  set_coeff_s_p(Teuchos::rcp(coeff_s_p,false));
  UnsetParameterVector unsetParameterVector(Teuchos::rcp(&coeff_s_p_,false));
  Teuchos::RCP<Epetra_MultiVector>
    dxds_star_ptr = Teuchos::rcp(new Epetra_MultiVector(*epetra_map_,np_,false));
  Epetra_MultiVector& dxds_star = *dxds_star_ptr;
  dxds_star.PutScalar(0.0);
  Epetra_Vector& gamma = *gamma_;
  int myN = dxds_star.MyLength();
  for ( int i=0 ; i<myN ; ++i ) {
    const int coeff_s_idx_i = this->coeff_s_idx(i);
    (*dxds_star(coeff_s_idx_i))[i] = dxds_exact( t, gamma[i], coeff_s(i) );
    // Note: Above, at least the column access will be validated in debug mode
    // but the row index i will not ever be.  Perhaps we can augment Epetra to
    // fix this?
  }
  return (dxds_star_ptr);
}


// Overridden from ParameterListAcceptor


void DiagonalTransientModel::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& paramList
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
  gamma_min_ = get<double>(*paramList_,Gamma_min_name);
  gamma_max_ = get<double>(*paramList_,Gamma_max_name);
  coeff_s_ = getArrayFromStringParameter<double>(*paramList_,Coeff_s_name);
  gamma_fit_ = getIntegralValue<EGammaFit>(*paramList_,Gamma_fit_name);
  x0_ = get<double>(*paramList_,x0_name);
  exactSolutionAsResponse_ = get<bool>(*paramList_,ExactSolutionAsResponse_name);
  initialize();
}


Teuchos::RCP<Teuchos::ParameterList> 
DiagonalTransientModel::getNonconstParameterList()
{
  return paramList_;
}


Teuchos::RCP<Teuchos::ParameterList>
DiagonalTransientModel::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}


Teuchos::RCP<const Teuchos::ParameterList>
DiagonalTransientModel::getParameterList() const
{
  return paramList_;
}


Teuchos::RCP<const Teuchos::ParameterList>
DiagonalTransientModel::getValidParameters() const
{
  using Teuchos::RCP; using Teuchos::ParameterList;
  using Teuchos::tuple;
  using Teuchos::setIntParameter; using Teuchos::setDoubleParameter;
  using Teuchos::setStringToIntegralParameter;
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set(Implicit_name,true);
    setDoubleParameter(
      Gamma_min_name, Gamma_min_default, "",
      &*pl
      );
    setDoubleParameter(
      Gamma_max_name, Gamma_max_default, "",
      &*pl
      );
    setStringToIntegralParameter<EGammaFit>(
      Gamma_fit_name, Gamma_fit_default, "",
      tuple<std::string>("Linear","Random"),
      tuple<EGammaFit>(GAMMA_FIT_LINEAR,GAMMA_FIT_RANDOM),
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
    pl->set( ExactSolutionAsResponse_name, ExactSolutionAsResponse_default );
    validPL = pl;
  }
  return validPL;
}


// Overridden from EpetraExt::ModelEvaluator


Teuchos::RCP<const Epetra_Map>
DiagonalTransientModel::get_x_map() const
{
  return epetra_map_;
}


Teuchos::RCP<const Epetra_Map>
DiagonalTransientModel::get_f_map() const
{
  return epetra_map_;
}


Teuchos::RCP<const Epetra_Map>
DiagonalTransientModel::get_p_map(int l) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( l, 0, Np_ );
#endif
  return map_p_[l];
}


Teuchos::RCP<const Teuchos::Array<std::string> >
DiagonalTransientModel::get_p_names(int l) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( l, 0, Np_ );
#endif
  return names_p_[l];
}


Teuchos::RCP<const Epetra_Map>
DiagonalTransientModel::get_g_map(int j) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( j, 0, Ng_ );
#endif
  return map_g_[j];
}


Teuchos::RCP<const Epetra_Vector>
DiagonalTransientModel::get_x_init() const
{
  return x_init_;
}


Teuchos::RCP<const Epetra_Vector>
DiagonalTransientModel::get_x_dot_init() const
{
  return x_dot_init_;
}


Teuchos::RCP<const Epetra_Vector>
DiagonalTransientModel::get_p_init(int l) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT( l == 0 );
#endif
  return p_init_[l];
}


Teuchos::RCP<Epetra_Operator>
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
  outArgs.set_Np_Ng(Np_,Ng_);
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
      DERIV_LINEARITY_NONCONST,
      DERIV_RANK_DEFICIENT,
      true // supportsAdjoint
      )
    );
  if (exactSolutionAsResponse_) {
    outArgs.setSupports(OUT_ARG_DgDp,0,0,DERIV_MV_BY_COL);
    outArgs.set_DgDp_properties(
      0,0,DerivativeProperties(
        DERIV_LINEARITY_NONCONST,
        DERIV_RANK_DEFICIENT,
        true // supportsAdjoint
        )
      );
  }
  return outArgs;
}


void DiagonalTransientModel::evalModel(
  const InArgs& inArgs, const OutArgs& outArgs
  ) const
{

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::null;
  using Teuchos::dyn_cast;

  const Epetra_Vector &x = *(inArgs.get_x());
  const double t = inArgs.get_t();
  if (Np_) set_coeff_s_p(inArgs.get_p(0)); // Sets for coeff_s(...) function!
  UnsetParameterVector unsetParameterVector(rcp(&coeff_s_p_,false));
  // Note: Above, the destructor for unsetParameterVector will ensure that the
  // RCP to the parameter vector will be unset no matter if an exception is
  // thrown or not.

  Epetra_Operator *W_out = ( implicit_ ? outArgs.get_W().get() : 0 );
  Epetra_Vector *f_out = outArgs.get_f().get();
  Epetra_MultiVector *DfDp_out = 0;
  if (Np_) DfDp_out = get_DfDp_mv(0,outArgs).get();
  Epetra_Vector *g_out = 0;
  Epetra_MultiVector *DgDp_out = 0;
  if (exactSolutionAsResponse_) {
    g_out = outArgs.get_g(0).get();
    DgDp_out = get_DgDp_mv(0,0,outArgs,DERIV_MV_BY_COL).get();
  }

  const Epetra_Vector &gamma = *gamma_;

  int localNumElements = x.MyLength();

  if (f_out) {
    Epetra_Vector &f = *f_out;
    if (implicit_) {
      const Epetra_Vector *x_dot_in = inArgs.get_x_dot().get();
      if (x_dot_in) {
        const Epetra_Vector &x_dot = *x_dot_in;
        for ( int i=0 ; i<localNumElements ; ++i )
          f[i] = x_dot[i] - f_func(x[i],t,gamma[i],coeff_s(i));
      }
      else {
        for ( int i=0 ; i<localNumElements ; ++i )
          f[i] = - f_func(x[i],t,gamma[i],coeff_s(i));
      }
    }
    else {
      for ( int i=0 ; i<localNumElements ; ++i )
        f[i] = f_func(x[i],t,gamma[i],coeff_s(i));
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
      values[0] = alpha - beta*gamma[i];
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
          = - d_evalR_d_s(t,gamma[i],coeff_s(i));
      }
    }
    else {
      for( int i = 0; i < localNumElements; ++i ) {
        DfDp[coeff_s_idx(i)][i]
          = + d_evalR_d_s(t,gamma[i],coeff_s(i));
      }
    }
  }

  if (g_out) {
    *g_out = *getExactSolution(t,&*coeff_s_p_);
    // Note: Above will wipe out coeff_s_p_ as a side effect!
  }

  if (DgDp_out) {
    *DgDp_out = *getExactSensSolution(t,&*coeff_s_p_);
    // Note: Above will wipe out coeff_s_p_ as a side effect!
  }

  // Warning: From here on out coeff_s_p_ is unset!
  
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
  // Create gamma
  //

  gamma_ = Teuchos::rcp(new Epetra_Vector(*epetra_map_));
  Epetra_Vector &gamma = *gamma_;
  switch(gamma_fit_) {
    case GAMMA_FIT_LINEAR: {
      const int N = gamma.GlobalLength();
      const double slope = (gamma_max_ - gamma_min_)/(N-1);
      const int MyLength = gamma.MyLength();
      if (1==MyLength) {
        gamma[0] = gamma_min_;
      }
      else {
        for ( int i=0 ; i<MyLength ; ++i )
        {
          gamma[i] = slope*(procRank*MyLength+i)+gamma_min_;
        }
      }
      break;
    }
    case GAMMA_FIT_RANDOM: {
      unsigned int seed = Teuchos::ScalarTraits<int>::random()+10*procRank; 
      seed *= seed;
      gamma.SetSeed(seed);
      gamma.Random(); // fill with random numbers in (-1,1)
      // Scale random numbers to (gamma_min_,gamma_max_)
      const double slope = (gamma_min_ - gamma_max_)/2.0;
      const double tmp = (gamma_max_ + gamma_min_)/2.0;
      int MyLength = gamma.MyLength();
      for (int i=0 ; i<MyLength ; ++i)
      {
        gamma[i] = slope*gamma[i] + tmp;
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
  names_p_.clear();
  {
    Teuchos::RCP<Teuchos::Array<std::string> >
      names_p_0 = Teuchos::rcp(new Teuchos::Array<std::string>());
    for ( int i = 0; i < np_; ++i )
      names_p_0->push_back("coeff_s("+Teuchos::toString(i)+")");
    names_p_.push_back(names_p_0);
  }

  coeff_s_idx_.clear();
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
  // Setup exact solution as response function
  //

  if (exactSolutionAsResponse_) {
    Ng_ = 1;
    map_g_.clear();
    map_g_.push_back(
      rcp( new Epetra_LocalMap(1,0,*epetra_comm_) )
      );
  }
  else {
    Ng_ = 0;
  }

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
  p_init_.clear();
  p_init_.push_back(
    rcp( new Epetra_Vector( ::Copy, *map_p_[0], &coeff_s_[0] ) )
    );

}


void DiagonalTransientModel::set_coeff_s_p( 
  const Teuchos::RCP<const Epetra_Vector> &coeff_s_p
  ) const
{
  if (!is_null(coeff_s_p))
    coeff_s_p_ = coeff_s_p;
  else
    unset_coeff_s_p();
}


void DiagonalTransientModel::unset_coeff_s_p() const
{
  using Teuchos::as;
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT(
    as<int>(get_p_map(0)->NumGlobalElements()) == as<int>(coeff_s_.size()) );
#endif
  coeff_s_p_ = Teuchos::rcp(
    new Epetra_Vector(
      ::View,
      *get_p_map(0),
      const_cast<double*>(&coeff_s_[0])
      )
    );
  // Note: The above const cast is okay since the coeff_s_p_ RCP is to a const
  // Epetr_Vector!
}



} // namespace EpetraExt


// Nonmembers


Teuchos::RCP<EpetraExt::DiagonalTransientModel>
EpetraExt::diagonalTransientModel(
  Teuchos::RCP<Epetra_Comm> const& epetra_comm,
  Teuchos::RCP<Teuchos::ParameterList> const& paramList
  )
{
  Teuchos::RCP<DiagonalTransientModel> diagonalTransientModel =
    Teuchos::rcp(new DiagonalTransientModel(epetra_comm));
  if (!is_null(paramList))
    diagonalTransientModel->setParameterList(paramList);
  return diagonalTransientModel;
}
