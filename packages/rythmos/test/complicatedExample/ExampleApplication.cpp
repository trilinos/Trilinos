//@HEADER

// ***********************************************************************
//
//                     Rythmos Package
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

#include "ExampleApplication.hpp"
#include "Teuchos_ScalarTraits.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif // HAVE_MPI

#include "Epetra_CrsMatrix.h"

#ifdef EXAMPLEAPPLICATION_DEBUG
#include <iostream>
#endif

ExampleApplication::ExampleApplication(Teuchos::RCP<Epetra_Comm> &epetra_comm_ptr, Teuchos::ParameterList &params)
{
  implicit_ = params.get<bool>( "Implicit" );
  lambda_min_ = params.get<double>( "Lambda_min" );
  lambda_max_ = params.get<double>( "Lambda_max" );
  lambda_fit_ = params.get<std::string>( "Lambda_fit" );
  numElements_ = params.get<int>( "NumElements" );
  x0_ = params.get<double>( "x0" );
  coeff_s_ = params.get<double>( "Coeff_s" );

  // Construct a Map with NumElements and index base of 0
  epetra_map_ptr_ = Teuchos::rcp( new Epetra_Map(numElements_, 0, *epetra_comm_ptr) );

  lambda_ptr_ = Teuchos::rcp(new Epetra_Vector(*epetra_map_ptr_));
  Epetra_Vector &lambda = *lambda_ptr_;
  if ( lambda_fit_ == "linear" )
  {
    int N = lambda.GlobalLength();
    double tmp;
    if (N == 1)
      tmp = 0; // lambda[0] = lambda_min_
    else
      tmp = (lambda_max_ - lambda_min_)/(N-1);
    int MyLength = lambda.MyLength();
    int MyPID = lambda.Comm().MyPID();
    for (int i=0 ; i<MyLength ; ++i)
    {
      lambda[i] = tmp*(MyPID*MyLength+i)+lambda_min_;
    }
  }
  else // if ( lambda_fit_ == "random" )
  {
    int MyPID = lambda.Comm().MyPID();
    unsigned int seed = Teuchos::ScalarTraits<int>::random()+10*MyPID; 
    seed *= seed;
    lambda.SetSeed(seed);
    lambda.Random(); // fill with random numbers in (-1,1)
    // Scale random numbers to (lambda_min_,lambda_max_)
    double slope = (lambda_min_ - lambda_max_)/2.0;
    double tmp = (lambda_max_ + lambda_min_)/2.0;
    int MyLength = lambda.MyLength();
    for (int i=0 ; i<MyLength ; ++i)
    {
      lambda[i] = slope*lambda[i] + tmp;
    }
  }

  if(implicit_) {
    int localNumElements = lambda.MyLength();
    W_graph_ = Teuchos::rcp(new Epetra_CrsGraph(::Copy,*epetra_map_ptr_,1));
    int indices[1];
    const int IB = epetra_map_ptr_->IndexBase();
    for( int i = 0; i < localNumElements; ++i ) {
      indices[0] = i + IB;  // global column
      W_graph_->InsertGlobalIndices(
        i + IB              // GlobalRow
        ,1                  // NumEntries
        ,indices            // Indices
        );
    }
    W_graph_->FillComplete();
  }
 
}

Teuchos::RCP<const Epetra_Vector> ExampleApplication::get_coeff() const
{
  return(lambda_ptr_);
}

Teuchos::RCP<const Epetra_Vector> ExampleApplication::getExactSolution(double t) const
{
  Teuchos::RCP<Epetra_Vector> x_star_ptr = Teuchos::rcp(new Epetra_Vector(*epetra_map_ptr_));
  Epetra_Vector& x_star = *x_star_ptr;
  Epetra_Vector& lambda = *lambda_ptr_;
  x_star.PutScalar(0.0);
  x_star.Scale(t,lambda);
  int myN = x_star.MyLength();
  if (coeff_s_ == 0.0)
  {
    for ( int i=0 ; i<myN ; ++i )
      x_star[i] = x0_*exp(x_star[i]);
  }
  else
  {
    for ( int i=0 ; i<myN ; ++i )
      x_star[i] = (x0_+(1.0/coeff_s_))*exp(x_star[i]) - exp(x_star[i])*cos(t*coeff_s_)/coeff_s_;
  }
  return(x_star_ptr);
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
ExampleApplication::get_x_map() const
{
  return epetra_map_ptr_;
}

Teuchos::RCP<const Epetra_Map>
ExampleApplication::get_f_map() const
{
  return epetra_map_ptr_;
}

Teuchos::RCP<const Epetra_Vector>
ExampleApplication::get_x_init() const
{
  Teuchos::RCP<Epetra_Vector> x_init = Teuchos::rcp(new Epetra_Vector(*epetra_map_ptr_));
  x_init->PutScalar(x0_);
  return x_init;
}

Teuchos::RCP<const Epetra_Vector>
ExampleApplication::get_x_dot_init() const
{
  Teuchos::RCP<Epetra_Vector> x_dot_init = Teuchos::rcp(new Epetra_Vector(*epetra_map_ptr_));
  x_dot_init->PutScalar(0.0);
  return x_dot_init;
}

Teuchos::RCP<Epetra_Operator>
ExampleApplication::create_W() const
{
  if(implicit_)
    return Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*W_graph_));
  return Teuchos::null;
}

EpetraExt::ModelEvaluator::InArgs
ExampleApplication::createInArgs() const
{
  InArgsSetup inArgs;
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
ExampleApplication::createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setSupports(OUT_ARG_f,true);
  if(implicit_) {
    outArgs.setSupports(OUT_ARG_W,true);
  }
  return outArgs;
}

void ExampleApplication::evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const
{
  const Epetra_Vector &x = *(inArgs.get_x());
  const double t = inArgs.get_t();
  const Epetra_Vector &lambda = *lambda_ptr_;
#ifdef EXAMPLEAPPLICATION_DEBUG
      std::cout << "----------------------------------------------------------------------" << std::endl;
      std::cout << "ExampleApplication::evalModel x = " << std::endl;
      x.Print(std::cout);
      std::cout << "ExampleApplication::evalModel lambda = " << std::endl;
      lambda.Print(std::cout);
#endif
  int localNumElements = x.MyLength();
  if(implicit_) 
  {
    const Epetra_Vector &x_dot = *inArgs.get_x_dot();
    if(outArgs.get_f().get()) 
    {
      Epetra_Vector &f = *outArgs.get_f();
      for (int i=0 ; i<localNumElements ; ++i)
      {
        f[i] = x_dot[i] - lambda[i]*x[i] - evalR(t,lambda[i],coeff_s_);
      }
#ifdef EXAMPLEAPPLICATION_DEBUG
      std::cout << "ExampleApplication::evalModel (implicit) x_dot = " << std::endl;
      x_dot.Print(std::cout);
      std::cout << "ExampleApplication::evalModel (implicit) f = " << std::endl;
      f.Print(std::cout);
#endif
    }
    Teuchos::RCP<Epetra_Operator> W;
    if( (W = outArgs.get_W()).get() ) {
      const double alpha = inArgs.get_alpha();
      const double beta = inArgs.get_beta();
      Epetra_CrsMatrix &crsW = Teuchos::dyn_cast<Epetra_CrsMatrix>(*W);
      double values[1];
      int indices[1];
      const int IB = epetra_map_ptr_->IndexBase();
      for( int i = 0; i < localNumElements; ++i ) 
      {
        values[0] = alpha - beta*lambda[i];
        indices[0] = i + IB;  // global column
        crsW.ReplaceGlobalValues(i + IB              // GlobalRow
                                 ,1                  // NumEntries
                                 ,values             // Values
                                 ,indices            // Indices
                                          );
      }
#ifdef EXAMPLEAPPLICATION_DEBUG
      std::cout << "ExampleApplication::evalModel (implicit) alpha, beta = " << std::endl;
      std::cout << "alpha = " << alpha << std::endl;
      std::cout << "beta = "  << beta  << std::endl;
      std::cout << "ExampleApplication::evalModel (implicit) W = " << std::endl;
      crsW.Print(std::cout);
#endif
    }
  }
  else 
  {
    Epetra_Vector &f = *outArgs.get_f();
    for (int i=0 ; i<localNumElements ; ++i)
    {
      f[i] = lambda[i]*x[i]+evalR(t,lambda[i],coeff_s_);
    }
#ifdef EXAMPLEAPPLICATION_DEBUG
    std::cout << "ExampleApplication::evalModel (explicit) f = " << std::endl;
    f.Print(std::cout);
#endif
  }
#ifdef EXAMPLEAPPLICATION_DEBUG
  std::cout << "----------------------------------------------------------------------" << std::endl;
#endif
}

double ExampleApplication::evalR(const double& t, const double& lambda, const double& s) const
{
  return(exp(lambda*t)*sin(s*t));
}

