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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "VanDerPolOscillator.hpp"

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#  include "mpi.h"
#else
#  include "Epetra_SerialComm.h"
#endif // HAVE_MPI

#include "Epetra_LocalMap.h"
#include "Epetra_CrsMatrix.h"

#include "Sacado.hpp"

VanDerPolOscillator::VanDerPolOscillator(Teuchos::RCP<Epetra_Comm> &epetra_comm_ptr, Teuchos::ParameterList &params)
{
  implicit_ = params.get<bool>( "implicit" );
  omega_ = params.get<double>( "omega" );
  string outfile_name = params.get<string>( "Output File Name" );
  outfile_.open(outfile_name.c_str());

  // Construct a replicated local map with 2 elements
  epetra_comm_ptr_ = epetra_comm_ptr;
  epetra_map_ptr_ = 
    Teuchos::rcp( new Epetra_LocalMap(2, 0, *epetra_comm_ptr_) );

  // Construct initial solution vector
  x0_ = Teuchos::rcp(new Epetra_Vector(*epetra_map_ptr_));
  (*x0_)[0] = params.get<double>( "x0_1" );
  (*x0_)[1] = params.get<double>( "x0_2" );

  // Construct map for Jacobian
  if (implicit_) {
    const int dim=2;
    W_graph_ = Teuchos::rcp(new Epetra_CrsGraph(::Copy,*epetra_map_ptr_,dim));
    int indices[] = {0, 1};
    W_graph_->InsertGlobalIndices(0, dim, indices);
    W_graph_->InsertGlobalIndices(1, dim, indices);
    W_graph_->FillComplete();
  }
}

VanDerPolOscillator::~VanDerPolOscillator()
{
  outfile_.close();
}

void
VanDerPolOscillator::saveSolution(const Epetra_Vector& x, double t)
{
  outfile_ << t << " " << x[0] << " " << x[1] << std::endl;
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
VanDerPolOscillator::get_x_map() const
{
  return epetra_map_ptr_;
}

Teuchos::RCP<const Epetra_Map>
VanDerPolOscillator::get_f_map() const
{
  return epetra_map_ptr_;
}

Teuchos::RCP<const Epetra_Vector>
VanDerPolOscillator::get_x_init() const
{
  return x0_;
}

Teuchos::RCP<Epetra_Operator>
VanDerPolOscillator::create_W() const
{
  if (implicit_)
    return Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*W_graph_));
  else
    return Teuchos::null;
}

EpetraExt::ModelEvaluator::InArgs
VanDerPolOscillator::createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.setSupports(IN_ARG_x_poly,true);
  if(implicit_) {
    inArgs.setSupports(IN_ARG_x_dot,true);
    inArgs.setSupports(IN_ARG_alpha,true);
    inArgs.setSupports(IN_ARG_beta,true);
  }
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
VanDerPolOscillator::createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_f_poly,true);
  if(implicit_) {
    outArgs.setSupports(OUT_ARG_W,true);
  }
  return outArgs;
}

void VanDerPolOscillator::evalModel( const InArgs& inArgs, 
				     const OutArgs& outArgs ) const
{
  // compute f(x)
  Teuchos::RCP<const Epetra_Vector> x = inArgs.get_x();
  Teuchos::RCP<Epetra_Vector> f = outArgs.get_f();
  if ( (x != Teuchos::null) && (f != Teuchos::null) ) {
    evalVField((*x)[0],(*x)[1],(*f)[0],(*f)[1]);
  }

  // compute f([x])
  Teuchos::RCP<const Teuchos::Polynomial<Epetra_Vector> > x_poly = 
    inArgs.get_x_poly();
  Teuchos::RCP<Teuchos::Polynomial<Epetra_Vector> > f_poly = 
    outArgs.get_f_poly();
  if ( (x_poly != Teuchos::null) && (f_poly != Teuchos::null) ) {
    unsigned int d = x_poly->degree();
    Sacado::Tay::Taylor<double> x1(d,0.0);
    Sacado::Tay::Taylor<double> x2(d,0.0);
    Sacado::Tay::Taylor<double> f1(d,0.0);
    Sacado::Tay::Taylor<double> f2(d,0.0);

    for (unsigned int i=0; i<=d; i++) {
      x1.fastAccessCoeff(i) = (*(x_poly->getCoefficient(i)))[0];
      x2.fastAccessCoeff(i) = (*(x_poly->getCoefficient(i)))[1];
    }

    evalVField(x1,x2,f1,f2);

    for (unsigned int i=0; i<=d; i++) {
      (*(f_poly->getCoefficient(i)))[0] = f1.coeff(i);
      (*(f_poly->getCoefficient(i)))[1] = f2.coeff(i);
    }
  }

  // compute W
  Teuchos::RCP<Epetra_Operator> W = outArgs.get_W();
  if (W != Teuchos::null) {
    const double alpha = inArgs.get_alpha();
    const double beta = inArgs.get_beta();
    Epetra_CrsMatrix &crsW = Teuchos::dyn_cast<Epetra_CrsMatrix>(*W);
    const int dim = 2;
    double values_1[2];
    double values_2[2];
    int indices[] = {0,1};

    Sacado::Fad::DFad<double> x1(dim,0,(*x)[0]);
    Sacado::Fad::DFad<double> x2(dim,1,(*x)[1]);
    Sacado::Fad::DFad<double> f1;
    Sacado::Fad::DFad<double> f2;

    evalVField(x1,x2,f1,f2);

    values_1[0] = alpha * f1.fastAccessDx(0) - beta;
    values_1[1] = alpha * f1.fastAccessDx(1);
    values_2[0] = alpha * f2.fastAccessDx(0);
    values_2[0] = alpha * f2.fastAccessDx(1) - beta;
   
    crsW.ReplaceGlobalValues(0,dim,values_1,indices);
    crsW.ReplaceGlobalValues(1,dim,values_2,indices);
  }
}

