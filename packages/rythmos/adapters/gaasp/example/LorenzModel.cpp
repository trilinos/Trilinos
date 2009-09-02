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

#include "LorenzModel.hpp"
#include "Teuchos_ScalarTraits.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif // HAVE_MPI

#include "Epetra_CrsMatrix.h"

#ifdef LORENZMODEL_DEBUG
#include <iostream>
#endif

LorenzModel::LorenzModel(Teuchos::RCP<Epetra_Comm> &epetra_comm_ptr_, Teuchos::ParameterList &params)
{
  param0 = params.get<double>( "Parameter 0", 10.0 );
  param1 = params.get<double>( "Parameter 1", 28.0 );
  param2 = params.get<double>( "Parameter 2", 8.0/3.0 );
  ic0 = params.get<double>( "IC 0", -6.9742 );
  ic1 = params.get<double>( "IC 1", -7.008 );
  ic2 = params.get<double>( "IC 2", 25.1377 );
  numElements_ = 3;

  // Construct a Map with NumElements and index base of 0
  epetra_map_ptr_ = Teuchos::rcp( new Epetra_Map(numElements_, 0, *epetra_comm_ptr_) );

}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
LorenzModel::get_x_map() const
{
  return epetra_map_ptr_;
}

Teuchos::RCP<const Epetra_Map>
LorenzModel::get_f_map() const
{
  return epetra_map_ptr_;
}

Teuchos::RCP<const Epetra_Vector>
LorenzModel::get_x_init() const
{
  Teuchos::RCP<Epetra_Vector> x_init = Teuchos::rcp(new Epetra_Vector(*epetra_map_ptr_));
  Epetra_Vector& x = *x_init;
  x[0] = ic0;
  x[1] = ic1;
  x[2] = ic2;
  return x_init;
}

Teuchos::RCP<const Epetra_Vector>
LorenzModel::get_x_dot_init() const
{
  return(Teuchos::null);
}

Teuchos::RCP<Epetra_Operator>
LorenzModel::create_W() const
{
  return(Teuchos::null);
}

EpetraExt::ModelEvaluator::InArgs
LorenzModel::createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setSupports(IN_ARG_t,true);
  inArgs.setSupports(IN_ARG_x,true);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
LorenzModel::createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setSupports(OUT_ARG_f,true);
  return outArgs;
}

void LorenzModel::evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const
{
  const Epetra_Vector &yin = *(inArgs.get_x());
  const double t = inArgs.get_t(); // ignored
#ifdef LORENZMODEL_DEBUG
      std::cout << "----------------------------------------------------------------------" << std::endl;
      std::cout << "LorenzModel::evalModel yin = " << std::endl;
      yin.Print(std::cout);
#endif
  Epetra_Vector &yout = *outArgs.get_f();

  yout[0] = -param0 * yin[0] + param0 * yin[1];
  yout[1] = param1 * yin[0] - yin[1] - yin[0]*yin[2];
  yout[2] = -param2*yin[2] + yin[0]*yin[1];

#ifdef LORENZMODEL_DEBUG
  std::cout << "LorenzModel::evalModel (explicit) f = " << std::endl;
  yout.Print(std::cout);
#endif
#ifdef LORENZMODEL_DEBUG
  std::cout << "----------------------------------------------------------------------" << std::endl;
#endif
}


