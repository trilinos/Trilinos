// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// 
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include "MockModelEval_B.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

MockModelEval_B::MockModelEval_B(const MPI_Comm appComm) 
{

#ifdef HAVE_MPI
    Comm = rcp(new Epetra_MpiComm(appComm));
#else
    Comm = rcp(new Epetra_SerialComm);
#endif

    //set up map and initial guess for solution vector
    const int vecLength = 1;
    x_map = rcp(new Epetra_Map(vecLength, 0, *Comm));
    x_vec = rcp(new Epetra_Vector(*x_map));
    x_vec->PutScalar(0.0);
    x_dot_vec = rcp(new Epetra_Vector(*x_map));
    x_dot_vec->PutScalar(1.0);

    //set up responses
    const int numResponses = 1;
    g_map = rcp(new Epetra_LocalMap(numResponses, 0, *Comm));

    //set up parameters
    const int numParameters= 1;
    p_map = rcp(new Epetra_LocalMap(numParameters, 0, *Comm));
    p_init = rcp(new Epetra_Vector(*p_map));
    for (int i=0; i<numParameters; i++) (*p_init)[i]= 1.0;

}

MockModelEval_B::~MockModelEval_B()
{
}

RCP<const Epetra_Map> MockModelEval_B::get_x_map() const
{
  return x_map;
}

RCP<const Epetra_Map> MockModelEval_B::get_f_map() const
{
  return x_map;
}

RCP<Epetra_Operator> MockModelEval_B::create_W() const
{
  return Teuchos::null;
}

RCP<const Epetra_Map> MockModelEval_B::get_p_map(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);
  return p_map;
}

RCP<const  Teuchos::Array<std::string> > MockModelEval_B::get_p_names(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_names() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);

  Teuchos::Ordinal num_p = p_init->MyLength();
  RCP<Teuchos::Array<std::string> > p_names = 
      rcp(new Teuchos::Array<std::string>(num_p) );
  for (int i=0; i<num_p; i++) {
    std::stringstream ss;
    ss << "Parameter " << i;
    const std::string name = ss.str();
    (*p_names)[i] = name;
  }
  return p_names;
}


RCP<const Epetra_Map> MockModelEval_B::get_g_map(int j) const
{
  TEST_FOR_EXCEPTION(j != 0, std::logic_error,
                     std::endl <<
                     "Error!  MockModelEval_B::get_g_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     j << std::endl);
  return g_map;
}

RCP<const Epetra_Vector> MockModelEval_B::get_x_init() const
{
  return x_vec;
}

RCP<const Epetra_Vector> MockModelEval_B::get_x_dot_init() const
{
  return x_dot_vec;
}

RCP<const Epetra_Vector> MockModelEval_B::get_p_init(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, std::logic_error,
                     std::endl <<
                     "Error!  App::ModelEval::get_p_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " <<
                     l << std::endl);
  return p_init;
}

EpetraExt::ModelEvaluator::InArgs MockModelEval_B::createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(1);
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.setSupports(IN_ARG_t,true);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs MockModelEval_B::createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(1, 1);

  outArgs.setSupports(OUT_ARG_f,true);
  return outArgs;
}

void MockModelEval_B::evalModel( const InArgs& inArgs,
                              const OutArgs& outArgs ) const
{

  // Parse InArgs
  RCP<const Epetra_Vector> p_in = inArgs.get_p(0);
  if (!p_in.get()) cout << "ERROR: MockModelEval_B requires p as inargs" << endl;

  RCP<const Epetra_Vector> x_in = inArgs.get_x();
  if (!x_in.get()) cout << "ERROR: MockModelEval_B requires x as inargs" << endl;
  int myVecLength = x_in->MyLength();

  // Parse OutArgs

  RCP<Epetra_Vector> f_out = outArgs.get_f(); 
  RCP<Epetra_Vector> g_out = outArgs.get_g(0); 

  if (f_out != Teuchos::null) {
    for (int i=0; i<myVecLength; i++) {
      int gid = x_in->Map().GID(i);
       (*f_out)[i] = 0.0*( (*x_in)[i] - gid ) - (*p_in)[0];
    }
  }

  if (g_out != Teuchos::null) x_in->MeanValue(&(*g_out)[0]);
} 
