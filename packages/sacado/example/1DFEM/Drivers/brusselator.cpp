// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "CZeroDiscretization.hpp"
#include "BrusselatorPDE.hpp"
#include "ConstantDirichletBC.hpp"
#include "GaussianQuadrature2.hpp"
#include "Application.hpp"

#include <iostream>

#include "Sacado_Fad_DFad.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

typedef Sacado::Fad::DFad<double> FadType;

int main(int argc, char *argv[]) {
  unsigned int nelem = 100;
  double h = 1.0/nelem;
  double alpha = 0.25;
  double beta = 0.0;
  double D1 = 1.0/40.0;
  double D2 = 1.0/40.0;

  try {

    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init(&argc,&argv);
#endif

    Teuchos::RefCountPtr<Epetra_Comm> Comm;

    // Create a communicator for Epetra objects
#ifdef HAVE_MPI
    Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    Comm = Teuchos::rcp(new Epetra_SerialComm);
#endif

    Teuchos::RefCountPtr< AbstractPDE<double> > real_pde = 
      Teuchos::rcp(new BrusselatorPDE<double>(alpha, beta, D1, D2));

    Teuchos::RefCountPtr< AbstractPDE<FadType> > fad_pde = 
      Teuchos::rcp(new BrusselatorPDE<FadType>(alpha, beta, D1, D2));

    Teuchos::RefCountPtr<AbstractQuadrature> quad = 
      Teuchos::rcp(new GaussianQuadrature2);

    real_pde->init(quad->numPoints(), 2);
    fad_pde->init(quad->numPoints(), 2);

    vector<double> x(nelem+1);
    for (unsigned int i=0; i<=nelem; i++)
      x[i] = h*i;

    CZeroDiscretization disc(x, real_pde->numEquations(), Comm);
    disc.createMesh();
    disc.createMaps();
    disc.createJacobianGraphs();

    Teuchos::RefCountPtr<Epetra_Map> map = disc.getMap();

    std::vector< Teuchos::RefCountPtr<const AbstractBC> > bc(4);
    bc[0] = Teuchos::rcp(new ConstantDirichletBC(map->MinAllGID(),
						 alpha));
    bc[1] = Teuchos::rcp(new ConstantDirichletBC(map->MinAllGID()+1,
						 beta/alpha));
    bc[2] = Teuchos::rcp(new ConstantDirichletBC(map->MaxAllGID()-1,
						 alpha));
    bc[3] = Teuchos::rcp(new ConstantDirichletBC(map->MaxAllGID(),
						 beta/alpha));

    Application app(Teuchos::rcp(&disc, false), bc, quad);

    Epetra_Vector u(*map);
    u.PutScalar(1.0);
    Epetra_Vector f(*map);
    Epetra_CrsMatrix jac(Copy, (*disc.getJacobianGraph()));

    app.computeGlobalJacobian(*fad_pde, u, f, jac);

    f.Print(std::cout);
    jac.Print(std::cout);

#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif

  }
  
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
  catch (string& s) {
    std::cout << s << std::endl;
  }
  catch (char *s) {
    std::cout << s << std::endl;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" <<std:: endl;
  }

}
