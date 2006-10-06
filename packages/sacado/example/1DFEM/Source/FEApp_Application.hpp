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

#ifndef FEAPP_APPLICATION_HPP
#define FEAPP_APPLICATION_HPP

#include <vector>

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"

#include "FEApp_AbstractBC.hpp"
#include "FEApp_AbstractPDE.hpp"
#include "FEApp_QuadratureFactory.hpp"
#include "FEApp_DiscretizationFactory.hpp"

#include "Sacado_Fad_DFad.hpp"
#include "Sacado_TemplateManager.hpp"
#include "Sacado_MPL_vector.hpp"

typedef double RealType;
typedef Sacado::Fad::DFad<double> FadType;
typedef Sacado::mpl::vector<RealType, FadType> ValidTypes;

namespace FEApp {

  class Application {
  public:

    //! Constructor 
    template <typename BuilderT>
    Application(
		const BuilderT& pdeBuilder,
		const std::vector<double>& coords,
		unsigned int num_equations,
		const Teuchos::RefCountPtr<const Epetra_Comm>& comm,
		const Teuchos::RefCountPtr<Teuchos::ParameterList>& params);

    //! Destructor
    ~Application();

    //! Get DOF map
    Teuchos::RefCountPtr<Epetra_Map> getMap();

    //! Get Jacobian graph
    Teuchos::RefCountPtr<Epetra_CrsGraph> getJacobianGraph();

    //! Set boundary conditions
    void setBCs(
       const std::vector< Teuchos::RefCountPtr<const FEApp::AbstractBC> > bcs);

    //! Compute global residual
    void computeGlobalResidual(const Epetra_Vector& x,
			       Epetra_Vector& f);

    //! Compute global Jacobian
    void computeGlobalJacobian(const Epetra_Vector& x,
			       Epetra_Vector& f,
			       Epetra_CrsMatrix& jac);

  private:
    
    //! Private to prohibit copying
    Application(const Application&);

    //! Private to prohibit copying
    Application& operator=(const Application&);

  protected:
    
    //! Element discretization
    Teuchos::RefCountPtr<FEApp::AbstractDiscretization> disc;

    //! Boundary conditions
    std::vector< Teuchos::RefCountPtr<const FEApp::AbstractBC> > bc;

    //! Quadrature rule
    Teuchos::RefCountPtr<const FEApp::AbstractQuadrature> quad;

    //! PDE equations
    FEApp::AbstractPDE_TemplateManager<ValidTypes> pdeTM;

    //! Importer for overlapped data
    Teuchos::RefCountPtr<Epetra_Import> importer;

    //! Exporter for overlapped data
    Teuchos::RefCountPtr<Epetra_Export> exporter;

    //! Overlapped solution vector
    Teuchos::RefCountPtr<Epetra_Vector> overlapped_x;

    //! Overlapped residual vector
    Teuchos::RefCountPtr<Epetra_Vector> overlapped_f;

    //! Overlapped Jacobian matrix
    Teuchos::RefCountPtr<Epetra_CrsMatrix> overlapped_jac;

  };

}

template <typename BuilderT>
FEApp::Application::Application(
		   const BuilderT& pdeBuilder,
		   const std::vector<double>& coords,
		   unsigned int num_equations,
		   const Teuchos::RefCountPtr<const Epetra_Comm>& comm,
		   const Teuchos::RefCountPtr<Teuchos::ParameterList>& params) 
{
  // Create quadrature object
  Teuchos::RefCountPtr<Teuchos::ParameterList> quadParams = 
    Teuchos::rcp(&(params->sublist("Quadrature")),false);
  FEApp::QuadratureFactory quadFactory(quadParams);
  quad = quadFactory.create();

  // Create discretization object
  Teuchos::RefCountPtr<Teuchos::ParameterList> discParams = 
    Teuchos::rcp(&(params->sublist("Discretization")),false);
  FEApp::DiscretizationFactory discFactory(discParams);
  disc = discFactory.create(coords, num_equations, comm);
  disc->createMesh();
  disc->createMaps();
  disc->createJacobianGraphs();

  // Create Epetra objects
  importer = Teuchos::rcp(new Epetra_Import(*(disc->getOverlapMap()), 
					    *(disc->getMap())));
  exporter = Teuchos::rcp(new Epetra_Export(*(disc->getOverlapMap()), 
					    *(disc->getMap())));
  overlapped_x = Teuchos::rcp(new Epetra_Vector(*(disc->getOverlapMap())));
  overlapped_f = Teuchos::rcp(new Epetra_Vector(*(disc->getOverlapMap())));
  overlapped_jac = 
    Teuchos::rcp(new Epetra_CrsMatrix(Copy, 
				      *(disc->getOverlapJacobianGraph())));

  // Initialize PDEs
  typedef FEApp::AbstractPDE_TemplateManager<ValidTypes>::iterator iterator;
  pdeTM.buildObjects(pdeBuilder);
  int nqp = quad->numPoints();
  int nn = disc->getNumNodesPerElement();
  for (iterator it = pdeTM.begin(); it != pdeTM.end(); ++it)
    it->init(nqp, nn);
}

#endif // FEAPP_APPLICATION_HPP
