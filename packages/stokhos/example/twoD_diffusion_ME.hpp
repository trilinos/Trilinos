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

#ifndef TWOD_DIFFUSION_ME_HPP
#define TWOD_DIFFUSION_ME_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

#include "EpetraExt_ModelEvaluator.h"

#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Import.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

#include "Stokhos.hpp"
#include "Stokhos_Epetra.hpp"

//! ModelEvaluator for a linear 2-D diffusion problem
class twoD_diffusion_ME : public EpetraExt::ModelEvaluator {
public:

  //! Constructor
  twoD_diffusion_ME(
    const Teuchos::RCP<Epetra_Comm>& comm, int n, int d,
    double s = 0.1, double mu = 0.2, 
    const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& basis = 
    Teuchos::null,
    bool log_normal = false);

  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  //! Return solution vector map
  Teuchos::RCP<const Epetra_Map> get_x_map() const;
  
  //! Return residual vector map
  Teuchos::RCP<const Epetra_Map> get_f_map() const;
  
  //! Return parameter vector map
  Teuchos::RCP<const Epetra_Map> get_p_map(int l) const;
  
  //! Return parameter vector map
  Teuchos::RCP<const Epetra_Map> get_p_sg_map(int l) const;

  //! Return response function map
  Teuchos::RCP<const Epetra_Map> get_g_map(int j) const;

  //! Return response function map
  Teuchos::RCP<const Epetra_Map> get_g_sg_map(int j) const;
  
  //! Return array of parameter names
  Teuchos::RCP<const Teuchos::Array<std::string> > 
  get_p_names(int l) const;
  
  //! Return array of parameter names
  Teuchos::RCP<const Teuchos::Array<std::string> > 
  get_p_sg_names(int l) const;
  
  //! Return initial solution
  Teuchos::RCP<const Epetra_Vector> get_x_init() const;
  
  //! Return initial parameters
  Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;

  //! Create W = alpha*M + beta*J matrix
  Teuchos::RCP<Epetra_Operator> create_W() const;
  
  //! Create InArgs
  InArgs createInArgs() const;
  
  //! Create OutArgs
  OutArgs createOutArgs() const;
  
  //! Evaluate model on InArgs
  void evalModel(const InArgs& inArgs, const OutArgs& outArgs) const;
  
  //@}

protected:

  //! Fill coefficient matrix given function to evaluate diffusion coefficient
  template <typename FuncT>
  void fillMatrices(const FuncT& func, int sz);

protected:

  double h;
  Teuchos::Array<double> mesh;
  Teuchos::Array<int> bcIndices;

  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis;
  bool log_normal;

  //! Solution vector map
  Teuchos::RCP<Epetra_Map> x_map;

  //! Overlapped solution vector map
  Teuchos::RCP<Epetra_Map> x_overlapped_map;

  //! Importer to overlapped distribution
  Teuchos::RCP<Epetra_Import> importer;

  //! Initial guess
  Teuchos::RCP<Epetra_Vector> x_init;

  //! Overlapped solution vector
  Teuchos::RCP<Epetra_Vector> x_overlapped;

  //! Parameter vector map
  Teuchos::RCP<Epetra_Map> p_map;

  //! Response vector map
  Teuchos::RCP<Epetra_Map> g_map;

  //! Initial parameters
  Teuchos::RCP<Epetra_Vector> p_init;
  
  //! Parameter names
  Teuchos::RCP< Teuchos::Array<std::string> > p_names;

  //! Jacobian graph
  Teuchos::RCP<Epetra_CrsGraph> graph;

  //! KL coefficients of operator  
  Teuchos::Array<Teuchos::RCP<Epetra_CrsMatrix> > A_k;

  //! Deterministic RHS
  Teuchos::RCP<Epetra_Vector> b;

  //! Vectors to store matrix-vector products in SG residual calculation
  mutable Teuchos::Array< Teuchos::RCP<Epetra_Vector> > sg_kx_vec_all;

  //! Matrix to store deterministic operator
  Teuchos::RCP<Epetra_CrsMatrix> A;
  
};

template <typename FuncT>
void
twoD_diffusion_ME::
fillMatrices(const FuncT& func, int sz)
{
  int NumMyElements = x_map->NumMyElements();
  int *MyGlobalElements = x_map->MyGlobalElements();
  double two;
  int n = mesh.size();
  int n2 = n*n;
  double h2 = h*h;
  int Indices[4];
  double Values[4];
  int NumEntries;

  A_k.resize(sz);
  for (int k=0; k<sz; k++) {
    A_k[k] = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *graph));
    for( int i=0 ; i<NumMyElements; ++i ) {
      if (bcIndices[i] == 1 && k == 0) two = 1; //Enforce BC in mean matrix
      if (bcIndices[i] == 0) {
	Indices[0] = MyGlobalElements[i]-n; //Down
	Indices[1] = MyGlobalElements[i]-1; //left
	Indices[2] = MyGlobalElements[i]+1; //right
	Indices[3] = MyGlobalElements[i]+n; //up
	NumEntries = 4;
	Values[0] = 
	  -(1/h2)*func(mesh[(MyGlobalElements[i]%n2)%n], 
		       mesh[(MyGlobalElements[i]%n2)/n]-(h/2), k);
	Values[1] = 
	  -(1/h2)*func(mesh[(MyGlobalElements[i]%n2)%n]-(h/2), 
		       mesh[(MyGlobalElements[i]%n2)/n], k);
	Values[2] = 
	  -(1/h2)*func(mesh[(MyGlobalElements[i]%n2)%n]+(h/2), 
		       mesh[(MyGlobalElements[i]%n2)/n], k);
	Values[3] = 
	  -(1/h2)*func(mesh[(MyGlobalElements[i]%n2)%n], 
		       mesh[(MyGlobalElements[i]%n2)/n]+(h/2), k);
	
	two = 
	  ( func(mesh[(MyGlobalElements[i]%n2)%n]-(h/2), 
		 mesh[(MyGlobalElements[i]%n2)/n], k) + 
	    func(mesh[(MyGlobalElements[i]%n2)%n]+(h/2), 
		 mesh[(MyGlobalElements[i]%n2)/n], k) + 
	    func(mesh[(MyGlobalElements[i]%n2)%n], 
		 mesh[(MyGlobalElements[i]%n2)/n]-(h/2), k) +
	    func(mesh[(MyGlobalElements[i]%n2)%n], 
		 mesh[(MyGlobalElements[i]%n2)/n]+(h/2), k) ) / h2;
      }
      if(bcIndices[i] == 0) 
	A_k[k]->ReplaceGlobalValues(MyGlobalElements[i], NumEntries, Values, 
			      Indices);
      if (bcIndices[i]==0 || k == 0) 
	A_k[k]->ReplaceGlobalValues(MyGlobalElements[i], 1, &two, 
			      MyGlobalElements+i);
    }
    A_k[k]->FillComplete();
  }
}

#endif // TWOD_DIFFUSION_ME_HPP
