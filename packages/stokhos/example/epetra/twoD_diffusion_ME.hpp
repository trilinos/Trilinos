// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TWOD_DIFFUSION_ME_HPP
#define TWOD_DIFFUSION_ME_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"

#include "EpetraExt_ModelEvaluator.h"

#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Import.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

#include "Stokhos.hpp"
#include "Stokhos_Epetra.hpp"
#include "Stokhos_AbstractPreconditionerFactory.hpp"

//! ModelEvaluator for a linear 2-D diffusion problem
class twoD_diffusion_ME : public EpetraExt::ModelEvaluator {
public:

  //! Constructor
  twoD_diffusion_ME(
    const Teuchos::RCP<const Epetra_Comm>& comm, int n, int d,
    double s = 0.1, double mu = 0.2, 
    const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& basis = 
    Teuchos::null,
    bool log_normal = false,
    bool eliminate_bcs = false,
    const Teuchos::RCP<Teuchos::ParameterList>& precParams = Teuchos::null);

  //! Destructor
  ~twoD_diffusion_ME();

  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  //! Return solution vector map
  Teuchos::RCP<const Epetra_Map> get_x_map() const;
  
  //! Return residual vector map
  Teuchos::RCP<const Epetra_Map> get_f_map() const;
  
  //! Return parameter vector map
  Teuchos::RCP<const Epetra_Map> get_p_map(int l) const; 
  
  //! Return response function map
  Teuchos::RCP<const Epetra_Map> get_g_map(int j) const;
  
  //! Return array of parameter names
  Teuchos::RCP<const Teuchos::Array<std::string> > 
  get_p_names(int l) const;
    
  //! Return initial solution
  Teuchos::RCP<const Epetra_Vector> get_x_init() const;
  
  //! Return initial parameters
  Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;

  //! Create W = alpha*M + beta*J matrix
  Teuchos::RCP<Epetra_Operator> create_W() const;

  //! Create preconditioner for W
  Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner> create_WPrec() const;
  
  //! Create InArgs
  InArgs createInArgs() const;
  
  //! Create OutArgs
  OutArgs createOutArgs() const;
  
  //! Evaluate model on InArgs
  void evalModel(const InArgs& inArgs, const OutArgs& outArgs) const;
  
  //@}

  //! Get mean matrix
  Teuchos::RCP<Epetra_CrsMatrix> get_mean() const { return A_k[0]; }

protected:

  //! Fill coefficient matrix given function to evaluate diffusion coefficient
  template <typename FuncT>
  void fillMatrices(const FuncT& func, int sz);

protected:

  double h;
  struct MeshPoint {
    MeshPoint() : up(-1), down(-1), left(-1), right(-1), boundary(false) {}
    double x, y;
    int up, down, left, right;
    bool boundary;
  };
  Teuchos::Array<MeshPoint> mesh;
  Teuchos::Array<int> bcIndices;

  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis;
  bool log_normal;
  bool eliminate_bcs;

  Teuchos::RCP<Teuchos::ParameterList> precParams;
  Teuchos::RCP<Stokhos::AbstractPreconditionerFactory> precFactory;

  //! Solution vector map
  Teuchos::RCP<Epetra_Map> x_map;

  //! Importer to overlapped distribution
  Teuchos::RCP<Epetra_Import> importer;

  //! Initial guess
  Teuchos::RCP<Epetra_Vector> x_init;

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

  //! Array to store a point for basis evaluation
  mutable Teuchos::Array<double> point;

  //! Array to store values of basis at a point
  mutable Teuchos::Array<double> basis_vals;
  
};

template <typename FuncT>
void
twoD_diffusion_ME::
fillMatrices(const FuncT& func, int sz)
{
  int NumMyElements = x_map->NumMyElements();
  int *MyGlobalElements = x_map->MyGlobalElements();
  double h2 = h*h;
  double val;

  A_k.resize(sz);
  for (int k=0; k<sz; k++) {
    A_k[k] = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *graph));
    for( int i=0 ; i<NumMyElements; ++i ) {

      // Center
      int global_idx = MyGlobalElements[i];
      if (mesh[global_idx].boundary) {
	if (k == 0)
	  val = 1.0; //Enforce BC in mean matrix
	else
	  val = 0.0;
	A_k[k]->ReplaceGlobalValues(global_idx, 1, &val, &global_idx);
      }
      else {
	double a_down = 
	  -func(mesh[global_idx].x, mesh[global_idx].y-h/2.0, k)/h2;
	double a_left = 
	  -func(mesh[global_idx].x-h/2.0, mesh[global_idx].y, k)/h2;
	double a_right = 
	  -func(mesh[global_idx].x+h/2.0, mesh[global_idx].y, k)/h2;
	double a_up = 
	  -func(mesh[global_idx].x, mesh[global_idx].y+h/2.0, k)/h2;

	// Center
	val = -(a_down + a_left + a_right + a_up);
	A_k[k]->ReplaceGlobalValues(global_idx, 1, &val, &global_idx);

	// Down
	if (!(eliminate_bcs && mesh[mesh[global_idx].down].boundary))
	  A_k[k]->ReplaceGlobalValues(global_idx, 1, &a_down, 
				      &mesh[global_idx].down);

	// Left
	if (!(eliminate_bcs && mesh[mesh[global_idx].left].boundary))
	  A_k[k]->ReplaceGlobalValues(global_idx, 1, &a_left, 
				      &mesh[global_idx].left);

	// Right
	if (!(eliminate_bcs && mesh[mesh[global_idx].right].boundary))
	  A_k[k]->ReplaceGlobalValues(global_idx, 1, &a_right, 
				      &mesh[global_idx].right);

	// Up
	if (!(eliminate_bcs && mesh[mesh[global_idx].up].boundary))
	  A_k[k]->ReplaceGlobalValues(global_idx, 1, &a_up, 
				      &mesh[global_idx].up);
      }
    }
    A_k[k]->FillComplete();
    A_k[k]->OptimizeStorage();
  }
}

#endif // TWOD_DIFFUSION_ME_HPP
