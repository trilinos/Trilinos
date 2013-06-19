// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef TWOD_DIFFUSION_PROBLEM_HPP
#define TWOD_DIFFUSION_PROBLEM_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Import.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

#include "Stokhos.hpp"
#include "Stokhos_Epetra.hpp"

//! A linear 2-D diffusion problem
class twoD_diffusion_problem  {
public:

  //! Constructor
  twoD_diffusion_problem(
    const Teuchos::RCP<const Epetra_Comm>& comm, int n, int d,
    double s = 0.1, double mu = 0.2, 
    const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& basis = 
    Teuchos::null,
    bool log_normal = false,
    bool eliminate_bcs = false);

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
  
  //! Compute residual
  void computeResidual(const Epetra_Vector& x,
		       const Epetra_Vector& p,
		       Epetra_Vector& f);

  //! Compute Jacobian
  void computeJacobian(const Epetra_Vector& x,
		       const Epetra_Vector& p,
		       Epetra_Operator& J);

  //! Compute response
  void computeResponse(const Epetra_Vector& x,
		       const Epetra_Vector& p,
		       Epetra_Vector& g);

  //! Compute SG residual
  void computeSGResidual(const Stokhos::EpetraVectorOrthogPoly& x_sg,
			 const Stokhos::EpetraVectorOrthogPoly& p_sg,
			 Stokhos::OrthogPolyExpansion<int,double>& expn,
			 Stokhos::EpetraVectorOrthogPoly& f_sg);

  //! Compute Jacobian
  void computeSGJacobian(const Stokhos::EpetraVectorOrthogPoly& x_sg,
			 const Stokhos::EpetraVectorOrthogPoly& p_sg,
			 Stokhos::EpetraOperatorOrthogPoly& J_sg);

  //! Compute SG response
  void computeSGResponse(const Stokhos::EpetraVectorOrthogPoly& x_sg,
			 const Stokhos::EpetraVectorOrthogPoly& p_sg,
			 Stokhos::EpetraVectorOrthogPoly& g_sg);
  
  //@}

  //! Get mean matrix
  Teuchos::RCP<Epetra_CrsMatrix> get_mean() const { return A_k[0]; }

protected:

  //! Fill coefficient matrix given function to evaluate diffusion coefficient
  template <typename FuncT>
  void fillMatrices(const FuncT& func, int sz);

  //! Compute A matrix
  void compute_A(const Epetra_Vector& p);

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

#endif // TWOD_DIFFUSION_ME_HPP
