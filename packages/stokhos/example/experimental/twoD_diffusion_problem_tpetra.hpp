// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TWOD_DIFFUSION_PROBLEM_TPETRA_HPP
#define TWOD_DIFFUSION_PROBLEM_TPETRA_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Comm.hpp"

#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "Stokhos.hpp"
#include "Stokhos_KL_ExponentialRandomField.hpp"
//#include "Stokhos_Epetra.hpp"

//! A linear 2-D diffusion problem
template <typename Scalar,
	  typename MeshScalar,
	  typename BasisScalar,
	  typename LocalOrdinal,
	  typename GlobalOrdinal,
	  typename Node>
class twoD_diffusion_problem  {
public:

  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Vector;
  typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_Operator;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsMatrix;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> Tpetra_CrsGraph;
  typedef Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Import;

  //! Constructor
  twoD_diffusion_problem(
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm, 
    LocalOrdinal n, LocalOrdinal d,
    BasisScalar s = 0.1, BasisScalar mu = 0.2, 
    bool log_normal = false,
    bool eliminate_bcs = false);

  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  //! Return solution vector map
  Teuchos::RCP<const Tpetra_Map> get_x_map() const;
  
  //! Return residual vector map
  Teuchos::RCP<const Tpetra_Map> get_f_map() const;
  
  //! Return parameter vector map
  Teuchos::RCP<const Tpetra_Map> get_p_map(LocalOrdinal l) const; 
  
  //! Return response function map
  Teuchos::RCP<const Tpetra_Map> get_g_map(LocalOrdinal j) const;
  
  //! Return array of parameter names
  Teuchos::RCP<const Teuchos::Array<std::string> > 
  get_p_names(LocalOrdinal l) const;
    
  //! Return initial solution
  Teuchos::RCP<const Tpetra_Vector> get_x_init() const;
  
  //! Return initial parameters
  Teuchos::RCP<const Tpetra_Vector> get_p_init(LocalOrdinal l) const;

  //! Create W = alpha*M + beta*J matrix
  Teuchos::RCP<Tpetra_CrsMatrix> create_W() const;
  
  //! Compute residual
  void computeResidual(const Tpetra_Vector& x,
		       const Tpetra_Vector& p,
		       Tpetra_Vector& f);

  //! Compute Jacobian
  void computeJacobian(const Tpetra_Vector& x,
		       const Tpetra_Vector& p,
		       Tpetra_CrsMatrix& J);

  //! Compute response
  void computeResponse(const Tpetra_Vector& x,
		       const Tpetra_Vector& p,
		       Tpetra_Vector& g);
  
  //@}

protected:

  //! Fill coefficient matrix given function to evaluate diffusion coefficient
  template <typename FuncT>
  void computeA(const FuncT& func, const Tpetra_Vector& p, 
		Tpetra_CrsMatrix& jac);

protected:

  // Mesh
  MeshScalar h;
  struct MeshPoint {
    MeshPoint() : up(-1), down(-1), left(-1), right(-1), boundary(false) {}
    MeshScalar x, y;
    LocalOrdinal up, down, left, right;
    bool boundary;
  };
  Teuchos::Array<MeshPoint> mesh;
  Teuchos::Array<GlobalOrdinal> bcIndices;
  bool log_normal;
  bool eliminate_bcs;

  //! Solution vector map
  Teuchos::RCP<const Tpetra_Map> x_map;

  //! Importer to overlapped distribution
  Teuchos::RCP<Tpetra_Import> importer;

  //! Initial guess
  Teuchos::RCP<Tpetra_Vector> x_init;

  //! Parameter vector map
  Teuchos::RCP<const Tpetra_Map> p_map;

  //! Response vector map
  Teuchos::RCP<const Tpetra_Map> g_map;

  //! Initial parameters
  Teuchos::RCP<Tpetra_Vector> p_init;
  
  //! Parameter names
  Teuchos::RCP< Teuchos::Array<std::string> > p_names;

  //! Jacobian graph
  Teuchos::RCP<Tpetra_CrsGraph> graph;

  //! RHS
  Teuchos::RCP<Tpetra_Vector> b;

  //! Matrix to store operator
  Teuchos::RCP<Tpetra_CrsMatrix> A;

  // Functor representing a diffusion function given by a KL expansion
  struct KL_Diffusion_Func {
    mutable Teuchos::Array<MeshScalar> point;
    Teuchos::RCP< Stokhos::KL::ExponentialRandomField<MeshScalar> > rf;
    KL_Diffusion_Func(MeshScalar xyLeft, MeshScalar xyRight, 
		      BasisScalar mean, BasisScalar std_dev, 
		      MeshScalar L, LocalOrdinal num_KL);
    Scalar operator() (MeshScalar x, MeshScalar y, 
		       const Teuchos::Array<Scalar>& rv) const;
  };

  // Functor representing a diffusion function given by a log-normal 
  // PC expansion
  template <typename DiffusionFunc>
  struct LogNormal_Diffusion_Func {
    const DiffusionFunc& func;
    LogNormal_Diffusion_Func(const DiffusionFunc& func_) : func(func_) {}
    Scalar operator() (MeshScalar x, MeshScalar y, 
		       const Teuchos::Array<Scalar>& rv) const {
      return std::exp(func(x,y,rv));
    }
  };

  Teuchos::RCP<KL_Diffusion_Func> klFunc;
  Teuchos::RCP<LogNormal_Diffusion_Func<KL_Diffusion_Func> > lnFunc;

};

#include "twoD_diffusion_problem_tpetra_def.hpp"

#endif // TWOD_DIFFUSION_ME_HPP
