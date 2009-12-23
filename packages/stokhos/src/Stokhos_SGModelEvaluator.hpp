// $Id$ 
// $Source$ 
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

#ifndef STOKHOS_SGMODELEVALUATOR_HPP
#define STOKHOS_SGMODELEVALUATOR_HPP

#include <vector>

#include "EpetraExt_ModelEvaluator.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "EpetraExt_BlockVector.h"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_VectorOrthogPoly.hpp"
#include "Stokhos_VectorOrthogPolyTraitsEpetra.hpp"
#include "Stokhos_PreconditionerFactory.hpp"

namespace Stokhos {

  //! Nonlinear, stochastic Galerkin ModelEvaluator
  /*!
   * SGModelEvaluator is an implementation of EpetraExt::ModelEvaluator that
   * generates a nonlinear problem from a stochastic Galerkin
   * expansion.  It wraps a supplied ModelEvaluator that supports the SG
   * versions of p, x, and possibly x_dot InArgs, and f and W OutArgs, and 
   * translates those into a new nonlinear problem.  It does so by 
   * concatenating all of the SG components of p, x, x_dot, and f into extended
   * block vectors that form the parameters, solution vector, time derivative
   * vector and residual for the new nonlinear problem.  For dealing with the
   * W matrix two methods are supported:  forming a fully-assembled SG matrix
   * and a "matrix free" method.  The choice is selected by setting the
   * "Jacobian Method" parameter of the parameter list supplied to the
   * constructor, which can be either "Fully Assembled" or "Matrix Free".  In
   * the first case, the W operator of the underlying model evaluator must be
   * an Epetra_CrsMatrix.  In the second case, a preconditioner for the mean
   * block must also be supplied via the "Preconditioner Factory" parameter
   * of this list.  This preconditioner factory must implement the
   * Stokhos::PreconditionerFactory interface also supplied in this file.
   * Currently using a preconditioner for the mean is the only option
   * available for preconditioning the SG system when using the matrix-free
   * method.
   */
  class SGModelEvaluator : public EpetraExt::ModelEvaluator {
  public:

    // Constructor
    SGModelEvaluator(
      const Teuchos::RCP<EpetraExt::ModelEvaluator>& me,
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis,
      const Teuchos::Array<int>& sg_p_index,
      const Teuchos::Array<int>& sg_g_index,
      const Teuchos::RCP<Teuchos::ParameterList>& params,
      const Teuchos::RCP<const Epetra_Comm>& comm,
      const Teuchos::Array< Stokhos::VectorOrthogPoly<Epetra_Vector> >& initial_p_sg,
      const Stokhos::VectorOrthogPoly<Epetra_Vector>* initial_x_sg = NULL);

    /** \name Overridden from EpetraExt::ModelEvaluator . */
    //@{

    //! Return solution vector map
    Teuchos::RCP<const Epetra_Map> get_x_map() const;

    //! Return residual vector map
    Teuchos::RCP<const Epetra_Map> get_f_map() const;

    //! Return parameter vector map
    Teuchos::RCP<const Epetra_Map> get_p_map(int l) const;

    //! Return response map
    Teuchos::RCP<const Epetra_Map> get_g_map(int l) const;

    //! Return array of parameter names
    Teuchos::RCP<const Teuchos::Array<std::string> > 
    get_p_names(int l) const;

    //! Return initial solution
    Teuchos::RCP<const Epetra_Vector> get_x_init() const;

    //! Return initial parameters
    Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;

    //! Create W = alpha*M + beta*J matrix
    Teuchos::RCP<Epetra_Operator> create_W() const;

    //! Create preconditioner operator
    Teuchos::RCP<Epetra_Operator> create_M() const;

    //! Create InArgs
    InArgs createInArgs() const;

    //! Create OutArgs
    OutArgs createOutArgs() const;

    //! Evaluate model on InArgs
    void evalModel(const InArgs& inArgs, const OutArgs& outArgs) const;

    //@}

    //! Set initial solution vector
    /*! 
     * This is NOT a virtual ModelEvaluator method, and is just a convenience
     * for users of this class.
     */
    void set_x_init(const Epetra_Vector& x_in);

  protected:

    //! Underlying model evaluator
    Teuchos::RCP<EpetraExt::ModelEvaluator> me;

    //! Stochastic Galerking basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> > sg_basis;

    //! Index of stochastic parameters
    Teuchos::Array<int> sg_p_index;

    //! Index of stochastic responses
    Teuchos::Array<int> sg_g_index;

    //! Algorithmic parameters
    Teuchos::RCP<Teuchos::ParameterList> params;

    //! Number of stochastic blocks
    unsigned int num_sg_blocks;

    //! Whether we support x (and thus f and W)
    bool supports_x;

    //! Underlying unknown map
    Teuchos::RCP<const Epetra_Map> x_map;

    //! Underlying residual map
    Teuchos::RCP<const Epetra_Map> f_map;

    //! Parallel SG communicator
    Teuchos::RCP<const Epetra_Comm> sg_comm;

    //! Block SG unknown map
    Teuchos::RCP<const Epetra_Map> sg_x_map;

    //! Block SG residual map
    Teuchos::RCP<const Epetra_Map> sg_f_map;

    //! Number of stochastic parameter vectors
    int num_p;

    //! Block SG parameter map
    Teuchos::Array< Teuchos::RCP<const Epetra_Map> > sg_p_map;

    //! SG coefficient parameter names
    Teuchos::Array< Teuchos::RCP< Teuchos::Array<std::string> > > sg_p_names;

    //! Number of stochastic response vectors
    int num_g;

    //! Block SG response map
    Teuchos::Array< Teuchos::RCP<const Epetra_Map> > sg_g_map;

    //! Triple product tensor
    Teuchos::RCP< const Stokhos::Sparse3Tensor<int, double> > Cijk;

    //! x_dot stochastic Galerkin components
    Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_Vector> > x_dot_sg_blocks;

    //! x stochastic Galerkin components
    Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_Vector> > x_sg_blocks;

    //! p stochastic Galerkin components
    Teuchos::Array< Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_Vector> > > p_sg_blocks;

    //! f stochastic Galerkin components
    mutable Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_Vector> > f_sg_blocks;

    //! W stochastic Galerkin components
    mutable Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_Operator> > W_sg_blocks;

    //! g stochastic Galerkin components
    mutable Teuchos::Array< Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_Vector> > > g_sg_blocks;

    //! dg/dp stochastic Galerkin components
    mutable Teuchos::Array< Teuchos::Array< Teuchos::RCP< Stokhos::VectorOrthogPoly<Derivative> > > > dgdp_sg_blocks;

    //! Method for creating block Jacobian
    enum EJacobianMethod {
      MATRIX_FREE,
      KL_REDUCED_MATRIX_FREE,
      FULLY_ASSEMBLED
    };

    //! Method for creating block Jacobian
    EJacobianMethod jacobianMethod;

    std::vector< std::vector<int> > rowStencil;
    std::vector<int> rowIndex;

    //! SG initial x
    Teuchos::RCP<EpetraExt::BlockVector> sg_x_init;

    //! SG initial p
    Teuchos::Array< Teuchos::RCP<EpetraExt::BlockVector> > sg_p_init;

    //! SG Preconditioner factory
    Teuchos::RCP<Stokhos::PreconditionerFactory> precFactory;

    //! Whether to always evaluate W with f
    bool eval_W_with_f;

    //! W pointer for evaluating W with f
    mutable Teuchos::RCP<Epetra_Operator> my_W;

  };

}

#endif // FEAPP_MODELEVALUATOR_HPP
