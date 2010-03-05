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

#ifndef STOKHOS_SGINVERSEMODELEVALUATOR_HPP
#define STOKHOS_SGINVERSEMODELEVALUATOR_HPP

#include "EpetraExt_ModelEvaluator.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

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
  class SGInverseModelEvaluator : public EpetraExt::ModelEvaluator {
  public:

    // Constructor
    SGInverseModelEvaluator(
      const Teuchos::RCP<EpetraExt::ModelEvaluator>& me,
      const Teuchos::Array<int>& sg_p_index,
      const Teuchos::Array<int>& sg_g_index);

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

    //! Return initial parameters
    Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;

    //! Create InArgs
    InArgs createInArgs() const;

    //! Create OutArgs
    OutArgs createOutArgs() const;

    //! Evaluate model on InArgs
    void evalModel(const InArgs& inArgs, const OutArgs& outArgs) const;

    //@}

  protected:

    //! Underlying model evaluator
    Teuchos::RCP<EpetraExt::ModelEvaluator> me;

    //! Index of stochastic parameters
    Teuchos::Array<int> sg_p_index;

    //! Index of stochastic responses
    Teuchos::Array<int> sg_g_index;

    //! Number of stochastic blocks
    int num_sg_blocks;

    //! Number of stochastic parameter vectors
    int num_p_sg;

    //! Number of stochastic response vectors
    int num_g_sg;

    //! Block SG p vectors
    Teuchos::Array< Teuchos::RCP< Epetra_Vector > > block_p;

    //! Block SG g vectors
    mutable Teuchos::Array< Teuchos::RCP< Epetra_Vector > > block_g;

    //! Block SG dg/dp vectors
    mutable Teuchos::Array< Teuchos::Array< Teuchos::RCP< Epetra_MultiVector > > > block_dgdp;

  };

}

#endif // STOKHOS_SGMODELEVALUATOR_HPP
