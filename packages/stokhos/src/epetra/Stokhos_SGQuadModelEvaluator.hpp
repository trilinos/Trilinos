// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_SGQUADMODELEVALUATOR_HPP
#define STOKHOS_SGQUADMODELEVALUATOR_HPP

#include "EpetraExt_ModelEvaluator.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

namespace Stokhos {

  /*!
   * \brief ModelEvaluator adaptor that implements the stochastic Galerkin
   * residual and Jacobian computations using quadrature.
   */
  /*!
   * This class provides a ModelEvaluator implementation to adapt a non-SG
   * capable ModelEvaluator to one that can be used by 
   * Stokhos::SGModelEvaluator.  It does so be implementing the SG residual
   * and Jacobian calculations by sampling a deterministic ModelEvaluator
   * at a set of quadrature points.
   */
  class SGQuadModelEvaluator : public EpetraExt::ModelEvaluator {
  public:

    // Constructor
    SGQuadModelEvaluator(
      const Teuchos::RCP<EpetraExt::ModelEvaluator>& me);
    
    /** \name Overridden from EpetraExt::ModelEvaluator . */
    //@{

    //! Return solution vector map
    Teuchos::RCP<const Epetra_Map> get_x_map() const;

    //! Return residual vector map
    Teuchos::RCP<const Epetra_Map> get_f_map() const;

    //! Return parameter vector map
    Teuchos::RCP<const Epetra_Map> get_p_map(int l) const;

    //! Return observation vector map
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

    //! Number of parameter vectors
    int num_p;

    //! Number of response vectors
    int num_g;

    //! Time derivative vector
    Teuchos::RCP<Epetra_Vector> x_dot_qp;

    //! Solution vector
    Teuchos::RCP<Epetra_Vector> x_qp;

    //! Parameter vectors
    Teuchos::Array< Teuchos::RCP<Epetra_Vector> > p_qp;

    //! Residual vector
    Teuchos::RCP<Epetra_Vector> f_qp;

    //! W operator
    Teuchos::RCP<Epetra_Operator> W_qp;

    //! Residual derivatives
    Teuchos::Array<EpetraExt::ModelEvaluator::Derivative> dfdp_qp;

    //! Response vectors
    Teuchos::Array< Teuchos::RCP<Epetra_Vector> > g_qp;

    //! Response derivative
    Teuchos::Array<EpetraExt::ModelEvaluator::Derivative> dgdx_qp;

    //! Response derivative
    Teuchos::Array<EpetraExt::ModelEvaluator::Derivative> dgdx_dot_qp;

    //! Response sensitivities
    Teuchos::Array< Teuchos::Array<EpetraExt::ModelEvaluator::Derivative> > dgdp_qp;

  };

}

#endif // STOKHOS_SGQUADMODELEVALUATOR_HPP
