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

#ifndef STOKHOS_SGQUADMPMODELEVALUATOR_HPP
#define STOKHOS_SGQUADMPMODELEVALUATOR_HPP

#include "EpetraExt_ModelEvaluator.h"
#include "EpetraExt_MultiComm.h"

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
   * and Jacobian calculations by sampling a multi-point ModelEvaluator
   * at a set of quadrature points.
   */
  class SGQuadMPModelEvaluator : public EpetraExt::ModelEvaluator {
  public:

    // Constructor
    SGQuadMPModelEvaluator(
      const Teuchos::RCP<EpetraExt::ModelEvaluator>& me,
      const Teuchos::RCP<const EpetraExt::MultiComm>& mp_comm,
      const Teuchos::RCP<const Epetra_Map>& mp_block_map);
    
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

    //! Parallel MP communicator
    Teuchos::RCP<const EpetraExt::MultiComm> mp_comm;

    //! Map for layout of parallel MP blocks
    Teuchos::RCP<const Epetra_Map> mp_block_map;

    //! SG quadrature
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > sg_quad;

    //! Number of parameter vectors
    int num_p;

    //! Number of response vectors
    int num_g;

    //! Number of multipoint parameter vectors
    int num_p_mp;

    //! Number of multipoint response vectors
    int num_g_mp;

    //! Index map between block-p and p_mp maps
    Teuchos::Array<int> mp_p_index_map;

    //! Index map between block-g and g_mp maps
    Teuchos::Array<int> mp_g_index_map;

    //! Time derivative vector
    mp_vector_t x_dot_mp;

    //! Solution vector
    mp_vector_t x_mp;

    //! Parameter vectors
    Teuchos::Array<mp_vector_t> p_mp;

    //! Residual vector
    mp_vector_t f_mp;

    //! W operator
    mp_operator_t W_mp;

    //! Residual derivatives
    Teuchos::Array<MPDerivative> dfdp_mp;

    //! Response vectors
    Teuchos::Array<mp_vector_t> g_mp;

    //! Response derivative
    Teuchos::Array<MPDerivative> dgdx_mp;

    //! Response derivative
    Teuchos::Array<MPDerivative> dgdx_dot_mp;

    //! Response sensitivities
    Teuchos::Array< Teuchos::Array<MPDerivative> > dgdp_mp;

  };

}

#endif // STOKHOS_SGQUADMPMODELEVALUATOR_HPP
