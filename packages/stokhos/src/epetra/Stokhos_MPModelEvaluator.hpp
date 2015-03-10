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

#ifndef STOKHOS_MPMODELEVALUATOR_HPP
#define STOKHOS_MPMODELEVALUATOR_HPP

#include "EpetraExt_ModelEvaluator.h"
#include "EpetraExt_MultiComm.h"
#include "EpetraExt_BlockVector.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Stokhos_ParallelData.hpp"

#include "Stokhos_ProductContainer.hpp"
#include "Stokhos_VectorOrthogPolyTraitsEpetra.hpp"
#include "Stokhos_ProductEpetraVector.hpp"
#include "Stokhos_ProductEpetraMultiVector.hpp"
#include "Stokhos_BlockDiagonalOperator.hpp"

namespace Stokhos {

  //! Multi-point model evaluator
  /*!
   * Transforms a multi-point nonlinear problem to single nonlinear problem.
   */
  class MPModelEvaluator : public EpetraExt::ModelEvaluator {
  public:

    // Constructor
    MPModelEvaluator(
      const Teuchos::RCP<EpetraExt::ModelEvaluator>& me,
      const Teuchos::RCP<const EpetraExt::MultiComm>& mp_comm,
      const Teuchos::RCP<const Epetra_Map>& mp_block_map,
      const Teuchos::RCP<Teuchos::ParameterList>& params);

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
    Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner> create_WPrec() const;

    //! Create MP operator representing dg/dxdot
    Teuchos::RCP<Epetra_Operator> create_DgDx_dot_op(int j) const;

    //! Create MP operator representing dg/dx
    Teuchos::RCP<Epetra_Operator> create_DgDx_op(int j) const;

    //! Create MP operator representing dg/dp
    Teuchos::RCP<Epetra_Operator> create_DgDp_op(int j, int i) const;

    //! Create MP operator representing df/dp
    Teuchos::RCP<Epetra_Operator> create_DfDp_op(int i) const;

    //! Create InArgs
    InArgs createInArgs() const;

    //! Create OutArgs
    OutArgs createOutArgs() const;

    //! Evaluate model on InArgs
    void evalModel(const InArgs& inArgs, const OutArgs& outArgs) const;

    //@}

    //! Set initial multi-point solution
    void set_x_mp_init(const Stokhos::ProductEpetraVector& x_mp_in);

    //! Return initial SG x
    Teuchos::RCP<const Stokhos::ProductEpetraVector> get_x_mp_init() const;

    //! Set initial multi-point parameter
    void set_p_mp_init(int i, const Stokhos::ProductEpetraVector& p_mp_in);

    //! Return initial SG parameters
    Teuchos::RCP<const Stokhos::ProductEpetraVector> get_p_mp_init(int l) const;

    //! Get indices of MP parameters
    /*!
     * These indices determine which parameter vectors support MP
     */
    Teuchos::Array<int> get_p_mp_map_indices() const;

    //! Get indices of MP responses
    /*!
     * These indices determine which response vectors support MP
     */
    Teuchos::Array<int> get_g_mp_map_indices() const;

    //! Get base maps of MP responses
    Teuchos::Array< Teuchos::RCP<const Epetra_Map> > get_g_mp_base_maps() const;

    //! Create multi-point vector using x map and owned mp map
    Teuchos::RCP<Stokhos::ProductEpetraVector> 
    create_x_mp(Epetra_DataAccess CV = Copy, 
		const Epetra_Vector* v = NULL) const;

    //! Create multi-point vector using x map
    Teuchos::RCP<Stokhos::ProductEpetraMultiVector> 
    create_x_mv_mp(int num_vecs,
		Epetra_DataAccess CV = Copy, 
		const Epetra_MultiVector* v = NULL) const;

    //! Create multi-point vector using p map
    Teuchos::RCP<Stokhos::ProductEpetraVector> 
    create_p_mp(int l, Epetra_DataAccess CV = Copy, 
		const Epetra_Vector* v = NULL) const;

    //! Create multi-point vector using p map
    Teuchos::RCP<Stokhos::ProductEpetraMultiVector> 
    create_p_mv_mp(int l, int num_vecs, Epetra_DataAccess CV = Copy, 
		   const Epetra_MultiVector* v = NULL) const;

    //! Create multi-point vector using f map
    Teuchos::RCP<Stokhos::ProductEpetraVector> 
    create_f_mp(Epetra_DataAccess CV = Copy, 
		const Epetra_Vector* v = NULL) const;

    //! Create multi-point multi-vector using f map
    Teuchos::RCP<Stokhos::ProductEpetraMultiVector> 
    create_f_mv_mp(int num_vecs, Epetra_DataAccess CV = Copy, 
		   const Epetra_MultiVector* v = NULL) const;

    //! Create multi-point vector using g map
    Teuchos::RCP<Stokhos::ProductEpetraVector> 
    create_g_mp(int l, Epetra_DataAccess CV = Copy, 
		const Epetra_Vector* v = NULL) const;

    //! Create multi-point multi-vector using g map
    Teuchos::RCP<Stokhos::ProductEpetraMultiVector> 
    create_g_mv_mp(int l, int num_vecs, Epetra_DataAccess CV = Copy, 
		const Epetra_MultiVector* v = NULL) const;

  protected:

    //! Underlying model evaluator
    Teuchos::RCP<EpetraExt::ModelEvaluator> me;

    //! Number of blocks
    unsigned int num_mp_blocks;

    //! Parallel MP communicator
    Teuchos::RCP<const EpetraExt::MultiComm> mp_comm;

    //! Map for layout of parallel MP blocks
    Teuchos::RCP<const Epetra_Map> mp_block_map;

    //! Algorithmic parameters
    Teuchos::RCP<Teuchos::ParameterList> params;

    //! Whether we support x (and thus f and W)
    bool supports_x;

    //! Underlying unknown map
    Teuchos::RCP<const Epetra_Map> x_map;

    //! Underlying residual map
    Teuchos::RCP<const Epetra_Map> f_map;

    //! Block MP unknown map
    Teuchos::RCP<const Epetra_Map> mp_x_map;

    //! Block MP residual map
    Teuchos::RCP<const Epetra_Map> mp_f_map;

    //! Number of parameter vectors of underlying model evaluator
    int num_p;

    //! Number of multi-point parameter vectors
    int num_p_mp;

    //! Index map between block-p and p_mp maps
    Teuchos::Array<int> mp_p_index_map;

    //! Block MP parameter map
    Teuchos::Array< Teuchos::RCP<const Epetra_Map> > mp_p_map;

    //! MP coefficient parameter names
    Teuchos::Array< Teuchos::RCP< Teuchos::Array<std::string> > > mp_p_names;

    //! Number of response vectors of underlying model evaluator
    int num_g;

    //! Number of multi-point response vectors
    int num_g_mp;

    //! Index map between block-g and g_mp maps
    Teuchos::Array<int> mp_g_index_map;

    //! Block MP response map
    Teuchos::Array< Teuchos::RCP<const Epetra_Map> > mp_g_map;

    //! W multi-point components
    mutable Teuchos::RCP< Stokhos::ProductEpetraOperator> W_mp_blocks;

    //! MP initial x
    Teuchos::RCP<Stokhos::ProductEpetraVector> mp_x_init;

    //! MP initial p
    Teuchos::Array< Teuchos::RCP<ProductEpetraVector> > mp_p_init;

    //! W pointer for evaluating preconditioner
    mutable Teuchos::RCP<Stokhos::BlockDiagonalOperator> my_W;

    //! x pointer for evaluating preconditioner
    mutable Teuchos::RCP<Epetra_Vector> my_x;

  };

}

#endif // FEAPP_MODELEVALUATOR_HPP
