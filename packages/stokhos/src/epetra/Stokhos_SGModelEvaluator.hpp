// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_SGMODELEVALUATOR_HPP
#define STOKHOS_SGMODELEVALUATOR_HPP

#include <vector>

#include "EpetraExt_ModelEvaluator.h"
#include "EpetraExt_MultiComm.h"
#include "EpetraExt_BlockVector.h"

#include "Stokhos_SGModelEvaluatorBase.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Stokhos_ParallelData.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_OrthogPolyExpansion.hpp"
#include "Stokhos_EpetraVectorOrthogPoly.hpp"
#include "Stokhos_EpetraMultiVectorOrthogPoly.hpp"
#include "Stokhos_EpetraOperatorOrthogPoly.hpp"
#include "Stokhos_SGOperator.hpp"
#include "Stokhos_EpetraSparse3Tensor.hpp"
#include "Stokhos_SGPreconditionerFactory.hpp"

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
  class SGModelEvaluator : public Stokhos::SGModelEvaluatorBase {
  public:

    // Constructor
    SGModelEvaluator(
      const Teuchos::RCP<EpetraExt::ModelEvaluator>& me,
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis,
      const Teuchos::RCP<const Stokhos::Quadrature<int,double> >& sg_quad,
      const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> >& sg_exp,
      const Teuchos::RCP<const Stokhos::ParallelData>& sg_parallel_data,
      const Teuchos::RCP<Teuchos::ParameterList>& params,
      bool scaleOP = true);

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

    //! Create SG operator representing dg/dxdot
    Teuchos::RCP<Epetra_Operator> create_DgDx_dot_op(int j) const;

    //! Create SG operator representing dg/dx
    Teuchos::RCP<Epetra_Operator> create_DgDx_op(int j) const;

    //! Create SG operator representing dg/dp
    Teuchos::RCP<Epetra_Operator> create_DgDp_op(int j, int i) const;

    //! Create SG operator representing df/dp
    Teuchos::RCP<Epetra_Operator> create_DfDp_op(int i) const;

    //! Create InArgs
    InArgs createInArgs() const;

    //! Create OutArgs
    OutArgs createOutArgs() const;

    //! Evaluate model on InArgs
    void evalModel(const InArgs& inArgs, const OutArgs& outArgs) const;

    //@}

    /** \name Overridden from Stokhos::SGModelEvaluatorBase . */
    //@{

    //! Set initial solution polynomial
    void set_x_sg_init(const Stokhos::EpetraVectorOrthogPoly& x_sg_in);

    //! Return initial SG x
    Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly> get_x_sg_init() const;

    //! Set initial parameter polynomial
    void set_p_sg_init(int i, const Stokhos::EpetraVectorOrthogPoly& p_sg_in);

    //! Return initial SG parameters
    Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly> get_p_sg_init(int l) const;

    //! Get indices of SG parameters
    /*!
     * These indices determine which parameter vectors support SG
     */
    Teuchos::Array<int> get_p_sg_map_indices() const;

    //! Get indices of SG responses
    /*!
     * These indices determine which response vectors support SG
     */
    Teuchos::Array<int> get_g_sg_map_indices() const;

    //! Get base maps of SG responses
    Teuchos::Array< Teuchos::RCP<const Epetra_Map> > get_g_sg_base_maps() const;

    //! Return overlap stochastic map
    Teuchos::RCP<const Epetra_BlockMap> get_overlap_stochastic_map() const;

    //! Return x sg overlap map
    Teuchos::RCP<const Epetra_BlockMap> get_x_sg_overlap_map() const;

    //! Return x sg importer
    Teuchos::RCP<const Epetra_Import> get_x_sg_importer() const;

    //! Create vector orthog poly using x map and owned sg map
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly>
    create_x_sg(Epetra_DataAccess CV = Copy,
                const Epetra_Vector* v = NULL) const;

    //! Create vector orthog poly using x map and overlap sg map
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly>
    create_x_sg_overlap(Epetra_DataAccess CV = Copy,
                        const Epetra_Vector* v = NULL) const;

    //! Create vector orthog poly using x map and owned sg map
    Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly>
    create_x_mv_sg(int num_vecs,
                   Epetra_DataAccess CV = Copy,
                   const Epetra_MultiVector* v = NULL) const;

    //! Create vector orthog poly using x map and overlap sg map
    Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly>
    create_x_mv_sg_overlap(int num_vecs,
                           Epetra_DataAccess CV = Copy,
                           const Epetra_MultiVector* v = NULL) const;

    //! Create vector orthog poly using p map
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly>
    create_p_sg(int l, Epetra_DataAccess CV = Copy,
                const Epetra_Vector* v = NULL) const;

    //! Create vector orthog poly using f map and owned sg map
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly>
    create_f_sg(Epetra_DataAccess CV = Copy,
                const Epetra_Vector* v = NULL) const;

    //! Create vector orthog poly using f map and overlap sg map
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly>
    create_f_sg_overlap(Epetra_DataAccess CV = Copy,
                        const Epetra_Vector* v = NULL) const;

    //! Create multi-vector orthog poly using f map and owned sg map
    Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly>
    create_f_mv_sg(int num_vecs, Epetra_DataAccess CV = Copy,
                   const Epetra_MultiVector* v = NULL) const;

    //! Create multi-vector orthog poly using f map and overlap sg map
    Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly>
    create_f_mv_sg_overlap(int num_vecs, Epetra_DataAccess CV = Copy,
                           const Epetra_MultiVector* v = NULL) const;

    //! Create vector orthog poly using g map
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly>
    create_g_sg(int l, Epetra_DataAccess CV = Copy,
                const Epetra_Vector* v = NULL) const;

    //! Create multi-vector orthog poly using g map
    Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly>
    create_g_mv_sg(int l, int num_vecs, Epetra_DataAccess CV = Copy,
                   const Epetra_MultiVector* v = NULL) const;

    //@}

    //! Import parallel solution vector
    Teuchos::RCP<EpetraExt::BlockVector>
    import_solution(const Epetra_Vector& x) const;

    //! Import parallel solution vector
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly>
    import_solution_poly(const Epetra_Vector& x) const;

    //! Export parallel solution vector
    Teuchos::RCP<EpetraExt::BlockVector>
    export_solution(const Epetra_Vector& x_overlapped) const;

    //! Import parallel residual vector
    Teuchos::RCP<EpetraExt::BlockVector>
    import_residual(const Epetra_Vector& f) const;

    //! Export parallel residual vector
    Teuchos::RCP<EpetraExt::BlockVector>
    export_residual(const Epetra_Vector& f_overlapped) const;

  protected:

    //! Underlying model evaluator
    Teuchos::RCP<EpetraExt::ModelEvaluator> me;

    //! Stochastic Galerkin basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int, double> > sg_basis;

    //! Stochastic Galerkin quadrature
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > sg_quad;

    //! Stochastic Galerkin expansion
    Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > sg_exp;

    //! Algorithmic parameters
    Teuchos::RCP<Teuchos::ParameterList> params;

    //! Number of stochastic blocks
    unsigned int num_sg_blocks;

    //! Number of W stochastic blocks (may be smaller than num_sg_blocks)
    unsigned int num_W_blocks;

    //! Number of p stochastic blocks (may be smaller than num_sg_blocks)
    unsigned int num_p_blocks;

    //! Whether we support x (and thus f and W)
    bool supports_x;

    //! Underlying unknown map
    Teuchos::RCP<const Epetra_Map> x_map;

    //! Underlying residual map
    Teuchos::RCP<const Epetra_Map> f_map;

    //! Parallel SG data
    Teuchos::RCP<const Stokhos::ParallelData> sg_parallel_data;

    //! Parallel SG communicator
    Teuchos::RCP<const EpetraExt::MultiComm> sg_comm;

    //! Epetra Cijk
    Teuchos::RCP<const Stokhos::EpetraSparse3Tensor> epetraCijk;

    //! Serial Epetra Cijk for dgdx*
    Teuchos::RCP<const Stokhos::EpetraSparse3Tensor> serialCijk;

    //! Map for stochastic blocks
    Teuchos::RCP<const Epetra_BlockMap> stoch_row_map;

    //! Overlapped map for stochastic blocks (local map)
    Teuchos::RCP<const Epetra_BlockMap> overlapped_stoch_row_map;

    //! Overlapped map for p stochastic blocks (local map)
    Teuchos::RCP<const Epetra_BlockMap> overlapped_stoch_p_map;

    //! Block SG unknown map
    Teuchos::RCP<const Epetra_Map> sg_x_map;

    //! Block SG overlapped unknown map
    Teuchos::RCP<const Epetra_Map> sg_overlapped_x_map;

    //! Block SG residual map
    Teuchos::RCP<const Epetra_Map> sg_f_map;

    //! Block SG overlapped residual map
    Teuchos::RCP<const Epetra_Map> sg_overlapped_f_map;

    //! Importer from SG to SG-overlapped maps
    Teuchos::RCP<Epetra_Import> sg_overlapped_x_importer;

    //! Exporter from SG-overlapped to SG maps
    Teuchos::RCP<Epetra_Export> sg_overlapped_f_exporter;

    //! Number of parameter vectors of underlying model evaluator
    int num_p;

    //! Number of stochastic parameter vectors
    int num_p_sg;

    //! Index map between block-p and p_sg maps
    Teuchos::Array<int> sg_p_index_map;

    //! Block SG parameter map
    Teuchos::Array< Teuchos::RCP<const Epetra_Map> > sg_p_map;

    //! SG coefficient parameter names
    Teuchos::Array< Teuchos::RCP< Teuchos::Array<std::string> > > sg_p_names;

    //! Number of response vectors of underlying model evaluator
    int num_g;

    //! Number of stochastic response vectors
    int num_g_sg;

    //! Index map between block-g and g_sg maps
    Teuchos::Array<int> sg_g_index_map;

    //! Block SG response map
    Teuchos::Array< Teuchos::RCP<const Epetra_Map> > sg_g_map;

    //! x_dot stochastic Galerkin components
    Teuchos::RCP< Stokhos::EpetraVectorOrthogPoly > x_dot_sg_blocks;

    //! x stochastic Galerkin components
    Teuchos::RCP< Stokhos::EpetraVectorOrthogPoly > x_sg_blocks;

    //! f stochastic Galerkin components
    mutable Teuchos::RCP< Stokhos::EpetraVectorOrthogPoly > f_sg_blocks;

    //! W stochastic Galerkin components
    mutable Teuchos::RCP< Stokhos::EpetraOperatorOrthogPoly > W_sg_blocks;

    //! SG initial x
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_x_init;

    //! SG initial p
    Teuchos::Array< Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> > sg_p_init;

    //! Whether to always evaluate W with f
    bool eval_W_with_f;

    //! W pointer for evaluating W with f
    mutable Teuchos::RCP<Stokhos::SGOperator> my_W;

    //! x pointer for evaluating preconditioner
    mutable Teuchos::RCP<Epetra_Vector> my_x;

    bool scaleOP;

    //! Preconditioner factory
    Teuchos::RCP< Stokhos::SGPreconditionerFactory > sg_prec_factory;

  };

}

#endif // FEAPP_MODELEVALUATOR_HPP
