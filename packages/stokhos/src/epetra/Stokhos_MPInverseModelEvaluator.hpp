// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_MPINVERSEMODELEVALUATOR_HPP
#define STOKHOS_MPINVERSEMODELEVALUATOR_HPP

#include "EpetraExt_ModelEvaluator.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

namespace Stokhos {

  //! Nonlinear, inverse multi-point ModelEvaluator
  /*!
   * MPInverseModelEvaluator is an implementation of EpetraExt::ModelEvaluator
   * that does the inverse of MPModelEvalutor, namely it takes MP versions of
   * the p InArgs and g and dg/dp OutArgs, and converts them to block vectors
   * that are passed to the underlying model evaluator.  This allows block
   * nonlinear problems to appear as MP problems.
   */
  class MPInverseModelEvaluator : public EpetraExt::ModelEvaluator {
  public:

    // Constructor
    MPInverseModelEvaluator(
      const Teuchos::RCP<EpetraExt::ModelEvaluator>& me,
      const Teuchos::Array<int>& mp_p_index_map,
      const Teuchos::Array<int>& mp_g_index_map,
      const Teuchos::Array< Teuchos::RCP<const Epetra_Map> >& base_g_maps);

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

    //! Mapping between multipoint block parameters and mp parameters
    Teuchos::Array<int> mp_p_index_map;

    //! Mapping between stochastic block responses and sg responses
    Teuchos::Array<int> mp_g_index_map;

    //! Base maps of block g vectors
    Teuchos::Array< Teuchos::RCP<const Epetra_Map> > base_g_maps;

    //! Number of parameters
    int num_p;

    //! Number of responses
    int num_g;

    //! Number of multi-point parameter vectors
    int num_p_mp;

    //! Number of multi-point response vectors
    int num_g_mp;

  };

}

#endif // STOKHOS_MPMODELEVALUATOR_HPP
