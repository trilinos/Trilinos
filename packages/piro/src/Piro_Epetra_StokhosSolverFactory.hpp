// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_EPETRA_STOKHOS_SOLVER_FACTORY_H
#define PIRO_EPETRA_STOKHOS_SOLVER_FACTORY_H

#include "EpetraExt_ModelEvaluator.h"
#include "Teuchos_VerboseObject.hpp"
#include "Piro_Epetra_StokhosNOXObserver.hpp"

#include "Stokhos_SGModelEvaluatorBase.hpp"
#include "Stokhos_SGInverseModelEvaluator.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_Quadrature.hpp"
#include "Stokhos_OrthogPolyExpansion.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#include "Stokhos_ParallelData.hpp"
#include "EpetraExt_MultiComm.h"

namespace Piro {
namespace Epetra {

  class StokhosSolverFactory :
    public Teuchos::VerboseObject<StokhosSolverFactory> {
  public:

    //! SG method
    enum SG_METHOD {
      SG_AD,
      SG_GLOBAL,
      SG_NI,
      SG_MPNI
    };

    //! SG ModelEvaluator method
    enum SG_ME_METHOD {
      SG_ME_DEFAULT,
      SG_ME_INTERLACED,
      SG_ME_ADAPTIVE
    };

    //! Constructor
    StokhosSolverFactory(const Teuchos::RCP<Teuchos::ParameterList>& piroParams,
                         const Teuchos::RCP<const Epetra_Comm>& globalComm);

    //! Reset Stokhos solver parameters
    void resetSolverParameters(const Teuchos::ParameterList& new_solver_params);

    /** \name Factory methods */
    //@{

    //! Create stochastic model evaluator
    Teuchos::RCP<Stokhos::SGModelEvaluatorBase> createSGModel(
      const Teuchos::RCP<EpetraExt::ModelEvaluator>& model);

    //! Create stochastic observer
    Teuchos::RCP<NOX::Epetra::Observer> createSGObserver(
        const Teuchos::RCP<NOX::Epetra::Observer>& noxObserver);

    //! Create stochastic solver
    Teuchos::RCP<EpetraExt::ModelEvaluator> createSGSolver(
      const Teuchos::RCP<EpetraExt::ModelEvaluator>& sg_model,
      const Teuchos::RCP<NOX::Epetra::Observer>& sg_observer = Teuchos::null);

    //! Create stochastic solver adapter
    Teuchos::RCP<Stokhos::SGInverseModelEvaluator> createSGSolverAdapter(
      const Teuchos::RCP<EpetraExt::ModelEvaluator>& sg_solver);

    //! Create response statistic model evaluator
    Teuchos::RCP<EpetraExt::ModelEvaluator> createRSModel(
      const Teuchos::RCP<EpetraExt::ModelEvaluator>& sg_model);

    //@}

    /** \name Accessors */
    //@{

    //! Get spatial comm
    Teuchos::RCP<const Epetra_Comm> getSpatialComm() const;

    //! Get stochastic comm
    Teuchos::RCP<const Epetra_Comm> getStochasticComm() const;

    //! Get global multi-comm
    Teuchos::RCP<const EpetraExt::MultiComm> getGlobalMultiComm() const;

    //! Get stochastic basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >
    getBasis() const { return basis; }

    //! Get quadrature rule
    Teuchos::RCP<const Stokhos::Quadrature<int,double> >
    getQuad() const { return quad; }

    //! Get SG method
    SG_METHOD getSGMethod() const { return sg_method; }

    //! Get SG ME method
    SG_ME_METHOD getSGMEMethod() const { return sg_me_method; }

    Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> >
    getExpansion() const { return expansion; }

    Teuchos::RCP<Stokhos::ParallelData> getParallelData() const
    { return sg_parallel_data; }

    //@}

  private:

    //! Get valid parameters
    Teuchos::RCP<const Teuchos::ParameterList>
     getValidSGParameters() const;

  private:

    enum SG_SOLVER {
      SG_KRYLOV,
      SG_GS,
      SG_JACOBI
    };

    Teuchos::RCP<Teuchos::ParameterList> piroParams;
    Teuchos::RCP<Teuchos::ParameterList> sgSolverParams;

    SG_METHOD sg_method;
    SG_ME_METHOD sg_me_method;
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis;
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad;
    Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion;
    Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk;
    Teuchos::RCP<const EpetraExt::MultiComm> sg_comm;
    Teuchos::RCP<Stokhos::ParallelData> sg_parallel_data;

    Teuchos::RCP<EpetraExt::ModelEvaluator> model;
    Teuchos::RCP<Stokhos::SGModelEvaluatorBase> sg_nonlin_model;

  };

}
}
#endif
