// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#ifndef PIRO_EPETRA_STOKHOS_SOLVER_FACTORY_H
#define PIRO_EPETRA_STOKHOS_SOLVER_FACTORY_H

#include "EpetraExt_ModelEvaluator.h"
#include "Teuchos_VerboseObject.hpp"
#include "Piro_Epetra_StokhosNOXObserver.hpp"

#include "Stokhos_SGModelEvaluator.hpp"
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

    //! Constructor
    StokhosSolverFactory(const Teuchos::RCP<Teuchos::ParameterList>& piroParams,
			 const Teuchos::RCP<const Epetra_Comm>& globalComm);

    //! Reset Stokhos solver parameters
    void resetSolverParameters(const Teuchos::ParameterList& new_solver_params);

    /** \name Factory methods */
    //@{

    //! Create stochastic model evaluator
    Teuchos::RCP<Stokhos::SGModelEvaluator> createSGModel(
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
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis;
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad;
    Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion;
    Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk;
    Teuchos::RCP<const EpetraExt::MultiComm> sg_comm;
    Teuchos::RCP<Stokhos::ParallelData> sg_parallel_data;

    Teuchos::RCP<EpetraExt::ModelEvaluator> model;
    Teuchos::RCP<Stokhos::SGModelEvaluator> sg_nonlin_model;
    
  };
  
}
}
#endif
