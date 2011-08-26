// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// 
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#ifndef PIRO_EPETRA_STOKHOS_SOLVER_FACTORY_H
#define PIRO_EPETRA_STOKHOS_SOLVER_FACTORY_H

#include "EpetraExt_ModelEvaluator.h"
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

  class StokhosSolverFactory {
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

    //! Destructor
    ~StokhosSolverFactory();

    /** \name Factory methods */
    //@{

    //! Create stochastic model evaluator
    Teuchos::RCP<Stokhos::SGModelEvaluator> createSGModel(
      const Teuchos::RCP<EpetraExt::ModelEvaluator>& model,
      const Teuchos::RCP<NOX::Epetra::Observer>& noxObserver = Teuchos::null);

    //! Create stochastic solver
    Teuchos::RCP<EpetraExt::ModelEvaluator> createSGSolver(
      const Teuchos::RCP<EpetraExt::ModelEvaluator>& sg_model);

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
