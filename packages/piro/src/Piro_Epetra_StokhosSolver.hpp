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

#ifndef PIRO_EPETRA_STOKHOSSOLVER_H
#define PIRO_EPETRA_STOKHOSSOLVER_H

#include "EpetraExt_ModelEvaluator.h"
#include "Piro_Epetra_StokhosNOXObserver.hpp"

#include "Stokhos_SGModelEvaluator.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_Quadrature.hpp"
#include "Stokhos_OrthogPolyExpansion.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#include "Stokhos_ParallelData.hpp"
#include "EpetraExt_MultiComm.h"

namespace Piro {
namespace Epetra {

  class StokhosSolver : public EpetraExt::ModelEvaluator {
  public:

    /** \name Constructors/initializers */
    //@{

    //! Constructor
    StokhosSolver(const Teuchos::RCP<Teuchos::ParameterList>& piroParams,
		  const Teuchos::RCP<const Epetra_Comm>& globalComm);

    //! Get spatial comm
    Teuchos::RCP<const Epetra_Comm> getSpatialComm() const;
    
    //! Get stochastic comm
    Teuchos::RCP<const Epetra_Comm> getStochasticComm() const;

    //! Get global multi-comm
    Teuchos::RCP<const EpetraExt::MultiComm> getGlobalMultiComm() const;

    //! Setup rest of model evaluator
    void setup(const Teuchos::RCP<EpetraExt::ModelEvaluator>& model,
	       Teuchos::RCP<NOX::Epetra::Observer> noxObserver = Teuchos::null);


    //@}

    ~StokhosSolver();


    /** \name Overridden from EpetraExt::ModelEvaluator . */
    //@{
    
    Teuchos::RCP<const Epetra_Map> get_g_map(int j) const;
    Teuchos::RCP<const Epetra_Map> get_g_sg_map(int j) const;
    /** \brief . */
    Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;
    Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly> get_p_sg_init(int l) const;

    /** \brief . */
    //  Teuchos::RCP<Epetra_Operator> create_W() const;
    /** \brief . */
    EpetraExt::ModelEvaluator::InArgs createInArgs() const;
    /** \brief . */
    EpetraExt::ModelEvaluator::OutArgs createOutArgs() const;
    /** \brief . */
    void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >
    getBasis() const { return basis; }
    Teuchos::RCP<const Stokhos::Quadrature<int,double> >
    getQuad() const { return quad; }
    Teuchos::RCP<Stokhos::SGModelEvaluator>
    get_sg_model() const { return sg_nonlin_model; }

    //! Set initial solution polynomial
    void set_x_sg_init(const Stokhos::EpetraVectorOrthogPoly& x_sg_in) {
      sg_nonlin_model->set_x_sg_init(x_sg_in);
    }

    //! Return initial SG x
    Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly> get_x_sg_init() const {
      return sg_nonlin_model->get_x_sg_init();
    }

    //! Set initial parameter polynomial
    void set_p_sg_init(int i, const Stokhos::EpetraVectorOrthogPoly& p_sg_in) {
      sg_nonlin_model->set_p_sg_init(i, p_sg_in);
    }

    //! Create vector orthog poly using x map and owned sg map
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> 
    create_x_sg(Epetra_DataAccess CV = Copy, 
		const Epetra_Vector* v = NULL) const {
      return sg_nonlin_model->create_x_sg(CV, v);
    }

    //! Create vector orthog poly using p map
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> 
    create_p_sg(int l, Epetra_DataAccess CV = Copy, 
		const Epetra_Vector* v = NULL) const {
      return sg_nonlin_model->create_p_sg(l, CV, v);
    }

    //! Create vector orthog poly using g map
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> 
    create_g_sg(int l, Epetra_DataAccess CV = Copy, 
		const Epetra_Vector* v = NULL) const {
      return sg_nonlin_model->create_g_sg(l, CV, v);
    }

    //! Create multi-vector orthog poly using g map
    Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> 
    create_g_mv_sg(int l, int num_vecs, Epetra_DataAccess CV = Copy, 
		   const Epetra_MultiVector* v = NULL) const {
      return sg_nonlin_model->create_g_mv_sg(l, num_vecs, CV, v);
    }
    
  private:
    /** \brief . */
    Teuchos::RCP<const Epetra_Map> get_x_map() const;
    /** \brief . */
    Teuchos::RCP<const Epetra_Map> get_f_map() const;
    /** \brief . */
    Teuchos::RCP<const Epetra_Vector> get_x_init() const;
    /** \brief . */
    Teuchos::RCP<const Epetra_Map> get_p_map(int l) const;
    /** \brief . */
    Teuchos::RCP<const Epetra_Map> get_p_sg_map(int l) const;
    /** \brief . */
    void setProblemParamDefaults(Teuchos::ParameterList* appParams_);
    /** \brief . */
    void setSolverParamDefaults(Teuchos::ParameterList* appParams_);

    Teuchos::RCP<const Teuchos::ParameterList>
     getValidSGParameters() const;

    //@}
    
  private:
    
    enum SG_METHOD {
      SG_AD,
      SG_GLOBAL,
      SG_NI,
      SG_MPNI
    };

    enum SG_SOLVER { 
      SG_KRYLOV, 
      SG_GS, 
      SG_JACOBI 
    };

    Teuchos::RCP<Teuchos::ParameterList> piroParams;

    SG_METHOD sg_method;
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis;
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad;
    Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion;
    Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk;
    Teuchos::RCP<const EpetraExt::MultiComm> sg_comm;
    Teuchos::RCP<Stokhos::ParallelData> sg_parallel_data;

    //These are set in the constructor and used in evalModel
    Teuchos::RCP<EpetraExt::ModelEvaluator> sg_solver;
    Teuchos::RCP<Stokhos::SGModelEvaluator> sg_nonlin_model;
    
  };
  
}
}
#endif
