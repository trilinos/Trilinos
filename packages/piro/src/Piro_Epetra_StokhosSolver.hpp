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

#ifndef PIRO_EPETRA_STOKHOS_SOLVER_H
#define PIRO_EPETRA_STOKHOS_SOLVER_H

#include "EpetraExt_ModelEvaluator.h"
#include "Piro_Epetra_StokhosSolverFactory.hpp"

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
	       const Teuchos::RCP<NOX::Epetra::Observer>& noxObserver = Teuchos::null);

    //! Reset Stokhos solver parameters
    void resetSolverParameters(const Teuchos::ParameterList& new_solver_params);


    //@}

    ~StokhosSolver();


    /** \name Overridden from EpetraExt::ModelEvaluator . */
    //@{
    
    /** \brief . */
    Teuchos::RCP<const Epetra_Map> get_p_map(int l) const;

    /** \brief . */
    Teuchos::RCP<const Epetra_Map> get_g_map(int j) const;
    
    /** \brief . */
    Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;

    /** \brief . */
    //  Teuchos::RCP<Epetra_Operator> create_W() const;
    /** \brief . */
    EpetraExt::ModelEvaluator::InArgs createInArgs() const;
    /** \brief . */
    EpetraExt::ModelEvaluator::OutArgs createOutArgs() const;
    /** \brief . */
    void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

    //@}

    /** \name Accessors */
    //@{

    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >
    getBasis() const { return sg_solver_factory.getBasis(); }

    Teuchos::RCP<const Stokhos::Quadrature<int,double> >
    getQuad() const { return sg_solver_factory.getQuad(); }

    Teuchos::RCP<Stokhos::SGModelEvaluator>
    get_sg_model() const { return sg_nonlin_model; }

    //! Set initial solution polynomial
    void set_x_sg_init(const Stokhos::EpetraVectorOrthogPoly& x_sg_in) {
      sg_nonlin_model->set_x_sg_init(x_sg_in);
    }

    //! Return initial SG x
    Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly> 
    get_x_sg_init() const {
      return sg_nonlin_model->get_x_sg_init();
    }

    //! Set initial parameter polynomial
    void set_p_sg_init(int i, const Stokhos::EpetraVectorOrthogPoly& p_sg_in) {
      sg_nonlin_model->set_p_sg_init(i, p_sg_in);
    }

    //! Get initial parameter polynomial
    Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly> 
    get_p_sg_init(int l) const {
      return sg_nonlin_model->get_p_sg_init(l);
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
		const Epetra_Vector* v = NULL) const;

    //! Create multi-vector orthog poly using g map
    Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> 
    create_g_mv_sg(int l, int num_vecs, Epetra_DataAccess CV = Copy, 
		   const Epetra_MultiVector* v = NULL) const;

    //@}
    
  private:

    /** \brief . */
    Teuchos::RCP<const Epetra_Map> get_x_map() const;
    /** \brief . */
    Teuchos::RCP<const Epetra_Map> get_f_map() const;
    /** \brief . */
    Teuchos::RCP<const Epetra_Vector> get_x_init() const;
    
  private:
    
    Teuchos::RCP<Teuchos::ParameterList> piroParams;
    Piro::Epetra::StokhosSolverFactory sg_solver_factory;
    Teuchos::RCP<Stokhos::SGModelEvaluator> sg_nonlin_model;
    Teuchos::RCP<Stokhos::SGInverseModelEvaluator> sg_solver;
    
  };
  
}
}
#endif
