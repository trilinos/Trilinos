/*
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
*/

#ifndef PIRO_EPETRA_NE_COUPLED_MODEL_EVALUATOR_HPP
#define PIRO_EPETRA_NE_COUPLED_MODEL_EVALUATOR_HPP

#include "EpetraExt_ModelEvaluator.h"
#include "Teuchos_VerboseObject.hpp"
#include "Piro_Epetra_StokhosSolver.hpp"

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Stokhos_EpetraVectorOrthogPoly.hpp"
#include "Stokhos_EpetraMultiVectorOrthogPoly.hpp"

namespace Piro {
  namespace Epetra {

    class AbstractNetworkModel {
    public:

      //! Constructor
      AbstractNetworkModel() {}

      //! Destructor
      virtual ~AbstractNetworkModel() {}

      //! evaluate model
      virtual void evalModel(
	const Teuchos::Array<EpetraExt::ModelEvaluator::InArgs>& model_inargs, 
	const Teuchos::Array<EpetraExt::ModelEvaluator::OutArgs>& model_outargs,
	const EpetraExt::ModelEvaluator::InArgs& network_inargs, 
	const EpetraExt::ModelEvaluator::OutArgs& network_outargs,
	const Teuchos::Array<int>& n_p,
	const Teuchos::Array<int>& n_g,
	const Teuchos::Array< Teuchos::RCP<Epetra_Vector> >& p,
	const Teuchos::Array< Teuchos::RCP<Epetra_Vector> >& g,
	const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dgdp,
	const Teuchos::Array<EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation>& dgdp_layout,
	const Teuchos::Array<EpetraExt::ModelEvaluator::OutArgs::sg_vector_t>& p_sg,
	const Teuchos::Array<EpetraExt::ModelEvaluator::OutArgs::sg_vector_t>& g_sg,
	const Teuchos::Array<Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> >& dgdp_sg,
	const Teuchos::Array<EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation>& dgdp_sg_layout) const = 0;

    };

    class ParamToResponseNetworkModel : 
      public AbstractNetworkModel {

    public:

      //! Constructor
      ParamToResponseNetworkModel() {}

      //! Destructor
      virtual ~ParamToResponseNetworkModel() {}

      //! evaluate model
      virtual void evalModel(
	const Teuchos::Array<EpetraExt::ModelEvaluator::InArgs>& model_inargs, 
	const Teuchos::Array<EpetraExt::ModelEvaluator::OutArgs>& model_outargs,
	const EpetraExt::ModelEvaluator::InArgs& network_inargs, 
	const EpetraExt::ModelEvaluator::OutArgs& network_outargs,
	const Teuchos::Array<int>& n_p,
	const Teuchos::Array<int>& n_g,
	const Teuchos::Array< Teuchos::RCP<Epetra_Vector> >& p,
	const Teuchos::Array< Teuchos::RCP<Epetra_Vector> >& g,
	const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dgdp,
	const Teuchos::Array<EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation>& dgdp_layout,
	const Teuchos::Array<EpetraExt::ModelEvaluator::OutArgs::sg_vector_t>& p_sg,
	const Teuchos::Array<EpetraExt::ModelEvaluator::OutArgs::sg_vector_t>& g_sg,
	const Teuchos::Array<Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> >& dgdp_sg,
	const Teuchos::Array<EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation>& dgdp_sg_layout) const;

    };

    class NECoupledModelEvaluator : 
      public EpetraExt::ModelEvaluator,
      public Teuchos::VerboseObject<NECoupledModelEvaluator> {
    public:

      /** \brief . */
      NECoupledModelEvaluator(
	const Teuchos::Array<Teuchos::RCP<EpetraExt::ModelEvaluator> >& models,
	const Teuchos::Array<Teuchos::RCP<Teuchos::ParameterList> >& piroParams,
	const Teuchos::RCP<AbstractNetworkModel>& network_model,
	const Teuchos::RCP<Teuchos::ParameterList>& params,
	const Teuchos::RCP<const Epetra_Comm>& comm,
	const Teuchos::Array< Teuchos::RCP<NOX::Epetra::Observer> >& observers =
	Teuchos::Array<Teuchos::RCP<NOX::Epetra::Observer> >());

      /** \name Overridden from EpetraExt::ModelEvaluator . */
      //@{

      /** \brief . */
      Teuchos::RefCountPtr<const Epetra_Map> get_x_map() const;
      /** \brief . */
      Teuchos::RefCountPtr<const Epetra_Map> get_f_map() const;
      /** \brief . */
      Teuchos::RefCountPtr<const Epetra_Vector> get_x_init() const;
      /** \brief . */
      Teuchos::RefCountPtr<const Epetra_Map> get_p_map(int l) const;
      /** \brief . */
      Teuchos::RefCountPtr<const Epetra_Map> get_g_map(int j) const;
      //! Return array of parameter names
      Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
      /** \brief . */
      Teuchos::RefCountPtr<const Epetra_Vector> get_p_init(int l) const;
      /** \brief . */
      Teuchos::RefCountPtr<Epetra_Operator> create_W() const;
      /** \brief . */
      InArgs createInArgs() const;
      /** \brief . */
      OutArgs createOutArgs() const;
      /** \brief . */
      void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

      //@}

    protected:

      void do_dimension_reduction(
	int model_index,
	const InArgs& inArgs,
	const InArgs& solver_inargs, 
	const OutArgs& solver_outargs,
	const Teuchos::RCP<EpetraExt::ModelEvaluator>& model,
	const Teuchos::RCP<EpetraExt::ModelEvaluator>& solver,
	const Teuchos::RCP<Teuchos::ParameterList>& solver_params,
	InArgs& reduced_inargs, 
	OutArgs& reduced_outargs,
	Teuchos::RCP<EpetraExt::ModelEvaluator>& reduced_solver,
	Teuchos::RCP<Teuchos::ParameterList>& reduced_params) const;

      void do_dimension_projection(
	int model_index,
	const InArgs& inArgs, 
	const InArgs& reduced_inargs, 
	const OutArgs& reduced_outargs,
	OutArgs& solver_outargs) const;

    private:

      // /////////////////////////////////////
      // Private member data
      
      typedef Stokhos::StandardStorage<int,double> StorageType;

      Teuchos::Array<Teuchos::RCP<EpetraExt::ModelEvaluator> > models;
      Teuchos::Array< Teuchos::RCP<Teuchos::ParameterList> > piroParams;
      Teuchos::RCP<AbstractNetworkModel> network_model;
      Teuchos::RCP<Teuchos::ParameterList> params;
      Teuchos::RCP<const Epetra_Comm> comm;
      Teuchos::Array< Teuchos::RCP<NOX::Epetra::Observer> > observers;
      
      Teuchos::Array< Teuchos::RCP<EpetraExt::ModelEvaluator> > solvers;
      Teuchos::Array< Teuchos::RCP<Piro::Epetra::StokhosSolver> > sgSolvers;
      int n_models;
      Teuchos::Array<int> p_indices;
      Teuchos::Array<int> g_indices;     
      Teuchos::Array<int> n_p;
      Teuchos::Array<int> n_g;
      Teuchos::Array<int> num_params;
      Teuchos::Array<int> num_responses;
      int num_params_total;
      int num_responses_total;
      bool supports_W;
      Teuchos::Array< std::pair<int,int> > param_map;
      Teuchos::Array< std::pair<int,int> > response_map;

      mutable Teuchos::Array<EpetraExt::ModelEvaluator::InArgs> solver_inargs; 
      mutable Teuchos::Array<EpetraExt::ModelEvaluator::OutArgs> solver_outargs;

      Teuchos::Array< Teuchos::RCP<const Epetra_Map> > p_maps;
      Teuchos::Array< Teuchos::RCP<const Epetra_Map> > g_maps;
      
      Teuchos::RCP<Epetra_Map> x_map;
      Teuchos::RCP<Epetra_Map> f_map;
      Teuchos::RCP<Epetra_Map> x_overlap_map;
      Teuchos::RCP<Epetra_Map> f_overlap_map;
      Teuchos::RCP<Epetra_Import> x_importer;
      Teuchos::RCP<Epetra_Export> f_exporter;
      Teuchos::RCP<Epetra_Vector> x_overlap;
      Teuchos::RCP<Epetra_Vector> f_overlap;
      Teuchos::RCP<Epetra_CrsGraph> W_graph;
      Teuchos::RCP<Epetra_CrsGraph> W_overlap_graph;
      Teuchos::RCP<Epetra_CrsMatrix> W_overlap;
      Teuchos::RCP<Epetra_Vector> x_init;

      Teuchos::Array< Teuchos::RCP<Epetra_Vector> > p;
      Teuchos::Array< Teuchos::RCP<Epetra_Vector> > g;
      Teuchos::Array< EDerivativeMultiVectorOrientation > dgdp_layout;
      Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> > dgdp;
     
      // Stochastic Galerkin data
      bool supports_x_sg;
      bool supports_f_sg;
      bool supports_W_sg;
      mutable Teuchos::RCP<const Epetra_BlockMap> sg_overlap_map;
      mutable OutArgs::sg_vector_t x_sg_overlap;
      mutable OutArgs::sg_vector_t f_sg_overlap;
      mutable OutArgs::sg_operator_t W_sg_overlap;
      mutable Teuchos::Array<OutArgs::sg_vector_t> p_sg;
      mutable Teuchos::Array<OutArgs::sg_vector_t> g_sg;
      mutable Teuchos::Array<EDerivativeMultiVectorOrientation> dgdp_sg_layout;
      mutable Teuchos::Array<Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> > dgdp_sg;

      Teuchos::Array<int> reduce_dimension;
      mutable Teuchos::RCP<const Stokhos::Quadrature<int,double> > st_quad;
    };

  }

}

#endif
