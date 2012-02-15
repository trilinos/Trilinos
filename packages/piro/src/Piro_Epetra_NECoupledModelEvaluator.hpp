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

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Stokhos_EpetraVectorOrthogPoly.hpp"
#include "Stokhos_EpetraMultiVectorOrthogPoly.hpp"

namespace Piro {
  namespace Epetra {

    class NECoupledModelEvaluator : 
      public EpetraExt::ModelEvaluator,
      public Teuchos::VerboseObject<NECoupledModelEvaluator> {
    public:

      /** \brief . */
      NECoupledModelEvaluator(
	const Teuchos::RCP<EpetraExt::ModelEvaluator>& modelA, 
	const Teuchos::RCP<EpetraExt::ModelEvaluator>& modelB,
	const Teuchos::RCP<Teuchos::ParameterList>& piroParamsA,
	const Teuchos::RCP<Teuchos::ParameterList>& piroParamsB,
	const Teuchos::RCP<Teuchos::ParameterList>& params,
	const Teuchos::RCP<const Epetra_Comm>& comm);

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
	const InArgs& inArgs,
	const InArgs& solver_inargs, 
	const OutArgs& solver_outargs,
	const Teuchos::RCP<EpetraExt::ModelEvaluator>& model,
	const Teuchos::RCP<EpetraExt::ModelEvaluator>& solver,
	const Teuchos::RCP<Teuchos::ParameterList>& solver_params,
	InArgs& reduced_inargs, 
	OutArgs& reduced_outargs,
	Teuchos::RCP<EpetraExt::ModelEvaluator>& reduced_solver,
	Teuchos::RCP<Teuchos::ParameterList>& reduced_params,
	Teuchos::RCP<const Teuchos::Array< Teuchos::Array<double> > >& red_basis_vals) const;

      void do_dimension_projection(
	const InArgs& inArgs, 
	const InArgs& reduced_inargs, 
	const OutArgs& reduced_outargs,
	const Teuchos::RCP<const Teuchos::Array< Teuchos::Array<double> > >& red_basis_vals,
	OutArgs& solver_outargs) const;

    private:

      // /////////////////////////////////////
      // Private member data
      
      typedef Stokhos::StandardStorage<int,double> StorageType;

      Teuchos::RCP<EpetraExt::ModelEvaluator> modelA;
      Teuchos::RCP<EpetraExt::ModelEvaluator> modelB;
      Teuchos::RCP<Teuchos::ParameterList> piroParamsA;
      Teuchos::RCP<Teuchos::ParameterList> piroParamsB;
      Teuchos::RCP<Teuchos::ParameterList> params;
      Teuchos::RCP<const Epetra_Comm> comm;
      Teuchos::RCP<EpetraExt::ModelEvaluator> solverA;
      Teuchos::RCP<EpetraExt::ModelEvaluator> solverB;
      Teuchos::RCP<Piro::Epetra::StokhosSolver> sgSolverA;
      Teuchos::RCP<Piro::Epetra::StokhosSolver> sgSolverB;
      int pIndexA;
      int pIndexB;
      int gIndexA;
      int gIndexB;
      int n_p_A;
      int n_p_B;
      int n_g_A;
      int n_g_B;
      int num_params_A;
      int num_params_B;
      int num_params_total;
      int num_responses_A;
      int num_responses_B;
      int num_responses_total;
      bool supports_W;
      Teuchos::Array<int> param_map;
      Teuchos::Array<int> response_map;
      Teuchos::RCP<const Epetra_Map> p_map_A;
      Teuchos::RCP<const Epetra_Map> p_map_B;
      Teuchos::RCP<const Epetra_Map> g_map_A;
      Teuchos::RCP<const Epetra_Map> g_map_B;
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
      Teuchos::RCP<Epetra_Vector> p_A;
      Teuchos::RCP<Epetra_Vector> p_B;
      Teuchos::RCP<Epetra_Vector> g_A;
      Teuchos::RCP<Epetra_Vector> g_B;
      EDerivativeMultiVectorOrientation dgdp_A_layout;
      EDerivativeMultiVectorOrientation dgdp_B_layout;
      Teuchos::RCP<Epetra_MultiVector> dgdp_A;
      Teuchos::RCP<Epetra_MultiVector> dgdp_B;
     
      // Stochastic Galerkin data
      bool supports_x_sg;
      bool supports_f_sg;
      bool supports_W_sg;
      mutable Teuchos::RCP<const Epetra_BlockMap> sg_overlap_map;
      mutable OutArgs::sg_vector_t p_A_sg;
      mutable OutArgs::sg_vector_t p_B_sg;
      mutable OutArgs::sg_vector_t g_A_sg;
      mutable OutArgs::sg_vector_t g_B_sg;
      mutable EDerivativeMultiVectorOrientation dgdp_A_sg_layout;
      mutable EDerivativeMultiVectorOrientation dgdp_B_sg_layout;
      mutable Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > dgdp_A_sg;
      mutable Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > dgdp_B_sg;

      bool reduce_dimension;
      mutable Teuchos::RCP<const Stokhos::Quadrature<int,double> > st_quad;
      mutable Teuchos::RCP<const Teuchos::Array< Teuchos::Array<double> > > red_basis_vals_A;
      mutable Teuchos::RCP<const Teuchos::Array< Teuchos::Array<double> > > red_basis_vals_B;
    };

  }

}

#endif
