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

#ifndef PIRO_NOXSOLVER_H
#define PIRO_NOXSOLVER_H

#include <iostream>
// NOX Objects
#include "NOX.H"
#include "NOX_Thyra.H"

#include "Thyra_ModelEvaluatorDefaultBase.hpp"
#include "Thyra_ResponseOnlyModelEvaluatorBase.hpp"
#include "Teuchos_RCP.hpp"

/** \brief Thyra-based Model Evaluator for NOX solves */

namespace Piro {

template <typename Scalar>
class NOXSolver
    : public Thyra::ResponseOnlyModelEvaluatorBase<Scalar>
{

  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  NOXSolver(Teuchos::RCP<Teuchos::ParameterList> appParams,
            Teuchos::RCP< Thyra::ModelEvaluatorDefaultBase<Scalar> > model
            );

  //@}

  ~NOXSolver();

  /** \name Overridden from Thyra::ModelEvaluatorBase . */
  //@{

  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int i) const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int i) const;

  /** \brief . */
  void evalModelImpl( const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
                  const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs ) const;

  private:

  //@}

  private:

   //These are set in the constructor and used in evalModel
   mutable Teuchos::RCP<Teuchos::ParameterList> appParams;
   Teuchos::RCP< Thyra::ModelEvaluatorDefaultBase<Scalar> > model;
   Teuchos::RCP<Teuchos::FancyOStream> out;

   Teuchos::RCP<NOX::Solver::Generic> solver;

   int num_p;
   int num_g;
};

}

#include "Piro_NOXSolver_Def.hpp"
#endif
