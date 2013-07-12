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

#ifndef PIRO_EPETRA_TRAPEZOIDRULESOLVER_H
#define PIRO_EPETRA_TRAPEZOIDRULESOLVER_H

#include <iostream>

#include "Epetra_Vector.h"
#include "EpetraExt_ModelEvaluator.h"
#include "NOX_Epetra_Observer.H"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Piro_Epetra_NOXSolver.hpp"

namespace Piro {
namespace Epetra {

/** \ingroup Piro_Epetra_solver_grp */
class TrapezoidDecorator
    : public EpetraExt::ModelEvaluator
{

  typedef double Scalar;

  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  TrapezoidDecorator( 
                Teuchos::RCP<EpetraExt::ModelEvaluator>& model
                );

  //@}

  ~TrapezoidDecorator();


  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  Teuchos::RCP<const Epetra_Map> get_g_map(int j) const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;
  /** \brief . */
  Teuchos::RCP<Epetra_Operator> create_W() const;
  /** \brief . */
  EpetraExt::ModelEvaluator::InArgs createInArgs() const;
  /** \brief . */
  EpetraExt::ModelEvaluator::OutArgs createOutArgs() const;
  /** \brief . */
  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  //! Method to give info to compute xDotDot(x), so that the
  // NOX solver can treat the time dep problem as steady 
  void injectData(const Teuchos::RCP<Epetra_Vector>& x_, 
                  const Teuchos::RCP<Epetra_Vector>& x_pred_, 
                  double fdt2_, double time_ );

  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_x_map() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_f_map() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_x_init() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_x_dot_init() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_p_map(int l) const;

  //@}

  private:

   //These are set in the constructor and used in evalModel
   Teuchos::RCP<EpetraExt::ModelEvaluator> model;
   Teuchos::RCP<Epetra_Vector> xDotDot;
   Teuchos::RCP<Epetra_Vector> x_pred;
   Teuchos::RCP<Epetra_Vector> x_save;
   double fdt2;
   double time; 

};


class TrapezoidRuleSolver
    : public EpetraExt::ModelEvaluator
{

  typedef double Scalar;

  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  TrapezoidRuleSolver(Teuchos::RCP<Teuchos::ParameterList> appParams,
                Teuchos::RCP<EpetraExt::ModelEvaluator> model,
                Teuchos::RCP<NOX::Epetra::Observer> observer = Teuchos::null
                );

  //@}

  ~TrapezoidRuleSolver();


  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

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
  void setProblemParamDefaults(Teuchos::ParameterList* appParams_);
  /** \brief . */
  void setSolverParamDefaults(Teuchos::ParameterList* appParams_);
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList>
    getValidTrapezoidRuleParameters() const;

  //@}

  private:
   //These are set in the constructor and used in evalModel
   mutable Teuchos::RCP<Teuchos::ParameterList> appParams;
   Teuchos::RCP<Piro::Epetra::TrapezoidDecorator> model;
   Teuchos::RCP<Piro::Epetra::NOXSolver> noxSolver;
   Teuchos::RCP<NOX::Epetra::Observer> observer;
   Teuchos::RCP<Teuchos::FancyOStream> out;
   Teuchos::EVerbosityLevel solnVerbLevel;

   int num_p;
   int num_g;

   int numTimeSteps;
   double t_init, t_final, delta_t;
};




}
}
#endif
