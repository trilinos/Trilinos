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

#include "NOX.H"
#include "NOX_Epetra.H"
#include "Epetra_Vector.h"
#include "Epetra_LocalMap.h"
#include "NOX_Epetra_ModelEvaluatorInterface.H"
#include <NOX_Epetra_MultiVector.H>
#include <NOX_Epetra_Observer.H>

#include "LOCA_GlobalData.H"
#include "LOCA_Epetra_TransposeLinearSystem_AbstractStrategy.H"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "EpetraExt_ModelEvaluator.h"

/** \brief Epetra-based NOX Solver
 *
 * This class will 
 * 
 * ToDo: Finish documentation!
 */

namespace Piro {
namespace Epetra {
class NOXSolver
    : public EpetraExt::ModelEvaluator
{
  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  NOXSolver(Teuchos::RCP<Teuchos::ParameterList> piroParams,
            Teuchos::RCP<EpetraExt::ModelEvaluator> model,
            Teuchos::RCP<NOX::Epetra::Observer> observer = Teuchos::null,
	    Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> interface = 
	    Teuchos::null,
	    Teuchos::RCP<NOX::Epetra::LinearSystem> linsys = Teuchos::null
            );


  //@}

  ~NOXSolver();


  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  Teuchos::RCP<const Epetra_Map> get_g_map(int j) const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;
  /** \brief . */
  Teuchos::RCP<Epetra_Operator> create_DgDp_op( int j, int l ) const;
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
  void setProblemParamDefaults(Teuchos::ParameterList* piroParams_);
  /** \brief . */
  void setSolverParamDefaults(Teuchos::ParameterList* piroParams_);

  //@}

  private:

   //These are set in the constructor and used in evalModel
   mutable Teuchos::RCP<Teuchos::ParameterList> piroParams;
   Teuchos::RCP<EpetraExt::ModelEvaluator> model;
   Teuchos::RCP<NOX::Epetra::Observer> observer;
   NOX::Utils utils;

   Teuchos::RCP<NOX::Solver::Generic> solver;
   Teuchos::RCP<NOX::Epetra::Vector> currentSolution;
   Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> interface;
   int num_p;
   int num_g;

   Teuchos::RCP<NOX::Epetra::Group> grp;
   Teuchos::RCP<LOCA::GlobalData> globalData;
   Teuchos::RCP<LOCA::Epetra::TransposeLinearSystem::AbstractStrategy> tls_strategy;

   enum DerivativeLayout { OP, COL, ROW };
};

}
}
#endif
