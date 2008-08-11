//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef EXAMPLE_APPLICATION_HPP
#define EXAMPLE_APPLICATION_HPP

#include "EpetraExt_ModelEvaluator.h"
#include "Teuchos_ParameterList.hpp"
#include "Rythmos_ConfigDefs.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"

#include "1DfemTransient.H"

/** \brief Base interface for evaluating a stateless "model".
 *
 * ToDo: Finish Documentation!
 */
class ExampleApplication1Dfem : public EpetraExt::ModelEvaluator {
public:

  // Constructor
  ExampleApplication1Dfem(Teuchos::RCP<Epetra_Comm> &epetra_comm_ptr, Teuchos::ParameterList &params);

  // Initialization
  void initialize(Teuchos::RCP<Epetra_Comm> &epetra_comm_ptr, Teuchos::ParameterList &params);

  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_x_map() const;

  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_f_map() const;

  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_x_init() const;

  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_x_dot_init() const;

  /** \brief . */
  Teuchos::RCP<Epetra_Operator> create_W() const;

  /** \brief . */
  InArgs createInArgs() const;

  /** \brief . */
  OutArgs createOutArgs() const;

  /** \brief . */
  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  /** \brief . */
  Teuchos::RCP<Epetra_Vector> get_exact_solution( double t ) const;

  //@}

private:

    // Epetra Comm:
    Teuchos::RCP<Epetra_Comm> epetra_comm_ptr_;
    // Epetra Map:
    Teuchos::RCP<const Epetra_Map> epetra_map_ptr_;
    
    // Global number of unknowns:
    int numElements_;

    Teuchos::RCP<Epetra_CrsGraph> W_graph_;

    // This ModelEvaluator object will own a TransientInterface object
    Teuchos::RCP<TransientInterface> problemInterfacePtr_;
};

#endif // EXAMPLE_APPLICATION_HPP
