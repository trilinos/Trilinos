// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef FEAPP_MODELEVALUATOR_HPP
#define FEAPP_MODELEVALUATOR_HPP

#include "Teuchos_RefCountPtr.hpp"
#include "EpetraExt_ModelEvaluator.h"

#include "FEApp_Application.hpp"

namespace FEApp {

  class ModelEvaluator : public EpetraExt::ModelEvaluator {
  public:

    // Constructor
    ModelEvaluator(const Teuchos::RefCountPtr<FEApp::Application>& app);

    /** \name Overridden from EpetraExt::ModelEvaluator . */
    //@{

    //! Return solution vector map
    Teuchos::RefCountPtr<const Epetra_Map> get_x_map() const;

    //! Return residual vector map
    Teuchos::RefCountPtr<const Epetra_Map> get_f_map() const;

    //! Return initial solution
    Teuchos::RefCountPtr<const Epetra_Vector> get_x_init() const;

    //! Create W = alpha*M + beta*J matrix
    Teuchos::RefCountPtr<Epetra_Operator> create_W() const;

    //! Create InArgs
    InArgs createInArgs() const;

    //! Create OutArgs
    OutArgs createOutArgs() const;

    //! Evaluate model on InArgs
    void evalModel(const InArgs& inArgs, const OutArgs& outArgs) const;

    //@}

  protected:

    //! Application object
    Teuchos::RefCountPtr<FEApp::Application> app;

  };

}

#endif // FEAPP_MODELEVALUATOR_HPP
