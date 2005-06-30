// @HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef EXAMPLE_APPLICATION_HPP
#define EXAMPLE_APPLICATION_HPP

#include "EpetraExt_ModelEvaluator.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Rythmos_ConfigDefs.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"

/** \brief Base interface for evaluating a stateless "model".
 *
 * ToDo: Finish Documentation!
 */
class ExampleApplication2 : public EpetraExt::ModelEvaluator {
public:

  // Constructor
  ExampleApplication2(Teuchos::ParameterList &params);

  // return ODE decay coefficient
  Teuchos::RefCountPtr<const Epetra_Vector> get_coeff() const;

  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Map> get_x_map() const;

  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Map> get_f_map() const;

  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_x_init() const;

  /** \brief . */
  InArgs createInArgs() const;

  /** \brief . */
  OutArgs createOutArgs() const;

  /** \brief . */
  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  //@}

private:

    // Epetra Comm:
    Teuchos::RefCountPtr<Epetra_Comm> epetra_comm_ptr_;
    // Epetra Map:
    Teuchos::RefCountPtr<Epetra_Map> epetra_map_ptr_;
    
    // Global number of unknowns:
    int numElements_;
    // Coefficients for ODE
    double lambda_min_;
    double lambda_max_;
    std::string lambda_fit_;
    Teuchos::RefCountPtr<Epetra_Vector> lambda_ptr_;
    // Constant initial condition for the problem:
    double x0_;

};

#endif // EXAMPLE_APPLICATION_HPP
