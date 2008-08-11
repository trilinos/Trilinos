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

#ifndef VANDERPOL_OSCILLATOR_HPP
#define VANDERPOL_OSCILLATOR_HPP

#include "EpetraExt_ModelEvaluator.h"
#include "Teuchos_ParameterList.hpp"
#include "Rythmos_ConfigDefs.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"

#include <fstream>

/** \brief Base interface for evaluating a stateless "model".
 *
 * ToDo: Finish Documentation!
 */
class VanDerPolOscillator : public EpetraExt::ModelEvaluator {
public:

  // Constructor
  VanDerPolOscillator(Teuchos::RCP<Epetra_Comm> &epetra_comm_ptr, 
		      Teuchos::ParameterList &params);

  // Destructor
  ~VanDerPolOscillator();

  void saveSolution(const Epetra_Vector& x, double t);

  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_x_map() const;

  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_f_map() const;

  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_x_init() const;

  /** \brief . */
  Teuchos::RCP<Epetra_Operator> create_W() const;

  /** \brief . */
  InArgs createInArgs() const;

  /** \brief . */
  OutArgs createOutArgs() const;

  /** \brief . */
  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  //@}

protected:

  /** \brief . */
  template <typename Scalar>
  void evalVField( const Scalar& x1, const Scalar& x2,
		   Scalar& f1, Scalar& f2) const {
    f1 = x2;
    f2 = ( (1.0 - x1*x1)*x2-x1) / omega_;
  }

protected:

  // Epetra Comm:
  Teuchos::RCP<Epetra_Comm> epetra_comm_ptr_;

  // Epetra Map:
  Teuchos::RCP<Epetra_Map> epetra_map_ptr_;
    
  // Coefficients for ODE
  double omega_;
  
  // Initial condition for the problem:
  Teuchos::RCP<Epetra_Vector> x0_;

  // File to output solutions to
  std::ofstream outfile_;

  // Do we need Jacobian's?
  bool implicit_;

  // Graph of Jacobian
  Teuchos::RCP<Epetra_CrsGraph> W_graph_;
};

#endif // EXAMPLE_APPLICATION_HPP
