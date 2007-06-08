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

#ifndef EPETRA_EXT_DIAGONAL_TRANSIENT_MODEL_HPP
#define EPETRA_EXT_DIAGONAL_TRANSIENT_MODEL_HPP


#include "EpetraExt_ModelEvaluator.h"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_Array.hpp"


class Epetra_Comm;
class Epetra_CrsGraph;


namespace EpetraExt {


/** \brief Simple transient diagonal model for an implicit or explicit ODE.
 *
 * ToDo: Finish Documentation!
 */
class DiagonalTransientModel
  : public ::EpetraExt::ModelEvaluator,
    public Teuchos::VerboseObject<DiagonalTransientModel>,
    public Teuchos::ParameterListAcceptor
{
public:

  /** \name Constructors, Initializers, Misc. */
  //@{

  /** \brief . */
  DiagonalTransientModel(
    Teuchos::RefCountPtr<Epetra_Comm> const& epetra_comm
    );

  /** \brief Return the exact solution as a function of time. */
  Teuchos::RefCountPtr<const Epetra_Vector>
  getExactSolution(
    const double t, const Epetra_Vector *coeff_s_p = 0
    ) const;

  /** \brief Return the exact sensitivity of x as a function of time. */
  Teuchos::RefCountPtr<const Epetra_MultiVector>
  getExactSensSolution(
    const double t, const Epetra_Vector *coeff_s_p = 0
    ) const;

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> getParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Map> get_x_map() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Map> get_f_map() const;
  /** \breif . */
  Teuchos::RefCountPtr<const Epetra_Map> get_p_map(int l) const;
  /** \breif . */
  Teuchos::RefCountPtr<const Epetra_Map> get_g_map(int j) const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_x_init() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_x_dot_init() const;
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

public:

  enum ELambdaFit { LAMBDA_FIT_LINEAR, LAMBDA_FIT_RANDOM };

private:

  // /////////////////////////////////////
  // Private types

  typedef Teuchos::Array<double> coeff_s_t;
  typedef Teuchos::Array<int> coeff_s_idx_t;
  typedef Teuchos::Array<Teuchos::RefCountPtr<const Epetra_Map> >  RCP_Eptra_Map_Array_t;
  typedef Teuchos::Array<Teuchos::RefCountPtr<Epetra_Vector> > RCP_Eptra_Vector_Array_t;

  // /////////////////////////////////////
  // Private member data

  Teuchos::RefCountPtr<Teuchos::ParameterList> paramList_;
  Teuchos::RefCountPtr<Epetra_Comm> epetra_comm_;
  Teuchos::RefCountPtr<Epetra_Map> epetra_map_;
  bool implicit_;
  int numElements_;
  double lambda_min_;
  double lambda_max_;
  coeff_s_t coeff_s_;
  coeff_s_idx_t coeff_s_idx_;
  ELambdaFit lambda_fit_;
  double x0_;
  bool exactSolutionAsResponse_;
  Teuchos::RefCountPtr<Epetra_Vector> lambda_;
  Teuchos::RefCountPtr<Epetra_CrsGraph> W_graph_;
  int Np_;
  int np_;
  int Ng_;
  RCP_Eptra_Map_Array_t map_p_;
  RCP_Eptra_Map_Array_t map_g_;
  RCP_Eptra_Vector_Array_t p_init_;
  Teuchos::RefCountPtr<Epetra_Vector> x_init_;
  Teuchos::RefCountPtr<Epetra_Vector> x_dot_init_;

  mutable Teuchos::RefCountPtr<const Epetra_Vector> coeff_s_p_;

  bool isIntialized_;

  // /////////////////////////////////////
  // Private member functions

  void initialize();

  void set_coeff_s_p( 
    const Teuchos::RefCountPtr<const Epetra_Vector> &coeff_s_p
    ) const;

  void unset_coeff_s_p() const;

  int coeff_s_idx(int i) const
    {
      return coeff_s_idx_[i];
    }

  double coeff_s(int i) const
    {
      return (*coeff_s_p_)[coeff_s_idx(i)];
    }

};


/** \brief Nonmember constructor.
 *
 * \relates DiagonalTransientModel.
 */
Teuchos::RefCountPtr<DiagonalTransientModel>
diagonalTransientModel(
  Teuchos::RefCountPtr<Epetra_Comm> const& epetra_comm,
  Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList = Teuchos::null
  );


} // namespace EpetraExt


// RAB: Note, I wrapped this example code in a namespace mainly to make the
// later definition of the nonmember functions safe (see the Thyra
// Coding Guildelines document).


#endif // EPETRA_EXT_DIAGONAL_TRANSIENT_MODEL_HPP
