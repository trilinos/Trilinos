//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER

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
 * The explicit ODE form of the model is:

 \verbatim

  x_dot(i) = f_hat(x(i), gamma(i), s(i), t), for i = 0...n-1, on t in [0,t_f]

 \endverbatim

 * where:

 \verbatim

  f_hat(x(i), gamma(i), s(i), t) = gama(i)*x(i) + exp(gamma(i)*t)*sin(s(i),t)

 \endverbatim

 * The implicit ODE form of the model i:


 \verbatim


  f(i)(x_dot(i), x(i), t) = x_dot(i) - f_hat(x(i), gamma(i), s(i), t),
  
    for i = 0...n-1, on t in [0,t_f]

 \endverbatim

 * This is a diagonal problem so it does not make the greatest test problem
 * but it does make it easy to derive tests for as a starter.
 *
 * The coefficients <tt>s</tt> can be exposed as model parameters and are
 * called <tt>coeff_s_p</tt> in the code.  The selection of the coefficients
 * is handled through the
 
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
    Teuchos::RCP<Epetra_Comm> const& epetra_comm
    );
  
  /** \brief Return the model vector <tt>gamma</tt>, */
  Teuchos::RCP<const Epetra_Vector> get_gamma() const;

  /** \brief Return the exact solution as a function of time. */
  Teuchos::RCP<const Epetra_Vector>
  getExactSolution(
    const double t, const Epetra_Vector *coeff_s_p = 0
    ) const;

  /** \brief Return the exact sensitivity of x as a function of time. */
  Teuchos::RCP<const Epetra_MultiVector>
  getExactSensSolution(
    const double t, const Epetra_Vector *coeff_s_p = 0
    ) const;

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_x_map() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Map> get_f_map() const;
  /** \breif . */
  Teuchos::RCP<const Epetra_Map> get_p_map(int l) const;
  /** \breif . */
  Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
  /** \breif . */
  Teuchos::RCP<const Epetra_Map> get_g_map(int j) const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_x_init() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_x_dot_init() const;
  /** \brief . */
  Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;
  /** \brief . */
  Teuchos::RCP<Epetra_Operator> create_W() const;
  /** \brief . */
  InArgs createInArgs() const;
  /** \brief . */
  OutArgs createOutArgs() const;
  /** \brief . */
  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  //@}

public:

  enum EGammaFit { GAMMA_FIT_LINEAR, GAMMA_FIT_RANDOM };

private:

  // /////////////////////////////////////
  // Private types

  typedef Teuchos::Array<double> coeff_s_t;
  typedef Teuchos::Array<int> coeff_s_idx_t;
  typedef Teuchos::Array<Teuchos::RCP<const Epetra_Map> >  RCP_Eptra_Map_Array_t;
  typedef Teuchos::Array<Teuchos::RCP<Epetra_Vector> > RCP_Eptra_Vector_Array_t;
  typedef Teuchos::Array<Teuchos::RCP<Teuchos::Array<std::string> > > RCP_Array_String_Array_t;
  

  // /////////////////////////////////////
  // Private member data

  Teuchos::RCP<Teuchos::ParameterList> paramList_;
  Teuchos::RCP<Epetra_Comm> epetra_comm_;
  Teuchos::RCP<Epetra_Map> epetra_map_;
  bool implicit_;
  int numElements_;
  double gamma_min_;
  double gamma_max_;
  coeff_s_t coeff_s_;
  coeff_s_idx_t coeff_s_idx_;
  EGammaFit gamma_fit_;
  double x0_;
  bool exactSolutionAsResponse_;
  Teuchos::RCP<Epetra_Vector> gamma_;
  Teuchos::RCP<Epetra_CrsGraph> W_graph_;
  int Np_;
  int np_;
  int Ng_;
  RCP_Eptra_Map_Array_t map_p_;
  RCP_Array_String_Array_t names_p_;
  RCP_Eptra_Map_Array_t map_g_;
  RCP_Eptra_Vector_Array_t p_init_;
  Teuchos::RCP<Epetra_Vector> x_init_;
  Teuchos::RCP<Epetra_Vector> x_dot_init_;

  mutable Teuchos::RCP<const Epetra_Vector> coeff_s_p_;

  bool isIntialized_;

  // /////////////////////////////////////
  // Private member functions

  void initialize();

  void set_coeff_s_p( 
    const Teuchos::RCP<const Epetra_Vector> &coeff_s_p
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
Teuchos::RCP<DiagonalTransientModel>
diagonalTransientModel(
  Teuchos::RCP<Epetra_Comm> const& epetra_comm,
  Teuchos::RCP<Teuchos::ParameterList> const& paramList = Teuchos::null
  );


} // namespace EpetraExt


// RAB: Note, I wrapped this example code in a namespace mainly to make the
// later definition of the nonmember functions safe (see the Thyra
// Coding Guildelines document).


#endif // EPETRA_EXT_DIAGONAL_TRANSIENT_MODEL_HPP
