/*
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
*/

#ifndef GLP_APP_ADV_DIFF_REACT_OPT_MODEL_HPP
#define GLP_APP_ADV_DIFF_REACT_OPT_MODEL_HPP

#include "EpetraExt_ModelEvaluator.h"
#include "GLpApp_GLpYUEpetraDataPool.hpp"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Array.hpp"

namespace GLpApp {

/** \brief PDE-constrained inverse problem based on a 2D discretization of a
 * diffusion/reaction system.
 *
 * The model evaluator subclass is used to represent the
 * simulation-constrained optimization problem:

 \verbatim

  min   g(x,p)
  s.t.  f(x,p) = 0;

 \endverbatim

 * where:<ul>
 *
 * <li><tt>x</tt> is the vector of discretized concentations of the species in
 * the 2D domain.
 *
 * <li><tt>p</tt> is the global vector of coefficients of a sine series basis
 * (see <tt>B_bar</tt> below).
 *
 * <li><tt>f(x,p) = A*x + reationRate*Ny(x) + B*(B_bar*p)</tt> is the
 * discretized 2D diffusion/reaction PDE.
 *
 * <li><tt>g(x,p) = 0.5 * (x-q)'*H*(x-q) +
 * 0.5*regBeta*(B_bar*p)'*R*(B_bar*p)</tt> is the least squares objective
 * function.
 *
 * <li><tt>A</tt> is the discretized Laplacian operator for the diffusion part
 * of the PDE state equation.  This matrix is constant, square and
 * singular</tt>.
 *
 * <li><tt>B</tt> is the sensitivity of the flux boundary conditions.  This is
 * a constant rectangular matrix.
 *
 * <li><tt>B_bar</tt> are the sine series coefficients with a column dimension
 * of <tt>np</tt>.
 *
 * <li><tt>Ny(x)</tt> is the nonlinear terms for the discretized reaction over
 * the 2D domain.
 *
 * <li><tt>reactionRate</tt> is the relative reaction rate which must take on
 * a non-zero value to form a solvable problem.
 *
 * <li><tt>H</tt> is the symmetric positive definite mass matrix for the
 * problem (i.e. the discretization of the inner product operator over the 2D
 * domain).
 *
 * <li><tt>q</tt> is a matching or target vector for the state <tt>x</tt> over
 * the 2D domain of the problem.
 *
 * <li><tt>R</tt> is the symmetric positive definite discretization of the
 * inner product of the flux function over the boundary of the 2D domain.
 *
 * <li><tt>regBeta</tt> is a regularization parameter that must be greater
 * than zero.
 *
 * </ul>
 *
 * The nuts and bolts of the implementation for this problem are contained in
 * the C++ class <tt>GLpApp::GLpYUEpetraDataPool</tt> that was originally
 * implemented by Denis Ridzal while a student at Rice University.  The class
 * <tt>GLpApp::GLpYUEpetraDataPool</tt> implements the basic operators and
 * nonlinear functions but this class puts them together to form a valid the
 * model in terms of a model evaluator interface.
 *
 * This example problem demonstrates a few different aspects of the
 * <tt>EpetraExt::ModelEvaluator</tt> interface:</ul>
 *
 * <li>How to manage parallel vector data.  The state variables in <tt>x</tt>
 * are managed as fully distributed parallel data while the flux sine-series
 * parameter coefficients <tt>p</tt> are managed as locally replicated data.
 *
 * <tt>Demonstrates shared compuation between the objective function
 * <tt>g(x,p)</tt> and the simulation equality constraints <tt>f(x,p)</tt> and
 * their derivatives.  The intermediate vector <tt>B_bar*p</tt> is computed
 * only once and is shared with the computation of <tt>g</tt> and <tt>f</tt>.
 * The intermediate vector <tt>R*(B_bar*p)</tt> is computed once and shared
 * between the computation of <tt>g</tt> and <tt>DgDp</tt>.
 *
 * </ul>
 *
 * The functions <tt>AdvDiffReactOptModel()</tt> <tt>createInArgs()</tt>,
 * <tt>createOutArgs()</tt> and <tt>evalModel()</tt> are fairly cleanly
 * written and are appropriate to be studied in order to show how to implement
 * other parallel simulation-constrained problems based on Epetra objects.
 *
 * The mesh for the 2D domain can either be read in as a mesh data file give
 * the files name or can be generated automatically on a square 2D domain.
 *
 * The program <tt>triangle</tt> can be used to generate meshes for arbitary
 * 2D geometries and then <tt>metis</tt> can be used to partition the mesh to
 * multiple domains.  Instructions for how to use <tt>triangle</tt> and
 * <tt>metis</tt> to generate meshes is described ???here???.
 *
 * Instead of reading in a mesh file, a square 2D mesh can be automatically
 * generated given just the length in the <tt>x</tt> and <tt>y</tt> directions
 * and the number of local elements in each direction.  Currently, the square
 * mesh is only partitioned in the <tt>x</tt> direction and therefore will not
 * demonstrate great parallel scalability for large numbers of processors due
 * to excessive amounts of shared boundary between processes.
 *
 * ToDo: Finish Documentation!
 */
class AdvDiffReactOptModel
  : public EpetraExt::ModelEvaluator
  , public Teuchos::VerboseObject<AdvDiffReactOptModel>
{
public:

  /** \brief Constructor. */
  AdvDiffReactOptModel(
    const Teuchos::RefCountPtr<const Epetra_Comm>  &comm
    ,const double                                  beta
    ,const double                                  len_x     // Ignored if meshFile is *not* empty
    ,const double                                  len_y     // Ignored if meshFile is *not* empty
    ,const int                                     local_nx  // Ignored if meshFile is *not* empty
    ,const int                                     local_ny  // Ignored if meshFile is *not* empty
    ,const char                                    meshFile[]
    ,const int                                     np
    ,const double                                  x0
    ,const double                                  p0
    ,const double                                  reactionRate
    ,const bool                                    normalizeBasis
    ,const bool                                    supportDerivatives
    );

  /** \brief . */
  void set_q( Teuchos::RefCountPtr<const Epetra_Vector> const& q );

  /** \brief . */
  Teuchos::RefCountPtr<GLpApp::GLpYUEpetraDataPool> getDataPool();

  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_MultiVector> get_B_bar() const;

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
  Teuchos::RefCountPtr<const Epetra_Vector> get_p_init(int l) const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_x_lower_bounds() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_x_upper_bounds() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_p_lower_bounds(int l) const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_p_upper_bounds(int l) const;
  /** \brief . */
  Teuchos::RefCountPtr<Epetra_Operator> create_W() const;
  /** \brief . */
  Teuchos::RefCountPtr<Epetra_Operator> create_DfDp_op(int l) const;
  /** \brief . */
  InArgs createInArgs() const;
  /** \brief . */
  OutArgs createOutArgs() const;
  /** \brief . */
  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  //@}

private:

  // /////////////////////////////////////
  // Private types

  typedef Teuchos::Array<Teuchos::RefCountPtr<const Epetra_Map> >  RCP_Eptra_Map_Array_t;
  typedef Teuchos::Array<Teuchos::RefCountPtr<Epetra_Vector> >     RCP_Eptra_Vector_Array_t;

  // /////////////////////////////////////
  // Private member data

  static const int Np_         = 2; // Number of axiliary parameters
  static const int p_bndy_idx  = 0; // index for boundary flux parameters
  static const int p_rx_idx    = 1; // index for reaction rate parameter

  bool      supportDerivatives_;

  bool      isInitialized_;

  Teuchos::RefCountPtr<GLpApp::GLpYUEpetraDataPool>   dat_;
  int                                                 np_;
  Teuchos::RefCountPtr<const Epetra_Vector>           q_;

  Teuchos::RefCountPtr<const Epetra_Map>              map_p_bar_;
  Teuchos::RefCountPtr<Epetra_MultiVector>            B_bar_;

  Teuchos::RefCountPtr<const Epetra_Comm>  epetra_comm_;
  Teuchos::RefCountPtr<const Epetra_Map>   map_x_;
  RCP_Eptra_Map_Array_t                    map_p_;
  Teuchos::RefCountPtr<const Epetra_Map>   map_f_;
  Teuchos::RefCountPtr<const Epetra_Map>   map_g_;

  Teuchos::RefCountPtr<Epetra_Vector> x0_;
  Teuchos::RefCountPtr<Epetra_Vector> xL_;
  Teuchos::RefCountPtr<Epetra_Vector> xU_;
  RCP_Eptra_Vector_Array_t            p0_;
  RCP_Eptra_Vector_Array_t            pL_;
  RCP_Eptra_Vector_Array_t            pU_;
  Teuchos::RefCountPtr<Epetra_Vector> gL_;
  Teuchos::RefCountPtr<Epetra_Vector> gU_;

  Teuchos::RefCountPtr<Epetra_CrsGraph>  W_graph_;

};

} // namespace GLpApp

#endif // GLP_APP_ADV_DIFF_REACT_OPT_MODEL_HPP
