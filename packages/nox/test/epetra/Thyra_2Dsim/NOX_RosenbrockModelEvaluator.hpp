/*
//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
//@HEADER
*/

#ifndef NOX_THYRA_ROSENBROCK_MODEL_EVALUATOR_HPP
#define NOX_THYRA_ROSENBROCK_MODEL_EVALUATOR_HPP

#include "EpetraExt_ModelEvaluator.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"

namespace NOX {

/** \brief ModelEvaluator for the Roesenbrock problem

   This example is the "%Rosenbrock function" from Jorge J. More',
   Burton S. Garbow, and Kenneth E. Hillstrom, Testing Unconstrained
   Optimization Software, ACM TOMS, Vol. 7, No. 1, March 1981,
   pp. 14-41.

   It comes originally from H. H. %Rosenbrock, An Automatic Method for
   Finding the Greatest or Least Value of a Function, J. Comput.
   3(1960):175-184.

   The function is defined as
   \f[
   F(x) = \left[
   \begin{array}{c}
   10 (x[2] - x[1]^2) \\
   1 - x[1]
   \end{array}
   \right]
   \f]

   The initial guess is given by
   \f[
   x_0 = \left[
   \begin{array}{c}
   -1.2\\
   1
   \end{array}
   \right]
   \f]

   The solution is
   \f[
   x_* = \left[
   \begin{array}{c}
   1\\
   1
   \end{array}
   \right]
   \f]
 */
class RosenbrockModelEvaluator : public EpetraExt::ModelEvaluator {
public:

  RosenbrockModelEvaluator(const Teuchos::RCP<const Epetra_Comm>& comm);

  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  Teuchos::RCP<const Epetra_Map> get_x_map() const;

  Teuchos::RCP<const Epetra_Map> get_f_map() const;

  Teuchos::RCP<const Epetra_Vector> get_x_init() const;

  Teuchos::RCP<Epetra_Operator> create_W() const;

  InArgs createInArgs() const;

  OutArgs createOutArgs() const;

  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  //@}

  Teuchos::RCP<const Epetra_Vector> get_analytic_solution() const;

private:

  double    d_;

  bool      isInitialized_;
  
  Teuchos::RCP<const Epetra_Comm>  epetra_comm_;
  Teuchos::RCP<const Epetra_Map>   map_x_;
  
  Teuchos::RCP<Epetra_Vector> x0_;
  Teuchos::RCP<Epetra_Vector> x_analytic_;
  
  Teuchos::RCP<Epetra_CrsGraph>  W_graph_;
  
};

}

#endif
