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

#ifndef EPETRA_MODEL_EVAL_2D_SIM_HPP
#define EPETRA_MODEL_EVAL_2D_SIM_HPP

#include "EpetraExt_ModelEvaluator.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"

/** \brief Simple example ModelEvaluator subclass for a 2x2 set of
 * parameterized nonlinear equations.
 *
 * The equations modeled are:
 \verbatim

    f[0] =       x[0]      + x[1]*x[1] - p[0];
    f[1] = d * ( x[0]*x[0] - x[1]      - p[1] );

 \endverbatim
 */
class EpetraModelEval2DSim : public EpetraExt::ModelEvaluator {
public:

  /** \brief . */
  EpetraModelEval2DSim(
		const double         d                   = 10.0
		,const double        p0                  = 2.0
		,const double        p1                  = 0.0
		,const double        x00                 = 1.0
		,const double        x01                 = 1.0
    ,const bool          showGetInvalidArg   = false
    );

  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Map> get_x_map() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Map> get_f_map() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_x_init() const;
  /** \brief . */
  Teuchos::RefCountPtr<Epetra_Operator> create_W() const;
  /** \brief . */
  InArgs createInArgs() const;
  /** \brief . */
  OutArgs createOutArgs() const;
  /** \brief . */
  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  //@}

private:

  // /////////////////////////////////////
  // Private member data

  double    d_;
  bool      showGetInvalidArg_;

  bool      isInitialized_;
  
  Teuchos::RefCountPtr<const Epetra_Comm>  epetra_comm_;
  Teuchos::RefCountPtr<const Epetra_Map>   map_x_;
  
  Teuchos::RefCountPtr<Epetra_Vector> x0_;
  Teuchos::RefCountPtr<Epetra_Vector> p_;
  
  Teuchos::RefCountPtr<Epetra_CrsGraph>  W_graph_;
  
};

#endif // EPETRA_MODEL_EVAL_2D_SIM_HPP
