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

#include "EpetraModelEval2DSim.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_CrsMatrix.h"

EpetraModelEval2DSim::EpetraModelEval2DSim(
  const double         d
  ,const double        p0
  ,const double        p1
  ,const double        x00
  ,const double        x01
  ,const bool          showGetInvalidArg
  )
  :d_(d),showGetInvalidArg_(showGetInvalidArg)
{
  using Teuchos::rcp;

	epetra_comm_ = rcp(new Epetra_SerialComm());

  const int nx = 2;

  map_x_ = rcp(new Epetra_Map(nx,0,*epetra_comm_));
  x0_ = rcp(new Epetra_Vector(*map_x_));  (*x0_)[0] = x00; (*x0_)[1] = x01;
  p_  = rcp(new Epetra_Vector(*map_x_));  (*p_)[0] = p0; (*p_)[1] = p1;

  // Initialize the graph for W CrsMatrix object
  W_graph_ = rcp(new Epetra_CrsGraph(::Copy,*map_x_,nx));
  {
    int indices[nx] = { 0, 1 };
    for( int i = 0; i < nx; ++i )
      W_graph_->InsertGlobalIndices(i,nx,indices);
  }
  W_graph_->FillComplete();

  isInitialized_ = true;

}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RefCountPtr<const Epetra_Map>
EpetraModelEval2DSim::get_x_map() const
{
  return map_x_;
}

Teuchos::RefCountPtr<const Epetra_Map>
EpetraModelEval2DSim::get_f_map() const
{
  return map_x_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
EpetraModelEval2DSim::get_x_init() const
{
  return x0_;
}

Teuchos::RefCountPtr<Epetra_Operator>
EpetraModelEval2DSim::create_W() const
{
  return Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*W_graph_));
}

EpetraExt::ModelEvaluator::InArgs
EpetraModelEval2DSim::createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(IN_ARG_x,true);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
EpetraModelEval2DSim::createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W,true);
  outArgs.set_W_properties(
    DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_FULL
      ,true // supportsAdjoint
      )
    );
  return outArgs;
}

void EpetraModelEval2DSim::evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const
{
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;
  //
  // Get the input arguments
  //
  const Epetra_Vector &x = *inArgs.get_x();
  //
  // Get the output arguments
  //
  Epetra_Vector       *f_out = outArgs.get_f().get();
  Epetra_Operator     *W_out = outArgs.get_W().get();
  if(showGetInvalidArg_) {
    Epetra_Vector *g_out = outArgs.get_g(0).get();
  }
  //
  // Compute the functions
  //
  const Epetra_Vector &p = *p_;
  if(f_out) {
    Epetra_Vector &f = *f_out;
    f[0] =        x[0]      + x[1]*x[1] - p[0];
    f[1] = d_ * ( x[0]*x[0] - x[1]      - p[1] );
  }
  if(W_out) {
    Epetra_CrsMatrix &DfDx = dyn_cast<Epetra_CrsMatrix>(*W_out);
    DfDx.PutScalar(0.0);
    //
    // Fill W = DfDx
    //
    // W = DfDx = [      1.0,  2*x[1] ]
    //            [ 2*d*x[0],     -d  ]
    //
    double values[2];
    int indexes[2];
    // Row [0]
    values[0] = 1.0;           indexes[0] = 0;
    values[1] = 2.0*x[1];      indexes[1] = 1;
    DfDx.SumIntoGlobalValues( 0, 2, values, indexes );
    // Row [1]
    values[0] = 2.0*d_*x[0];   indexes[0] = 0;
    values[1] = -d_;           indexes[1] = 1;
    DfDx.SumIntoGlobalValues( 1, 2, values, indexes );
  }
}
