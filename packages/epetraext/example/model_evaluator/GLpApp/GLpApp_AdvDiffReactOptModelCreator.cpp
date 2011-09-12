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

#include "GLpApp_AdvDiffReactOptModelCreator.hpp"

namespace GLpApp {

AdvDiffReactOptModelCreator::AdvDiffReactOptModelCreator()
  :len_x_(1.0)
  ,len_y_(1.0)
  ,local_nx_(3)
  ,local_ny_(4)
  ,geomFileBase_("")
  ,np_(1)
  ,normalizeBasis_(false)
  ,beta_(0.0)
  ,reactionRate_(1.0)
  ,x0_(0.0)
  ,p0_(1.0)
  ,supportDerivatives_(true)
{}

void AdvDiffReactOptModelCreator::setupCLP(
  Teuchos::CommandLineProcessor *clp
  )
{
  clp->setOption( "len-x", &len_x_, "Mesh dimension in the x direction (Overridden by --geom-file-base)." );
  clp->setOption( "len-y", &len_y_, "Mesh dimension in the y direction (Overridden by --geom-file-base)." );
  clp->setOption( "local-nx", &local_nx_, "Number of local discretization segments in the x direction (Overridden by --geom-file-base)." );
  clp->setOption( "local-ny", &local_ny_, "Number of local discretization segments in the y direction (Overridden by --geom-file-base)." );
  clp->setOption( "geom-file-base", &geomFileBase_, "Base name of geometry file to read the mesh from." );
  clp->setOption( "np", &np_, "The number of optimization parameters p (If < 0 then all of boundary is used)" );
  clp->setOption( "normalize-basis", "no-normalize-basis", &normalizeBasis_, "Normalize the basis for the parameters p or not." );
  clp->setOption( "beta", &beta_, "Regularization." );
  clp->setOption( "reaction-rate", &reactionRate_, "The rate of the reaction" );
  clp->setOption( "x0", &x0_, "Initial guess for the state." );
  clp->setOption( "p0", &p0_, "Initial guess or nonminal value for optimization parameters." );
  clp->setOption( "support-derivatives","no-support-derivatives",&supportDerivatives_,"Support derivatives or not." );
}

Teuchos::RefCountPtr<AdvDiffReactOptModel>
AdvDiffReactOptModelCreator::createModel(
  const Teuchos::RefCountPtr<const Epetra_Comm>     &comm
  ,std::ostream                                     *out
  ) const
{
  return Teuchos::rcp(
    new GLpApp::AdvDiffReactOptModel(
      comm,beta_,len_x_,len_y_,local_nx_,local_ny_,geomFileBase_.c_str()
      ,np_,x0_,p0_,reactionRate_,normalizeBasis_,supportDerivatives_
      )
    );
}

} // namespace GLpApp
