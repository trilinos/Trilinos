// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_STK_SetupLOWSFactory_impl_hpp__
#define __Panzer_STK_SetupLOWSFactory_impl_hpp__

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_ParameterListCallback.hpp"
#include "Panzer_STK_ParameterListCallbackBlocked.hpp"

#include "Teuchos_AbstractFactoryStd.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "EpetraExt_VectorOut.h"

#include "ml_rbm.h"

#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"

#ifdef PANZER_HAVE_TEKO
#include "Teko_StratimikosFactory.hpp"
#endif

#ifdef PANZER_HAVE_MUELU
#include <Thyra_MueLuPreconditionerFactory.hpp>
#include <Thyra_MueLuRefMaxwellPreconditionerFactory.hpp>
#include "Stratimikos_MueLuHelpers.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#endif

#ifdef PANZER_HAVE_IFPACK2
#include <Thyra_Ifpack2PreconditionerFactory.hpp>
#endif

namespace panzer_stk {

}

#endif
