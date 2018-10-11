// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
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
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef PYTRILINOS_ANASAZI_HEADERS_
#define PYTRILINOS_ANASAZI_HEADERS_

// Configuration
#include "PyTrilinos_config.h"
#include "AnasaziConfigDefs.hpp"

// Anasazi include files
#ifdef HAVE_PYTRILINOS_EPETRA
#include "Anasaziepetra_DLLExportMacro.h"
#endif
#include "AnasaziVersion.cpp"
#include "AnasaziTypes.hpp"
#include "AnasaziOutputManager.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziMultiVec.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziOperator.hpp"
#include "AnasaziEigenproblem.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziStatusTest.hpp"
#include "AnasaziStatusTestCombo.hpp"
#include "AnasaziStatusTestMaxIters.hpp"
#include "AnasaziStatusTestOutput.hpp"
#include "AnasaziStatusTestResNorm.hpp"
#include "AnasaziOrthoManager.hpp"
#include "AnasaziMatOrthoManager.hpp"
#include "AnasaziBasicOrthoManager.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziEigensolverDecl.hpp"
#include "AnasaziSolverManager.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#ifdef HAVE_PYTRILINOS_EPETRA
#include "AnasaziEpetraAdapter.hpp"
#endif

#endif
