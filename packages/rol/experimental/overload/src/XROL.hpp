// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#pragma once

#include "cxxstd.hpp"

// Teuchos Includes
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_ParameterList.hpp"

// ROL Includes
#include "ROL_Types.hpp"


// Depends only on cxxstd.hpp, 
// Independent of other XROL headers
#include "XROL_Defines.hpp"
#include "XROL_ElementwiseFunction.hpp"
#include "XROL_Exception.hpp"
#include "XROL_Output.hpp"
#include "XROL_VectorTraits.hpp"
#include "XROL_Vector.hpp"
#include "XROL_TypeCheck.hpp"



// Depends on XROL_VectorTraits.hpp and XROL_Vector.hpp
#include "XROL_VectorCheck.hpp"
#include "XROL_StdVector.hpp"
#include "XROL_ArrayVector.hpp"
#include "XROL_Vector_SimOpt.hpp"
#include "XROL_PartitionedVector.hpp"
#include "XROL_Objective.hpp"


// Depends on XROL_Objective
#include "XROL_ObjectiveImpl.hpp"
#include "XROL_ObjectiveCheck.hpp"
#include "XROL_ObjectiveVectors.hpp"

#include "XROL_Objective_ExtendedInterface.hpp"

//#include "XROL_TestObjective.hpp"
