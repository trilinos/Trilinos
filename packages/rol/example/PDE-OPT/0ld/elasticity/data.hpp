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

#ifndef ROL_PDEOPT_ELASTICITY_DATA_H
#define ROL_PDEOPT_ELASTICITY_DATA_H

#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "../TOOLS/elasticity.hpp"

template<class Real>
class ElasticityData : public Elasticity <Real> {

private:

public:

ElasticityData(const ROL::Ptr<const Teuchos::Comm<int> > &comm,
             const Teuchos::RCP<Teuchos::ParameterList> &parlist,
             const ROL::Ptr<std::ostream> &outStream) 
{
	this->Initialize (comm, parlist, outStream);
	this->SetParallelStructure();
	this->SetUpLocalIntrepidArrays();
	this->ComputeLocalSystemMats();
	this->AssembleSystemGraph();
	this->AssembleSystemMats();
	//Setup DBC information, do not specify any bc sides, use coordinates to determine the BC instead
	std::vector<int> dbc_side {};
	this->SetUpMyDBCInfo(true, dbc_side);
	//Setup all loads 
    	this->process_loading_information(parlist);
	this->ComputeLocalForceVec();
	this->AssembleRHSVector();
	//
	this->EnforceDBC();
	this->ConstructSolvers();
}

}; // class ElasticityData

#endif
