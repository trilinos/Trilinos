// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PDEOPT_ELASTICITY_DATA_H
#define ROL_PDEOPT_ELASTICITY_DATA_H

#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "../TOOLS/elasticity.hpp"

template<class Real>
class ElasticityData : public Elasticity <Real> {

private:
  using GO = typename Tpetra::Map<>::global_ordinal_type;

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
	std::vector<GO> dbc_side {};
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
