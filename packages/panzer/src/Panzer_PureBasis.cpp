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

#include "Panzer_PureBasis.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"
#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include <sstream>

panzer::PureBasis::
PureBasis(const std::string & basis_type,
	  const int basis_order,
	  const int num_cells,
	  const Teuchos::RCP<const shards::CellTopology> & cell_topology) :
  topology_(cell_topology),
  num_cells_(num_cells)
{
  initialize(basis_type,basis_order);
}

panzer::PureBasis::
PureBasis(const std::string & basis_type,const int basis_order,const CellData & cell_data) :
  topology_(cell_data.getCellTopology()),
  num_cells_(cell_data.numCells())
{
  initialize(basis_type,basis_order);
}

void panzer::PureBasis::initialize(const std::string & basis_type,const int basis_order)
{
  intrepid_basis_ = panzer::createIntrepidBasis<double,Intrepid::FieldContainer<double> >(basis_type, basis_order, topology_);

  basis_type_ = basis_type;

  std::ostringstream os;
  os << basis_type_ << ":" << basis_order;
  basis_name_ = os.str();

  // For deprecated basis descriptions, we have to patch the names by
  // not tacking on the basis order (this knowledge is already
  // embedded in the name)
  if (basis_type_ == "Q1" ||
      basis_type_ == "Q2" ||
      basis_type_ == "T1" ||
      basis_type_ == "T2" ||
      basis_type_ == "TEdge1" ||
      basis_type_ == "QEdge1")
    basis_name_ = basis_type_;

  field_basis_name_ = "Basis: " + basis_name_;
  field_basis_name_D1_ = "Grad Basis: " + basis_name_;
  field_basis_name_D2_ = "D2 Basis: " + basis_name_;

  
  if(  basis_type_ == "HGrad" || basis_name_=="Q1" || basis_name_=="Q2" || basis_name_=="T1" || basis_name_=="T2")
    element_space_ = HGRAD;
  else if(basis_type_=="HCurl" || basis_name_=="TEdge1" || basis_name_=="QEdge1")
    element_space_ = HCURL;
  else if(basis_type_=="HDiv")
    element_space_ = HDIV;
  else { TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
				    "PureBasis::initializeIntrospection - Invalid basis name \"" 
				    << basis_type_ << "\""); }
  
  switch(getElementSpace()) {
  case HGRAD:
     basis_rank_ = 0;
     break;
  case HCURL:
     basis_rank_ = 1;
     break;
  case HDIV:
     basis_rank_ = 1;
     break;
  default:
     TEUCHOS_ASSERT(false);
     break;
  };

  using Teuchos::rcp;
  using PHX::MDALayout;

  functional = rcp(new MDALayout<Cell,BASIS>(numCells(), cardinality()));

  functional_grad = rcp(new MDALayout<Cell,BASIS,Dim>(numCells(),
						      cardinality(),
						      dimension()));

  coordinates = rcp(new MDALayout<Cell,BASIS,Dim>(numCells(),
		   			          cardinality(),
						  dimension()));

  functional_D2 = rcp(new MDALayout<Cell,BASIS,Dim,Dim>(numCells(),
							cardinality(),
							dimension(),
							dimension()));
}

int panzer::PureBasis::cardinality() const
{
  return intrepid_basis_->getCardinality();
}

int panzer::PureBasis::numCells() const
{
  return num_cells_;
}

int panzer::PureBasis::dimension() const
{
  return topology_->getDimension();
}

std::string panzer::PureBasis::type() const
{
  return basis_type_;
}

int panzer::PureBasis::order() const
{
  return intrepid_basis_->getDegree();
}

std::string panzer::PureBasis::name() const
{
  return basis_name_;
}

std::string panzer::PureBasis::fieldName() const
{
  return field_basis_name_;
}

std::string panzer::PureBasis::fieldNameD1() const
{
  return field_basis_name_D1_;
}    
 
std::string panzer::PureBasis::fieldNameD2() const
{
  return field_basis_name_D2_;
}    

Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > 
panzer::PureBasis::getIntrepidBasis() const
{
   return intrepid_basis_;
}
