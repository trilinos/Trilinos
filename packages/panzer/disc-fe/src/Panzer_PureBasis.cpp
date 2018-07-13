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
#include "Panzer_BasisDescriptor.hpp"
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

panzer::PureBasis::
PureBasis(const panzer::BasisDescriptor & description,
          const Teuchos::RCP<const shards::CellTopology> & cell_topology,
          const int num_cells):
  topology_(cell_topology),
  num_cells_(num_cells)
{
  initialize(description.getType(), description.getOrder());
}

void panzer::PureBasis::initialize(const std::string & in_basis_type,const int in_basis_order)
{
  // Support for deprecated basis descriptions
  std::string basis_type = in_basis_type;
  int basis_order = in_basis_order;

  if (basis_type=="Q1" || basis_type=="T1") {
    basis_type = "HGrad";
    basis_order = 1;
  }
  else if (basis_type == "Q2" || basis_type=="T2") {
    basis_type = "HGrad";
    basis_order = 2;
  }
  else if (basis_type == "TEdge1" || basis_type=="QEdge1") {
    basis_type = "HCurl";
    basis_order = 1;
  }
  else if(basis_type == "Const") {
    basis_type = "Const";
    basis_order = 0;
  }
  // End deprecated basis support

  intrepid_basis_ = panzer::createIntrepid2Basis<PHX::Device::execution_space,double,double>(basis_type, basis_order, *topology_);

  basis_type_ = basis_type;

  std::ostringstream os;
  os << basis_type_ << ":" << basis_order;
  basis_name_ = os.str();

  field_basis_name_ = "Basis: " + basis_name_;
  field_basis_name_D1_ = "Grad Basis: " + basis_name_;
  field_basis_name_D2_ = "D2 Basis: " + basis_name_;

  if(  basis_type_ == "HGrad")
    element_space_ = HGRAD;
  else if(basis_type_=="HCurl")
    element_space_ = HCURL;
  else if(basis_type_=="HDiv")
    element_space_ = HDIV;
  else if(basis_type_=="Const")
    element_space_ = CONST;
  else { TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
				    "PureBasis::initializeIntrospection - Invalid basis name \"" 
				    << basis_type_ << "\""); }
  
  switch(getElementSpace()) {
  case CONST:
     basis_rank_ = 0;
     break;
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

  cell_data = rcp(new MDALayout<Cell>(numCells()));

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
  
  local_mat_layout = Teuchos::rcp(new PHX::MDALayout<panzer::Cell, panzer::BASIS, panzer::BASIS>(
                     this->numCells(), this->cardinality(), this->cardinality()));  
  
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

Teuchos::RCP< Intrepid2::Basis<PHX::Device::execution_space,double,double> > 
panzer::PureBasis::getIntrepid2Basis() const
{
   return intrepid_basis_;
}

bool 
panzer::PureBasis::supportsBasisCoordinates() const
{
  // typedef Kokkos::DynRankView<double,PHX::Device> Array;
  // Teuchos::RCP<const Intrepid2::DofCoordsInterface<Array> > coord_interface 
  //     = Teuchos::rcp_dynamic_cast<const Intrepid2::DofCoordsInterface<Array> >(getIntrepid2Basis());

  // return !Teuchos::is_null(coord_interface);

  return true;
}
