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

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_PointRule.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Panzer_CellTopologyInfo.hpp"

panzer::BasisIRLayout::
BasisIRLayout(std::string basis_type, const int basis_order, const panzer::PointRule& int_rule) :
  basis_name(basis_type),
  field_basis_name("Basis: " + basis_name),
  field_basis_name_D1("Grad Basis: " + basis_name),
  field_basis_name_D2("D2 Basis: " + basis_name)
{
  basis_data = Teuchos::rcp(new PureBasis(basis_type,basis_order,int_rule.workset_size,int_rule.topology));

  setup(basis_data->getIntrepidBasis(),int_rule);
}

panzer::BasisIRLayout::
BasisIRLayout(const Teuchos::RCP<const panzer::PureBasis> & b, const panzer::PointRule& int_rule) :
  basis_name(b->name()),
  field_basis_name(b->fieldName()),
  field_basis_name_D1(b->fieldNameD1()),
  field_basis_name_D2(b->fieldNameD2())
{
  basis_data = b;

  setup(basis_data->getIntrepidBasis(),int_rule);
}

void panzer::BasisIRLayout::
setup(const Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > & iBasis,
      const panzer::PointRule & int_rule)
{
  intrepid_basis = iBasis;

  cardinality = intrepid_basis->getCardinality();
  num_cells = int_rule.dl_vector->dimension(0);
  num_ip = int_rule.dl_vector->dimension(1);
  dimension = int_rule.dl_vector->dimension(2);
  
  using Teuchos::rcp;
  using PHX::MDALayout;
  
  basis_ref = rcp(new MDALayout<BASIS,IP>(cardinality, num_ip));
  
  basis = 
    rcp(new MDALayout<Cell,BASIS,IP>(num_cells, cardinality, num_ip));
  
  basis_grad_ref = 
    rcp(new MDALayout<BASIS,IP,Dim>(cardinality, num_ip, dimension));
  
  basis_grad = rcp(new MDALayout<Cell,BASIS,IP,Dim>(num_cells,
						    cardinality,
						    num_ip,
						    dimension));

  basis_D2_ref =  rcp(new MDALayout<BASIS,IP,Dim,Dim>(cardinality, 
						      num_ip, 
						      dimension, 
						      dimension));
  
  basis_D2 = rcp(new MDALayout<Cell,BASIS,IP,Dim,Dim>(num_cells,
						      cardinality,
						      num_ip,
						      dimension,
						      dimension));

  functional = rcp(new MDALayout<Cell,BASIS>(num_cells, cardinality));

  functional_grad = rcp(new MDALayout<Cell,BASIS,Dim>(num_cells,
						      cardinality,
						      dimension));

  functional_D2 = rcp(new MDALayout<Cell,BASIS,Dim,Dim>(num_cells,
							cardinality,
							dimension,
							dimension));

  const Teuchos::RCP<const shards::CellTopology>& topology = basis_data->getCellTopology();
  cell_topo_info = rcp(new panzer::CellTopologyInfo(num_cells, topology) );
  
}

int panzer::BasisIRLayout::getCardinality() const
{
  return cardinality;
}

int panzer::BasisIRLayout::getNumCells() const
{
  return num_cells;
}

int panzer::BasisIRLayout::getNumPoints() const
{
  return num_ip;
}

int panzer::BasisIRLayout::getDimension() const
{
  return dimension;
}

std::string panzer::BasisIRLayout::name() const
{
  return basis_name;
}

std::string panzer::BasisIRLayout::fieldName() const
{
  return field_basis_name;
}

std::string panzer::BasisIRLayout::fieldNameD1() const
{
  return field_basis_name_D1;
}    
 
std::string panzer::BasisIRLayout::fieldNameD2() const
{
  return field_basis_name_D2;
}    

Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > 
panzer::BasisIRLayout::getIntrepidBasis() const
{
   return intrepid_basis;
}

Teuchos::RCP< const panzer::PureBasis>
panzer::BasisIRLayout::getBasis() const
{
   return basis_data;
}

void panzer::BasisIRLayout::print(std::ostream & os) const
{
   os << "Name = " << name() 
      << ", Dimension = " << getDimension()
      << ", Cells = " << getNumCells() 
      << ", Num Points = " << getNumPoints();
}
