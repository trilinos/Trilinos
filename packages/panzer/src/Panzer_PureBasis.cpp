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
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

panzer::PureBasis::
PureBasis(const std::string & basis_type,int numCells,const Teuchos::RCP<const shards::CellTopology> & cellTopo) :
  basis_name(basis_type),
  field_basis_name("Basis: " + basis_type),
  field_basis_name_D1("Grad Basis: " + basis_type),
  field_basis_name_D2("D2 Basis: " + basis_type)
{
  dimension = cellTopo->getDimension();

  topology = cellTopo;
  intrepid_basis = panzer::createIntrepidBasis<double,Intrepid::FieldContainer<double> >(basis_type, dimension, cellTopo);

  cardinality = intrepid_basis->getCardinality();
  num_cells = numCells;

  // is this a HGRAD, HCURL, etc basis
  initializeIntrospection(basis_type);

  initializeDataLayouts();
}

panzer::PureBasis::
PureBasis(const std::string & basis_type,const CellData & cell_data) :
  basis_name(basis_type),
  field_basis_name("Basis: " + basis_type),
  field_basis_name_D1("Grad Basis: " + basis_type),
  field_basis_name_D2("D2 Basis: " + basis_type)
{
  dimension = cell_data.baseCellDimension();

  topology = cell_data.getCellTopology();
  intrepid_basis = panzer::createIntrepidBasis<double,Intrepid::FieldContainer<double> >(basis_type, dimension, cell_data.getCellTopology());

  cardinality = intrepid_basis->getCardinality();
  num_cells = cell_data.numCells();

  // is this a HGRAD, HCURL, etc basis
  initializeIntrospection(basis_type);

  initializeDataLayouts();
}

int panzer::PureBasis::getCardinality() const
{
  return cardinality;
}

int panzer::PureBasis::getNumCells() const
{
  return num_cells;
}

int panzer::PureBasis::getDimension() const
{
  return dimension;
}

std::string panzer::PureBasis::name() const
{
  return basis_name;
}

std::string panzer::PureBasis::fieldName() const
{
  return field_basis_name;
}

std::string panzer::PureBasis::fieldNameD1() const
{
  return field_basis_name_D1;
}    
 
std::string panzer::PureBasis::fieldNameD2() const
{
  return field_basis_name_D2;
}    

Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > 
panzer::PureBasis::getIntrepidBasis() const
{
   return intrepid_basis;
}

void panzer::PureBasis::initializeIntrospection(const std::string & name)
{
   if(  name=="Q1" || name=="Q2"
     || name=="T1" || name=="T2")
      elementSpace = HGRAD;
   else if(name=="TEdge1")
      elementSpace = HCURL;
   else if(name=="QEdge1")
      elementSpace = HCURL;
   else { TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                                     "PureBasis::initializeIntrospection - Invalid basis name \"" 
                                     << name << "\""); }

  switch(getElementSpace()) {
  case HGRAD:
     basisRank = 0;
     break;
  case HCURL:
     basisRank = 1;
     break;
  case HDIV:
     basisRank = 1;
     break;
  default:
     TEUCHOS_ASSERT(false);
     break;
  };

}

void panzer::PureBasis::initializeDataLayouts()
{
  using Teuchos::rcp;
  using PHX::MDALayout;

  functional = rcp(new MDALayout<Cell,BASIS>(num_cells, cardinality));

  functional_grad = rcp(new MDALayout<Cell,BASIS,Dim>(num_cells,
						      cardinality,
						      dimension));

  coordinates = rcp(new MDALayout<Cell,BASIS,Dim>(num_cells,
		   			          cardinality,
						  dimension));

  functional_D2 = rcp(new MDALayout<Cell,BASIS,Dim,Dim>(num_cells,
							cardinality,
							dimension,
							dimension));
}
