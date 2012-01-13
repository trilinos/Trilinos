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
  using Teuchos::rcp;
  using PHX::MDALayout;

  dimension = cellTopo->getDimension();

  intrepid_basis = panzer::createIntrepidBasis<double,Intrepid::FieldContainer<double> >(basis_type, dimension, cellTopo);

  cardinality = intrepid_basis->getCardinality();
  num_cells = numCells;

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

panzer::PureBasis::
PureBasis(const std::string & basis_type,const CellData & cell_data) :
  basis_name(basis_type),
  field_basis_name("Basis: " + basis_type),
  field_basis_name_D1("Grad Basis: " + basis_type),
  field_basis_name_D2("D2 Basis: " + basis_type)
{
  using Teuchos::rcp;
  using PHX::MDALayout;

  dimension = cell_data.baseCellDimension();

  intrepid_basis = panzer::createIntrepidBasis<double,Intrepid::FieldContainer<double> >(basis_type, dimension, cell_data.getCellTopology());

  cardinality = intrepid_basis->getCardinality();
  num_cells = cell_data.numCells();

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
