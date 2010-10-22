
#include "Panzer_Basis.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"

#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

panzer::Basis::
Basis(std::string basis_type, const panzer::IntegrationRule& int_rule) :
  basis_name(basis_type),
  field_basis_name("Basis: " + basis_type),
  field_basis_name_D1("Grad Basis: " + basis_type),
  field_basis_name_D2("D2 Basis: " + basis_type)
{
  //Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > >
  intrepid_basis = panzer::createIntrepidBasis<double,Intrepid::FieldContainer<double> >(basis_type, int_rule.spatial_dimension);

  cardinality = intrepid_basis->getCardinality();
  num_cells = int_rule.dl_vector->dimension(0);
  num_ip = int_rule.dl_vector->dimension(1);
  dimension = int_rule.dl_vector->dimension(2);
  int_rule_degree = int_rule.cubature_degree;

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
}

int panzer::Basis::getCardinality() const
{
  return cardinality;
}

int panzer::Basis::integrationRuleDegree() const
{
  return int_rule_degree;
}

int panzer::Basis::getNumCells() const
{
  return num_cells;
}

int panzer::Basis::getNumIntPoints() const
{
  return num_ip;
}

int panzer::Basis::getDimension() const
{
  return dimension;
}

std::string panzer::Basis::name() const
{
  return basis_name;
}

std::string panzer::Basis::fieldName() const
{
  return field_basis_name;
}

std::string panzer::Basis::fieldNameD1() const
{
  return field_basis_name_D1;
}    
 
std::string panzer::Basis::fieldNameD2() const
{
  return field_basis_name_D2;
}    

Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > 
panzer::Basis::getIntrepidBasis() const
{
   return intrepid_basis;
}
