#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

panzer::BasisIRLayout::
BasisIRLayout(std::string basis_type, const panzer::IntegrationRule& int_rule) :
  basis_name(basis_type),
  field_basis_name("Basis: " + basis_type),
  field_basis_name_D1("Grad Basis: " + basis_type),
  field_basis_name_D2("D2 Basis: " + basis_type)
{
  basis_data = Teuchos::rcp(new PureBasis(basis_type,int_rule.workset_size,int_rule.topology));

  setup(basis_data->getIntrepidBasis(),int_rule);
}

panzer::BasisIRLayout::
BasisIRLayout(const Teuchos::RCP<const panzer::PureBasis> & b, const panzer::IntegrationRule& int_rule) :
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
      const panzer::IntegrationRule & int_rule)
{
  intrepid_basis = iBasis;

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

int panzer::BasisIRLayout::getCardinality() const
{
  return cardinality;
}

int panzer::BasisIRLayout::integrationRuleDegree() const
{
  return int_rule_degree;
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
      << ", Quad Degree = " << integrationRuleDegree() 
      << ", Quad Points = " << getNumPoints();
}
