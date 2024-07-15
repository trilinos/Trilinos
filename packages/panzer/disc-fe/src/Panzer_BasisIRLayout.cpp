// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_PointRule.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"
#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Panzer_CellTopologyInfo.hpp"


// ***********************************************************************
// Nonmember ctors

Teuchos::RCP<panzer::BasisIRLayout> 
panzer::basisIRLayout(std::string basis_type, const int basis_order, const panzer::PointRule& pt_rule)
{
  return Teuchos::rcp(new panzer::BasisIRLayout(basis_type,basis_order,pt_rule),true);
}

Teuchos::RCP<panzer::BasisIRLayout> 
panzer::basisIRLayout(const Teuchos::RCP<const PureBasis> & b, const PointRule& pt_rule)
{
  return Teuchos::rcp(new panzer::BasisIRLayout(b,pt_rule),true);
}


// ***********************************************************************
// Class implementation

panzer::BasisIRLayout::
BasisIRLayout(std::string basis_type, const int basis_order, const panzer::PointRule& point_rule)
{
  basis_data_ = Teuchos::rcp(new PureBasis(basis_type,basis_order,point_rule.workset_size,point_rule.topology));

  setup(point_rule);
}

panzer::BasisIRLayout::
BasisIRLayout(const Teuchos::RCP<const panzer::PureBasis> & b, const panzer::PointRule& point_rule) :
  basis_data_(b)
{
  setup(point_rule);
}

void panzer::BasisIRLayout::
setup(const panzer::PointRule & point_rule)
{
  basis_name_ = basis_data_->name() + ":" + point_rule.getName();
  num_cells_ = point_rule.dl_vector->extent(0);
  num_points_ = point_rule.dl_vector->extent(1);
  dimension_ = point_rule.dl_vector->extent(2);
  
  using Teuchos::rcp;
  using PHX::MDALayout;
  
  basis_ref = rcp(new MDALayout<BASIS,IP>(cardinality(), numPoints()));
  
  basis = 
    rcp(new MDALayout<Cell,BASIS,IP>(numCells(), cardinality(), numPoints()));
  
  basis_grad_ref = 
    rcp(new MDALayout<BASIS,IP,Dim>(cardinality(), numPoints(), dimension()));
  
  basis_grad = rcp(new MDALayout<Cell,BASIS,IP,Dim>(numCells(),
						    cardinality(),
						    numPoints(),
						    dimension()));

  basis_D2_ref =  rcp(new MDALayout<BASIS,IP,Dim,Dim>(cardinality(), 
						      numPoints(), 
						      dimension(), 
						      dimension()));
  
  basis_D2 = rcp(new MDALayout<Cell,BASIS,IP,Dim,Dim>(numCells(),
						      cardinality(),
						      numPoints(),
						      dimension(),
						      dimension()));

  functional = rcp(new MDALayout<Cell,BASIS>(numCells(), cardinality()));

  functional_grad = rcp(new MDALayout<Cell,BASIS,Dim>(numCells(),
						      cardinality(),
						      dimension()));

  functional_D2 = rcp(new MDALayout<Cell,BASIS,Dim,Dim>(numCells(),
							cardinality(),
							dimension(),
							dimension()));

  const Teuchos::RCP<const shards::CellTopology>& topology = basis_data_->getCellTopology();
  cell_topo_info = rcp(new panzer::CellTopologyInfo(numCells(), topology) );
  
}

int panzer::BasisIRLayout::cardinality() const
{
  return basis_data_->cardinality();
}

int panzer::BasisIRLayout::numCells() const
{
  return num_cells_;
}

int panzer::BasisIRLayout::numPoints() const
{
  return num_points_;
}

int panzer::BasisIRLayout::dimension() const
{
  return dimension_;
}

std::string panzer::BasisIRLayout::name() const
{
  return basis_name_;
}

std::string panzer::BasisIRLayout::fieldName() const
{
  return basis_data_->fieldName();
}

std::string panzer::BasisIRLayout::fieldNameD1() const
{
  return basis_data_->fieldNameD1();
}    
 
std::string panzer::BasisIRLayout::fieldNameD2() const
{
  return basis_data_->fieldNameD2();
}    

Teuchos::RCP< Intrepid2::Basis<PHX::Device::execution_space,double,double> > 
panzer::BasisIRLayout::getIntrepid2Basis() const
{
   return basis_data_->getIntrepid2Basis();
}

Teuchos::RCP< const panzer::PureBasis>
panzer::BasisIRLayout::getBasis() const
{
   return basis_data_;
}

void panzer::BasisIRLayout::print(std::ostream & os) const
{
   os << "Name = " << name() 
      << ", Dimension = " << dimension()
      << ", Cells = " << numCells() 
      << ", Num Points = " << numPoints();
}
