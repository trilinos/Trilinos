// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_BASIS_HPP
#define PANZER_BASIS_HPP

#include <string>
#include "Teuchos_RCP.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Intrepid2_Basis.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_CellTopologyInfo.hpp"

namespace panzer {

  class PointRule;
  class CellTopologyInfo;
  class BasisIRLayout;

  //! Nonmember constructor
  Teuchos::RCP<panzer::BasisIRLayout> 
  basisIRLayout(std::string basis_type, const int basis_order, const PointRule& pt_rule);

  //! Nonmember constructor
  Teuchos::RCP<panzer::BasisIRLayout>
  basisIRLayout(const Teuchos::RCP<const PureBasis> & b, const PointRule& pt_rule);

  class BasisIRLayout { 

  public:
    
    BasisIRLayout(std::string basis_type, const int basis_order, const PointRule& int_rule);
    BasisIRLayout(const Teuchos::RCP<const PureBasis> & b, const PointRule& int_rule);

    void setup(const panzer::PointRule & int_rule);

    int cardinality() const;
    
    int numCells() const;
    
    int numPoints() const;
    
    int dimension() const;

    //! Unique key for workset indexing composed of basis name and point rule name
    std::string name() const;

    std::string fieldName() const;
    
    std::string fieldNameD1() const;
    
    std::string fieldNameD2() const;

    Teuchos::RCP< Intrepid2::Basis<PHX::Device::execution_space,double,double> > 
    getIntrepid2Basis() const;

    Teuchos::RCP<const PureBasis> getBasis() const;

    void print(std::ostream & os) const;

    Teuchos::RCP<const CellTopologyInfo> getCellTopologyInfo() const
    { return cell_topo_info; }
    

  public:
    
    //! <BASIS,IP>
    Teuchos::RCP<PHX::DataLayout> basis_ref;
    //! <Cell,BASIS,IP>
    Teuchos::RCP<PHX::DataLayout> basis;
    //! <BASIS,IP,Dim>
    Teuchos::RCP<PHX::DataLayout> basis_grad_ref;
    //! <Cell,BASIS,IP,Dim>
    Teuchos::RCP<PHX::DataLayout> basis_grad;
    //! <BASIS,IP,Dim,Dim>
    Teuchos::RCP<PHX::DataLayout> basis_D2_ref;
    //! <Cell,BASIS,IP,Dim,Dim>
    Teuchos::RCP<PHX::DataLayout> basis_D2;

    //! <Cell,Basis>
    Teuchos::RCP<PHX::DataLayout> functional;
    //! <Cell,Basis,Dim>
    Teuchos::RCP<PHX::DataLayout> functional_grad;
    //! <Cell,Basis,Dim,Dim>
    Teuchos::RCP<PHX::DataLayout> functional_D2;

  private:
    std::string basis_name_;
    int num_cells_;
    int num_points_;
    int dimension_;
    
    Teuchos::RCP<const PureBasis> basis_data_;
    
    Teuchos::RCP<const CellTopologyInfo> cell_topo_info;
  };

  typedef std::pair<std::string,Teuchos::RCP<panzer::BasisIRLayout> > StrBasisPair;

  //! Simple binary comparison class to help with sorting
  struct StrBasisComp {
    bool operator() (const StrBasisPair & lhs, const StrBasisPair & rhs) const
    {return lhs.first<rhs.first;}
  };

}

#endif
