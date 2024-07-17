// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_PureBasis_HPP
#define PANZER_PureBasis_HPP

#include <string>
#include "Teuchos_RCP.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Kokkos_DynRankView.hpp"
#include "Intrepid2_Basis.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"

namespace panzer {

  class CellData;
  class BasisDescriptor;

  //! Description and data layouts associated with a particular basis
  class PureBasis { 

  public:
    typedef enum { HGRAD=0, HCURL=1, HDIV=2, HVOL=3, CONST=4 } EElementSpace;
    
    /** Build a basis given a type, order and CellData object
      \param[in] basis_type String name that describes the type of basis ("HGrad", "HDiv", "HCurl", or "HVol")
      \param[in] basis_order Order of the basis
      \param[in] cell_data Description of the basis
    */
    PureBasis(const std::string & basis_type,const int basis_order,const CellData & cell_data);

    /** Build a basis given a type, order, number of cells (for data layouts) and shards topology
      \param[in] basis_type String name that describes the type of basis ("HGrad", "HDiv", "HCurl", or "HVol")
      \param[in] basis_order Order of the basis
      \param[in] num_cells Number of cells used in the data layouts for this basis
      \param[in] cell_topo A shards topology description
    */
    PureBasis(const std::string & basis_type,const int basis_order,const int num_cells,const Teuchos::RCP<const shards::CellTopology> & cell_topo);


    /** Build a basis given a type, order, number of cells (for data layouts) and shards topology
      \param[in] description Description of basis
      \param[in] cell_topo A shards topology description
      \param[in] num_cells Number of cells used in the data layouts for this basis
    */
    PureBasis(const panzer::BasisDescriptor & description, const Teuchos::RCP<const shards::CellTopology> & cell_topology, const int num_cells);

    //! Returns the number of basis coefficients
    int cardinality() const;

    //! Returns the number of cells in the data layouts
    int numCells() const;
    
    //! Returns the dimension of the basis from the topology
    int dimension() const;

    //! Returns the basis type
    std::string type() const;
    
    //! Returns the polynomial order of the basis
    int order() const;
    
    //! A unique key that is the combination of the basis type and basis order
    std::string name() const;
    
    std::string fieldName() const;
    
    std::string fieldNameD1() const;
    
    std::string fieldNameD2() const;

    Teuchos::RCP< Intrepid2::Basis<PHX::Device::execution_space,double,double> > 
    getIntrepid2Basis() const;

    template <typename ExecutionSpace,typename OutputValueType, typename PointValueType>
    Teuchos::RCP< Intrepid2::Basis<ExecutionSpace,OutputValueType,PointValueType> > 
    getIntrepid2Basis() const
    { return panzer::createIntrepid2Basis<ExecutionSpace,OutputValueType,PointValueType>(type(), order(), *(getCellTopology())); }

    EElementSpace getElementSpace() const
    { return element_space_; }

    bool requiresOrientations() const
    { 
      return intrepid_basis_->requireOrientation(); 
    }

    bool supportsGrad() const
    { return getElementSpace()==HGRAD; }

    bool supportsCurl() const
    { return getElementSpace()==HCURL; }

    bool supportsDiv() const
    { return getElementSpace()==HDIV; }

    bool isVectorBasis() const
    { return getElementSpace()==HCURL || getElementSpace()==HDIV; }

    bool isScalarBasis() const
    { return getElementSpace()==HGRAD || getElementSpace()==CONST || getElementSpace()==HVOL; }

    int getBasisRank() const
    { return basis_rank_; }

    bool supportsBasisCoordinates() const;

    Teuchos::RCP<const shards::CellTopology> getCellTopology() const
    { return topology_; }

  public:
    //! <Cell> 
    Teuchos::RCP<PHX::DataLayout> cell_data;
    //! <Cell,Basis> or <Cell,Basis>
    Teuchos::RCP<PHX::DataLayout> functional;
    //! <Cell,Basis,Dim>
    Teuchos::RCP<PHX::DataLayout> functional_grad;
    //! <Cell,Basis,Dim>
    Teuchos::RCP<PHX::DataLayout> functional_D2;
    //! <Cell,Basis,Dim>
    Teuchos::RCP<PHX::DataLayout> coordinates;
    //! <Cell,Basis,Basis>
    Teuchos::RCP<PHX::DataLayout> local_mat_layout;

  private:
    
    //! Initialize the basis object
    void initialize(const std::string & basis_type,const int basis_order);

  private:

    Teuchos::RCP<const shards::CellTopology> topology_;
    Teuchos::RCP< Intrepid2::Basis<PHX::Device::execution_space,double,double> > intrepid_basis_;

    std::string basis_type_;
    std::string basis_name_;
    std::string field_basis_name_;
    std::string field_basis_name_D1_;
    std::string field_basis_name_D2_;

    int num_cells_;

    EElementSpace element_space_;
    int basis_rank_;
  };

  typedef std::pair<std::string,Teuchos::RCP<panzer::PureBasis> > StrPureBasisPair;

  //! Simple binary comparison class to help with sorting
  struct StrPureBasisComp {
    bool operator() (const StrPureBasisPair & lhs, const StrPureBasisPair & rhs) const
    {return lhs.first<rhs.first;}
  };

}

#endif
