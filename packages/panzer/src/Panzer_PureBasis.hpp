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

#ifndef PANZER_PureBasis_HPP
#define PANZER_PureBasis_HPP

#include <string>
#include "Teuchos_RCP.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_Basis.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"

namespace panzer {

  class CellData;

  //! Description and data layouts associated with a particular basis
  class PureBasis { 

  public:
    typedef enum { HGRAD=0, HCURL=1, HDIV=2 } EElementSpace;
    
    /** Build a basis given a type, order and CellData object
      \param[in] basis_type String name that describes the type of basis
      \param[in] basis_order Order of the basis
      \param[in] cell_data Description of the basis
    */
    PureBasis(const std::string & basis_type,const int basis_order,const CellData & cell_data);

    /** Build a basis given a type, order, number of cells (for data layouts) and shards topology
      \param[in] basis_type String name that describes the type of basis
      \param[in] basis_order Order of the basis
      \param[in] num_cells Number of cells used in the data layouts for this basis
      \param[in] cell_topo A shards topology description
    */
    PureBasis(const std::string & basis_type,const int basis_order,const int num_cells,const Teuchos::RCP<const shards::CellTopology> & cell_topo);

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

    Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > 
    getIntrepidBasis() const;

    template <typename ScalarT,typename ArrayT>
    Teuchos::RCP< Intrepid::Basis<ScalarT,ArrayT> > 
    getIntrepidBasis() const
    { return panzer::createIntrepidBasis<ScalarT,ArrayT>(type(), order(), getCellTopology()); }

    EElementSpace getElementSpace() const
    { return element_space_; }

    bool requiresOrientations() const
    { return getElementSpace()==HCURL || getElementSpace()==HDIV; }

    bool supportsGrad() const
    { return getElementSpace()==HGRAD; }

    bool supportsCurl() const
    { return getElementSpace()==HCURL; }

    bool supportsDiv() const
    { return getElementSpace()==HDIV; }

    bool isVectorBasis() const
    { return getElementSpace()==HCURL || getElementSpace()==HDIV; }

    bool isScalarBasis() const
    { return getElementSpace()==HGRAD; }

    int getBasisRank() const
    { return basis_rank_; }

    Teuchos::RCP<const shards::CellTopology> getCellTopology() const
    { return topology_; }

  public:
    //! <Cell,Basis> or <Cell,Basis>
    Teuchos::RCP<PHX::DataLayout> functional;
    //! <Cell,Basis,Dim>
    Teuchos::RCP<PHX::DataLayout> functional_grad;
    //! <Cell,Basis,Dim>
    Teuchos::RCP<PHX::DataLayout> functional_D2;
    //! <Cell,Basis,Dim>
    Teuchos::RCP<PHX::DataLayout> coordinates;

  private:
    
    //! Initialize the basis object
    void initialize(const std::string & basis_type,const int basis_order);

  private:

    Teuchos::RCP<const shards::CellTopology> topology_;
    Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > intrepid_basis_;

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
