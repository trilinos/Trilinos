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

#include "Panzer_Dimension.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"

namespace panzer {

  class PureBasis { 

  public:
    typedef enum { HGRAD=0, HCURL=1, HDIV=2 } EElementSpace;
    
    PureBasis(const std::string & basis_type,const int basis_order,const CellData & cell_data);
    PureBasis(const std::string & basis_type,const int basis_order,const int numCells,const Teuchos::RCP<const shards::CellTopology> & cellTopo);

    int getCardinality() const;

    int getNumCells() const;
    
    int getDimension() const;

    std::string type() const;
    
    int order() const;
    
    //! A unique key that is the combination of the basis type and basis order
    std::string key() const;
    
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
    { return elementSpace; }

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
    { return basisRank; }

    Teuchos::RCP<const shards::CellTopology> getCellTopology() const
    { return topology; }

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

    //! Initialize all introspection data for this basis type: Intrepid does not provide this
    void initializeIntrospection(const std::string & name);

    //! Initialize data layouts for this basis: uses num_cells, cardinality, dimension
    void initializeDataLayouts();

    Teuchos::RCP<const shards::CellTopology> topology;
    Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > intrepid_basis;

    const std::string basis_type;
    std::string basis_key;
    const std::string field_basis_name;
    const std::string field_basis_name_D1;
    const std::string field_basis_name_D2;

    int cardinality;
    int num_cells;
    int dimension;

    EElementSpace elementSpace;
    int basisRank;
  };

  typedef std::pair<std::string,Teuchos::RCP<panzer::PureBasis> > StrPureBasisPair;

  //! Simple binary comparison class to help with sorting
  struct StrPureBasisComp {
    bool operator() (const StrPureBasisPair & lhs, const StrPureBasisPair & rhs) const
    {return lhs.first<rhs.first;}
  };

}

#endif
