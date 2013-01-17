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


#ifndef PANZER_BASIS_HPP
#define PANZER_BASIS_HPP

#include <string>
#include "Teuchos_RCP.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_Basis.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_CellTopologyInfo.hpp"

namespace panzer {

  class PointRule;
  class CellTopologyInfo;   

  class BasisIRLayout { 

  public:
    
    BasisIRLayout(std::string basis_type, const int basis_order, const PointRule& int_rule);
    BasisIRLayout(const Teuchos::RCP<const PureBasis> & b, const PointRule& int_rule);

    void setup(const Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > & iBasis,
               const panzer::PointRule & int_rule);

    int getCardinality() const;
    
    int getNumCells() const;
    
    int getNumPoints() const;
    
    int getDimension() const;

    std::string name() const;
    
    std::string fieldName() const;
    
    std::string fieldNameD1() const;
    
    std::string fieldNameD2() const;

    Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > 
    getIntrepidBasis() const;

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
    Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > intrepid_basis;
    const std::string basis_name;
    const std::string field_basis_name;
    const std::string field_basis_name_D1;
    const std::string field_basis_name_D2;

    int cardinality;
    int num_cells;
    int num_ip;
    int dimension;
    // int int_rule_degree;
    
    Teuchos::RCP<const PureBasis> basis_data;
    
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
