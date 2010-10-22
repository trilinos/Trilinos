
#ifndef PANZER_BASIS_HPP
#define PANZER_BASIS_HPP

#include <string>
#include "Teuchos_RCP.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_Basis.hpp"
#include "Panzer_Dimension.hpp"

namespace panzer {

  class IntegrationRule;

  class Basis { 

  public:
    
    Basis(std::string basis_type, const panzer::IntegrationRule& int_rule);

    int getCardinality() const;
    
    int getNumCells() const;
    
    int getNumIntPoints() const;
    
    int getDimension() const;

    int integrationRuleDegree() const;
    
    std::string name() const;
    
    std::string fieldName() const;
    
    std::string fieldNameD1() const;
    
    std::string fieldNameD2() const;

    Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > 
    getIntrepidBasis() const;

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
    int int_rule_degree;
  };

}

#endif

