
#ifndef PANZER_BASIS_HPP
#define PANZER_BASIS_HPP

#include <string>
#include "Teuchos_RCP.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_Basis.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_PureBasis.hpp"

namespace panzer {

  class PointRule;

  class BasisIRLayout { 

  public:
    
    BasisIRLayout(std::string basis_type, const PointRule& int_rule);
    BasisIRLayout(const Teuchos::RCP<const PureBasis> & b, const PointRule& int_rule);

    void setup(const Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > & iBasis,
               const panzer::PointRule & int_rule);

    int getCardinality() const;
    
    int getNumCells() const;
    
    int getNumPoints() const;
    
    int getDimension() const;

    // int integrationRuleDegree() const;
    
    std::string name() const;
    
    std::string fieldName() const;
    
    std::string fieldNameD1() const;
    
    std::string fieldNameD2() const;

    Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > 
    getIntrepidBasis() const;

    Teuchos::RCP<const PureBasis> getBasis() const;

    void print(std::ostream & os) const;

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
  };

  typedef std::pair<std::string,Teuchos::RCP<panzer::BasisIRLayout> > StrBasisPair;

  //! Simple binary comparison class to help with sorting
  struct StrBasisComp {
    bool operator() (const StrBasisPair & lhs, const StrBasisPair & rhs) const
    {return lhs.first<rhs.first;}
  };

}

#endif
