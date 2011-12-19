#ifndef PANZER_PureBasis_HPP
#define PANZER_PureBasis_HPP

#include <string>
#include "Teuchos_RCP.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_Basis.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_CellData.hpp"

namespace panzer {

  class PureBasis { 

  public:
    
    PureBasis(std::string basis_type,const CellData & cell_data);

    int getCardinality() const;

    int getNumCells() const;
    
    int getDimension() const;

    std::string name() const;
    
    std::string fieldName() const;
    
    std::string fieldNameD1() const;
    
    std::string fieldNameD2() const;

    Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > 
    getIntrepidBasis() const;

  public:
    //! <Cell,Basis>
    Teuchos::RCP<PHX::DataLayout> functional;
    //! <Cell,Basis,Dim>
    Teuchos::RCP<PHX::DataLayout> functional_grad;
    //! <Cell,Basis,Dim,Dim>
    Teuchos::RCP<PHX::DataLayout> functional_D2;

  private:
    Teuchos::RCP<const shards::CellTopology> topology;
    Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > intrepid_basis;

    const std::string basis_name;
    const std::string field_basis_name;
    const std::string field_basis_name_D1;
    const std::string field_basis_name_D2;

    int cardinality;
    int num_cells;
    int dimension;
  };

  typedef std::pair<std::string,Teuchos::RCP<panzer::PureBasis> > StrPureBasisPair;

  //! Simple binary comparison class to help with sorting
  struct StrPureBasisComp {
    bool operator() (const StrPureBasisPair & lhs, const StrPureBasisPair & rhs) const
    {return lhs.first<rhs.first;}
  };

}

#endif
