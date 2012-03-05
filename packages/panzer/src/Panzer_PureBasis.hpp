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
    typedef enum { HGRAD=0, HCURL=1, HDIV=2 } EElementSpace;
    
    PureBasis(const std::string & basis_type,const CellData & cell_data);
    PureBasis(const std::string & basis_type,int numCells,const Teuchos::RCP<const shards::CellTopology> & cellTopo);

    int getCardinality() const;

    int getNumCells() const;
    
    int getDimension() const;

    std::string name() const;
    
    std::string fieldName() const;
    
    std::string fieldNameD1() const;
    
    std::string fieldNameD2() const;

    Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > 
    getIntrepidBasis() const;

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

    const std::string basis_name;
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
