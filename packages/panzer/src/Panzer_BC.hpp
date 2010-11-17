
#ifndef PANZER_BC_H
#define PANZER_BC_H

#include <string>
#include <iostream>
#include <functional>
#include <cstddef>
#include "Teuchos_RCP.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace panzer {

  //! Type of boundary condition.
  enum BCType {
    BCT_Dirichlet,
    BCT_Neumann
  };

  //! Stores input information for a boundary condition.
  class BC {
    
  public:

    //! \brief Ctor.
    BC(std::size_t bc_id,
       BCType bc_type,
       std::string sideset_id,
       std::string element_block_id,
       std::string equation_set_name,
       std::string strategy);

    //! \brief Ctor with Teuchos::ParameterList for extra params.
    BC(std::size_t bc_id,
       BCType bc_type,
       std::string sideset_id,
       std::string element_block_id,
       std::string equation_set_name,
       std::string strategy,
       const Teuchos::ParameterList& p);

    //! Dtor.
    ~BC();

    //! Returns a unique identifier for this bc - needed for unique parameter setting in LOCA and for map key comparisons (strict weak ordering).
    std::size_t bcID() const;

    //! Returns the boundary condition type (Dirichlet or Neumann).
    BCType bcType() const;

    //! Returns the set id.
    std::string sidesetID() const;

    //! Returns the element block id associated with this sideset.
    std::string elementBlockID() const;

    //! Returns the unknown name/keyword.
    std::string equationSetName() const;

    //! Returns the keyword used to construct a bc strategy.
    std::string strategy() const;

    //! Returns a parameter list with user defined parameters for bc.
    Teuchos::RCP<const Teuchos::ParameterList> params() const;

    //! A unique string identifier for this boundary condition.
    std::string identifier() const;

    //! Print object using an ostream.
    void print(std::ostream& os) const;

  private:

    std::size_t m_bc_id;

    BCType m_bc_type;

    std::string m_sideset_id;

    std::string m_element_block_id;

    std::string m_equation_set_name;

    std::string m_strategy;

    Teuchos::RCP<Teuchos::ParameterList> m_params;

  };

  std::ostream& 
    operator<<(std::ostream & os, const panzer::BC& bc);

  struct LessBC {
    
    bool operator()(const panzer::BC& left, 
		    const panzer::BC& right) const
    {
      return left.bcID() < right.bcID();
    }
  };

}

#endif
