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


#ifndef PANZER_BC_HPP
#define PANZER_BC_HPP

#include <string>
#include <iostream>
#include <functional>
#include <cstddef>
#include <vector>

#include <unordered_map>

#include "Teuchos_RCP.hpp"

#include "Panzer_GlobalData.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace panzer {

  class BC;
  class WorksetDescriptor;

  /** \brief Nonmember constructor to build BC objects from a ParameterList
      \relates panzer::BC
  */
  void buildBCs(std::vector<panzer::BC>& bcs, const Teuchos::ParameterList& p, const Teuchos::RCP<panzer::GlobalData> global_data);

  //! Type of boundary condition.
  enum BCType {
    BCT_Dirichlet,
    BCT_Neumann,
    BCT_Interface
  };

  //! Stores input information for a boundary condition.
  class BC {
  public:
    // types supporting hashing
    struct BCHash {
      std::hash<std::string> hash;
      std::size_t operator()(const BC & bc) const
      { return this->hash(bc.elementBlockID() + "_" + bc.sidesetID());}
    };

    struct BCEquality {
      bool operator()(const BC & bc1,const BC & bc2) const
      { return bc1.elementBlockID()==bc2.elementBlockID() && bc1.sidesetID()==bc2.sidesetID(); }
    };
    
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

    //! \brief Ctor based on ParameterList
    BC(std::size_t bc_id,
       const Teuchos::ParameterList& p);

    //! \brief Ctor based on ParameterList
    BC(std::size_t bc_id,
       const Teuchos::ParameterList& p,
       const Teuchos::RCP<panzer::GlobalData> gd);

    //! Dtor.
    ~BC();

    //! Returns a unique identifier for this bc - needed for unique parameter setting in LOCA and for map key comparisons (strict weak ordering).
    std::size_t bcID() const;

    //! Returns the boundary condition type (Dirichlet or Neumann or Interface).
    BCType bcType() const;

    //! Returns the set id.
    std::string sidesetID() const;

    //! Returns the element block id associated with this sideset.
    std::string elementBlockID() const;

    //! Returns the second element block id associated with this sideset.
    std::string elementBlockID2() const;

    //! Returns the unknown name/keyword.
    std::string equationSetName() const;

    //! Returns the second unknown name/keyword.
    std::string equationSetName2() const;

    //! Returns the keyword used to construct a bc strategy.
    std::string strategy() const;

    //! Returns a parameter list with user defined parameters for bc.
    Teuchos::RCP<const Teuchos::ParameterList> params() const;

    //! Returns the RCP to the global data.
    Teuchos::RCP<panzer::GlobalData> global_data() const;

    //! Returns a nonconst parameter list with user defined parameters for bc.  Nonconst is meant to be used for parameter list validation.
    Teuchos::RCP<Teuchos::ParameterList> nonconstParams() const;

    //! A unique string identifier for this boundary condition.
    std::string identifier() const;

    //! Print object using an ostream.
    void print(std::ostream& os) const;

  private:

    void validateParameters(Teuchos::ParameterList& p) const;

  private:

    std::size_t m_bc_id;

    BCType m_bc_type;

    std::string m_sideset_id;

    std::string m_element_block_id;

    std::string m_element_block_id2;

    std::string m_equation_set_name;

    std::string m_equation_set_name2;

    std::string m_strategy;

    Teuchos::RCP<Teuchos::ParameterList> m_params;

    Teuchos::RCP<panzer::GlobalData> m_gd;
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

  WorksetDescriptor bcDescriptor(const panzer::BC & bc);

}

#endif
