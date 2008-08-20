// @HEADER
// ************************************************************************
// 
//            Phalanx: A Partial Differential Equation Assembly 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#ifndef PHX_DATA_LAYOUT_MDALAYOUT
#define PHX_DATA_LAYOUT_MDALAYOUT

#include "Phalanx_DataLayout.hpp"
#include <string>
#include <iostream>
#include <vector>
#include "Phalanx_DimTag.hpp"

namespace PHX {

  // ******************************************************************
  template<std::size_t, std::size_t> struct
  Check_num_ctor_arguments_equal_to_num_template_arguments;

  // ******************************************************************
  template<typename Tag0, typename Tag1, typename Tag2, typename Tag,
	   typename Tag4, typename Tag5, typename Tag6, typename Tag7>
  struct DLTagList;

  // ******************************************************************
  /*! \brief A concrete implementation of the DataLayout class for compile time multidimensional arrays.
  */
  template<typename Tag0 = void, typename Tag1 = void, 
	   typename Tag2 = void, typename Tag3 = void,
	   typename Tag4 = void, typename Tag5 = void,
	   typename Tag6 = void, typename Tag7 = void>
  class MDALayout;

  // ******************************************************************
  template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	   typename Tag4, typename Tag5, typename Tag6, typename Tag7>
  class MDALayout : public DataLayout {

  public:

    typedef std::size_t size_type;

    typedef PHX::DLTagList<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> tag_list;

    enum { MaxRank = 8 };

    enum { Rank = tag_list::Rank };

    MDALayout(std::size_t size1, std::size_t size2, std::size_t size3, 
	      std::size_t size4, std::size_t size5, std::size_t size6,
	      std::size_t size7, std::size_t size8);

    MDALayout(std::size_t size1, std::size_t size2, std::size_t size3, 
	      std::size_t size4, std::size_t size5, std::size_t size6,
	      std::size_t size7);

    MDALayout(std::size_t size1, std::size_t size2, std::size_t size3, 
	      std::size_t size4, std::size_t size5, std::size_t size6);

    MDALayout(std::size_t size1, std::size_t size2, std::size_t size3, 
	      std::size_t size4, std::size_t size5);

    MDALayout(std::size_t size1, std::size_t size2, std::size_t size3, 
	      std::size_t size4);

    MDALayout(std::size_t size1, std::size_t size2, std::size_t size3);

    MDALayout(std::size_t size1, std::size_t size2);

    MDALayout(std::size_t size1);

    ~MDALayout() {}

    virtual bool operator==(const DataLayout& right) const;

    virtual std::size_t rank() const; 

    virtual void dimensions(std::vector<std::size_t>& dim) const; 

    virtual std::size_t size() const;

    virtual const std::string identifier() const;

    virtual const std::size_t dimension(std::size_t ordinal) const;

    virtual void print(std::ostream& os, int offset) const;

  private:
    
    std::string createIdentifier();

  private:

    std::vector<const char*> m_dim_name;

    std::size_t m_dim_size[Rank];

    std::size_t m_size;

    std::string m_identifier;

  };

  template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	   typename Tag4, typename Tag5, typename Tag6, typename Tag7>
  std::ostream& operator<<(std::ostream& os, 
			   const PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,
			   Tag5,Tag6,Tag7>& t);

}

#include "Phalanx_DataLayout_MDALayout_Def.hpp"

#endif
