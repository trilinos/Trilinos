// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
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

namespace PHX {

  // ******************************************************************
  template<std::size_t, std::size_t> struct
  Check_num_ctor_arguments_equal_to_num_template_arguments;

  // ******************************************************************
  template<typename Tag0, typename Tag1, typename Tag2, typename Tag,
	   typename Tag4, typename Tag5, typename Tag6, typename Tag7>
  struct DLTagList;

  // ******************************************************************
  /*! \brief A concrete implementation of the DataLayout class for compile time checked multidimensional arrays.
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

    typedef PHX::DLTagList<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> tag_list;

    enum { MaxRank = 8 };

    enum { Rank = tag_list::Rank };

    MDALayout(size_type size1, size_type size2, size_type size3, 
	      size_type size4, size_type size5, size_type size6,
	      size_type size7, size_type size8);

    MDALayout(size_type size1, size_type size2, size_type size3, 
	      size_type size4, size_type size5, size_type size6,
	      size_type size7);

    MDALayout(size_type size1, size_type size2, size_type size3, 
	      size_type size4, size_type size5, size_type size6);

    MDALayout(size_type size1, size_type size2, size_type size3, 
	      size_type size4, size_type size5);

    MDALayout(size_type size1, size_type size2, size_type size3, 
	      size_type size4);

    MDALayout(size_type size1, size_type size2, size_type size3);

    MDALayout(size_type size1, size_type size2);

    MDALayout(size_type size1);

    ~MDALayout() {}

    virtual bool operator==(const DataLayout& right) const;

    virtual size_type rank() const; 

    virtual size_type dimension(size_type ordinal) const;

    virtual void dimensions(std::vector<size_type>& dim) const; 

    virtual std::string name(size_type ordinal) const;

    virtual void names(std::vector<std::string>& names) const; 

    virtual size_type size() const;

    virtual std::string identifier() const;

    virtual void print(std::ostream& os, int offset) const;

  private:
    
    std::string createIdentifier();

  private:

    std::vector<const char*> m_dim_name;

    size_type m_dim_size[Rank];

    size_type m_size;

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
