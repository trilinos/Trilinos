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

#ifndef PHX_DATA_LAYOUT_MDALAYOUT_DEF
#define PHX_DATA_LAYOUT_MDALAYOUT_DEF

#include <iostream>
#include <sstream>
#include <typeinfo>

//**********************************************************************
namespace PHX{

  template<> struct 
  Check_num_ctor_arguments_equal_to_num_template_arguments<8,8> {};
  template<> struct 
  Check_num_ctor_arguments_equal_to_num_template_arguments<7,7> {};
  template<> struct 
  Check_num_ctor_arguments_equal_to_num_template_arguments<6,6> {};
  template<> struct 
  Check_num_ctor_arguments_equal_to_num_template_arguments<5,5> {};
  template<> struct 
  Check_num_ctor_arguments_equal_to_num_template_arguments<4,4> {};
  template<> struct 
  Check_num_ctor_arguments_equal_to_num_template_arguments<3,3> {};
  template<> struct 
  Check_num_ctor_arguments_equal_to_num_template_arguments<2,2> {};
  template<> struct 
  Check_num_ctor_arguments_equal_to_num_template_arguments<1,1> {};
}
//**********************************************************************
namespace PHX{
  
template<typename T0, typename T1, typename T2, typename T3,
	 typename T4, typename T5, typename T6, typename T7>
struct DLTagList
{ enum { Rank = 8 }; };

template<typename T0, typename T1, typename T2, typename T3,
	 typename T4, typename T5, typename T6>
struct DLTagList<T0,T1,T2,T3,T4,T5,T6,void>
{ enum { Rank = 7 }; };

template<typename T0, typename T1, typename T2, typename T3,
	 typename T4, typename T5>
struct DLTagList<T0,T1,T2,T3,T4,T5,void,void>
{ enum { Rank = 6 }; };

template<typename T0, typename T1, typename T2, typename T3,
	 typename T4>
struct DLTagList<T0,T1,T2,T3,T4,void,void,void>
{ enum { Rank = 5 }; };

template<typename T0, typename T1, typename T2, typename T3>
struct DLTagList<T0,T1,T2,T3,void,void,void,void>
{ enum { Rank = 4 }; };

template<typename T0, typename T1, typename T2>
struct DLTagList<T0,T1,T2,void,void,void,void,void>
{ enum { Rank = 3 }; };

template<typename T0, typename T1>
struct DLTagList<T0,T1,void,void,void,void,void,void>
{ enum { Rank = 2 }; };

template<typename T0>
struct DLTagList<T0,void,void,void,void,void,void,void>
{ enum { Rank = 1 }; };
}
//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(std::size_t size1, std::size_t size2, std::size_t size3, 
	  std::size_t size4, std::size_t size5, std::size_t size6,
	  std::size_t size7, std::size_t size8)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,8>();

  m_dim_size[0] = size1;
  m_dim_size[1] = size2;
  m_dim_size[2] = size3;
  m_dim_size[3] = size4;
  m_dim_size[4] = size5;
  m_dim_size[5] = size6;
  m_dim_size[6] = size7;
  m_dim_size[7] = size8;

  m_dim_name.push_back(Tag0::descriptor().name());
  m_dim_name.push_back(Tag1::descriptor().name());
  m_dim_name.push_back(Tag2::descriptor().name());
  m_dim_name.push_back(Tag3::descriptor().name());
  m_dim_name.push_back(Tag4::descriptor().name());
  m_dim_name.push_back(Tag5::descriptor().name());
  m_dim_name.push_back(Tag6::descriptor().name());
  m_dim_name.push_back(Tag7::descriptor().name());

  m_size = size1 * size2 * size3 * size4 * size5 * size6 * size7 * size8;

  m_identifier = this->createIdentifier();
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(std::size_t size1, std::size_t size2, std::size_t size3, 
	  std::size_t size4, std::size_t size5, std::size_t size6,
	  std::size_t size7)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,7>();

  m_dim_size[0] = size1;
  m_dim_size[1] = size2;
  m_dim_size[2] = size3;
  m_dim_size[3] = size4;
  m_dim_size[4] = size5;
  m_dim_size[5] = size6;
  m_dim_size[6] = size7;

  m_dim_name.push_back(Tag0::descriptor().name());
  m_dim_name.push_back(Tag1::descriptor().name());
  m_dim_name.push_back(Tag2::descriptor().name());
  m_dim_name.push_back(Tag3::descriptor().name());
  m_dim_name.push_back(Tag4::descriptor().name());
  m_dim_name.push_back(Tag5::descriptor().name());
  m_dim_name.push_back(Tag6::descriptor().name());

  m_size = size1 * size2 * size3 * size4 * size5 * size6 * size7;

  m_identifier = this->createIdentifier();
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(std::size_t size1, std::size_t size2, std::size_t size3, 
	  std::size_t size4, std::size_t size5, std::size_t size6)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,6>();

  m_dim_size[0] = size1;
  m_dim_size[1] = size2;
  m_dim_size[2] = size3;
  m_dim_size[3] = size4;
  m_dim_size[4] = size5;
  m_dim_size[5] = size6;

  m_dim_name.push_back(Tag0::descriptor().name());
  m_dim_name.push_back(Tag1::descriptor().name());
  m_dim_name.push_back(Tag2::descriptor().name());
  m_dim_name.push_back(Tag3::descriptor().name());
  m_dim_name.push_back(Tag4::descriptor().name());
  m_dim_name.push_back(Tag5::descriptor().name());

  m_size = size1 * size2 * size3 * size4 * size5 * size6;

  m_identifier = this->createIdentifier();
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(std::size_t size1, std::size_t size2, std::size_t size3, 
	  std::size_t size4, std::size_t size5)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,5>();

  m_dim_size[0] = size1;
  m_dim_size[1] = size2;
  m_dim_size[2] = size3;
  m_dim_size[3] = size4;
  m_dim_size[4] = size5;

  m_dim_name.push_back(Tag0::descriptor().name());
  m_dim_name.push_back(Tag1::descriptor().name());
  m_dim_name.push_back(Tag2::descriptor().name());
  m_dim_name.push_back(Tag3::descriptor().name());
  m_dim_name.push_back(Tag4::descriptor().name());

  m_size = size1 * size2 * size3 * size4 * size5;

  m_identifier = this->createIdentifier();
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(std::size_t size1, std::size_t size2, std::size_t size3, 
	  std::size_t size4)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,4>();

  m_dim_size[0] = size1;
  m_dim_size[1] = size2;
  m_dim_size[2] = size3;
  m_dim_size[3] = size4;

  m_dim_name.push_back(Tag0::descriptor().name());
  m_dim_name.push_back(Tag1::descriptor().name());
  m_dim_name.push_back(Tag2::descriptor().name());
  m_dim_name.push_back(Tag3::descriptor().name());

  m_size = size1 * size2 * size3 * size4;

  m_identifier = this->createIdentifier();
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(std::size_t size1, std::size_t size2, std::size_t size3)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,3>();

  m_dim_size[0] = size1;
  m_dim_size[1] = size2;
  m_dim_size[2] = size3;

  m_dim_name.push_back(Tag0::descriptor().name());
  m_dim_name.push_back(Tag1::descriptor().name());
  m_dim_name.push_back(Tag2::descriptor().name());

  m_size = size1 * size2 * size3;

  m_identifier = this->createIdentifier();
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(std::size_t size1, std::size_t size2)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,2>();

  m_dim_size[0] = size1;
  m_dim_size[1] = size2;

  m_dim_name.push_back(Tag0::descriptor().name());
  m_dim_name.push_back(Tag1::descriptor().name());

  m_size = size1 * size2;

  m_identifier = this->createIdentifier();
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(std::size_t size1)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,1>();

  m_dim_size[0] = size1;

  m_dim_name.push_back(Tag0::descriptor().name());

  m_size = size1;

  m_identifier = this->createIdentifier();
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
bool PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
operator==(const PHX::DataLayout& right) const
{
  const PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>* tmp = 0;
  tmp = dynamic_cast< const PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>* >(&right);

  if (tmp == 0)
    return false;

  for (std::size_t i=0; i < Rank; ++i)
    if (m_dim_size[i] != tmp->m_dim_size[i])
      return false;

  return (this->size() == tmp->size());
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
std::size_t PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
size() const
{ return m_size; }

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
const std::string PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
identifier() const
{ 
  return m_identifier; 
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
void PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
print(std::ostream& os, int offset) const
{
  os << m_identifier;
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
std::string PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
createIdentifier()
{
  std::ostringstream os;
  os << "MDA<";
  for (std::size_t i=0; i < m_dim_name.size(); ++i) {
    if (i > 0)
      os << ",";
    os << std::string(m_dim_name[i]);
  }
  os << ">(";
  for (std::size_t i=0; i < Rank; ++i) {
    if (i > 0)
      os << ",";
    os << m_dim_size[i];
  }
  os << ")";
  
  return os.str();
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
std::ostream& PHX::operator<<(std::ostream& os, 
			      const PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,
			      Tag5,Tag6,Tag7>& v)
{
  v.print(os,0);
  return os;
}
 
//**********************************************************************

#endif
