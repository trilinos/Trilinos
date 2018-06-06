// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
#include "Teuchos_Assert.hpp"

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
MDALayout(size_type size1, size_type size2, size_type size3, 
	  size_type size4, size_type size5, size_type size6,
	  size_type size7, size_type size8)
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

  m_dim_name.push_back(Tag0::tag().name());
  m_dim_name.push_back(Tag1::tag().name());
  m_dim_name.push_back(Tag2::tag().name());
  m_dim_name.push_back(Tag3::tag().name());
  m_dim_name.push_back(Tag4::tag().name());
  m_dim_name.push_back(Tag5::tag().name());
  m_dim_name.push_back(Tag6::tag().name());
  m_dim_name.push_back(Tag7::tag().name());

  m_size = size1 * size2 * size3 * size4 * size5 * size6 * size7 * size8;

  m_identifier = this->createIdentifier();
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(size_type size1, size_type size2, size_type size3, 
	  size_type size4, size_type size5, size_type size6,
	  size_type size7)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,7>();

  m_dim_size[0] = size1;
  m_dim_size[1] = size2;
  m_dim_size[2] = size3;
  m_dim_size[3] = size4;
  m_dim_size[4] = size5;
  m_dim_size[5] = size6;
  m_dim_size[6] = size7;

  m_dim_name.push_back(Tag0::tag().name());
  m_dim_name.push_back(Tag1::tag().name());
  m_dim_name.push_back(Tag2::tag().name());
  m_dim_name.push_back(Tag3::tag().name());
  m_dim_name.push_back(Tag4::tag().name());
  m_dim_name.push_back(Tag5::tag().name());
  m_dim_name.push_back(Tag6::tag().name());

  m_size = size1 * size2 * size3 * size4 * size5 * size6 * size7;

  m_identifier = this->createIdentifier();
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(size_type size1, size_type size2, size_type size3, 
	  size_type size4, size_type size5, size_type size6)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,6>();

  m_dim_size[0] = size1;
  m_dim_size[1] = size2;
  m_dim_size[2] = size3;
  m_dim_size[3] = size4;
  m_dim_size[4] = size5;
  m_dim_size[5] = size6;

  m_dim_name.push_back(Tag0::tag().name());
  m_dim_name.push_back(Tag1::tag().name());
  m_dim_name.push_back(Tag2::tag().name());
  m_dim_name.push_back(Tag3::tag().name());
  m_dim_name.push_back(Tag4::tag().name());
  m_dim_name.push_back(Tag5::tag().name());

  m_size = size1 * size2 * size3 * size4 * size5 * size6;

  m_identifier = this->createIdentifier();
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(size_type size1, size_type size2, size_type size3, 
	  size_type size4, size_type size5)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,5>();

  m_dim_size[0] = size1;
  m_dim_size[1] = size2;
  m_dim_size[2] = size3;
  m_dim_size[3] = size4;
  m_dim_size[4] = size5;

  m_dim_name.push_back(Tag0::tag().name());
  m_dim_name.push_back(Tag1::tag().name());
  m_dim_name.push_back(Tag2::tag().name());
  m_dim_name.push_back(Tag3::tag().name());
  m_dim_name.push_back(Tag4::tag().name());

  m_size = size1 * size2 * size3 * size4 * size5;

  m_identifier = this->createIdentifier();
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(size_type size1, size_type size2, size_type size3, 
	  size_type size4)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,4>();

  m_dim_size[0] = size1;
  m_dim_size[1] = size2;
  m_dim_size[2] = size3;
  m_dim_size[3] = size4;

  m_dim_name.push_back(Tag0::tag().name());
  m_dim_name.push_back(Tag1::tag().name());
  m_dim_name.push_back(Tag2::tag().name());
  m_dim_name.push_back(Tag3::tag().name());

  m_size = size1 * size2 * size3 * size4;

  m_identifier = this->createIdentifier();
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(size_type size1, size_type size2, size_type size3)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,3>();

  m_dim_size[0] = size1;
  m_dim_size[1] = size2;
  m_dim_size[2] = size3;

  m_dim_name.push_back(Tag0::tag().name());
  m_dim_name.push_back(Tag1::tag().name());
  m_dim_name.push_back(Tag2::tag().name());

  m_size = size1 * size2 * size3;

  m_identifier = this->createIdentifier();
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(size_type size1, size_type size2)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,2>();

  m_dim_size[0] = size1;
  m_dim_size[1] = size2;

  m_dim_name.push_back(Tag0::tag().name());
  m_dim_name.push_back(Tag1::tag().name());

  m_size = size1 * size2;

  m_identifier = this->createIdentifier();
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(size_type size1)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,1>();

  m_dim_size[0] = size1;

  m_dim_name.push_back(Tag0::tag().name());

  m_size = size1;

  m_identifier = this->createIdentifier();
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(const std::string& prefix,
	  size_type size1, size_type size2, size_type size3, 
	  size_type size4, size_type size5, size_type size6,
	  size_type size7, size_type size8)
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

  m_dim_name.push_back(Tag0::tag().name());
  m_dim_name.push_back(Tag1::tag().name());
  m_dim_name.push_back(Tag2::tag().name());
  m_dim_name.push_back(Tag3::tag().name());
  m_dim_name.push_back(Tag4::tag().name());
  m_dim_name.push_back(Tag5::tag().name());
  m_dim_name.push_back(Tag6::tag().name());
  m_dim_name.push_back(Tag7::tag().name());

  m_size = size1 * size2 * size3 * size4 * size5 * size6 * size7 * size8;

  m_identifier = this->createIdentifier(prefix);
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(const std::string& prefix,
	  size_type size1, size_type size2, size_type size3, 
	  size_type size4, size_type size5, size_type size6,
	  size_type size7)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,7>();

  m_dim_size[0] = size1;
  m_dim_size[1] = size2;
  m_dim_size[2] = size3;
  m_dim_size[3] = size4;
  m_dim_size[4] = size5;
  m_dim_size[5] = size6;
  m_dim_size[6] = size7;

  m_dim_name.push_back(Tag0::tag().name());
  m_dim_name.push_back(Tag1::tag().name());
  m_dim_name.push_back(Tag2::tag().name());
  m_dim_name.push_back(Tag3::tag().name());
  m_dim_name.push_back(Tag4::tag().name());
  m_dim_name.push_back(Tag5::tag().name());
  m_dim_name.push_back(Tag6::tag().name());

  m_size = size1 * size2 * size3 * size4 * size5 * size6 * size7;

  m_identifier = this->createIdentifier(prefix);
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(const std::string& prefix,
	  size_type size1, size_type size2, size_type size3, 
	  size_type size4, size_type size5, size_type size6)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,6>();

  m_dim_size[0] = size1;
  m_dim_size[1] = size2;
  m_dim_size[2] = size3;
  m_dim_size[3] = size4;
  m_dim_size[4] = size5;
  m_dim_size[5] = size6;

  m_dim_name.push_back(Tag0::tag().name());
  m_dim_name.push_back(Tag1::tag().name());
  m_dim_name.push_back(Tag2::tag().name());
  m_dim_name.push_back(Tag3::tag().name());
  m_dim_name.push_back(Tag4::tag().name());
  m_dim_name.push_back(Tag5::tag().name());

  m_size = size1 * size2 * size3 * size4 * size5 * size6;

  m_identifier = this->createIdentifier(prefix);
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(const std::string& prefix,
	  size_type size1, size_type size2, size_type size3, 
	  size_type size4, size_type size5)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,5>();

  m_dim_size[0] = size1;
  m_dim_size[1] = size2;
  m_dim_size[2] = size3;
  m_dim_size[3] = size4;
  m_dim_size[4] = size5;

  m_dim_name.push_back(Tag0::tag().name());
  m_dim_name.push_back(Tag1::tag().name());
  m_dim_name.push_back(Tag2::tag().name());
  m_dim_name.push_back(Tag3::tag().name());
  m_dim_name.push_back(Tag4::tag().name());

  m_size = size1 * size2 * size3 * size4 * size5;

  m_identifier = this->createIdentifier(prefix);
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(const std::string& prefix,
	  size_type size1, size_type size2, size_type size3, 
	  size_type size4)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,4>();

  m_dim_size[0] = size1;
  m_dim_size[1] = size2;
  m_dim_size[2] = size3;
  m_dim_size[3] = size4;

  m_dim_name.push_back(Tag0::tag().name());
  m_dim_name.push_back(Tag1::tag().name());
  m_dim_name.push_back(Tag2::tag().name());
  m_dim_name.push_back(Tag3::tag().name());

  m_size = size1 * size2 * size3 * size4;

  m_identifier = this->createIdentifier(prefix);
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(const std::string& prefix,
	  size_type size1, size_type size2, size_type size3)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,3>();

  m_dim_size[0] = size1;
  m_dim_size[1] = size2;
  m_dim_size[2] = size3;

  m_dim_name.push_back(Tag0::tag().name());
  m_dim_name.push_back(Tag1::tag().name());
  m_dim_name.push_back(Tag2::tag().name());

  m_size = size1 * size2 * size3;

  m_identifier = this->createIdentifier(prefix);
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(const std::string& prefix,
	  size_type size1, size_type size2)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,2>();

  m_dim_size[0] = size1;
  m_dim_size[1] = size2;

  m_dim_name.push_back(Tag0::tag().name());
  m_dim_name.push_back(Tag1::tag().name());

  m_size = size1 * size2;

  m_identifier = this->createIdentifier(prefix);
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
MDALayout(const std::string& prefix,
	  size_type size1)
{ 
  PHX::Check_num_ctor_arguments_equal_to_num_template_arguments<Rank,1>();

  m_dim_size[0] = size1;

  m_dim_name.push_back(Tag0::tag().name());

  m_size = size1;

  m_identifier = this->createIdentifier(prefix);
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

  for (size_type i=0; i < Rank; ++i)
    if (m_dim_size[i] != tmp->m_dim_size[i])
      return false;

  return (this->size() == tmp->size());
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::Device::size_type 
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::rank() const
{ return Rank; }

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
void PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
dimensions(std::vector<PHX::Device::size_type>& dim) const
{ 
  dim.resize(Rank);
  for(std::size_t i=0; i < dim.size(); ++i)
    dim[i] = m_dim_size[i];
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
void PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
names(std::vector<std::string>& names) const
{ 
  names.resize(Rank);
  for(std::size_t i=0; i < names.size(); ++i)
    names[i] = m_dim_name[i];
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::Device::size_type PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
size() const
{ return m_size; }

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
std::string PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
identifier() const
{ 
  return m_identifier; 
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::Device::size_type 
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
dimension(size_type ordinal) const
{ 
  this->checkForValidRank(ordinal);
  return m_dim_size[ordinal];
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
PHX::Device::size_type 
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
extent(size_type ordinal) const
{ 
  this->checkForValidRank(ordinal);
  return m_dim_size[ordinal];
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
int 
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
extent_int(size_type ordinal) const
{
  this->checkForValidRank(ordinal);
  return static_cast<int>(m_dim_size[ordinal]);
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
std::string
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
name(size_type ordinal) const
{ 
  this->checkForValidRank(ordinal);
  return m_dim_name[ordinal];
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
void PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
print(std::ostream& os, int /* offset */) const
{
  os << m_identifier;
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
std::string PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
createIdentifier(const std::string& prefix)
{
  std::ostringstream os;
  os << prefix << "<";
  for (std::size_t i=0; i < m_dim_name.size(); ++i) {
    if (i > 0)
      os << ",";
    os << std::string(m_dim_name[i]);
  }
  os << ">(";
  for (size_type i=0; i < Rank; ++i) {
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
void
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
setExtentsOnDerivedClass(const std::vector<PHX::Device::size_type>& extents)
{
#ifdef PHX_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(extents.size() != Rank, std::runtime_error,
                             "ERROR - MDALayout::setExtentsOnDerivedClass()"
                             << " - Rank mismatch while setting extents!");
#endif
  for (std::size_t i=0; i < extents.size(); ++i)
    m_dim_size[i] = extents[i];
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
template<typename IndexType>
typename std::enable_if<std::is_signed<IndexType>::value>::type
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
checkForValidRank(const IndexType& ordinal) const
{
  if (ordinal > Rank-1 || ordinal < 0) {
    std::ostringstream os;
    os << "Requested Ordinal " << ordinal 
       << " is outside the valid range of 0 - " << Rank - 1
       << " in DataLayout object:\n"
       << m_identifier << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(ordinal > Rank-1 || ordinal < 0, 
			       std::runtime_error, os.str());
  }
}

//**********************************************************************
template<typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
template<typename IndexType>
typename std::enable_if<std::is_unsigned<IndexType>::value>::type
PHX::MDALayout<Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>::
checkForValidRank(const IndexType& ordinal) const
{
  if (ordinal > Rank-1) {
    std::ostringstream os;
    os << "Requested Ordinal " << ordinal 
       << " is outside the valid range of 0 - " << Rank - 1
       << " in DataLayout object:\n"
       << m_identifier << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(ordinal > Rank-1, 
			       std::runtime_error, os.str());
  }
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
