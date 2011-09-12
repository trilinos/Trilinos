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
