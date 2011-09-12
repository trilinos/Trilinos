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


#ifndef PHX_MD_FIELD_H
#define PHX_MD_FIELD_H

#include <iostream>
#include <string>
#include "Teuchos_ArrayRCP.hpp"
#include "Shards_Array.hpp"
#include "Phalanx_FieldTag_Tag.hpp"

namespace PHX {

  template<typename DataT,
	   typename Tag0 = void, typename Tag1 = void, typename Tag2 = void, 
	   typename Tag3 = void, typename Tag4 = void, typename Tag5 = void,
	   typename Tag6 = void, typename Tag7 = void>
  class MDField;

  // *************************************
  // Preferred compile time checked MDField
  // *************************************

  template<typename DataT,
	   typename Tag0, typename Tag1, typename Tag2, 
	   typename Tag3, typename Tag4, typename Tag5,
	   typename Tag6, typename Tag7>
  class MDField {
    
  public:

    typedef DataT value_type;

#ifdef PHX_USE_COMPILETIME_ARRAY
    typedef typename shards::Array<DataT,shards::NaturalOrder,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> array_type;
#else
    typedef typename shards::Array<DataT,shards::NaturalOrder,void,void,void,void,void,void,void,void> array_type;
#endif
    
    typedef typename array_type::size_type size_type;

    MDField(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& t);
    
    MDField(const PHX::Tag<DataT>& v);
    
    MDField();
    
    ~MDField();
    
    const PHX::FieldTag& fieldTag() const;

    DataT& operator()(size_type index1, size_type index2, size_type index3, 
		      size_type index4, size_type index5, size_type index6,
		      size_type index7, size_type index8);

    DataT& operator()(size_type index1, size_type index2, size_type index3, 
		      size_type index4, size_type index5, size_type index6,
		      size_type index7);

    DataT& operator()(size_type index1, size_type index2, size_type index3, 
		      size_type index4, size_type index5, size_type index6);
    
    DataT& operator()(size_type index1, size_type index2, size_type index3, 
		      size_type index4, size_type index5);
    
    DataT& operator()(size_type index1, size_type index2, size_type index3, 
		      size_type index4);
    
    DataT& operator()(size_type index1, size_type index2, size_type index3);
    
    DataT& operator()(size_type index1, size_type index2);
    
    DataT& operator()(size_type index1);
    
    DataT& operator[](size_type index);

    const DataT& operator()(size_type index1, size_type index2, 
			    size_type index3, size_type index4, 
			    size_type index5, size_type index6,
			    size_type index7, size_type index8) const;

    const DataT& operator()(size_type index1, size_type index2,
			    size_type index3, size_type index4,
			    size_type index5, size_type index6,
			    size_type index7) const;

    const DataT& operator()(size_type index1, size_type index2,
			    size_type index3, size_type index4,
			    size_type index5, size_type index6) const;
    
    const DataT& operator()(size_type index1, size_type index2,
			    size_type index3, size_type index4,
			    size_type index5) const;
    
    const DataT& operator()(size_type index1, size_type index2,
			    size_type index3, size_type index4) const;
    
    const DataT& operator()(size_type index1, size_type index2,
			    size_type index3) const;
    
    const DataT& operator()(size_type index1, size_type index2) const;
    
    const DataT& operator()(size_type index1) const;
    
    const DataT& operator[](size_type index) const;

    size_type rank() const;

    size_type dimension(size_type ord) const;

    void dimensions(std::vector<size_type>& dims);

    size_type size() const;

    void setFieldTag(const PHX::Tag<DataT>& t);
    
    void setFieldData(const Teuchos::ArrayRCP<DataT>& d);
    
    void print(std::ostream& os, bool printValues = false) const;

  private:
    
    PHX::Tag<DataT> m_tag;
    
    array_type m_field_data;

    Teuchos::ArrayRCP<DataT> m_array_rcp;

#ifdef PHX_DEBUG
    bool m_tag_set;
    bool m_data_set;
    static const std::string m_field_tag_error_msg;
    static const std::string m_field_data_error_msg;
#endif

  };
  
  template<typename DataT,
	   typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	   typename Tag4, typename Tag5, typename Tag6, typename Tag7>
  std::ostream& operator<<(std::ostream& os, 
			   const PHX::MDField<DataT, Tag0, Tag1, 
			   Tag2, Tag3, Tag4, Tag5, Tag6, Tag7>& h);
  
  // *************************************
  // Runtime time checked MDField
  // *************************************

  template<typename DataT>
  class MDField<DataT,void,void,void,void,void,void,void,void> {
    
  public:

    typedef DataT value_type;

    typedef typename shards::Array<DataT,shards::NaturalOrder,void,void,void,void,void,void,void,void> array_type;
    
    typedef typename array_type::size_type size_type;

    MDField(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& t);
    
    MDField(const PHX::Tag<DataT>& v);
    
    MDField();
    
    ~MDField();
    
    const PHX::FieldTag& fieldTag() const;

    DataT& operator()(size_type index1, size_type index2, size_type index3, 
		      size_type index4, size_type index5, size_type index6,
		      size_type index7, size_type index8);

    DataT& operator()(size_type index1, size_type index2, size_type index3, 
		      size_type index4, size_type index5, size_type index6,
		      size_type index7);

    DataT& operator()(size_type index1, size_type index2, size_type index3, 
		      size_type index4, size_type index5, size_type index6);
    
    DataT& operator()(size_type index1, size_type index2, size_type index3, 
		      size_type index4, size_type index5);
    
    DataT& operator()(size_type index1, size_type index2, size_type index3, 
		      size_type index4);
    
    DataT& operator()(size_type index1, size_type index2, size_type index3);
    
    DataT& operator()(size_type index1, size_type index2);
    
    DataT& operator()(size_type index1);
    
    DataT& operator[](size_type index);

    const DataT& operator()(size_type index1, size_type index2, 
			    size_type index3, size_type index4, 
			    size_type index5, size_type index6,
			    size_type index7, size_type index8) const;

    const DataT& operator()(size_type index1, size_type index2,
			    size_type index3, size_type index4,
			    size_type index5, size_type index6,
			    size_type index7) const;

    const DataT& operator()(size_type index1, size_type index2,
			    size_type index3, size_type index4,
			    size_type index5, size_type index6) const;
    
    const DataT& operator()(size_type index1, size_type index2,
			    size_type index3, size_type index4,
			    size_type index5) const;
    
    const DataT& operator()(size_type index1, size_type index2,
			    size_type index3, size_type index4) const;
    
    const DataT& operator()(size_type index1, size_type index2,
			    size_type index3) const;
    
    const DataT& operator()(size_type index1, size_type index2) const;
    
    const DataT& operator()(size_type index1) const;
    
    const DataT& operator[](size_type index) const;

    size_type rank() const;

    size_type dimension(size_type ord) const;

    void dimensions(std::vector<size_type>& dims);

    size_type size() const;

    void setFieldTag(const PHX::Tag<DataT>& t);
    
    void setFieldData(const Teuchos::ArrayRCP<DataT>& d);
    
    void print(std::ostream& os, bool printValues = false) const;

  private:
    
    PHX::Tag<DataT> m_tag;
    
    array_type m_field_data;

    Teuchos::ArrayRCP<DataT> m_array_rcp;

#ifdef PHX_DEBUG
    bool m_tag_set;
    bool m_data_set;
    static const std::string m_field_tag_error_msg;
    static const std::string m_field_data_error_msg;
#endif

  };
  
  template<typename DataT>
  std::ostream& operator<<(std::ostream& os, 
			   const PHX::MDField<DataT, void, void, 
			   void, void, void, void, void, void>& h);
  
} 

#include "Phalanx_MDField_Def.hpp"

#endif 
