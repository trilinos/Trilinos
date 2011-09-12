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


#ifndef PHX_DATA_CONTAINER_HPP
#define PHX_DATA_CONTAINER_HPP

#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_FieldTag_STL_Functors.hpp"
#include "Phalanx_DataContainer_Base.hpp"

namespace PHX {

  /*! \brief Container that holds all fields associated with a specific DataT.
    
      One DataContainer is instantiated for each data type in each
      evaluation type.

  */
  template <typename DataT, typename Traits>
  class DataContainer : public PHX::DataContainerBase<Traits> {
    
  public:

    DataContainer() {}
    
    ~DataContainer() {}
    
    Teuchos::ArrayRCP<DataT> getFieldData(const PHX::FieldTag& t);
    
    void allocateField(const Teuchos::RCP<PHX::FieldTag>& t,
		       typename Traits::Allocator& a);

    const std::type_info& dataTypeInfo() const;

    std::size_t getSizeOfDataType() const;

    void print(std::ostream& os) const;
    
  private:
    
    typedef std::map< Teuchos::RCP<const PHX::FieldTag>, 
                      Teuchos::ArrayRCP<DataT>, 
                      FTComp > m_data_t;

    m_data_t m_data;
    
  };
  
  template <typename DataT, typename Traits>
  std::ostream& operator<<(std::ostream& os, 
			   const PHX::DataContainer<DataT, Traits>& dc)
  {
    dc.print(os);
    return os;
  }

} 

#include "Phalanx_DataContainer_Def.hpp"

#endif 
