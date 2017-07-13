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

#ifndef PANZER_BASIS_DESCRIPTOR_HPP
#define PANZER_BASIS_DESCRIPTOR_HPP

#include <string>
#include <functional>

namespace panzer {

class BasisDescriptor
{
public:

  /** \brief Constructor for empty basis
   *
   */
  BasisDescriptor();

  /** \brief Destructor
   *
   */
  virtual ~BasisDescriptor() = default;

  /** \brief Constructor for basis description
   *
   * \param[in] basis_order Basis order (e.g. 1 could be piecewise linear)
   * \param[in] basis_type Basis type (e.g. HGrad, HDiv, HCurl, ...)
   */
  BasisDescriptor(const int basis_order, const std::string & basis_type);

  /** \brief Get type of basis
   *
   * \return Type of basis
   */
  const std::string & getType() const {return _basis_type;}

  /** \brief Get order of basis
   *
   * \return Order of basis
   */
  int getOrder() const {return _basis_order;}

  /** \brief Get unique key associated with basis of this order and type
   *  The key is used to sort through a map of BasisDescriptors.
   *
   * \return Unique basis key
   */
  std::size_t getKey() const {return _key;}

protected:

  /// Basis type (HGrad, HDiv, HCurl,...)
  std::string _basis_type;

  // Basis order (>0)
  int _basis_order;

  // Unique key associated with basis.
  std::size_t _key;

};

}


namespace std {

template <>
struct hash<panzer::BasisDescriptor>
{
  std::size_t operator()(const panzer::BasisDescriptor& desc) const;
};

}


#endif
