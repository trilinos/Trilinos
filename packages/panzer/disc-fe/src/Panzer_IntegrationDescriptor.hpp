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


#ifndef PANZER_INTEGRATION_DESCRIPTOR_HPP
#define PANZER_INTEGRATION_DESCRIPTOR_HPP

#include <functional>

namespace panzer {

class IntegrationDescriptor
{
public:

  /** \brief Possible integration types
   *
   */
  enum {
    NONE,       /// No integral specified - default state
    VOLUME,     /// Integral over volume
    SURFACE,    /// Integral over all sides of cells (closed surface integral)
    SIDE,       /// Integral over a specific side of cells (side must be set)
    CV_VOLUME,  /// Control volume integral
    CV_SIDE,    /// Control volume side integral
    CV_BOUNDARY /// Control volume boundary integral (side must be set)
  };

  /** \brief Constructor for empty integrator
   *
   */
  IntegrationDescriptor();

  /// Destructor
  virtual ~IntegrationDescriptor() = default;

  /** \brief Constructor for integrator description
   *
   * \param[in] cubature_order Order of polynomial to integrate (e.g. 2 could integrate a quadratic equation)
   * \param[in] integration_type Type of integration (e.g. Volume integral, Surface integral, Side integral...)
   * \param[in] side Side of cell to integrate over (default to -1 -> ignore sides)
   */
  IntegrationDescriptor(const int cubature_order, const int integration_type, const int side=-1);

  /** \brief Get type of integrator
   *
   * \return Type of integrator
   */
  const int & getType() const {return _integration_type;}

  /** \brief Get order of integrator
   *
   * \return Order of integrator
   */
  const int & getOrder() const {return _cubature_order;}

  /** \brief Get side associated with integration - this is for backward compatibility
   *
   * \return Side of cell (= Subcell index)
   */
  const int & getSide() const {return _side;}

  /** \brief Get unique key associated with integrator of this order and type
   *  The key is used to sort through a map of IntegrationDescriptors.
   *
   * \return Unique basis key
   */
  std::size_t getKey() const {return _key;}

protected:

  /** \brief Setup function
   *
   * \param[in] cubature_order Order of polynomial to integrate (e.g. 2 could integrate a quadratic equation)
   * \param[in] integration_type Type of integration (e.g. Volume integral, Surface integral, Side integral...)
   * \param[in] side Side of cell to integrate over (default to -1 -> ignore sides)
   */
  void setup(const int cubature_order, const int integration_type, const int side=-1);

  /// Type of integration
  int _integration_type;

  /// Order of integration (Order of polynomial this integrator is designed for)
  int _cubature_order;

  /// Side associated with integration - this is for backward compatibility
  int _side;

  /// Unique key associated with integrator
  std::size_t _key;

};

}


namespace std {

template <>
struct hash<panzer::IntegrationDescriptor>
{
  std::size_t operator()(const panzer::IntegrationDescriptor& desc) const;
};

}


#endif
