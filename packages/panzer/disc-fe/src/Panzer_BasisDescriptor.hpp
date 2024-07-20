// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_BASIS_DESCRIPTOR_HPP
#define PANZER_BASIS_DESCRIPTOR_HPP

#include <string>
#include <functional>

#include "Panzer_PointDescriptor.hpp"

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
   * \param[in] basis_order Basis order as defined by Intrepid2 (e.g. 1 could be piecewise linear)
   * \param[in] basis_type Basis type (a string: "HGrad", "HDiv", "HCurl", or "HVol")
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

  /** \brief Build a point descriptor that builds reference points for
   *  the DOF locations. This method throws if no points exist for this
   *  basis.
   */
  PointDescriptor getPointDescriptor() const;

protected:
 
  /// Basis type (HGrad, HDiv, HCurl, HVol)
  std::string _basis_type;

  // Basis order (>0)
  int _basis_order;

  // Unique key associated with basis.
  std::size_t _key;
};

bool operator==(const panzer::BasisDescriptor& left,
                const panzer::BasisDescriptor& right);

} // namespace panzer

namespace std {
  template <>
  struct hash<panzer::BasisDescriptor>
  {
    std::size_t operator()(const panzer::BasisDescriptor& desc) const;
  };
}


#endif
