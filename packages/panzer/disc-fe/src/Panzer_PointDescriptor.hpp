// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_PointDescriptor_hpp__
#define __Panzer_PointDescriptor_hpp__

#include <string>

#include "Teuchos_RCP.hpp"

// forward declarations
namespace panzer {
class PointGenerator;
}

// class declarations
namespace panzer {

class PointDescriptor
{
public:

  /// Default constructor, no version
  PointDescriptor() = delete;
 
  /** \brief Constructor for the point descriptor. 
    *
    * \param[in] type String that defines the "type" of this point descriptor,
    *                 used to generate unique hashes
    * \param[in] generator PointGenerator object for the points.
    */
  PointDescriptor(const std::string & type,const Teuchos::RCP<PointGenerator> & generator);

  /// Destructor
  virtual ~PointDescriptor() = default;

  /**
   * \brief Check if the point descriptor has a generator for generating point values
   */
  bool
  hasGenerator() const
  {return _generator != Teuchos::null;}

  /** Build a generator class that generates any reference points on 
    * a specified topology.
    *
    * \param[in] The cell topology to build the coordinates on
    */
  const PointGenerator & getGenerator() const { return *_generator; }

  /** \brief Get unique string associated with the type of point
   *         descriptor. This will be used generate a hash
   *        to sort through a map of PointDescriptors.
   *
   * \return  A string that uniquely describes this point descriptor
   */
  const std::string & getType() const { return _type; }

  /** \brief Get unique key associated with integrator of this order and type
   *  The key is used to sort through a map of IntegrationDescriptors.
   *
   * \return Unique basis key
   */
  std::size_t getKey() const {return _key;}

protected:

  /** \brief Setup the point descriptor. Protected and used by constructors.
    *
    * \param[in] type String that defines the "type" of this point descriptor,
    *                 used to generate unique hashes
    * \param[in] generator PointGenerator object for the points.
    */
  void setup(const std::string & type,const Teuchos::RCP<PointGenerator> & generator);

  /// Type string
  std::string _type;

  /// Unique key associated with integrator
  std::size_t _key;

  /// PointGenerator object to build the points
  Teuchos::RCP<PointGenerator> _generator;
};

}


namespace std {

template <>
struct hash<panzer::PointDescriptor>
{
  std::size_t operator()(const panzer::PointDescriptor& desc) const;
};

}


#endif
