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
