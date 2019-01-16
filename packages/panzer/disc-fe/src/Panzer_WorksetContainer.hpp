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

#ifndef __Panzer_WorksetContainer_hpp__
#define __Panzer_WorksetContainer_hpp__

// Trilinos includes
#include "Teuchos_RCP.hpp"

// Panzer includes
#include "Panzer_WorksetDescriptor.hpp"
#include "Panzer_WorksetFactoryBase.hpp"

// STL includes
#include <unordered_map>

namespace panzer
{

/** \brief Class that provides access to worksets on
  * each element block and side set.
  *
  * This class provides access to worksets on each 
  * element block and side set. This is done using an
  * optional lazy construction mechanism that builds the
  * worksets in a just in time fashion. Because the specifics
  * of a workset is constructed are based on the type of 
  * mesh database, each new implementation must inherit
  * from the <code>WorksetFactoryBase</code> class. This
  * class will then use that one to handle the lazy evaluation.
  */
class
WorksetContainer
{
public:

  // Default constructor with no factory
  WorksetContainer() = default;

  /**
   * \brief Base constructor that sets the workset factory
   *
   * \param[in] factory Workset factory with which to generate worksets
   */
  WorksetContainer(const Teuchos::RCP<const WorksetFactoryBase> & factory);

  /**
   * \brief Copy constructor
   *
   * \note This operation does not copy worksets - only the factory
   *
   * \param[in] wc Container to copy from
   */
  WorksetContainer(const WorksetContainer & wc);

  /// Default destructor
  ~WorksetContainer() = default;

  /**
   * \brief Set the workset factory
   *
   * \note A factory must be set before calling 'getWorksets'
   *
   * \param[in] factory Factory for generating worksets
   */
  void
  setFactory(const Teuchos::RCP<const WorksetFactoryBase> & factory);

  /**
   * \brief Get the factory associated with this container
   *
   * \return Workset factory - returns Teuchos::null if not set
   */
  Teuchos::RCP<const WorksetFactoryBase>
  getFactory() const;

  /**
   * \brief Get worksets associated with a given workset descriptor
   *
   * \throws If factory has not been set
   *
   * \note If you just want worksets and don't want to store then call getFactory()->getWorksets(description).
   *
   * \return Pointer to list of worksets associated with descriptor
   */
  Teuchos::RCP<std::vector<Workset> >
  getWorksets(const WorksetDescriptor & description);

  /// Clear worksets currently in container - does not clear factory
  void
  clearWorksets();

protected:

  /// Definition of worket map
  using WorksetMap = std::unordered_map<WorksetDescriptor,Teuchos::RCP<std::vector<Workset> > >;

  /// Factory used to generate worksets
  Teuchos::RCP<const WorksetFactoryBase> factory_;

  /// Workset storage
  WorksetMap worksets_;

};

}

#endif
