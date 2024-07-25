// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_WorksetFactoryBase_hpp__
#define __Panzer_WorksetFactoryBase_hpp__

#include <string>
#include <vector>
#include <map>

#include "Panzer_Workset.hpp"
#include "Panzer_WorksetDescriptor.hpp"
#include "Panzer_WorksetNeeds.hpp"

namespace panzer {

// Used to apply orientations to worksets constructed by factory
class OrientationsInterface;

/** Pure virtual base class used to construct 
  * worksets on volumes and side sets.
  */
class WorksetFactoryBase {
public:
   virtual ~WorksetFactoryBase() {}

   /** Build sets of boundary condition worksets for an interface case.
     */
   virtual
   Teuchos::RCP<std::map<unsigned,panzer::Workset> >
   getSideWorksets(const panzer::WorksetDescriptor & desc,
                   const panzer::WorksetNeeds & needs_a,
                   const panzer::WorksetNeeds & needs_b) const = 0;

   /** Build sets of boundary condition worksets
     */
   virtual
   Teuchos::RCP<std::map<unsigned,panzer::Workset> > 
   getSideWorksets(const panzer::WorksetDescriptor & desc,
		   const panzer::WorksetNeeds & needs) const = 0;

   /** Build workssets specified by the workset descriptor.
     */
   virtual
   Teuchos::RCP<std::vector<panzer::Workset> >
   getWorksets(const WorksetDescriptor & worksetDesc,
               const panzer::WorksetNeeds & needs) const = 0;

   /**
    * \brief Used to apply orientations to any bases added to the worksets
    *
    * \param[in] orientations Orientations object used to apply orientations to worksets
    */
   void
   setOrientationsInterface(const Teuchos::RCP<const panzer::OrientationsInterface> & orientations)
   {orientations_ = orientations;}

   /**
    * \brief Get the orientations associated with the worksets
    *
    * \return Orientations information
    */
   Teuchos::RCP<const OrientationsInterface>
   getOrientationsInterface() const
   {return orientations_;}

protected:

   /// Indexer used for applying orientations
   Teuchos::RCP<const OrientationsInterface> orientations_;
};

}

#endif
