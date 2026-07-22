// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_STK_WorksetFactory_hpp__
#define __Panzer_STK_WorksetFactory_hpp__

#include "Panzer_WorksetFactoryBase.hpp"
#include "Panzer_WorksetDescriptor.hpp"
#include "Panzer_WorksetNeeds.hpp"

#include "Panzer_STK_Interface.hpp"

namespace panzer
{
  struct LocalMeshInfo;
}

namespace panzer_stk {

/** Pure virtual base class used to construct 
  * worksets on volumes and side sets.
  */
class WorksetFactory : public panzer::WorksetFactoryBase {
public:
   WorksetFactory() {}

   WorksetFactory(const Teuchos::RCP<const STK_Interface> & mesh) : mesh_(mesh) {}

   virtual ~WorksetFactory() {}

  /** Set mesh
     */
   virtual
   void setMesh(const Teuchos::RCP<const panzer_stk::STK_Interface> & mesh);

   /** Build sets of boundary condition worksets
     */
   virtual
   Teuchos::RCP<std::map<unsigned,panzer::Workset> > 
   getSideWorksets(const panzer::WorksetDescriptor & desc,
                   const panzer::WorksetNeeds & needs) const;

   /** Build sets of boundary condition worksets for the BCT_Interface case.
     */
   virtual
   Teuchos::RCP<std::map<unsigned,panzer::Workset> >
   getSideWorksets(const panzer::WorksetDescriptor & desc,
                   const panzer::WorksetNeeds & needs_a,
                   const panzer::WorksetNeeds & needs_b) const;

   /** Build workssets specified by the workset descriptor.
     */
   virtual
   Teuchos::RCP<std::vector<panzer::Workset> >
   getWorksets(const panzer::WorksetDescriptor & worksetDesc,
               const panzer::WorksetNeeds & needs) const;

private:

   /// Mesh
   Teuchos::RCP<const STK_Interface> mesh_;

   // This needs to be set at the start, but is currently setup only if the
   // workset descriptor requiers it
   /// Alternative form of mesh
   mutable Teuchos::RCP<const panzer::LocalMeshInfo> mesh_info_;


};

}

#endif
