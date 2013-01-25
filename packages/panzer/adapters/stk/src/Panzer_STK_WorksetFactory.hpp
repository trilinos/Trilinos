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

#ifndef __Panzer_STK_WorksetFactory_hpp__
#define __Panzer_STK_WorksetFactory_hpp__

#include "Panzer_WorksetFactoryBase.hpp"

#include "Panzer_STK_Interface.hpp"

namespace panzer_stk {

/** Pure virtual base class used to construct 
  * worksets on volumes and side sets.
  */
class WorksetFactory : public panzer::WorksetFactoryBase {
public:
   WorksetFactory(const Teuchos::RCP<const STK_Interface> & mesh) : mesh_(mesh) {}
   virtual ~WorksetFactory() {}

   /** Build sets of volume worksets
     */ 
   virtual
   Teuchos::RCP<std::vector<panzer::Workset> >
   getVolumeWorksets(const std::string & eBlock,
                     const panzer::PhysicsBlock & pb) const;

   /** Build sets of boundary condition worksets
     */
   virtual
   Teuchos::RCP<std::map<unsigned,panzer::Workset> > 
   getSideWorksets(const panzer::BC & bc,
                 const panzer::PhysicsBlock & pb) const;

private:

   Teuchos::RCP<const STK_Interface> mesh_;

};

}

#endif
