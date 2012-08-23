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

#ifndef Panzer_STK_MeshFactory_hpp__
#define Panzer_STK_MeshFactory_hpp__

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterListAcceptorDefaultBase.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include "Panzer_STK_PeriodicBC_Parser.hpp"

namespace panzer_stk {

class STK_Interface;

/** Pure virtual interface that constructs a 
  * STK_Mesh interface object.
  */
class STK_MeshFactory : public Teuchos::ParameterListAcceptorDefaultBase {
public:
   /** Construct a STK_Inteface object described
     * by this factory.
     *
     * \param[in] parallelMach Descriptor for machine to build this mesh on.
     *
     * \returns Pointer to <code>STK_Interface</code> object with 
     *          <code>isModifiable()==false</code>.
     */ 
   virtual Teuchos::RCP<STK_Interface> buildMesh(stk::ParallelMachine parallelMach) const = 0;

   /** This builds all the meta data of the mesh. Does not call metaData->commit.
     * Allows user to add solution fields and other pieces. The mesh can be "completed"
     * by calling <code>completeMeshConstruction</code>.
     */
   virtual Teuchos::RCP<STK_Interface> buildUncommitedMesh(stk::ParallelMachine parallelMach) const = 0;

   /** Finishes building a mesh object started by <code>buildUncommitedMesh</code>.
     */
   virtual void completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const = 0;

   /** Parse the periodic boundary condition parameter list and build a vector of periodic boundary
     * conditions (a convenience function)
     */
   static void parsePeriodicBCList(const Teuchos::RCP<Teuchos::ParameterList> & pl,
                                   std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > & periodicBC)
   {
      panzer_stk::PeriodicBC_Parser parser;
      parser.setParameterList(pl);
      periodicBC = parser.getMatchers();
   }

protected:
   // vector of periodic boundary condition objects
   std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > periodicBCVec_; 
};

}

#endif
