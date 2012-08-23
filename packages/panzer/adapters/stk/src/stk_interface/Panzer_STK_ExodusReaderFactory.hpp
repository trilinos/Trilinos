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

#ifndef Panzer_STK_ExodusReaderFactory_hpp__
#define Panzer_STK_ExodusReaderFactory_hpp__

#include <string>

#include "Panzer_STK_config.hpp"
#include "Panzer_STK_MeshFactory.hpp"

#ifdef HAVE_IOSS

#include <stk_io/MeshReadWriteUtils.hpp>

namespace panzer_stk {

class STK_Interface;

/** Concrete mesh factory instantiation. This reads
  * a mesh from an exodus file and builds a STK_Interface object.
  *
  * Also, if a nonzero restart index (the Exodus indices are 1 based) is
  * specified then this will set the initial state time in the STK_Interface.
  * However, as prescribed by that interface the currentStateTime is only
  * set by the writeToExodus call, thus the currentStateTime of the created
  * STK_Interface object will be zero. It is up to the user to rectify this
  * when calling writeToExodus.
  */
class STK_ExodusReaderFactory : public STK_MeshFactory {
public:

   STK_ExodusReaderFactory();

   STK_ExodusReaderFactory(const std::string & fileName,int restartIndex=0);

   /** Construct a STK_Inteface object described
     * by this factory.
     *
     * \param[in] parallelMach Descriptor for machine to build this mesh on.
     *
     * \returns Pointer to <code>STK_Interface</code> object with 
     *          <code>isModifiable()==false</code>.
     */ 
   virtual Teuchos::RCP<STK_Interface> buildMesh(stk::ParallelMachine parallelMach) const;

   /** This builds all the meta data of the mesh. Does not call metaData->commit.
     * Allows user to add solution fields and other pieces. The mesh can be "completed"
     * by calling <code>completeMeshConstruction</code>.
     */
   virtual Teuchos::RCP<STK_Interface> buildUncommitedMesh(stk::ParallelMachine parallelMach) const;

   /** Finishes building a mesh object started by <code>buildUncommitedMesh</code>.
     */
   virtual void completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const;

   //! From ParameterListAcceptor
   void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList);

   //! From ParameterListAcceptor
   Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

   //! Get file name mean is read from.
   const std::string & getFileName() const
   { return fileName_; }

protected:

   void registerElementBlocks(STK_Interface & mesh,stk::io::MeshData & meshData) const;
   void registerSidesets(STK_Interface & mesh,stk::io::MeshData & meshData) const;
   void registerNodesets(STK_Interface & mesh,stk::io::MeshData & meshData) const;

   std::string fileName_;
   int restartIndex_;
};

}

#endif
#endif
