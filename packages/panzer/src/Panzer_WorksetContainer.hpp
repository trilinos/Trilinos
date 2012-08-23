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

#include "Teuchos_RCP.hpp"

#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_WorksetFactoryBase.hpp"

namespace panzer {

struct SideId {
   SideId(const BC & bc)
      : ss_id(bc.sidesetID()), eblk_id(bc.elementBlockID()) {}

   std::string ss_id;
   std::string eblk_id;
};

/** Required to distinguish between boundary conditions
  * and sides.
  */
struct LessSide {
   bool operator()(const SideId & left, 
                   const SideId  right) const
   { return   (left.ss_id+"_"+left.eblk_id 
            < right.ss_id+"_"+right.eblk_id); }
};

/** \brief Class that provides access to worksets on
  * each element block and side set.
  *
  * This class provides access to worksets on each 
  * element block and side set. This is done using an
  * optional lazy construction mechnism that builds the
  * worksets in a just in time fashion. Because the specifics
  * of a workset is constructed are based on the type of 
  * mesh database, each new implementation must inherit
  * from the <code>WorksetFactoryBase</code> class. This
  * class will then use that one to handle the lazy evaluation.
  */
class WorksetContainer {
public:
   //! Default contructor, starts with no workset factory objects
   WorksetContainer();

   /** Instantiate a workset object with a specified factory and input physics block
     * map.
     *
     * \param[in] factory Factory to be used for constructing worksets
     * \param[in] physicsBlocks Vector of physics blocks
     * \param[in] wkstSz Number of elements in a workset built by this container
     */ 
   WorksetContainer(const Teuchos::RCP<const WorksetFactoryBase> & factory,
                    const std::vector<Teuchos::RCP<PhysicsBlock> > & physicsBlocks,
                    std::size_t wkstSz);

   /** Copies the workset factory, the PhysicsBlock vector, and the workset size,
     * but not constructed worksets.
     */
   WorksetContainer(const WorksetContainer & wc);

   /** Set the workset factory, and as a consequence clear out all previously computed
     * worksets.
     */ 
   void setFactory(const Teuchos::RCP<const WorksetFactoryBase> & factory)
   { clear(); wkstFactory_ = factory; }

   //! Access the workset factory pointer.
   Teuchos::RCP<const WorksetFactoryBase> getFactory() const
   { return wkstFactory_; }

   //! The physics block vector
   void setPhysicsBlockVector(const std::vector<Teuchos::RCP<PhysicsBlock> > & physicsBlocks);

   //! set the workset size
   void setWorksetSize(std::size_t worksetSize)
   { worksetSize_ = worksetSize; }

   //! get the workset size
   std::size_t getWorksetSize() const
   { return worksetSize_; }

   /** Clear all allocated worksets, maintain the workset factory and element to physics
     * block map.
     */ 
   void clear();

   //! Look up an input physics block, throws an exception if it can be found.
   const PhysicsBlock & lookupPhysicsBlock(const std::string & eBlock) const;

   //! Access to volume worksets
   Teuchos::RCP<std::vector<Workset> > getVolumeWorksets(const std::string & eBlock);

   //! Access, and construction of volume worksets
   Teuchos::RCP<std::vector<Teuchos::RCP<std::vector<Workset> > > > getVolumeWorksets() const;
 
   //! Access, and construction of side worksets
   Teuchos::RCP<std::map<unsigned,Workset> > getSideWorksets(const BC & bc);

   //! Iterator access to volume worksets
   inline std::vector<Workset>::iterator begin(const std::string & eBlock)
   { return getVolumeWorksets(eBlock)->begin(); }

   //! Iterator access to volume worksets
   inline std::vector<Workset>::iterator end(const std::string & eBlock)
   { return getVolumeWorksets(eBlock)->end(); }

   //! Iterator access to side worksets
   inline std::map<unsigned,Workset>::iterator begin(const BC & bc)
   { return getSideWorksets(bc)->begin(); }

   //! Iterator access to side worksets
   inline std::map<unsigned,Workset>::iterator end(const BC & bc)
   { return getSideWorksets(bc)->end(); }

   /** Allocate worksets associated with the element blocks in a vector, this
     * will overwrite an previously constructed worksets.
     */
   void allocateVolumeWorksets(const std::vector<std::string> & eBlocks);

   /** Allocate worksets associated with the BC objects in a vector, this
     * will overwrite an previously constructed worksets.
     */
   void allocateSideWorksets(const std::vector<BC> & bcs);

private:
   typedef std::map<std::string,Teuchos::RCP<std::vector<Workset> > > VolumeMap;
   typedef std::map<SideId,Teuchos::RCP<std::map<unsigned,Workset> >,LessSide> SideMap;

   Teuchos::RCP<const WorksetFactoryBase> wkstFactory_;      //! How to construct worksets
   std::map<std::string,Teuchos::RCP<PhysicsBlock> > ebToPb_; //! Maps element blocks to input physics block objects

   VolumeMap volWorksets_;
   SideMap sideWorksets_;

   std::size_t worksetSize_;
};

/** Build a map of volume worksets from a list of element blocks. Note that this
  * may, if needed, allocate these worksets in the workset container object.
  *
  * \param[in] wc Workset container that contains, or will contain, the volume worksets.
  * \param[in] elementBlockNames Element blocks to build the worksets for.
  * \param[out] volumeWksts Map form element block ids to the volume worksets. This will
  *                         not be cleared but the worksets with shared element names 
  *                         will be overwritten.
  */
void getVolumeWorksetsFromContainer(WorksetContainer & wc,
                                    const std::vector<std::string> & elementBlockNames,
                                    std::map<std::string,Teuchos::RCP<std::vector<Workset> > > & volumeWksts);

/** Build a map of side worksets from a list of boundary conditions. Note that this
  * may, if needed, allocate these worksets in the workset container object.
  *
  * \param[in]  wc        Workset container that contains, or will contain, the volume worksets.
  * \param[in]  bcs       Boundary conditions to use.
  * \param[out] sideWksts Map form element block ids to the volume worksets. This will
  *                       not be cleared but the worksets with shared element names 
  *                       will be overwritten.
  */
void getSideWorksetsFromContainer(WorksetContainer & wc,
                                  const std::vector<BC> & bcs,
                                  std::map<BC,Teuchos::RCP<std::map<unsigned,Workset> >,LessBC> & sideWksts);

}

#endif
