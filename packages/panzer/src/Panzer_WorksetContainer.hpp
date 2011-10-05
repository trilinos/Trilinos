#ifndef __Panzer_WorksetContainer_hpp__
#define __Panzer_WorksetContainer_hpp__

#include "Teuchos_RCP.hpp"

#include "Panzer_InputPhysicsBlock.hpp"
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
     * \param[in] eb_to_ipb Map from "element blocks" to "input physics block" objects
     * \param[in] wkstSz Number of elements in a workset built by this container
     */ 
   WorksetContainer(const Teuchos::RCP<const WorksetFactoryBase> & factory,
                    const std::map<std::string,InputPhysicsBlock> & ebToIpb,
                    std::size_t wkstSz);

   /** Copies the workset factory, the InputPhysicsBlock map, and the workset size,
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

   //! The input physics block map
   void setInputPhysicsBlockMap(const std::map<std::string,InputPhysicsBlock> & ebToIpb)
   { ebToIpb_ = ebToIpb; }   
  
   //! get a reference to the input physics block map
   const std::map<std::string,InputPhysicsBlock> & getInputPhysicsBlockMap() const
   { return ebToIpb_; }

   /** Clear all allocated worksets, maintain the workset factory and element to physics
     * block map.
     */ 
   void clear();

   //! Look up an input physics block, throws an exception if it can be found.
   const InputPhysicsBlock & lookupInputPhysicsBlock(const std::string & eBlock) const;

   //! Access, and construction of volume worksets
   Teuchos::RCP<std::vector<Workset> > getVolumeWorksets(const std::string & eBlock);
 
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
   std::map<std::string,InputPhysicsBlock> ebToIpb_; //! Maps element blocks to input physics block objects

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
