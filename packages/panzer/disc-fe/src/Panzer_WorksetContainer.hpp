// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_WorksetContainer_hpp__
#define __Panzer_WorksetContainer_hpp__

#include "Teuchos_RCP.hpp"

#include "Intrepid2_Orientation.hpp"

#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_WorksetFactoryBase.hpp"
#include "Panzer_WorksetDescriptor.hpp" // what the workset is defined over
#include "Panzer_WorksetNeeds.hpp"      // whats in a workset basis/integration rules

namespace panzer {

// forward declaration
class GlobalIndexer;

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

   /** Instantiate a workset object with a specified factory and input workset needs
     * map.
     *
     * \param[in] factory Factory to be used for constructing worksets
     * \param[in] needs Workset needs mapped from the element blocks
     *                  (integration rules and basis values for each element block)
     */ 
   WorksetContainer(const Teuchos::RCP<const WorksetFactoryBase> & factory,
                    const std::map<std::string,WorksetNeeds> & needs = {});

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

   //! set the workset size
   void setWorksetSize(std::size_t worksetSize)
   { worksetSize_ = worksetSize; }

   //! get the workset size
   std::size_t getWorksetSize() const
   { return worksetSize_; }

   /** Set the needs for an element block. Note any old needs will
     * be erased, and all worksets will be deleted. This is equivalent
     * to calling <code>clear();</code> then <code>setNeeds()</code>. This clearing
     * is required in order to for the worksets to remain consistent.
     *
     * \param[in] eBlock Element block associated with the needs object
     * \param[in] needs Needs for the element block.
     */
   void setNeeds(const std::string & eBlock,const WorksetNeeds & needs);

   /** Clear all allocated worksets, maintain the workset factory and
     * element to physics block map.
     */ 
   void clear();

   /** Clear all allocated worksets for the volume assembly, maintain
     * the workset factory, element to physics block map and side
     * worksets. This is meant for cases where volume assembly is
     * needed during setup, but during a run only side worksets are
     * needed to adjust boundary conditions. This allows for
     * significant memory savings.
     */
   void clearVolumeWorksets();

   /** Clear all allocated worksets for the side assembly, maintain
     * the workset factory, element to physics block map and volume
     * worksets.
     */
   void clearSideWorksets();

   //! Look up an input physics block, throws an exception if it can not be found.
   const WorksetNeeds & lookupNeeds(const std::string & eBlock) const;

   //! Access to volume worksets
   Teuchos::RCP<std::vector<Workset> > getWorksets(const WorksetDescriptor & wd);

   //! Access, and construction of side worksets
   Teuchos::RCP<std::map<unsigned,Workset> > getSideWorksets(const WorksetDescriptor & desc);

   /** Set the global indexer. This is used solely for accessing the
     * orientations.
     */
   void setGlobalIndexer(const Teuchos::RCP<const panzer::GlobalIndexer> & ugi);

   /** Add a basis to the worksets (if required). If reuqired this will clear
     * the workset reconstructing all the arrays. Add to all element blocks.
     */
   void addBasis(const std::string & type,int order,const std::string & rep_field);

   /** Get the cell orientations used to build the basis values objects.
     */
   Teuchos::RCP<const std::vector<Intrepid2::Orientation> >  getOrientations() const
   { return orientations_; }

private:
   /** Set the orientations. Can only be called once, this also sets the internally stored
     * global indexer. If an exception is raised, saying it wasn't null then this method
     * has been previously called.
     */
   void applyOrientations(const Teuchos::RCP<const panzer::GlobalIndexer> & ugi);

   /** Using the stored global indexer, set the orientations for a volume workset on a
     * specified element block.
     */
   void applyOrientations(const std::vector<Intrepid2::Orientation> & orientations, 
                          const std::string & eBlock,
                          std::vector<Workset> & worksets) const;

   /** Using the stored global indexer, set the orientations for a side workset.
     */
   void applyOrientations(const std::vector<Intrepid2::Orientation> & orientations, 
                          const WorksetDescriptor & desc,
                          std::map<unsigned,Workset> & worksets) const;

   typedef std::unordered_map<WorksetDescriptor,Teuchos::RCP<std::vector<Workset> > > WorksetMap;
   typedef std::unordered_map<WorksetDescriptor,Teuchos::RCP<std::map<unsigned,Workset> > > SideMap;

   /** Using the stored global indexer, set the orientations for a volume workset on a
     * specified element block.
     */
   void applyOrientations(const std::string & eBlock,std::vector<Workset> & worksets) const;

   /** Using the stored global indexer, set the orientations for a side workset.
     */
   void applyOrientations(const WorksetDescriptor & desc,std::map<unsigned,Workset> & worksets) const;

   /** Set all the workset identifier in a vector.
     *
     * \param[in] wd       Workset descriptor, this defines the base point for the identifiers
     * \param[in] worksets Set the unique identifiers on these worksets
     */
   void setIdentifiers(const WorksetDescriptor & wd,std::vector<Workset> & worksets);

   /** Set all the workset identifier in a vector.
     *
     * \param[in] wd       Workset descriptor, this defines the base point for the identifiers
     * \param[in] worksets Set the unique identifiers on these worksets
     */
   void setIdentifiers(const WorksetDescriptor & wd,std::map<unsigned,Workset> & wkstMap);

   /** Check if WorksetNeeds have been registered
    *
    * \return True if there are element blocks with required needs objects
    */
   bool hasNeeds() const
   {return ebToNeeds_.size() > 0;}

   Teuchos::RCP<const WorksetFactoryBase> wkstFactory_;      //! How to construct worksets
   std::map<std::string,WorksetNeeds> ebToNeeds_; //! Maps element blocks to input physics block objects

   WorksetMap worksets_;
   SideMap sideWorksets_;

   std::size_t worksetSize_;

   Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer_;

   Teuchos::RCP<std::vector<Intrepid2::Orientation> >  orientations_;
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
