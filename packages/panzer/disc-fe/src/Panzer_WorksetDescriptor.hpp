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

#ifndef __Panzer_WorksetDescriptor_hpp__
#define __Panzer_WorksetDescriptor_hpp__

#include <string>
#include <ostream>
#include <functional>
#include "Panzer_HashUtils.hpp"

namespace panzer {

  //! Special values for the workset size. When the workset size is set on the WorksetDescriptor an interger value can be set, or one of these special enums can be used.
  enum WorksetSizeType : int {
    //! Workset size is set to the total number of local elements in the MPI process
    ALL_ELEMENTS=-1,
    //! Workset size is set to zero
    NO_ELEMENTS=0,
  };

/** Class provides a simple description of the types of worksets
 * that need to be constructed and used. This is a generalization of
 * using strings and pairs of string to represent the element blocks
 * and sidesets. It is primarily used in specifying the "domain" of the
 * assembly algorithm, that is which elements will be used in the assembly
 * process.
 *
 * There are currently two construction paths supported, a CLASSIC
 * path that is the legacy code that constructs worksets out of
 * locally owned elements. The new path uses the PARTITIONED
 * approach. The PARTITIONED approach supports creating worksets that
 * might require ghosted and/or virtual cells commonly used in DG
 * methods. For the CLASSIC mode, the workset size is determined by
 * the deprecated CellData object in the WorksetNeeds. This is the
 * current default for backwards compatibility. NOTE: For CLASSIC,
 * this means that the worksetSize in this object is ignored! For
 * PARTITIONED, used by DG codes, the worksetSize in the
 * WorksetDescriptor is always used.
 */
class
WorksetDescriptor
{
public:

  /// Default constructor
  WorksetDescriptor():
    worksetSize_(NO_ELEMENTS),
    include_ghosts_(false),
    group_by_subcell_(false),
    sideAssembly_(false),
    cascade_(false)
  {

  }

  /** Constructor specifying a lone element block with a requested workset size.
   *
   * \param[in] elementBlock Name of the element block
   * \param[in] worksetSize Requested workset size. This is an integer > 0 for a user specified size or can be of type WorksetSizeType for special cases
   */
  WorksetDescriptor(const std::string & elementBlock,
                    const int worksetSize=WorksetSizeType::ALL_ELEMENTS)
  : elementBlock_(elementBlock),
    worksetSize_(worksetSize),
    include_ghosts_(false),
    group_by_subcell_(false),
    sideAssembly_(false),
    cascade_(false)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(elementBlock_=="",std::runtime_error,
                                   "WorksetDescriptor constr: Element block name must be non-empty!");
  }

  /** Constructor that defines a side set. Note that the
   * specified sideset must be a non-empty string.
   *
   * Options for workset_size: ALL_ELEMENTS, NO_ELEMENTS, >0
   *   ALL_ELEMENTS -> workset size is set to largest possible value
   *   >0 -> workset size is set to this value (overwrites WorksetNeeds)
   *
   * \param[in] element_block Element block that includes the side
   * \param[in] sideset Side set that is being used
   * \param[in] worksetSize Requested workset size. This is an integer > 0 for a user specified size or can be of type WorksetSizeType for special cases
   */
  WorksetDescriptor(const std::string & elementBlock,
                    const std::string & sideset,
                    const int worksetSize=WorksetSizeType::ALL_ELEMENTS)
  : elementBlock_(elementBlock),
    sideset_(sideset),
    worksetSize_(worksetSize),
    include_ghosts_(false),
    group_by_subcell_(true),
    sideAssembly_(false),
    cascade_(false)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(elementBlock_=="",std::runtime_error,
                               "WorksetDescriptor constr: Element block name must be non-empty!");
    TEUCHOS_TEST_FOR_EXCEPTION(sideset_=="",std::runtime_error,
                               "WorksetDescriptor constr: Side set name must be non-empty!");
  }

  /** Constructor that defines a side set. Note that the
   * specified sideset must be a non-empty string.
   *
   * Options for workset_size: ALL_ELEMENTS, NO_ELEMENTS, >0
   *   ALL_ELEMENTS -> workset size is set to largest possible value
   *   >0 -> workset size is set to this value (overwrites WorksetNeeds)
   *
   * \param[in] element_block_0 Element block on one side of sideset_0
   * \param[in] element_block_1 Element block on other side of sideset_1
   * \param[in] sideset_0 Sideset of interest attached to element_block_0
   * \param[in] sideset_1 Sideset of interest attached to element_block_1
   * \param[in] worksetSize Requested workset size. This is an integer > 0 for a user specified size or can be of type WorksetSizeType for special cases
   *
   */
  WorksetDescriptor(const std::string & elementBlock_0,
                    const std::string & elementBlock_1,
                    const std::string & sideset_0,
                    const std::string & sideset_1,
                    const int worksetSize=WorksetSizeType::ALL_ELEMENTS)
  : elementBlock_(elementBlock_0),
    elementBlock_2_(elementBlock_1),
    sideset_(sideset_0),
    sideset_2_(sideset_1),
    worksetSize_(worksetSize),
    include_ghosts_(false),
    group_by_subcell_(true),
    sideAssembly_(true),
    cascade_(false)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(elementBlock_=="",std::runtime_error,
                               "WorksetDescriptor constr: Element block 0 name must be non-empty!");
    TEUCHOS_TEST_FOR_EXCEPTION(sideset_=="",std::runtime_error,
                               "WorksetDescriptor constr: Side set 0 name must be non-empty!");
    TEUCHOS_TEST_FOR_EXCEPTION(elementBlock_2_=="",std::runtime_error,
                               "WorksetDescriptor constr: Element block 1 name must be non-empty!");
    TEUCHOS_TEST_FOR_EXCEPTION(sideset_2_=="",std::runtime_error,
                               "WorksetDescriptor constr: Side set 1 name must be non-empty!");
  }

  void
  setElementBlock(const std::string & element_block)
  {
    elementBlock_ = element_block;
    elementBlock_2_ = "";
    TEUCHOS_TEST_FOR_EXCEPTION(elementBlock_=="",std::runtime_error,
                               "WorksetDescriptor::setElementBlock : Element block name must be non-empty!");
  }
  void
  setSideset(const std::string & sideset)
  {
    sideset_ = sideset;
    sideset_2_ = "";
  }
  void
  setWorksetSize(const int workset_size)
  {
    TEUCHOS_ASSERT(workset_size != NO_ELEMENTS);
    if(workset_size < 0){
      // Only one negative workset size is allowed
      TEUCHOS_ASSERT(workset_size == ALL_ELEMENTS);
    }
    worksetSize_ = workset_size;
  }
  void
  setInterface(const std::string & element_block,
               const std::string & sideset,
               const std::string & other_element_block,
               const std::string & other_sideset)
  {
    elementBlock_ = element_block;
    sideset_ = sideset;
    elementBlock_2_ = other_element_block;
    sideset_2_ = other_sideset;

    TEUCHOS_TEST_FOR_EXCEPTION(elementBlock_=="",std::runtime_error,
                               "WorksetDescriptor::setInterface : Element block name must be non-empty!");
    TEUCHOS_TEST_FOR_EXCEPTION(sideset_=="",std::runtime_error,
                               "WorksetDescriptor::setInterface : Sideset name must be non-empty!");
    TEUCHOS_TEST_FOR_EXCEPTION(elementBlock_2_=="",std::runtime_error,
                               "WorksetDescriptor::setInterface : Other element block name must be non-empty!");
    TEUCHOS_TEST_FOR_EXCEPTION(sideset_2_=="",std::runtime_error,
                               "WorksetDescriptor::setInterface : Other sideset name must be non-empty!");
  }
  void
  setIncludeGhosts(const bool include_ghosts)
  {
    include_ghosts_ = include_ghosts;
  }
  void
  setSideAssembly(const bool side_assembly)
  {
    sideAssembly_ = side_assembly;

    // Current limitation
    if(sideAssembly_){
      cascade_ = false;
      group_by_subcell_ = true;
    }
  }
  void
  setCascade(const bool cascade)
  {
    cascade_ = cascade;

    // Current limitation
    if(cascade_)
      sideAssembly_ = false;
  }
  void
  setGroupBySubcell(const bool group_by_subcell)
  {
    group_by_subcell_ = group_by_subcell;
  }

  /** \brief Get element block name
   *
   * Only set block if this is an interface between two element blocks
   *
   * \param[in] block Block id (0 or 1).
   * \return Name of element block
   */
  const std::string & getElementBlock(const int block=0) const
  { return (block==0) ? elementBlock_ : elementBlock_2_; }

  /** \brief Get sideset name
   *
   * Only set block if this is an interface between two element blocks
   *
   * \param[in] block Block id (0 or 1).
   * \return Name of sideset
   */
  const std::string & getSideset(const int block=0) const
  { return (block==0) ? sideset_ : sideset_2_; }

  //! Group worksets by subcell
  bool groupBySubcell() const
  {return group_by_subcell_;}

  //! Expects side set assembly on volume
  bool sideAssembly() const
  { return sideAssembly_; }

  /** \brief Identifies this workset as an interface between two element blocks
   *
   * \return True if element blocks exists on both sides of side
   */
  bool connectsElementBlocks() const
  {
    return useSideset() and elementBlock_2_ != "";
  }

  //! This descriptor is for a side set.
  bool useSideset() const
  //{ return useSideset_; }
  { return sideset_ != ""; }

  //! Build a cascade of worksets
  bool buildCascade() const
  { return cascade_; }

  //! Get the requested workset size (default -2 (workset size is set elsewhere), -1 (largest possible workset size), >0 (workset size))
  int getWorksetSize() const
  { return worksetSize_; }

  //! Check if workset should include ghost cells
  bool includeGhostCells() const
  {return include_ghosts_;}

private:

  //! Element block, required to be non-empty
  std::string elementBlock_;

  //! Element block on other side of side
  std::string elementBlock_2_;

  //! Side set, must be non-empty if <code>useSideset_</code> is true
  std::string sideset_;

  //! Side set on other side of side, must be non-empty if <code>useSideset_</code> is true and <code>elementBlock2_</code> is not empty
  std::string sideset_2_;

  //! Maximum number of owned cells in a workset
  int worksetSize_;

  //! Include ghost cells in worksets
  bool include_ghosts_;

  //! For sideset worksets, this call will group worksets by subcell dimension/index
  bool group_by_subcell_;

  /**
   * \brief This indicates if side quadrature rules are constructed
   * or volume rules are constructued. Ignored if useSideset_ is false.
   *
   * \note This is for backward compatability purposes - IntegrationDescriptor now handles this
   *
   */
  bool sideAssembly_;

  /**
   * \brief Request a cascade of worksets
   *
   * Cascade means that if you request a volume workset for a given block/sideset
   * you will get a list of worksets that are grouped by subcell dimension and subcell index
   *
   * i.e. If you have single 3D cell (think Hex) on the boundary of a sideset, a cascade will return 11 worksets:
   *
   *    Workset  | Subcell Dimension | Subcell Index
   *       0     |        0          |      0         (i.e. Node 0 of 8)
   *       1     |        0          |      1
   *       2     |        0          |      2
   *       3     |        0          |      3
   *       4     |        0          |      4
   *       5     |        1          |      0         (i.e. Edge 0 of 12)
   *       6     |        1          |      1
   *       7     |        1          |      2
   *       8     |        1          |      3
   *       9     |        1          |      4
   *      10     |        2          |      0         (i.e. Face 0 of 6)
   *
   * All of these above worksets could refer to the same cell.
   *
   * This is an advanced feature, and should only be used if it is understood.
   *
   */
  bool cascade_;

};

//! Equality operation for use with hash tables and maps
inline bool operator==(const WorksetDescriptor & a,const WorksetDescriptor & b)
{
  // Check all fields
  return (a.getElementBlock()==b.getElementBlock())
      && (a.getElementBlock(1)==b.getElementBlock(1))
      && (a.getSideset(1)==b.getSideset(1))
      && (a.getSideset()==b.getSideset())
      && (a.sideAssembly()==b.sideAssembly())
      && (a.useSideset()==b.useSideset())
      && (a.buildCascade()==b.buildCascade());
}

//! Hash function that satisifies the stl hash interface
inline std::size_t hash_value(const WorksetDescriptor & wd)
{
  std::size_t seed = 0;

  panzer::hash_combine(seed,wd.getElementBlock(0));
  if(wd.getElementBlock(1) != "")
    panzer::hash_combine(seed,wd.getElementBlock(1));
  if(wd.useSideset()) {

    panzer::hash_combine(seed,wd.getSideset(0));
    if(wd.getSideset(1) != "")
      panzer::hash_combine(seed,wd.getSideset(1));

    // optionally hash on side set and side assembly
    panzer::hash_combine(seed,wd.buildCascade());
    panzer::hash_combine(seed,wd.getSideset());
    panzer::hash_combine(seed,wd.sideAssembly());
  }

  return seed;
}

//! I/O utility
inline std::ostream & operator<<(std::ostream & os,const WorksetDescriptor & wd)
{
  if(wd.useSideset())
    os << "Side descriptor: "
       << "eblock = \"" << wd.getElementBlock() << "\", "
       << "ss = \"" << wd.getSideset() << "\", "
       << "side assembly = " << (wd.sideAssembly() ? "on" : "off")
       << "cascade = " << (wd.buildCascade() ? "on" : "off");
  else
    os << "Block descriptor: "
       << "eblock = \"" << wd.getElementBlock() << "\"";

  return os;
}

/**
 * \brief Builds a descriptor specifying an element block.
 *
 * \param[in] element_block Element block of cells
 * \param[in] workset_size Maximum number of owned cells in a workset (defaults to all)
 *
 * \return Descriptor
 */
inline
WorksetDescriptor
blockDescriptor(const std::string & element_block,
                const int workset_size = WorksetSizeType::ALL_ELEMENTS)
{
  return WorksetDescriptor(element_block,workset_size);
}

/**
 * \brief Builds a descriptor specifying an element block.
 *
 * \note These worksets will include ghost cells
 *
 * \param[in] element_block Element block of cells
 * \param[in] workset_size Maximum number of owned cells in a workset (defaults to all)
 *
 * \return Descriptor
 */
inline
WorksetDescriptor
blockGhostedDescriptor(const std::string & element_block,
                       const int workset_size = WorksetSizeType::ALL_ELEMENTS)
{
  WorksetDescriptor wd(element_block,workset_size);
  wd.setIncludeGhosts(true);
  return wd;
}

/**
 * \brief Builds a descriptor specifying cells on element block and sideset.
 *
 * \note This returns a volume assembled workset
 * \note This will only return cells that share a 'side' with the sideset (use sidesetCascadeDescriptor to get node+edge+side)
 *
 * \param[in] element_block Element block of cells
 * \param[in] sideset Sideset on which cells connect
 * \param[in] workset_size Maximum number of owned cells in a workset (defaults to all)
 *
 * \return Descriptor
 */
inline
WorksetDescriptor
sidesetDescriptor(const std::string & element_block,
                  const std::string & sideset,
                  const int workset_size = WorksetSizeType::ALL_ELEMENTS)
{
  WorksetDescriptor wd(element_block,sideset,workset_size);
  wd.setGroupBySubcell(true);
  wd.setSideAssembly(false);
  return wd;
}

/**
 * \brief Builds a descriptor specifying cells on element block and sideset.
 *
 * \note This returns a volume assembled workset
 * \note This will only return cells that share a 'side' with the sideset (use sidesetCascadeDescriptor to get node+edge+side)
 * \note These worksets will include ghost cells
 *
 * \param[in] element_block Element block of cells
 * \param[in] sideset Sideset on which cells connect
 * \param[in] workset_size Maximum number of owned cells in a workset (defaults to all)
 *
 * \return Descriptor
 */
inline
WorksetDescriptor
sidesetGhostedDescriptor(const std::string & element_block,
                         const std::string & sideset,
                         const int workset_size = WorksetSizeType::ALL_ELEMENTS)
{
  WorksetDescriptor wd(element_block,sideset,workset_size);
  wd.setGroupBySubcell(false); // <- this can be true or false
  wd.setSideAssembly(false);
  wd.setIncludeGhosts(true);
  return wd;
}

/**
 * \brief Builds a descriptor specifying cells sharing an interface between two element blocks
 *
 * \note This returns a volume assembled workset
 * \note This will only return cells that share a 'side' with the sideset
 *
 * \param[in] element_block Element block of cells (owned)
 * \param[in] other_element_block Element block of cells (owned and/or ghost)
 * \param[in] sideset Sideset on which cells connect
 * \param[in] workset_size Maximum number of owned cells in a workset (defaults to all)
 *
 * \return Descriptor
 */
inline
WorksetDescriptor
sidesetInterfaceDescriptor(const std::string & element_block,
                           const std::string & other_element_block,
                           const std::string & sideset,
                           const int workset_size = WorksetSizeType::ALL_ELEMENTS)
{
  WorksetDescriptor wd(element_block,other_element_block,sideset,sideset,workset_size);
  wd.setGroupBySubcell(true);
  wd.setSideAssembly(false);
  return wd;
}

/**
 * \brief Builds a descriptor specifying cells on element block and sideset separated by side.
 *
 * \note This is for backward compatibility only - use 'sidesetDescriptor' if possible
 *
 * \note This returns a side assembled workset
 * \note This will only return cells that share a 'side' with the sideset
 *
 * \param[in] element_block Element block of cells
 * \param[in] sideset Sideset on which cells connect
 * \param[in] workset_size Maximum number of owned cells in a workset (defaults to all)
 *
 * \return Descriptor
 */
inline
WorksetDescriptor
sidesetSideAssembledDescriptor(const std::string & element_block,
                               const std::string & sideset,
                               const int workset_size = WorksetSizeType::ALL_ELEMENTS)
{
  WorksetDescriptor wd(element_block,sideset,workset_size);
  wd.setGroupBySubcell(true);
  wd.setSideAssembly(true);
  return wd;
}

/**
 * \brief Builds a descriptor specifying cells sharing an interface between two element blocks
 *
 * \note This is for backward compatibility only - use 'sidesetInterfaceDescriptor' if possible
 *
 * \note This returns a side assembled workset
 * \note This will only return cells that share a 'side' with the sideset
 *
 * \param[in] element_block Element block of cells (owned)
 * \param[in] other_element_block Element block of cells (owned and/or ghost)
 * \param[in] sideset Sideset on which cells connect
 * \param[in] workset_size Maximum number of owned cells in a workset (defaults to all)
 *
 * \return Descriptor
 */
inline
WorksetDescriptor
sidesetInterfaceSideAssembledDescriptor(const std::string & element_block,
                                        const std::string & other_element_block,
                                        const std::string & sideset,
                                        const int workset_size = WorksetSizeType::ALL_ELEMENTS)
{
  WorksetDescriptor wd(element_block,other_element_block,sideset,sideset,workset_size);
  wd.setGroupBySubcell(true);
  wd.setSideAssembly(true);
  return wd;
}

/**
 * \brief Builds a cascade descriptor specifying an element block and a sideset.
 *
 * \note This returns a volume assembled workset descriptor
 * \note This will return cells that share a node+edge+side with the sideset
 *
 * \param[in] element_block Element block of cells
 * \param[in] sideset Sideset on which cells connect
 * \param[in] workset_size Maximum number of owned cells in a workset (defaults to all)
 *
 * \return Descriptor
 */
inline
WorksetDescriptor
sidesetCascadeDescriptor(const std::string & element_block,
                         const std::string & sideset,
                         const int workset_size = WorksetSizeType::ALL_ELEMENTS)
{
  WorksetDescriptor wd(element_block,sideset,workset_size);
  wd.setGroupBySubcell(true); // <- this probably doesn't do anything
  wd.setSideAssembly(false);
  wd.setCascade(true);
  return wd;
}

} // end panzer

namespace std {

template <>
struct hash<panzer::WorksetDescriptor>
{
  std::size_t operator()(const panzer::WorksetDescriptor& wd) const
  {
    return panzer::hash_value(wd);
  }
};

}

#endif
