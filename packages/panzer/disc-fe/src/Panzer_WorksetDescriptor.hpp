// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    //! Backwards compatibility mode that ignores the worksetSize in the WorksetDescriptor
    CLASSIC_MODE=-2,
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
class WorksetDescriptor {
public:

  /** Constructor specifying a lone element block with a requested workset size.
   *
   * \param[in] elementBlock Name of the element block
   * \param[in] worksetSize Requested workset size. This is an integer > 0 for a user specified size or can be of type WorksetSizeType for special cases
   * \param[in] requiresPartitioning If set to true, uses the new path for building worksets with partitioning
   * \param[in] applyOrientations If set to true, computes and applies orientations to relevant bases
   */
  WorksetDescriptor(const std::string & elementBlock,
                    const int worksetSize=WorksetSizeType::CLASSIC_MODE,
                    const bool requiresPartitioning=false,
                    const bool applyOrientations=true)
  : elementBlock_(elementBlock),
    worksetSize_(worksetSize),
    requiresPartitioning_(requiresPartitioning),
    applyOrientations_(applyOrientations),
    sideAssembly_(false)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(elementBlock_=="",std::runtime_error,
                                   "WorksetDescriptor constr: Element block name must be non-empty!");
  }

  /** Constructor that defines a side set. Note that the
   * specified sideset must be a non-empty string.
   *
   * \param[in] elementBlock Element block that includes the side
   * \param[in] sideset Side set that is being used
   * \param[in] sideAssembly Are integration rules and
   *                         basis functions evaluated on the
   *                         side or on the volume of the element.
   */
  WorksetDescriptor(const std::string & elementBlock,
                    const std::string & sideset,
                    const bool sideAssembly)
  : elementBlock_(elementBlock),
    sideset_(sideset),
    worksetSize_(CLASSIC_MODE),
    requiresPartitioning_(false),
    applyOrientations_(true),
    sideAssembly_(sideAssembly)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(elementBlock_=="",std::runtime_error,
                               "WorksetDescriptor constr: Element block name must be non-empty!");
    TEUCHOS_TEST_FOR_EXCEPTION(sideset_=="",std::runtime_error,
                               "WorksetDescriptor constr: Side set name must be non-empty!");
  }

  /** Constructor that defines a side set. Note that the
   * specified sideset must be a non-empty string.
   *
   * Options for workset_size: EMPTY, FULL, SPECIAL, >0
   *   EMPTY -> workset size is set by cellData in WorksetNeeds
   *   FULL -> workset size is set to largest possible value
   *   SPECIAL  -> Special case
   *   >0 -> workset size is set to this value (overwrites WorksetNeeds)
   *
   * \param[in] element_block Element block that includes the side
   * \param[in] sideset Side set that is being used
   * \param[in] worksetSize Requested workset size. This is an integer > 0 for a user specified size or can be of type WorksetSizeType for special cases
   * \param[in] requiresPartitioning If set to true, uses the new path for building worksets with partitioning
   * \param[in] applyOrientations If set to true, computes and applies orientations to relevant bases
   */
  WorksetDescriptor(const std::string & elementBlock,
                    const std::string & sideset,
                    const int worksetSize=WorksetSizeType::CLASSIC_MODE,
                    const bool requiresPartitioning=false,
                    const bool applyOrientations=true)
  : elementBlock_(elementBlock),
    sideset_(sideset),
    worksetSize_(worksetSize),
    requiresPartitioning_(requiresPartitioning),
    applyOrientations_(applyOrientations),
    sideAssembly_(false)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(elementBlock_=="",std::runtime_error,
                               "WorksetDescriptor constr: Element block name must be non-empty!");
    TEUCHOS_TEST_FOR_EXCEPTION(sideset_=="",std::runtime_error,
                               "WorksetDescriptor constr: Side set name must be non-empty!");
  }

  /** Constructor that defines a side set. Note that the
   * specified sideset must be a non-empty string.
   *
   * Options for workset_size: -2, -1, 0, >0
   *   -2 -> workset size is set by cellData in WorksetNeeds
   *   -1 -> workset size is set to largest possible value
   *   0  -> Special case
   *   >0 -> workset size is set to this value (overwrites WorksetNeeds)
   *
   * \param[in] element_block_0 Element block on one side of sideset_0
   * \param[in] element_block_1 Element block on other side of sideset_1
   * \param[in] sideset_0 Sideset of interest attached to element_block_0
   * \param[in] sideset_1 Sideset of interest attached to element_block_1
   * \param[in] worksetSize Requested workset size. This is an integer > 0 for a user specified size or can be of type WorksetSizeType for special cases
   * \param[in] requiresPartitioning If set to true, uses the new path for building worksets with partitioning
   * \param[in] applyOrientations If set to true, computes and applies orientations to relevant bases
   *
   */
  WorksetDescriptor(const std::string & elementBlock_0,
                    const std::string & elementBlock_1,
                    const std::string & sideset_0,
                    const std::string & sideset_1,
                    const int worksetSize=WorksetSizeType::CLASSIC_MODE,
                    const bool requiresPartitioning=false,
                    const bool applyOrientations=true)
  : elementBlock_(elementBlock_0),
    elementBlock_2_(elementBlock_1),
    sideset_(sideset_0),
    sideset_2_(sideset_1),
    worksetSize_(worksetSize),
    requiresPartitioning_(requiresPartitioning),
    applyOrientations_(applyOrientations),
    sideAssembly_(false)
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

  //! Copy constructor
  WorksetDescriptor(const WorksetDescriptor & src) = default;

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

  //! Expects side set assembly on volume
  //TEUCHOS_DEPRECATED
  bool sideAssembly() const
  { return sideAssembly_; }
//  { return useSideset(); }

  /** \brief Identifies this workset as an interface between two element blocks
   *
   * \return True if element blocks exists on both sides of side
   */
  bool connectsElementBlocks() const
  {
    return useSideset() and elementBlock_2_ != "";
  }

  /** \brief Do we need to partition the local mesh prior to generating worksets.
   *
   * Note that this is required if running surface integrals on a discontinuous discretization.
   *
   * \return True if partitioning is required
   */
  bool requiresPartitioning() const
  {
    return requiresPartitioning_;
  }

  //! This descriptor is for a side set.
  bool useSideset() const
  //{ return useSideset_; }
  { return sideset_ != ""; }

  //! Get the requested workset size (default -2 (workset size is set elsewhere), -1 (largest possible workset size), >0 (workset size))
  int getWorksetSize() const
  { return worksetSize_; }

  bool applyOrientations() const {return applyOrientations_;}

private:

  //! Element block, required to be non-empty
  std::string elementBlock_;

  //! Element block on other side of side
  std::string elementBlock_2_;

  //! Side set, must be non-empty if <code>useSideset_</code> is true
  std::string sideset_;

  //! Side set on other side of side, must be non-empty if <code>useSideset_</code> is true and <code>elementBlock2_</code> is not empty
  std::string sideset_2_;

  //! Requested workset size
  int worksetSize_;

  //! Marks if the mesh require partitioning before generating worksets
  bool requiresPartitioning_;

  //! Apply orientations - used for continuous discretizations with edge/face elements
  bool applyOrientations_;

  /** This indicates if side quadrature rules are constructed
   * or volume rules are constructued. Ignored if useSideset_
   * is false.
   */
  bool sideAssembly_;
};

//! Equality operation for use with hash tables and maps
inline bool operator==(const WorksetDescriptor & a,const WorksetDescriptor & b)
{
  if(a.useSideset())
    // if side set is in use, check all fields
    return    a.getElementBlock()==b.getElementBlock()
        && a.getSideset()==b.getSideset()
        && a.sideAssembly()==b.sideAssembly()
        && a.useSideset()==b.useSideset();
  else
    // otherwise check that both descriptor don't use side sets
    // and check the element block (the remaining fields are allowed
    // to be unset)
    return    a.getElementBlock()==b.getElementBlock()
        && a.useSideset()==b.useSideset();
}

//! Hash function that satisifies the stl hash interface
inline std::size_t hash_value(const WorksetDescriptor & wd)
{
  std::size_t seed = 0;

  panzer::hash_combine(seed,wd.getElementBlock());
  if(wd.useSideset()) {
    // optionally hash on side set and side assembly
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
    << "side assembly = " << (wd.sideAssembly() ? "on" : "off");
  else
    os << "Block descriptor: "
    << "eblock = \"" << wd.getElementBlock() << "\"";

  return os;
}

/** Builds a descriptor specifying an element block.
 */
//TEUCHOS_DEPRECATED
inline WorksetDescriptor blockDescriptor(const std::string & eBlock)
{ return WorksetDescriptor(eBlock); }

/** Builds a descriptor specifying a sideset, specify surface terms.
 */
//TEUCHOS_DEPRECATED
inline WorksetDescriptor sidesetDescriptor(const std::string & eBlock,const std::string & sideset)
{ return WorksetDescriptor(eBlock,sideset,false); }

/** Builds a descriptor specifying a sideset, however specify volumetric terms not
 * surface terms.
 */
//TEUCHOS_DEPRECATED
inline WorksetDescriptor sidesetVolumeDescriptor(const std::string & eBlock,const std::string & sideset)
{ return WorksetDescriptor(eBlock,sideset,true); }

}

namespace std {

template <>
struct hash<panzer::WorksetDescriptor>
{
  std::size_t operator()(const panzer::WorksetDescriptor& wd) const
  {
    std::size_t seed = 0;

    panzer::hash_combine(seed,wd.getElementBlock());
    panzer::hash_combine(seed,wd.requiresPartitioning());
    panzer::hash_combine(seed,wd.getWorksetSize());
    if(wd.useSideset()) {
      // optionally hash on side set and side assembly
      panzer::hash_combine(seed,wd.getSideset());
      panzer::hash_combine(seed,wd.sideAssembly());
    }

    return seed;
  }
};

}

#endif
