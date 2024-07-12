// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_OrientationContainerBase_hpp__
#define __Panzer_OrientationContainerBase_hpp__

namespace panzer {

/** This class is used to access orientations and 
  * provides a degree of seperation between the
  * BasisValues objects and the global indexer (which
  * computes and stores the orientation). The particular
  * thing that this is does is avoids the need for the
  * WorksetContainer/Factory to know anything about
  * the local or global ordinal types.
  */
template <typename Scalar,typename Array>
class OrientationContainerBase {
public:
  virtual ~OrientationContainerBase() {}

  /** Get the orientations for a number of cell local ids. This will
    * be dependent on a particular basis.
    *
    * \param[in] cell_local_ids Cells to build orientations for.
    * \param[in] orientations Array of orientations (previously allocated)
    *                         to be filled with the orientations of a
    *                         particular basis.
    */
  virtual void getOrientations(const std::string & blockId,
                               const std::vector<std::size_t> & cell_local_ids,
                               Array & orientations) const = 0;

};

}

#endif
