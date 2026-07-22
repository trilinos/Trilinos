// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_Workset_Builder_hpp__
#define __Panzer_Workset_Builder_hpp__

#include <vector>
#include <map>

#include "Teuchos_RCP.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_WorksetNeeds.hpp"

namespace panzer {
  
  class Workset;
  class WorksetDetails;

  template<typename ArrayT>
  Teuchos::RCP<std::vector<Workset> > 
  buildWorksets(const WorksetNeeds & needs,
                const std::string & elementBlock,
                const std::vector<std::size_t>& local_cell_ids,
                const ArrayT& node_coordinates);

  template<typename ArrayT>
  Teuchos::RCP<std::map<unsigned,Workset> >
  buildBCWorkset(const WorksetNeeds & needs,
                 const std::string & elementBlock,
                 const std::vector<std::size_t>& local_cell_ids,
                 const std::vector<std::size_t>& local_side_ids,
                 const ArrayT& node_coordinates,
                 const bool populate_value_arrays = true);

  template<typename ArrayT>
  Teuchos::RCP<std::map<unsigned,panzer::Workset> >
  buildBCWorkset(const WorksetNeeds & needs_a,
                 const std::string & blockid_a,
                 const std::vector<std::size_t>& local_cell_ids_a,
                 const std::vector<std::size_t>& local_side_ids_a,
                 const ArrayT& node_coordinates_a,
                 const panzer::WorksetNeeds & needs_b,
                 const std::string & blockid_b,
                 const std::vector<std::size_t>& local_cell_ids_b,
                 const std::vector<std::size_t>& local_side_ids_b,
                 const ArrayT& node_coordinates_b);

  /** This routine supports construction of worksets that are
    * more DG like. The elements are assumed to shared an
    * edge (or face) and the side id is specified accordingly.
    * Note that no checking of "sharing" is done when the workset
    * is constructed.
    */
  template<typename ArrayT>
  Teuchos::RCP<std::vector<Workset> > 
  buildEdgeWorksets(const WorksetNeeds & needs_a,
                   const std::string & eblock_a,
                   const std::vector<std::size_t>& local_cell_ids_a,
                   const std::vector<std::size_t>& local_side_ids_a,
                   const ArrayT& node_coordinates_a,
                   const WorksetNeeds & needs_b,
                   const std::string & eblock_b,
                   const std::vector<std::size_t>& local_cell_ids_b,
                   const std::vector<std::size_t>& local_side_ids_b,
                   const ArrayT& node_coordinates_b);

  template<typename ArrayT>
  std::vector<Workset>::iterator
  buildEdgeWorksets(const std::vector<std::size_t> & cell_indices,
                    const WorksetNeeds & needs_a,
                    const std::string & eblock_a,
                    const std::vector<std::size_t>& local_cell_ids_a,
                    const std::vector<std::size_t>& local_side_ids_a,
                    const ArrayT& node_coordinates_a,
                    const WorksetNeeds & needs_b,
                    const std::string & eblock_b,
                    const std::vector<std::size_t>& local_cell_ids_b,
                    const std::vector<std::size_t>& local_side_ids_b,
                    const ArrayT& node_coordinates_b,
                    std::vector<Workset>::iterator beg);

  /** Populate basis values and integration values data structures in
    * the WorksetDetails object being passed in. Note that this works for
    * both edge data structures and for volumetric data structures. Note
    * that for this function to work the WorksetDetails object must already
    * be populated with coordinates.
    *
    * \param[in] num_cells The number of cells contained in the workset
    * \param[in] isSide    Populate using a side subcell
    * \param[in] pb        Physics block that contains the integration rules
    *                      and basis functions
    * \param[in,out] details Object to be populated already containing the
    *                        node coordinates for the cells.
    *
    * \note This function is primarily for internal use. A user should not
    *       call it directly.
    */
  void populateValueArrays(std::size_t num_cells,
                           bool isSide,
                           const WorksetNeeds & pb,
                           WorksetDetails & details,
                           const Teuchos::RCP<WorksetDetails> other_details = Teuchos::null);
}

#endif
