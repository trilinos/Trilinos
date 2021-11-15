// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <adapt/RefinementInfoByType.hpp>
#include <stk_util/diag/PrintTable.hpp>
#include <percept/stk_mesh.hpp>
#include <adapt/Refiner.hpp>

  namespace percept {

    RefinementInfo::RefinementInfo(Refiner *ref) : m_refiner(ref)
    {
      m_totalNumElementsBeforeRefine = 0;
      m_totalNumSidesBeforeRefine    = 0;
      m_totalNumElementsAfterRefine  = 0;
      m_totalNumSidesAfterRefine     = 0;
      m_numMarkedSides               = 0;
      m_numMarkedElements            = 0;

      m_full_stats = false;
      if (ref->getMesh().getProperty("RefinementInfo.full_stats") == "true")
        m_full_stats = true;
      for (int i=0; i < 4; ++i)
        {
          m_numberOfMarkedSubDimEntities[i]   = 0;
          m_numberOfSubDimEntities[i]         = 0;
          m_totalNumEntitiesBeforeFilter[i]   = 0;
          m_totalNumEntitiesAfterFilter[i]    = 0;
          m_totalNumEntitiesBeforeEstimate[i] = 0;
          m_totalNumEntitiesAfterEstimate[i]  = 0;
        }
    }

    /** Estimate number of elements for each topology based on old/new saved in m_numOrigElems, m_numNewElems.
     *  Saves estimates in m_numOrigElemsLast, m_numNewElemsLast corresponding to iRefinePass.
     *  iRefinePass is used to predict number of new elements in query_only mode of UniformRefiner.
     *  Pass in 0 for iRefinePass unless you are using UniformRefiner in query_only mode, in which case
     *  pass in the current refinement pass
     */
    void RefinementInfo::estimateNew(int iRefinePass)
    {
      for (unsigned irank = 0; irank < m_refinementInfoByType.size(); irank++)
        {
          RefinementInfoCount numOrig = m_refinementInfoByType[irank].m_numOrigElems;
          RefinementInfoCount numNew = m_refinementInfoByType[irank].m_numNewElems;
          m_refinementInfoByType[irank].m_numOrigElemsLast = numOrig;
          m_refinementInfoByType[irank].m_numNewElemsLast = numNew;
          if (numOrig)
            {
              double refFactor = ((double)numNew)/((double)numOrig);
              double refFactorNew = std::pow(refFactor, ((double)(iRefinePass+1) ));
              double refFactorOld = std::pow(refFactor, ((double)(iRefinePass) ));
              numNew = (RefinementInfoCount)((double)numOrig * refFactorNew);
              numOrig = (RefinementInfoCount)((double)numOrig * refFactorOld);
              m_refinementInfoByType[irank].m_numOrigElemsLast = numOrig;
              m_refinementInfoByType[irank].m_numNewElemsLast = numNew;
            }
        }
    }

    /** iRefinePass is used to predict number of new elements in query_only mode of UniformRefiner.
     *  Pass in 0 for iRefinePass unless you are using UniformRefiner in query_only mode, in which case
     *  pass in the current refinement pass
     */
    void RefinementInfo::printTable(std::ostream& os, int iRefinePass, bool printAll, const std::string& extra_title)
    {
      if (m_refinementInfoByType.size() == 0) return;

      RefinementInfoCount numOrigTot = 0;
      RefinementInfoCount numNewTot = 0;
      estimateNew(iRefinePass);
      for (unsigned irank = 0; irank < m_refinementInfoByType.size(); irank++)
        {
          RefinementInfoCount numOrigLast = m_refinementInfoByType[irank].m_numOrigElemsLast;
          RefinementInfoCount numNewLast = m_refinementInfoByType[irank].m_numNewElemsLast;
          numOrigTot += numOrigLast;
          numNewTot += numNewLast;
        }

      RefinementInfoCount numOrigNodes = m_refinementInfoByType[0].m_numOrigNodes;
      RefinementInfoCount numNewNodes = m_refinementInfoByType[0].m_numNewNodes;

      stk::PrintTable table;
      table.setTitle("Refinement Info: "+extra_title+"\n");

      table << "|" << "                     " <<  "|" << stk::justify(stk::PrintTable::Cell::CENTER)
        << "Original" << stk::end_col << "     " << "|" << "New     " << stk::end_col << "     " << "|" << stk::end_header;
      table << stk::justify(stk::PrintTable::Cell::LEFT) ;
      table << "|" << "Element Topology Type" <<  "|" << "Elements" << stk::end_col << "Nodes" << "|" << "Elements" << stk::end_col << "Nodes" << "|" << stk::end_header;

      for (unsigned irank = 0; irank < m_refinementInfoByType.size(); irank++)
        {
          if (!printAll && m_refinementInfoByType[irank].m_numOrigElems == 0)
            continue;

          RefinementInfoCount numOrig = m_refinementInfoByType[irank].m_numOrigElemsLast;
          RefinementInfoCount numNew = m_refinementInfoByType[irank].m_numNewElemsLast;

          table << "|" << m_refinementInfoByType[irank].m_topology.getName() << "|"
                << numOrig << stk::end_col << " " << "|"
                << numNew  << stk::end_col << " " << "|"
                << stk::end_row;
        }

      table << "|" << "Totals" << "|"
            << numOrigTot << stk::end_col << numOrigNodes << "|"
            << numNewTot  << stk::end_col << numNewNodes << "|"
            << stk::end_row;

      os << "\n" << table;

      if (m_full_stats)
        {
          stk::PrintTable table2;
          table2.setTitle("Full Stats\n");

          table2 << "|" << " what " << "|" << "       total    " << "|" << "     refined or marked   " << "|" << "%" << stk::end_header;

          table2 << "|" << " #elems / marked           " << "|" << m_totalNumElementsBeforeRefine << "|" << m_numMarkedElements << "|"
                 << 100.0*(double(m_numMarkedElements)/double(m_totalNumElementsBeforeRefine?m_totalNumElementsBeforeRefine:1)) << stk::end_row;
          table2 << "|" << " #sides / marked           " << "|" << m_totalNumSidesBeforeRefine    << "|" << m_numMarkedSides    << "|"
                 << 100.0*(double(m_numMarkedSides)/double(m_totalNumSidesBeforeRefine?m_totalNumSidesBeforeRefine:1)) << stk::end_row;
          table2 << "|" << " marked_edges     " << "|" << m_numberOfSubDimEntities[1]    << "|" << m_numberOfMarkedSubDimEntities[1]    << "|"
                 << 100.0*(double(m_numberOfMarkedSubDimEntities[1])/double(m_numberOfSubDimEntities[1]?m_numberOfSubDimEntities[1]:1)) << stk::end_row;
          table2 << "|" << " marked_faces     " << "|" << m_numberOfSubDimEntities[2]    << "|" << m_numberOfMarkedSubDimEntities[2]    << "|"
                 << 100.0*(double(m_numberOfMarkedSubDimEntities[2])/double(m_numberOfSubDimEntities[2]?m_numberOfSubDimEntities[2]:1)) << stk::end_row;
          table2 << "|" << " marked_centroids " << "|" << m_numberOfSubDimEntities[3]    << "|" << m_numberOfMarkedSubDimEntities[3]    << "|"
                 << 100.0*(double(m_numberOfMarkedSubDimEntities[3])/double(m_numberOfSubDimEntities[3]?m_numberOfSubDimEntities[3]:1)) << stk::end_row;
          table2 << "|" << " filtered elems   " << "|" << m_totalNumEntitiesBeforeFilter[3]    << "|" << m_totalNumEntitiesAfterFilter[3]    << "|"
                 << 100.0*(double(m_totalNumEntitiesAfterFilter[3])/double(m_totalNumEntitiesBeforeFilter[3]?m_totalNumEntitiesBeforeFilter[3]:1)) << stk::end_row;
          table2 << "|" << " filtered sides   " << "|" << m_totalNumEntitiesBeforeFilter[2]    << "|" << m_totalNumEntitiesAfterFilter[2]    << "|"
                 << 100.0*(double(m_totalNumEntitiesAfterFilter[2])/double(m_totalNumEntitiesBeforeFilter[2]?m_totalNumEntitiesBeforeFilter[2]:1)) << stk::end_row;
          table2 << "|" << " estimated elems   " << "|" << m_totalNumEntitiesBeforeEstimate[3]    << "|" << m_totalNumEntitiesAfterEstimate[3]    << "|"
                 << 100.0*(double(m_totalNumEntitiesAfterEstimate[3])/double(m_totalNumEntitiesBeforeEstimate[3]?m_totalNumEntitiesBeforeEstimate[3]:1)) << stk::end_row;
          table2 << "|" << " estimated sides   " << "|" << m_totalNumEntitiesBeforeEstimate[2]    << "|" << m_totalNumEntitiesAfterEstimate[2]    << "|"
                 << 100.0*(double(m_totalNumEntitiesAfterEstimate[2])/double(m_totalNumEntitiesBeforeEstimate[2]?m_totalNumEntitiesBeforeEstimate[2]:1)) << stk::end_row;

          os << "\n" << table2;
        }
    }

    void RefinementInfo::countCurrentNodes(percept::PerceptMesh& eMesh)
    {
      stk::mesh::Selector selector(eMesh.get_fem_meta_data()->locally_owned_part());
      std::vector<size_t> count ;
      stk::mesh::count_entities( selector, *eMesh.get_bulk_data(), count );

      RefinementInfoCount nnodes = count[0];

      stk::ParallelMachine pm = eMesh.get_bulk_data()->parallel();
      stk::all_reduce( pm, stk::ReduceSum<1>( &nnodes ) );

      for (unsigned i = 0; i < m_refinementInfoByType.size(); i++)
        {
          m_refinementInfoByType[i].m_numNewNodes = nnodes;
        }
    }

    void RefinementInfo::full_stats_before_mark()
    {
      if (!m_full_stats) return;
      for (unsigned i=0; i < 4; ++i)
        {
          m_numberOfMarkedSubDimEntities[i] = 0;
          m_numberOfSubDimEntities[i] = 0;
        }
    }

    void RefinementInfo::full_stats_after_mark()
    {
      if (!m_full_stats) return;
      PerceptMesh& eMesh = m_refiner->getMesh();
      NodeRegistry& nr = m_refiner->getNodeRegistry();
      SubDimCellToDataMap::iterator
        iter,
        iter_begin = nr.getMap().begin(),
        iter_end = nr.getMap().end();

      const unsigned map_sde_size[5] = {0,
                                        (unsigned)eMesh.element_rank(), // only one means centroid
                                        (unsigned)eMesh.edge_rank(),    // 2 = edge
                                        (unsigned)eMesh.side_rank(),    // 3 or 4, face
                                        (unsigned)eMesh.side_rank() };
      for (iter = iter_begin; iter != iter_end; ++iter)
        {
          const SubDimCell_SDCEntityType& subDimEntity = (*iter).first;
          SubDimCellData& nodeId_elementOwnderId = (*iter).second;
          NodeIdsOnSubDimEntityType& nodeIds_onSE = std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnderId);

          unsigned idx = map_sde_size[subDimEntity.size()];
          if (nodeIds_onSE.size())
            ++m_numberOfMarkedSubDimEntities[idx];
          ++m_numberOfSubDimEntities[idx];
        }

      stk::ParallelMachine pm = eMesh.get_bulk_data()->parallel();
      stk::all_reduce( pm, stk::ReduceSum<4>( &m_numberOfSubDimEntities[0] ) );
      stk::all_reduce( pm, stk::ReduceSum<4>( &m_numberOfMarkedSubDimEntities[0] ) );
      //if (eMesh.get_rank() == 0)
    }

    void RefinementInfo::full_stats_before_refine()
    {
      if (!m_full_stats) return;
      PerceptMesh& eMesh = m_refiner->getMesh();
      countRefinedEntities(eMesh.side_rank(), m_totalNumSidesBeforeRefine, &m_numMarkedSides);
      countRefinedEntities(eMesh.element_rank(), m_totalNumElementsBeforeRefine, &m_numMarkedElements);
    }

    void RefinementInfo::full_stats_after_refine()
    {
      if (!m_full_stats) return;
      PerceptMesh& eMesh = m_refiner->getMesh();
      countRefinedEntities(eMesh.side_rank(), m_totalNumSidesAfterRefine, 0);
      countRefinedEntities(eMesh.element_rank(), m_totalNumElementsAfterRefine, 0);
    }

    void RefinementInfo::full_stats_before_filter(stk::mesh::EntityRank rank, RefinementInfoCount count)
    {
      stk::ParallelMachine pm = m_refiner->getMesh().get_bulk_data()->parallel();
      stk::all_reduce( pm, stk::ReduceSum<1>( &count ) );
      m_totalNumEntitiesBeforeFilter[rank] = count;
    }

    void RefinementInfo::full_stats_after_filter(stk::mesh::EntityRank rank, RefinementInfoCount count)
    {
      stk::ParallelMachine pm = m_refiner->getMesh().get_bulk_data()->parallel();
      stk::all_reduce( pm, stk::ReduceSum<1>( &count ) );
      m_totalNumEntitiesAfterFilter[rank] = count;
    }

    void RefinementInfo::full_stats_before_estimate(stk::mesh::EntityRank rank, RefinementInfoCount count)
    {
      stk::ParallelMachine pm = m_refiner->getMesh().get_bulk_data()->parallel();
      stk::all_reduce( pm, stk::ReduceSum<1>( &count ) );
      m_totalNumEntitiesBeforeEstimate[rank] = count;
    }

    void RefinementInfo::full_stats_after_estimate(stk::mesh::EntityRank rank, RefinementInfoCount count)
    {
      stk::ParallelMachine pm = m_refiner->getMesh().get_bulk_data()->parallel();
      stk::all_reduce( pm, stk::ReduceSum<1>( &count ) );
      m_totalNumEntitiesAfterEstimate[rank] = count;
    }

    void RefinementInfo::countRefinedEntities(stk::mesh::EntityRank rank, RefinementInfoCount& num_elem, RefinementInfoCount * num_elem_marked)
    {
      percept::PerceptMesh& eMesh = m_refiner->getMesh();
      NodeRegistry& nodeRegistry = m_refiner->getNodeRegistry();

      if (num_elem_marked)
        *num_elem_marked = 0;
      num_elem = 0;

      const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( rank );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          const unsigned num_elements_in_bucket = bucket.size();
          const CellTopologyData * const cell_topo_data = eMesh.get_cell_topology(bucket);
          if (cell_topo_data == 0)
            {
              continue;
            }
          if (!bucket.owned())
            continue;

          for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
            {
              stk::mesh::Entity element = bucket[iElement];

              if (eMesh.hasFamilyTree(element) && eMesh.isParentElement(element, true))
                continue;

              bool foundMarkedEdge = false;

              //const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);
              if (num_elem_marked)
                {
                  for (stk::mesh::EntityRank irank=stk::topology::EDGE_RANK; irank <= rank; ++irank)
                    {
                      unsigned numSubDimNeededEntities = 0;

                      if (irank == eMesh.edge_rank())
                        {
                          numSubDimNeededEntities = cell_topo_data->edge_count;
                        }
                      else if (irank == eMesh.side_rank())
                        {
                          numSubDimNeededEntities = cell_topo_data->side_count;
                        }
                      else if (irank == eMesh.element_rank())
                        {
                          numSubDimNeededEntities = 1;
                        }

                      for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                        {
                          NodeIdsOnSubDimEntityType* nodeIds_onSE_ptr = nodeRegistry.getNewNodesOnSubDimEntity(element, irank, iSubDimOrd);
                          if (nodeIds_onSE_ptr == 0)
                            {
                              continue;
                            }
                          NodeIdsOnSubDimEntityType& nodeIds_onSE = *nodeIds_onSE_ptr;

                          if (nodeIds_onSE.size() != 0)
                            {
                              foundMarkedEdge = true;
                              break;
                            }
                        }
                      if (foundMarkedEdge)
                        break;
                    }
                  if (foundMarkedEdge)
                    ++(*num_elem_marked);
                }
              ++num_elem;
            }
        }

      stk::ParallelMachine pm = eMesh.get_bulk_data()->parallel();
      stk::all_reduce( pm, stk::ReduceSum<1>( &num_elem ) );
      if (num_elem_marked)
        stk::all_reduce( pm, stk::ReduceSum<1>( num_elem_marked ) );


    }
  }
