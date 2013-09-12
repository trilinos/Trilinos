#include <exception>
#include <fstream>
#include <set>
#include <typeinfo>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <stk_adapt/RefinerUtil.hpp>
#include <stk_adapt/Refiner.hpp>
#include <stk_adapt/NodeRegistry.hpp>

#include <boost/unordered_set.hpp>

namespace stk {
  namespace adapt {

    using namespace std;
    using namespace percept;

#define EXTRA_PRINT_UR_GETBLOCKS 0

    // FIXME move this to a utils class
    BlockNamesType RefinerUtil::getBlockNames(std::string& block_name, unsigned proc_rank, percept::PerceptMesh& eMesh)
    {
      BlockNamesType blocks(stk::percept::EntityRankEnd+1u);
      if (block_name.length() == 0)
        return blocks;

      if (block_name.substr(0, 5) == "file:")
        {
          if (1) throw std::runtime_error("file: option Not implemented");
          std::string fileName = block_name.substr(5, block_name.length()-5);
          std::ifstream file(fileName.c_str());
          while(!file.eof())
            {
              std::string block;
              file >> block;
              if (block[0] != '#')
                {
                  if (block.substr(0,6) == "block_")
                    blocks[stk::mesh::MetaData::ELEMENT_RANK].push_back(block);
                  else if (block.substr(0,8) == "surface_")
                    blocks[eMesh.face_rank()].push_back(block);
                }

            }
        }
      else
        {
          std::string names = block_name;

          // pre-process to look for ".." range indicator

          std::string new_names = names;
          new_names = "";
          while(1)
            {
              if (!names.length())
                break;
              size_t ipos = names.find(',');
              bool last_one =  (ipos == std::string::npos);

              {
                std::string n1 = (last_one ? names : names.substr(0, ipos) );
                bool inc = true;
                //bool exc = false;
                if ('-' == n1[0])
                  {
                    //exc = true;
                    inc = false;
                  }
                else if ('+' == n1[0])
                  {
                  }
                else
                  {
                    n1 = "+" + n1;
                  }
                std::string plus_or_minus = (inc?"+":"-");
                std::string n2 = n1.substr(1, n1.length()-1);
                std::string id_string_start = "";
                std::string id_string_end = "";
                // leave open the possibility for other identifiers for range
                std::string dotdot = "..";
                int dotdot_len = dotdot.length();
                size_t pos_dotdot = n1.find(dotdot);
                if (pos_dotdot != std::string::npos)
                  {
                    if (EXTRA_PRINT_UR_GETBLOCKS) std::cout << "tmp with .., n1= " << n1 << " n2= " << n2 << std::endl;

                    if (n1.length() > 6 && n1.substr(1,6) == "block_")
                      {
                        // +block_123..block_125
                        // 0123456789^1234567890
                        id_string_start = n1.substr(7, pos_dotdot-7);
                        id_string_end = n1.substr(pos_dotdot+dotdot_len+6, n1.length()-(pos_dotdot+dotdot_len+6));
                      }
                    else if (n1.length() > 8 && n1.substr(1,8) == "surface_")
                      {
                        // error
                      }
                    else
                      {
                        // +12..45
                        // 012^456
                        //std::cout << "tmp pos_dotdot= " << pos_dotdot << std::endl;

                        id_string_start = n1.substr(1, pos_dotdot-1);
                        id_string_end = n1.substr(pos_dotdot+dotdot_len+0, n1.length()-(pos_dotdot+dotdot_len+0));
                      }

                    int id_start = 0;
                    int id_end = 0;
                    try {
                      id_start = boost::lexical_cast<int>(id_string_start);
                      id_end = boost::lexical_cast<int>(id_string_end);
                    }
                    catch (std::exception& X)
                      {
                        std::cout << "RefinerUtil::getBlockNames: exception: " << X.what() << std::endl;
                        std::cout << "RefinerUtil::getBlockNames: invalid range syntax in block_name: with .., id_string_start= "
                                  << id_string_start << " id_string_end= " << id_string_end << std::endl;
                        throw std::runtime_error("invalid input syntax");
                      }
                    catch ( const std::exception * X )
                      {
                        std::cout << "RefinerUtil::getBlockNames: exception: " << X->what() << std::endl;
                        std::cout << "RefinerUtil::getBlockNames: invalid range syntax in block_name: with .., id_string_start= "
                                  << id_string_start << " id_string_end= " << id_string_end << std::endl;
                        throw std::runtime_error("invalid input syntax");
                      }
                    catch( ... )
                      {
                        throw std::runtime_error("invalid input syntax");
                      }
                    if (EXTRA_PRINT_UR_GETBLOCKS)
                      {
                        std::cout << "tmp with .., id_string_start= " << id_string_start << " id_string_end= " << id_string_end << std::endl;
                        std::cout << "tmp with .., id_start= " << id_start << " id_end= " << id_end << std::endl;
                      }

                    for (int id=id_start; id <= id_end; id++)
                      {
                        new_names += plus_or_minus+boost::lexical_cast<std::string>(id)+(id == id_end ? "" : ",");
                      }
                    if (!last_one)
                      new_names += ",";
                    if (EXTRA_PRINT_UR_GETBLOCKS)
                      std::cout << "tmp new_names with .. = " << new_names << std::endl;
                  }
                else
                  {
                    new_names += n1 + (last_one? "":",");
                    if (EXTRA_PRINT_UR_GETBLOCKS)
                      std::cout << "tmp new_names without .. = " << new_names << std::endl;
                  }
                if (last_one)
                  {
                    break;
                  }
                else
                  {
                    names = names.substr(ipos+1, names.length()-(ipos+1));
                  }
              }
            }
          if (EXTRA_PRINT_UR_GETBLOCKS)
            std::cout << "tmp new_names after .. (range) processing = " << new_names << std::endl;

          names = new_names;
          std::string names_save = names;

          // post process to remove +name if -name exists
          new_names = "";
          while(1)
            {
              if (!names.length())
                break;
              size_t ipos = names.find(',');
              bool last_one =  (ipos == std::string::npos);

              std::string n1 = (last_one ? names : names.substr(0, ipos) );

              bool inc = true;
              //bool exc = false;
              if ('-' == n1[0])
                {
                  //exc = true;
                  inc = false;
                }
              else if ('+' == n1[0])
                {
                }
              else
                {
                  //error
                }
              std::string n2 = n1.substr(1, n1.length()-1);

              if (inc)
                {
                  size_t jpos = names_save.find("-"+n2);
                  if (jpos != std::string::npos)
                    {
                      // don't add it
                    }
                  else
                    {
                      new_names += n1 + (last_one? "":",");
                    }
                }
              else
                {
                  new_names += n1 + (last_one? "":",");
                }

              if (last_one)
                {
                  break;
                }
              else
                {
                  names = names.substr(ipos+1, names.length()-(ipos+1));
                }
            }


          if (EXTRA_PRINT_UR_GETBLOCKS)
            std::cout << "tmp new_names after post-proc to remove +name if -name exists= " << new_names << std::endl;
          if (new_names.length() && !proc_rank)
            {
              std::cout << "RefinerUtil:: --block_name option after processing for removing -name= " << new_names << std::endl;
            }

          // final step
          names = new_names;
          while(1)
            {
              if (!names.length())
                break;
              size_t ipos = names.find(',');
              bool last_one =  (ipos == std::string::npos);

              {
                std::string n1 = (last_one ? names : names.substr(0, ipos) );

                bool inc = true;
                //bool exc = false;
                if ('-' == n1[0])
                  {
                    //exc = true;
                    inc = false;
                  }
                else if ('+' == n1[0])
                  {
                  }
                else
                  {
                    n1 = "+" + n1;
                  }
                std::string n2 = n1.substr(1, n1.length()-1);

                //std::cout << "n1= " << n1 << " n2= " << n2 << std::endl;
                if (n1.length() > 6 && n1.substr(1,6) == "block_")
                  blocks[stk::mesh::MetaData::ELEMENT_RANK].push_back(n1);
                else if (n1.length() > 8 && n1.substr(1,8) == "surface_")
                  blocks[eMesh.face_rank()].push_back(n1);
                else
                  {
                    std::string pm = (inc?"+":"-");
                    blocks[stk::mesh::MetaData::ELEMENT_RANK].push_back(pm+"block_"+n2);
                  }
                if (last_one)
                  {
                    break;
                  }
                else
                  {
                    names = names.substr(ipos+1, names.length()-(ipos+1));
                  }
              }
            }
          if (EXTRA_PRINT_UR_GETBLOCKS)
            std::cout << "tmp RefinerUtil::getBlockNames: blocks = " << blocks << std::endl;
        }

      return blocks;
    }

    // FIXME move this to a utils class
    /**
     * This method looks for surfaces that share nodes with the blocks specified in @param blocks and if it finds
     * any surfaces (sidesets), they are added to the blocks so they get refined properly.
     * TODO: If a surface is shared by more than one block, an error is thrown.
     */

    BlockNamesType RefinerUtil::correctBlockNamesForPartPartConsistency(percept::PerceptMesh& eMesh, BlockNamesType& blocks)
    {
      if (EXTRA_PRINT_UR_GETBLOCKS) std::cout << "RefinerUtil::correctBlockNamesForPartPartConsistency..." << std::endl;

      if (blocks[stk::mesh::MetaData::ELEMENT_RANK].size() == 0)
        return blocks;

      stk::mesh::EntityRank subDimRank = (eMesh.get_spatial_dim() == 3 ? eMesh.face_rank() : eMesh.edge_rank());

      mesh::PartVector all_parts = eMesh.get_fem_meta_data()->get_parts();
      for (mesh::PartVector::iterator i_part = all_parts.begin(); i_part != all_parts.end(); ++i_part)
        {
          mesh::Part *  part = *i_part ;

          for (mesh::PartVector::iterator i_surfacePart = all_parts.begin(); i_surfacePart != all_parts.end(); ++i_surfacePart)
            {
              mesh::Part *  surfacePart = *i_surfacePart ;
              if ( stk::mesh::is_auto_declared_part(*surfacePart) )
                continue;

              const CellTopologyData * part_cell_topo_data = eMesh.get_cell_topology(*surfacePart);
              CellTopology surf_topo(part_cell_topo_data);
              //if (EXTRA_PRINT_UR_GETBLOCKS) std::cout << "tmp srk surfacePart= " << surfacePart->name() << " topo= " << (part_cell_topo_data?surf_topo.getName() : "NULL") << std::endl;

              if (part_cell_topo_data && part->primary_entity_rank() == stk::mesh::MetaData::ELEMENT_RANK && surfacePart->primary_entity_rank() == subDimRank)
                {
                  std::string partNamePlus = "+" + part->name();
                  std::vector<std::string>::iterator partInBlocks = std::find(blocks[stk::mesh::MetaData::ELEMENT_RANK].begin(), blocks[stk::mesh::MetaData::ELEMENT_RANK].end(), partNamePlus);
                  // if this part is not in the blocks list, skip it
                  if (partInBlocks == blocks[stk::mesh::MetaData::ELEMENT_RANK].end())
                    {
                      //if (EXTRA_PRINT_UR_GETBLOCKS) std::cout << "tmp srk skipping part= " << partNamePlus << std::endl;
                      continue;
                    }
                  std::string surfacePartNamePlus = "+" + surfacePart->name();
                  std::vector<std::string>::iterator surfacePartInBlocks = std::find(blocks[subDimRank].begin(), blocks[subDimRank].end(), surfacePartNamePlus);
                  // if this surface is already in the list, skip it
                  if (surfacePartInBlocks != blocks[subDimRank].end())
                    {
                      //if (EXTRA_PRINT_UR_GETBLOCKS) std::cout << "tmp srk skipping surf= " << surfacePartNamePlus << std::endl;
                      continue;
                    }
                  bool isBoundarySurface= eMesh.isBoundarySurface(*part, *surfacePart);

                  if (EXTRA_PRINT_UR_GETBLOCKS) std::cout << "tmp srk isBoundarySurface for part/surf= " << part->name() << " / " << surfacePart->name() << " = " << isBoundarySurface << std::endl;
                  if (isBoundarySurface)
                    {
                      if (EXTRA_PRINT_UR_GETBLOCKS) std::cout << "tmp part [" << part->name() << "] shares sideset [" << surfacePart->name() << "]" << std::endl;
                      blocks[subDimRank].push_back(std::string("+"+surfacePart->name()));
                    }
                  else
                    {
                      //std::cout << "tmp part [" << part->name() << "] doesn't shares sideset [" << surfacePart->name() << "]" << std::endl;
                    }
                }
            }
        }
      if (0) std::cout << "tmp RefinerUtil::correctBlockNamesForPartPartConsistency: blocks = " << blocks << std::endl;
      return blocks;
    }

    //static
    void RefinerUtil::
    addAncestorsToUnrefineList(percept::PerceptMesh& eMesh, int num_levels_to_add, ElementUnrefineCollection& elements_to_unref)
    {
      int num_levels_to_add_1 = num_levels_to_add;
      if (num_levels_to_add < 0) num_levels_to_add_1 = 1000;
      //std::cout << "here 1 elements_to_unref.size= " << elements_to_unref.size() << std::endl;
      for (int ilev=0; ilev < num_levels_to_add_1; ilev++)
        {
          ElementUnrefineCollection to_add(*eMesh.get_bulk_data());
          for (ElementUnrefineCollection::iterator iter = elements_to_unref.begin(); iter != elements_to_unref.end(); ++iter)
            {
              stk::mesh::Entity element = *iter;
              if (eMesh.hasFamilyTree(element))
                {
                  //std::cout << "here 2 elements_to_unref.size= " << elements_to_unref.size() << std::endl;
                  stk::mesh::Entity parent = eMesh.getParent(element, false);
                  if (elements_to_unref.find(parent) != elements_to_unref.end())
                    continue;
                  //std::cout << "here 3 elements_to_unref.size= " << elements_to_unref.size() << std::endl;
                  std::vector<stk::mesh::Entity> children;
                  bool hasChildren = eMesh.getChildren(parent, children, true, false);
                  if (hasChildren && children.size())
                    {
                      //std::cout << "here 4 elements_to_unref.size= " << elements_to_unref.size() << std::endl;
                      bool allChildrenInUnrefSet = true;
                      for (unsigned ichild=0; ichild < children.size(); ichild++)
                        {
                          if (elements_to_unref.find(children[ichild]) == elements_to_unref.end())
                            {
                              allChildrenInUnrefSet = false;
                              break;
                            }
                        }
                      if (allChildrenInUnrefSet)
                        {
                          //std::cout << "here 5 elements_to_unref.size= " << elements_to_unref.size() << std::endl;
                          to_add.insert(parent);
                        }
                    }
                }
            }
          //std::cout << "RefinerUtil::addAncestorsToUnrefineList ilev= " << ilev << " to_add.size= " << to_add.size() << std::endl;
          if (to_add.size())
            {
              elements_to_unref.insert(to_add.begin(), to_add.end());
            }
          else
            {
              break;
            }
        }

    }

    /// create missing edges after adapt - for edge-based simulators
    void RefinerUtil::
    create_missing_edges(percept::PerceptMesh& eMesh)
    {
      typedef MySubDimCell<SDCEntityType, 2, CompareSDCEntityType> SubDimCell2;
      typedef boost::unordered_set<SubDimCell2, my_fast_hash<SDCEntityType, 2>, my_fast_equal_to<SDCEntityType, 2> > SubDimCellSet;
      SubDimCellSet edge_set, new_edge_set;

      // put existing edges in the set
      const std::vector<stk::mesh::Bucket*> & edge_buckets = eMesh.get_bulk_data()->buckets( eMesh.edge_rank() );
      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = edge_buckets.begin() ; k != edge_buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          // add all edges since we want to avoid creating them even if they are shared from another proc
          //if (bucket.owned())
            {
              const unsigned num_edges_in_bucket = bucket.size();

              for (unsigned i_edge = 0; i_edge < num_edges_in_bucket; i_edge++)
                {
                  stk::mesh::Entity edge = bucket[i_edge];
                  const percept::MyPairIterRelation edge_nodes (eMesh, edge, eMesh.node_rank());
                  SubDimCell2 subDimEntity(eMesh);
                  subDimEntity.clear();
                  subDimEntity.insert(edge_nodes[0].entity());
                  subDimEntity.insert(edge_nodes[1].entity());
                  edge_set.insert(subDimEntity);
                }
            }
        }

      // visit elements (and sides if in 3D) and check if edges present; add as needed
      stk::mesh::EntityRank rank_start = eMesh.side_rank();
      if (eMesh.get_spatial_dim() == 2) rank_start = eMesh.element_rank();
      for (stk::mesh::EntityRank rank = rank_start; rank <= eMesh.element_rank(); ++rank)
        {
          const std::vector<stk::mesh::Bucket*> & buckets = eMesh.get_bulk_data()->buckets( rank );
          for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;
              if (bucket.owned())
                {
                  const CellTopologyData * cell_topo_data = eMesh.get_cell_topology(bucket);
                  CellTopology topo(cell_topo_data);
                  unsigned edge_count = cell_topo_data->edge_count;
                  const unsigned num_elements_in_bucket = bucket.size();

                  for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                    {
                      stk::mesh::Entity element = bucket[i_element];
                      const percept::MyPairIterRelation elem_nodes (eMesh, element, eMesh.node_rank());
                      for (unsigned iedgeOrd=0; iedgeOrd < edge_count; iedgeOrd++)
                        {
                          unsigned in0 = cell_topo_data->edge[iedgeOrd].node[0];
                          unsigned in1 = cell_topo_data->edge[iedgeOrd].node[1];
                          SubDimCell2 subDimEntity(eMesh);
                          subDimEntity.clear();
                          subDimEntity.insert(elem_nodes[in0].entity());
                          subDimEntity.insert(elem_nodes[in1].entity());
                          SubDimCellSet::iterator fnd = edge_set.find(subDimEntity);
                          if (fnd == edge_set.end())
                            {
                              edge_set.insert(subDimEntity);
                              new_edge_set.insert(subDimEntity);
                            }
                        }
                    }
                }
            }
        }

      // now create new/missing edges
      std::vector<stk::mesh::Entity> new_edges;
      eMesh.createEntities(eMesh.edge_rank(), new_edge_set.size(), new_edges);

      eMesh.get_bulk_data()->modification_begin();
      unsigned count=0;
      for (SubDimCellSet::iterator edge_it = new_edge_set.begin(); edge_it != new_edge_set.end(); ++edge_it, ++count)
        {
          const SubDimCell2& edge = *edge_it;
          eMesh.get_bulk_data()->declare_relation(new_edges[count], edge[0], 0);
          eMesh.get_bulk_data()->declare_relation(new_edges[count], edge[1], 1);
        }
      eMesh.get_bulk_data()->modification_end();
    }

  }
}
