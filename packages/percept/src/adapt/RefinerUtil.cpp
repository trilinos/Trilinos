// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <exception>
#include <fstream>
#include <set>
#include <typeinfo>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <adapt/SerializeNodeRegistry.hpp>
#include <adapt/RefinerUtil.hpp>
#include <adapt/Refiner.hpp>
#include <adapt/NodeRegistry.hpp>
#include <percept/PEnums.hpp>

#include <stk_mesh/base/MeshUtils.hpp>

namespace percept {

  using namespace percept;

#define EXTRA_PRINT_UR_GETBLOCKS 0

  static std::string strip_mult(const std::string& inp)
  {
    size_t pc = inp.find(":");
    if (pc == std::string::npos)
      return inp;
    else
      return inp.substr(0,pc);
  }

  BlockNamesType RefinerUtil::getBlockNames(const std::string& block_name_0, unsigned proc_rank, percept::PerceptMesh& eMesh, const std::string& geomFile)
  {
    std::string block_name = block_name_0;

    BlockNamesType blocks(percept::EntityRankEnd+1u);
    if (block_name.length() == 0)
      return blocks;

    if (block_name.substr(0, 5) == "file:")
      {
        //if (1) throw std::runtime_error("file: option Not implemented");
        std::string fileName = block_name.substr(5, block_name.length()-5);
        for (int iter=0; iter < 1000; ++iter)
          {
            if (fileName[0] == ' ')
              fileName = fileName.substr(1, fileName.length()-1);
            else
              break;
          }
        if (eMesh.get_rank() == 0)
          {
            std::cout << "block_names processing, found 'file:', will read from file = " << fileName << std::endl;
          }
        std::ifstream file(fileName.c_str());
        int line=0;
        block_name = "";
        std::string str_line;
        while(std::getline(file, str_line))
          {
            std::string block = str_line;
            if (eMesh.get_rank() == 0)
              std::cout << "file: line= " << line << " block= " << block << std::endl;
            if (block[0] != '#' && block.length())
              {
                if (line) block_name += ",";
                block_name += block;
                ++line;
              }
          }
        if (eMesh.get_rank() == 0)
          {
            std::cout << "block_names processing, after file read string is: " << block_name << std::endl;
          }
      }

    if (1)
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

                  std::string mult = "";
                  size_t NPOS = std::string::npos;
                  size_t cpos = id_string_end.find(":");
                  if (cpos != std::string::npos)
                    {
                      size_t xpos = id_string_end.find("x");
                      if (xpos == NPOS)
                        xpos = id_string_end.find("X");
                      VERIFY_OP_ON(xpos, !=, NPOS, "missing x or X in blocks specifier");
                      mult = id_string_end.substr(cpos);
                      id_string_end = id_string_end.substr(0, cpos);
                    }
                  if (EXTRA_PRINT_UR_GETBLOCKS) std::cout << "tmp with .., id_string_start= " << id_string_start << " id_string_end= " << id_string_end << " mult= " << mult << std::endl;

                  int id_start = 0;
                  int id_end = 0;
                  try {
                    id_start = std::stoi(id_string_start);
                    id_end = std::stoi(id_string_end);
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
                      new_names += plus_or_minus + std::to_string(id) + mult + (id == id_end ? "" : ",");
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
        bool has_plus = names.find('+') != std::string::npos;
        new_names = "";
        while(1)
          {
            if (!names.length())
              break;
            size_t pos_comma = names.find(',');
            bool last_one =  (pos_comma == std::string::npos);

            std::string n1 = (last_one ? names : names.substr(0, pos_comma) );

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
            if (EXTRA_PRINT_UR_GETBLOCKS) std::cout << "n1= " << n1 << " n2= " << n2 << " strip_mult(n2)= " << strip_mult(n2) << std::endl;
            n2 = strip_mult(n2);
            // if (n2 == "block_910")
            //   {
            //     std::cout << "n1= " << n1 << " n2= " << n2 << std::endl;
            //   }

            if (inc)
              {

                size_t jpos = names_save.find("-"+n2);
                bool add_it = true;
                if (jpos != std::string::npos)
                  {
                    // don't add it
                    add_it = false;
                    // check for full match
                    std::string ns;
                    for (size_t ipos=jpos; ipos < names_save.length(); ++ipos)
                      {
                        if (names_save[ipos] == ',')
                          break;
                        ns += names_save[ipos];
                      }
                    std::string ns2 = ns.substr(1, ns.length()-1);
                    if (ns2 != n2)
                      add_it = true;
                    if (EXTRA_PRINT_UR_GETBLOCKS)
                      {
                        std::cout << "ns= " << ns << " ns2= " << ns2 << " n2= " << n2 << std::endl;
                      }
                  }
                if (add_it)
                  {
                    new_names += n1 + (last_one? "":",");
                  }
              }
            else
              {
                if (!has_plus)
                  new_names += n1 + (last_one? "":",");
              }

            if (last_one)
              {
                break;
              }
            else
              {
                names = names.substr(pos_comma+1, names.length()-(pos_comma+1));
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
                blocks[stk::topology::ELEMENT_RANK].push_back(n1);
              else if (n1.length() > 8 && n1.substr(1,8) == "surface_")
                blocks[eMesh.face_rank()].push_back(n1);
              else
                {
                  std::string pm = (inc?"+":"-");
                  blocks[stk::topology::ELEMENT_RANK].push_back(pm+"block_"+n2);
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

    // change so all -block get removed and replaced with + only
    for (stk::mesh::EntityRank irank=eMesh.side_rank(); irank <= eMesh.element_rank(); irank++)
      {
        bool found_minus = false, found_plus = false;
        for (unsigned ib=0; ib < blocks[irank].size(); ++ib)
          {
            if (blocks[irank][ib][0] == '-')
              {
                found_minus = true;
              }
            if (blocks[irank][ib][0] == '+')
              {
                found_plus = true;
              }
          }
        if (found_minus && found_plus)
          {
            std::ostringstream errmsg;
            errmsg << "You've specified an inconsistent combination of block names, mixing -block_n... with +block_m...\n";
            errmsg << "You can only specify -block_n if there are no +block_m's or if the +block_m's contain the -block_n\n";
            errmsg << "   to be deleted, e.g. +1,2,3.7,-4\n";
            errmsg << "Curent processed blocks list = " << blocks;
            errmsg << "\nSpecified on input: " << block_name << std::endl;
            throw std::runtime_error(errmsg.str());
          }

        if (found_minus)
          {
            std::vector<std::string> new_names;
            stk::mesh::PartVector all_parts = eMesh.get_fem_meta_data()->get_parts();
            for (stk::mesh::PartVector::iterator i_part = all_parts.begin(); i_part != all_parts.end(); ++i_part)
              {
                stk::mesh::Part *  part = *i_part ;
                if (eMesh.is_auto_or_geom_part(geomFile, part))
                  continue;

                if (part->primary_entity_rank() == irank)
                  {
                    std::string partNamePlus = "+" + part->name();
                    std::string partNameMinus = "-" + part->name();
                    bool remove_it = false;
                    for (unsigned ib=0; ib < blocks[irank].size(); ++ib)
                      {
                        if (partNameMinus == blocks[irank][ib])
                          {
                            remove_it = true;
                            break;
                          }
                      }
                    if (!remove_it)
                      {
                        new_names.push_back(partNamePlus);
                      }
                  }
              }
            if (EXTRA_PRINT_UR_GETBLOCKS && !eMesh.get_rank())
              std::cout << "RefinerUtil::getBlockNames: irank= " << irank << " old blocks = \n" << blocks[irank] << "\n new= \n" << new_names << std::endl;
            blocks[irank] = new_names;
          }
      }

    return blocks;
  }

  static void make_parallel_consistent(PerceptMesh& eMesh, std::vector<std::string>& blocks)
  {
    char buf[1024];
    stk::CommSparse comm_all(eMesh.parallel());

    unsigned proc_size = comm_all.parallel_size();

    // pack
    for (unsigned pr=0; pr < proc_size; pr++)
      {
        comm_all.send_buffer( pr ).pack< unsigned > (blocks.size());
        for (unsigned ib=0; ib < blocks.size(); ib++)
          {
            comm_all.send_buffer( pr ).pack< unsigned > ( blocks[ib].length());
            comm_all.send_buffer( pr ).pack< char > (blocks[ib].c_str(), blocks[ib].length());
          }
      }
    // allocateBuffers
    {
      comm_all.allocate_buffers();
    }
    // pack
    for (unsigned pr=0; pr < proc_size; pr++)
      {
        comm_all.send_buffer( pr ).pack< unsigned > (blocks.size());
        for (unsigned ib=0; ib < blocks.size(); ib++)
          {
            comm_all.send_buffer( pr ).pack< unsigned > ( blocks[ib].length());
            comm_all.send_buffer( pr ).pack< char > (blocks[ib].c_str(), blocks[ib].length());
          }
      }
    // communicate
    comm_all.communicate();

    // unpack
    std::set<std::string> blocks_new(blocks.begin(), blocks.end());
    for (unsigned pr=0; pr < proc_size; pr++)
      {
        unsigned bsize=0;
        comm_all.recv_buffer( pr ).unpack< unsigned > (bsize);
        for (unsigned ib=0; ib < bsize; ib++)
          {
            unsigned len=0;
            comm_all.recv_buffer( pr ).unpack< unsigned > ( len );
            comm_all.recv_buffer( pr ).unpack< char > (buf, len);
            std::string str(buf, len);
            blocks_new.insert(str);
          }
      }
    blocks.resize(0);
    blocks.assign(blocks_new.begin(), blocks_new.end());
    std::sort(blocks.begin(), blocks.end());
    //std::cout << "make_parallel_consistent blocks= " << blocks << std::endl;
  }

  /**
   * This method looks for surfaces of blocks specified in @param blocks and if it finds
   * any surfaces (sidesets), they are added to the blocks so they get refined properly.
   * If a surface is shared by more than one block, and both blocks are not in the @param blocks list,
   *   an error is thrown.
   *
   * algorithm:
   *   for sides in side_part_map
   *     if all elem_part in side_part_map[side] are in blocks, ok, else throw
   *     if side not in blocks[side_rank], add it
   *
   */

#define DEBUG_CORRECT_BLOCKS_1 0

  BlockNamesType RefinerUtil::correctBlockNamesForPartPartConsistency(percept::PerceptMesh& eMesh, BlockNamesType& blocks, const std::string& geomFile)
  {
    if (blocks[stk::topology::ELEMENT_RANK].size() == 0)
      return blocks;

    SidePartMap side_part_map;
    Refiner::get_side_part_relations(eMesh, false, side_part_map);

    stk::mesh::EntityRank subDimRank = eMesh.side_rank();

    SidePartMap::iterator iter;
    for (iter = side_part_map.begin(); iter != side_part_map.end(); iter++)
      {
        stk::mesh::Part *side_part = iter->first;
        std::string side_part_name = side_part->name();
        const stk::mesh::PartVector *side_pv  = side_part->attribute<stk::mesh::PartVector>();
        stk::mesh::PartVector side_pv1;
        if (side_pv) side_pv1 = *side_pv;
        side_pv1.push_back(side_part);

        stk::mesh::PartVector elem_pv1;

        const stk::mesh::PartVector& epv = iter->second;
        for (unsigned iepv=0; iepv < epv.size(); iepv++)
          {
            stk::mesh::Part *elem_part = epv[iepv];
            std::string elem_part_name = elem_part->name();
            const stk::mesh::PartVector *elem_pv  = elem_part->attribute<stk::mesh::PartVector>();
            elem_pv1.push_back(elem_part);
            if (elem_pv) elem_pv1.insert(elem_pv1.end(), elem_pv->begin(), elem_pv->end());
          }

        for (unsigned iside_pv = 0; iside_pv < side_pv1.size(); iside_pv++)
          {
            stk::mesh::Part *surfacePart = side_pv1[iside_pv];
            const CellTopologyData * surfacePart_topo_data = eMesh.get_cell_topology(*surfacePart);
            CellTopology surf_topo(surfacePart_topo_data);
            if (!surfacePart_topo_data)
              continue;

            bool at_least_one_elem_part_in_block_names = false;
            for (unsigned ielem_pv = 0; ielem_pv < elem_pv1.size(); ielem_pv++)
              {
                stk::mesh::Part *elementPart = elem_pv1[ielem_pv];
                stk::mesh::Part *part = elementPart;
                std::string partNamePlus = "+" + part->name();

                std::vector<std::string>::iterator partInBlocks = std::find(blocks[eMesh.element_rank()].begin(), blocks[eMesh.element_rank()].end(), partNamePlus);
                if (partInBlocks != blocks[eMesh.element_rank()].end())
                  {
                    at_least_one_elem_part_in_block_names = true;
                    break;
                  }
              }

            if (!at_least_one_elem_part_in_block_names)
              continue;

            // if at_least_one_elem_part_in_block_names then all must be in there
            for (unsigned ielem_pv = 0; ielem_pv < elem_pv1.size(); ielem_pv++)
              {
                stk::mesh::Part *elementPart = elem_pv1[ielem_pv];
                stk::mesh::Part *part = elementPart;
                std::string partNamePlus = "+" + part->name();

                std::vector<std::string>::iterator partInBlocks = std::find(blocks[eMesh.element_rank()].begin(), blocks[eMesh.element_rank()].end(), partNamePlus);
                bool found_it = partInBlocks != blocks[eMesh.element_rank()].end();

                if (!found_it)
                  {
                    std::ostringstream errmsg;
                    errmsg << "ERROR: correctBlockNamesForPartPartConsistency: found a surface (" + surfacePart->name() + ") that shares a block (" + part->name() + ") that is not in the specified\n"
                      " list of blocks to be refined.  Re-run by adding the missing block(s), which are = \n";
                    for (unsigned jelem_pv = 0; jelem_pv < elem_pv1.size(); jelem_pv++)
                      {
                        stk::mesh::Part *jelementPart = elem_pv1[jelem_pv];
                        stk::mesh::Part *jpart = jelementPart;
                        std::string jpartNamePlus = "+" + jpart->name();
                        std::vector<std::string>::iterator jpartInBlocks = std::find(blocks[eMesh.element_rank()].begin(), blocks[eMesh.element_rank()].end(), jpartNamePlus);
                        if (jpartInBlocks == blocks[eMesh.element_rank()].end())
                          {
                            errmsg << " " << jpart->name();
                          }
                      }
                    errmsg << std::endl;
                  }

                if (surfacePart_topo_data && part->primary_entity_rank() == eMesh.element_rank() && surfacePart->primary_entity_rank() == subDimRank)
                  {
                    std::string surfacePartNamePlus = "+" + surfacePart->name();
                    std::vector<std::string>::iterator surfacePartInBlocks = std::find(blocks[subDimRank].begin(), blocks[subDimRank].end(), surfacePartNamePlus);
                    // if this surface is already in the list, skip it
                    if (surfacePartInBlocks != blocks[subDimRank].end())
                      {
                        continue;
                      }
                    bool isBoundarySurface= true; // by definition, the side part map is map of sides to shared element parts

                    if (isBoundarySurface)
                      {
                        blocks[subDimRank].push_back(std::string("+"+surfacePart->name()));
                      }
                  }
              }
          }
      }

    for (unsigned ibr=0; ibr < blocks.size(); ibr++)
      {
        make_parallel_consistent(eMesh,blocks[ibr]);
      }

    return blocks;
  }

  void RefinerUtil::remove_existing_nodes(PerceptMesh& eMesh)
  {
    const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( eMesh.node_rank() );
    
    stk::mesh::Part* new_nodes_part = eMesh.get_non_const_part("refine_new_nodes_part");
    VERIFY_OP_ON(new_nodes_part, != , 0, "new_nodes_part");

    std::vector<stk::mesh::Part*> remove_parts(1, new_nodes_part);
    std::vector<stk::mesh::Part*> add_parts;
    std::vector<stk::mesh::Entity> node_vec;
    
    stk::mesh::Selector removePartSelector(*new_nodes_part & eMesh.get_fem_meta_data()->locally_owned_part() );
    for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) {

      stk::mesh::Bucket & bucket = **k ;
      if (removePartSelector(bucket)) {
        
        const unsigned num_entity_in_bucket = bucket.size();
        for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++) {

          stk::mesh::Entity node = bucket[ientity];
          node_vec.push_back(node);
        }
      }
    }
    for (unsigned ii=0; ii < node_vec.size(); ii++) {
      eMesh.get_bulk_data()->change_entity_parts( node_vec[ii], add_parts, remove_parts );
    }
  }

  void RefinerUtil::add_new_nodes(PerceptMesh& eMesh)
  {
    const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( eMesh.node_rank() );
    
    stk::mesh::Part* new_nodes_part = eMesh.get_non_const_part("refine_new_nodes_part");
    VERIFY_OP_ON(new_nodes_part, != , 0, "new_nodes_part");

    std::vector<stk::mesh::Entity> new_nodes;
    for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) {

      stk::mesh::Bucket & bucket = **k ;
      if (bucket.owned()) {

        const unsigned num_entity_in_bucket = bucket.size();
        for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++) {

          stk::mesh::Entity node = bucket[ientity];
          NewNodesType::value_type *ndata = stk::mesh::field_data(*eMesh.m_new_nodes_field, node);
          if (ndata && ndata[0] != 0) {
            new_nodes.push_back(node);
          }
        }
      }
    }
    
    std::vector<stk::mesh::Part*> add_parts(1, new_nodes_part);
    std::vector<stk::mesh::Part*> remove_parts;
    for (unsigned ind = 0; ind < new_nodes.size(); ind++) {
      eMesh.get_bulk_data()->change_entity_parts( new_nodes[ind], add_parts, remove_parts );
    }
  }

  void RefinerUtil::collect_locally_owned_entities(PerceptMesh& eMesh,  stk::mesh::EntityRank rank,  std::vector<stk::mesh::Entity>& elements)
  {
    const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( rank );
    for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
      {
        stk::mesh::Bucket & bucket = **k ;
        const unsigned num_elements_in_bucket = bucket.size();
        if (bucket.owned())
          {
            for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
              {
                elements.push_back(bucket[iElement]);
              }
          }
      }
  }

  std::string rank_to_string(stk::mesh::EntityRank rank)
  {
    std::ostringstream orank;
    orank << rank;
    return orank.str();
  }
    
  void RefinerUtil::get_parent_entity_and_id(PerceptMesh& eMesh, stk::mesh::EntityRank rank, stk::mesh::Entity& element, stk::mesh::Entity& parent_elem, bool debug)
  {
    ParentElementType::value_type *parent_element_id_from_field = NULL;

    if (is_matching_rank(*eMesh.m_parent_element_field, element))
      {
        parent_element_id_from_field = stk::mesh::field_data( *eMesh.m_parent_element_field , element );
      }
    else if (eMesh.m_parent_element_field_side && is_matching_rank(*eMesh.m_parent_element_field_side, element))
      {
        parent_element_id_from_field = stk::mesh::field_data( *eMesh.m_parent_element_field_side , element );
      }
    VERIFY_OP_ON(parent_element_id_from_field, !=, 0, "error parent element field not set on element");
    
    stk::mesh::EntityId parent_element_id = static_cast<stk::mesh::EntityId>(parent_element_id_from_field[0]);

    if (parent_element_id) 
      parent_elem = eMesh.get_bulk_data()->get_entity(rank, parent_element_id);

    if (debug && rank == eMesh.side_rank())
      {
        std::cout << "RefinerUtil::rebuild_family_tree_child: side= " << eMesh.identifier(element) << " parent_element_id_from_field= "
                  << parent_element_id << " parent.valid= " << eMesh.is_valid(parent_elem);
        eMesh.print(element);
      }
    
    if (parent_element_id)
      {
        
        if (rank == eMesh.side_rank())
          {
            //VERIFY_OP_ON(eMesh.is_valid(parent_elem), ==, false, "root cause");
            parent_elem = stk::mesh::Entity();
            
            stk::mesh::EntityId id_new = 0;
            stk::mesh::ConnectivityOrdinal ord;
            eMesh.decipher_exodus_side_id(parent_element_id, id_new, ord);
            VERIFY_OP_ON(id_new, !=, 0, "bad id_new");
            if (id_new)
              {
                stk::mesh::Entity element_owner = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), id_new);
                
                VERIFY_OP_ON(eMesh.is_valid(element_owner), ==, true, "bad parent side elem");
                percept::MyPairIterRelation elem_to_side_rels (eMesh, element_owner, rank);
                for (unsigned jj=0; jj < elem_to_side_rels.size(); ++jj)
                  {
                    if (elem_to_side_rels[jj].relation_ordinal() == ord)
                      {
                        parent_elem = elem_to_side_rels[jj].entity();
                        break;
                      }
                  }
              }
          }
        
        VERIFY_OP_ON(eMesh.is_valid(parent_elem), ==, true, "bad parent found, rank= "+rank_to_string(rank));
        VERIFY_OP_ON(parent_elem, !=, element, "hmmm, parent=element");
      }
  }

  void RefinerUtil::find_or_set_parent_child_relation(PerceptMesh& eMesh, stk::mesh::EntityRank rank, stk::mesh::Entity& element, stk::mesh::Entity& parent_elem, size_t& i_ft, std::vector<stk::mesh::Entity>& ft_new_elements)
  {
    const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);

    stk::mesh::Entity family_tree = stk::mesh::Entity();
    bool found = false;
    unsigned ordinal = 0;
    
    if (eMesh.hasFamilyTree(parent_elem)) {
      
      //percept::MyPairIterRelation family_tree_relations (eMesh, element, FAMILY_TREE_RANK);
      percept::MyPairIterRelation parent_to_family_tree_relations (eMesh, parent_elem, FAMILY_TREE_RANK);
      VERIFY_OP_ON(parent_to_family_tree_relations.size(), > , 0, "hmmm");
      VERIFY_OP_ON(parent_to_family_tree_relations.size(), <=, 2, "bad family tree relations");
      
      bool isParentAlready = false;
      if (parent_to_family_tree_relations.size() == 1) {
        
        unsigned parent_elem_ft_level_0 = eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, parent_elem);
        family_tree = parent_to_family_tree_relations[parent_elem_ft_level_0].entity();
      }
      else if (parent_to_family_tree_relations.size() == 2) {
        
        unsigned parent_elem_ft_level_1 = eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_1, parent_elem);
        family_tree = parent_to_family_tree_relations[parent_elem_ft_level_1].entity();
        isParentAlready = true;
      }
      
      VERIFY_OP_ON(eMesh.is_valid(family_tree), ==, true, "hmmm");
      percept::MyPairIterRelation family_tree_relations (eMesh, family_tree, eMesh.entity_rank(parent_elem));
      VERIFY_OP_ON(family_tree_relations.size(), > , 0, "hmmm");
      if (family_tree_relations[0].entity() == parent_elem) {
        
        isParentAlready = true;
        VERIFY_OP_ON(eMesh.isParentElement(parent_elem), ==, true, "hmmm");
      }
      for (unsigned i = 1; i < family_tree_relations.size(); i++) {
        
        //if (family_tree_relations[i].relation_ordinal() == (ordinal + 1))
        if (family_tree_relations[i].entity() == element) {
          found = true;
          break;
        }
      }
      if (!found && isParentAlready) {
        VERIFY_OP_ON(family_tree_relations.size(), > , 0, "hmmm");
        ordinal = family_tree_relations.size() - 1;
      }
    }
    if (!found) {
      if (i_ft >= ft_new_elements.size()) {
        throw std::runtime_error("ran out of ft_new_elements");
      }
      stk::mesh::Entity familyTreeNewElement = ft_new_elements[i_ft++];
      UniformRefinerPatternBase::set_parent_child_relations(eMesh, parent_elem, element, familyTreeNewElement, ordinal);
    }
  }
  
  void RefinerUtil::rebuild_family_tree(PerceptMesh& eMesh, bool debug)
  {
    VERIFY_OP_ON(eMesh.m_parent_element_field_set, == , true, "parent_element_field_set");
    VERIFY_OP_ON(eMesh.m_new_nodes_field_set, == , true, "new_nodes_field_set");

    eMesh.get_bulk_data()->modification_begin();
    
    remove_existing_nodes(eMesh);    
    add_new_nodes(eMesh);

    stk::mesh::fixup_ghosted_to_shared_nodes(*eMesh.get_bulk_data());

    eMesh.get_bulk_data()->modification_end();

    for (stk::mesh::EntityRank rank = eMesh.side_rank(); rank <= eMesh.element_rank(); ++rank) {

      std::vector<stk::mesh::Entity> ft_new_elements;
      const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);
      
      std::vector<stk::mesh::Entity> elements;
      collect_locally_owned_entities(eMesh, rank, elements);
      
      size_t num_elem_needed = elements.size();
      
      eMesh.get_bulk_data()->modification_begin();
      
      eMesh.createEntities( FAMILY_TREE_RANK, 2*num_elem_needed, ft_new_elements);
      size_t i_ft = 0;
      
      for (size_t ii=0; ii < elements.size(); ++ii) {
        
        stk::mesh::Entity element = elements[ii];
        stk::mesh::Entity parent_elem = stk::mesh::Entity();
        
        get_parent_entity_and_id(eMesh, rank, element, parent_elem, debug);
        
        if (!eMesh.is_valid(parent_elem))
          continue;
        
        find_or_set_parent_child_relation(eMesh, rank, element, parent_elem, i_ft, ft_new_elements);
      }
      
      std::vector<stk::mesh::Entity> ft_delete;
      for (unsigned ii=i_ft; ii < ft_new_elements.size(); ++ii) {

        if (eMesh.get_bulk_data()->has_no_relations(ft_new_elements[ii]))
          ft_delete.push_back(ft_new_elements[ii]);
      }
      for (unsigned ii=0; ii < ft_delete.size(); ++ii) {
        
        if (!eMesh.get_bulk_data()->destroy_entity(ft_delete[ii]))
          throw std::runtime_error("bad destroy of family_tree");
      }
      stk::mesh::fixup_ghosted_to_shared_nodes(*eMesh.get_bulk_data());
      
      eMesh.get_bulk_data()->modification_end();
    }
    
    // get parts correct
    eMesh.get_bulk_data()->modification_begin();
    Refiner::set_active_part(eMesh);
    stk::mesh::fixup_ghosted_to_shared_nodes(*eMesh.get_bulk_data());
    eMesh.get_bulk_data()->modification_end();
  }

  void RefinerUtil::save_node_registry(PerceptMesh& eMesh, NodeRegistry& nodeRegistry, const std::string& msg, bool doComm)
  {
    if (eMesh.m_node_registry_field == 0)
      return;

    if (1)
      {
        std::vector<stk::mesh::Entity> vec;
        stk::mesh::Selector sel = eMesh.get_fem_meta_data()->universal_part();
        stk::mesh::get_selected_entities(sel , eMesh.get_bulk_data()->buckets(eMesh.node_rank()), vec);

        SubDimCell_SDCEntityType subDimEntity(&eMesh);

        // FIXME
        for (size_t ii=0; ii < vec.size(); ++ii)
          {
            stk::mesh::Entity node = vec[ii];

            VERIFY_OP_ON(eMesh.entity_rank(node), ==, eMesh.node_rank(), "bad rank");
            double *node_data = stk::mesh::field_data(*eMesh.m_node_registry_field, node);
            for (unsigned kk=0; kk < NUM_NR_FIELD_SLOTS; ++kk)
              {
                node_data[kk] = 0.0;
              }
          }
      }

    SubDimCellToDataMap& map = nodeRegistry.getMap();
    SubDimCellToDataMap::iterator iter;

    for (iter = map.begin(); iter != map.end(); ++iter)
      {
        const SubDimCell_SDCEntityType& subDimEntity = (*iter).first;
        SubDimCellData& nodeId_elementOwnerId = (*iter).second;

        // tuple storage: SDC_DATA_GLOBAL_NODE_IDS, SDC_DATA_OWNING_ELEMENT_KEY,  SDC_DATA_OWNING_SUBDIM_RANK, SDC_DATA_OWNING_SUBDIM_ORDINAL, SDC_DATA_SPACING
        NodeIdsOnSubDimEntityType& nodeIds_onSE      = std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId);

        stk::mesh::EntityKey       owningElementKey  = std::get<SDC_DATA_OWNING_ELEMENT_KEY>(nodeId_elementOwnerId);
        stk::mesh::EntityId        owningElementId   = owningElementKey.id();
        (void)owningElementId;
        stk::mesh::EntityRank      owningElementRank = owningElementKey.rank();
        //VERIFY_OP_ON(owningElementRank, ==, eMesh.element_rank(), "bad owningElementRank");
        unsigned               owningSubDimRank  = std::get<SDC_DATA_OWNING_SUBDIM_RANK>(nodeId_elementOwnerId);
        unsigned               owningSubDimOrd   = std::get<SDC_DATA_OWNING_SUBDIM_ORDINAL>(nodeId_elementOwnerId);
        VERIFY_OP_ON(owningSubDimOrd, >, 0, "hmm 2");
        --owningSubDimOrd;
        unsigned owningSubDimSize = subDimEntity.size();
        (void)owningSubDimSize;

        stk::mesh::Entity owningElement = eMesh.get_bulk_data()->get_entity(owningElementKey.rank(), owningElementKey.id());
        if (!eMesh.is_valid(owningElement))
          {
            continue;
          }

        for (unsigned jj = 0; jj < nodeIds_onSE.size(); ++jj)
          {
            stk::mesh::Entity node = nodeIds_onSE[jj];

            if (!eMesh.is_valid(node))
              {
                continue;
              }

            VERIFY_OP_ON(eMesh.entity_rank(node), ==, eMesh.node_rank(), "bad rank");
            NodeRegistryFieldType::value_type *node_data = stk::mesh::field_data(*eMesh.m_node_registry_field, node);
            if (!node_data)
              {
                continue;
              }

            node_data[NR_FIELD_OWNING_ELEMENT_ID]     = static_cast<NodeRegistryFieldType::value_type>(owningElementKey.id());
            node_data[NR_FIELD_OWNING_ELEMENT_RANK]   = static_cast<NodeRegistryFieldType::value_type>(owningElementRank);
            node_data[NR_FIELD_MARK]                  = static_cast<NodeRegistryFieldType::value_type>(nodeIds_onSE.m_mark);
            node_data[NR_FIELD_OWNING_SUBDIM_RANK]    = static_cast<NodeRegistryFieldType::value_type>(owningSubDimRank);
            node_data[NR_FIELD_OWNING_SUBDIM_ORDINAL] = static_cast<NodeRegistryFieldType::value_type>(owningSubDimOrd + 1);
            node_data[NR_FIELD_OWNING_SUBDIM_SIZE]    = static_cast<NodeRegistryFieldType::value_type>(owningSubDimSize);

          }
      }

    if (doComm)
      {
        std::vector< const stk::mesh::FieldBase *> fields;
        fields.push_back(eMesh.m_node_registry_field);
        stk::mesh::copy_owned_to_shared(*eMesh.get_bulk_data(), fields);
      }
  }

  void RefinerUtil::rebuild_node_registry(PerceptMesh& eMesh, NodeRegistry& nodeRegistry, bool initNR, PerceptMesh* eMeshNR, NodeRegistry *compareNR, bool skipEmpty)
  {
    if (eMesh.m_node_registry_field == 0)
      return;

    if (initNR)
      {
        nodeRegistry.initialize();
        nodeRegistry.getMap().clear();
        nodeRegistry.init_comm_all();
      }

    if (1)
      {
        std::vector< const stk::mesh::FieldBase *> fields;
        fields.push_back(eMesh.m_node_registry_field);
        stk::mesh::copy_owned_to_shared(*eMesh.get_bulk_data(), fields);
      }

    std::vector<stk::mesh::Entity> vec;
    stk::mesh::Selector sel = eMesh.get_fem_meta_data()->universal_part();
    stk::mesh::get_selected_entities(sel , eMesh.get_bulk_data()->buckets(eMesh.node_rank()), vec);

    SubDimCell_SDCEntityType subDimEntity(&eMesh);

    for (size_t inode=0; inode < vec.size(); ++inode)
      {
        stk::mesh::Entity node = vec[inode];

        double *node_data = stk::mesh::field_data(*eMesh.m_node_registry_field, node);

        stk::mesh::EntityId  owningElementId (  static_cast<stk::mesh::EntityId>(node_data[NR_FIELD_OWNING_ELEMENT_ID]) );
        stk::mesh::EntityRank  owningElementRank (  static_cast<stk::mesh::EntityRank>(node_data[NR_FIELD_OWNING_ELEMENT_RANK]) );
        stk::mesh::EntityKey owningElementKey(owningElementRank, owningElementId);
        stk::mesh::Entity owningElement = eMesh.get_bulk_data()->get_entity(owningElementRank, owningElementId);
        if (!eMesh.is_valid(owningElement))
          {
            continue;
          }
        VERIFY_OP_ON(eMesh.is_valid(owningElement), ==, true, "bad owningElement");
        stk::mesh::EntityRank owningSubDimRank = static_cast<stk::mesh::EntityRank>(node_data[NR_FIELD_OWNING_SUBDIM_RANK]);
        unsigned              owningSubDimOrd  = static_cast<unsigned>(node_data[NR_FIELD_OWNING_SUBDIM_ORDINAL]);
        VERIFY_OP_ON(owningSubDimOrd, >, 0, "hmmm 3");
        --owningSubDimOrd;
        unsigned              owningSubDimSize = static_cast<unsigned>(node_data[NR_FIELD_OWNING_SUBDIM_SIZE]);

        bool foundGhostNode = nodeRegistry.getSubDimEntity(subDimEntity, owningElement, owningSubDimRank, owningSubDimOrd);
        if (foundGhostNode && nodeRegistry.getCheckForGhostedNodes())
          {
            //std::cout << eMesh.rank() << " rebuild_node_registry foundGhostNode= " << foundGhostNode << std::endl;
            continue;
          }
        VERIFY_OP_ON(subDimEntity.size(), ==, owningSubDimSize, "bad owningSubDimSize");

        for (unsigned ii=0; ii < subDimEntity.size(); ++ii)
          {
            VERIFY_OP_ON(eMesh.is_valid(subDimEntity[ii]), ==, true, "bad node in rebuild_node_registry");
          }
        static SubDimCellData empty_SubDimCellData;
        SubDimCellData* nodeId_elementOwnerId_ptr = nodeRegistry.getFromMapPtr(subDimEntity);
        SubDimCellData& nodeId_elementOwnerId = (nodeId_elementOwnerId_ptr ? *nodeId_elementOwnerId_ptr : empty_SubDimCellData);
        bool is_empty = nodeId_elementOwnerId_ptr == 0;

        if (is_empty)
          {
            SubDimCellData data { NodeIdsOnSubDimEntityType(1, stk::mesh::Entity(), nodeRegistry.NR_MARK_NONE),
                owningElementKey, (unsigned char)owningSubDimRank, (unsigned char)(owningSubDimOrd + 1), {}};
            NodeIdsOnSubDimEntityType& nid_new = std::get<SDC_DATA_GLOBAL_NODE_IDS>(data);
            nid_new.resize(1);
            nid_new[0] = node;
            nid_new.m_entity_id_vector[0] = eMesh.id(node);
            nid_new.m_mark = static_cast<unsigned>(node_data[NR_FIELD_MARK]);

            nodeRegistry.putInMap(subDimEntity,  data);
          }
        else
          {
            NodeIdsOnSubDimEntityType& nodeIds_onSE = std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId);
            nodeIds_onSE.push_back(node);
            nodeIds_onSE.m_entity_id_vector.push_back(eMesh.id(node));
            nodeIds_onSE.m_mark = static_cast<unsigned>(node_data[NR_FIELD_MARK]);
            VERIFY_OP_ON(nodeIds_onSE.size(), ==, 1, "bad size on proc: "+eMesh.rank());
          }
      }
  }

}

