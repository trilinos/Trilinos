#include <exception>
#include <fstream>
#include <set>
#include <typeinfo>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <stk_adapt/RefinerUtil.hpp>

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
          blocks[eMesh.element_rank()].push_back(block);
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
        bool exc = false;
        if ('-' == n1[0]) 
        {
          exc = true;
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
      bool exc = false;
      if ('-' == n1[0]) 
      {
        exc = true;
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
        bool exc = false;
        if ('-' == n1[0]) 
        {
          exc = true;
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
          blocks[eMesh.element_rank()].push_back(n1);
        else if (n1.length() > 8 && n1.substr(1,8) == "surface_")
          blocks[eMesh.face_rank()].push_back(n1);
        else
        {
          std::string pm = (inc?"+":"-");
          blocks[eMesh.element_rank()].push_back(pm+"block_"+n2);
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
  if (blocks[eMesh.element_rank()].size() == 0)
    return blocks;

  stk::mesh::EntityRank subDimRank = (eMesh.getSpatialDim() == 3 ? eMesh.face_rank() : eMesh.edge_rank());

  mesh::PartVector all_parts = eMesh.getFEM_meta_data()->get_parts();
  for (mesh::PartVector::iterator i_part = all_parts.begin(); i_part != all_parts.end(); ++i_part)
  {
    mesh::Part *  part = *i_part ;

    for (mesh::PartVector::iterator i_surfacePart = all_parts.begin(); i_surfacePart != all_parts.end(); ++i_surfacePart)
    {
      mesh::Part *  surfacePart = *i_surfacePart ;
      const CellTopologyData * part_cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(*surfacePart);

      if (part_cell_topo_data && part->primary_entity_rank() == eMesh.element_rank() && surfacePart->primary_entity_rank() == subDimRank)
      {
        std::string partNamePlus = "+" + part->name();
        std::vector<std::string>::iterator partInBlocks = std::find(blocks[eMesh.element_rank()].begin(), blocks[eMesh.element_rank()].end(), partNamePlus);
        // if this part is not in the blocks list, skip it
        if (partInBlocks == blocks[eMesh.element_rank()].end())
        {
          continue;
        }
        std::string surfacePartNamePlus = "+" + surfacePart->name();
        std::vector<std::string>::iterator surfacePartInBlocks = std::find(blocks[subDimRank].begin(), blocks[subDimRank].end(), surfacePartNamePlus);
        // if this surface is already in the list, skip it
        if (surfacePartInBlocks != blocks[subDimRank].end())
        {
          continue;
        }
        if (eMesh.isBoundarySurface(*part, *surfacePart))
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

}
}
