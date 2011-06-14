/*
 * Copyright(C) 1999-2010
 * Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software.
 *         
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef IOSS_Ioex_Internals_h
#define IOSS_Ioex_Internals_h

#include <exodusII.h>
#include <vector>
#include <cstring>

namespace Ioss {
  class ElementBlock;
  class NodeSet;
  class SideBlock;
  class SideSet;
}
  /*!
   * This set of classes provides a thin wrapper around the exodusII
   * internals.  It supplants several of the exodusII API calls in
   * order to avoid ncredef calls which totally rewrite the existing
   * database and can be very expensive.  These routines provide all
   * required variable, dimension, and attribute definitions to the
   * underlying netcdf file with only a single ncredef call.
   *
   * To use the application must create an Ioex::Internals instance
   * and call the Ioex::Internals::write_meta_data() function.  This
   * function requires several classes as arguments including:
   * <ul>
   * <li> Mesh -- defines mesh global metadata
   * <li> Block -- defines metadata for each block
   * <li> NodeSet -- defines metadata for each nodeset
   * <li> SideSet -- defines metadata for each sideset
   * <li> CommunicationMetaData -- global metadata relating to
   * parallel info.
   * </ul>
   *
   * Calling Ioex::Internals::write_meta_data(), replaces the
   * following exodusII and nemesis API calls:
   * <ul>
   * <li> ex_put_init(),
   * <li> ex_put_elem_block(),
   * <li> ex_put_node_set_param(),
   * <li> ex_put_side_set_param(),
   * <li> ne_put_init_info(),
   * <li> ne_put_loadbal_param(),
   * <li> ne_put_cmap_params(),
   * </ul>
   */
namespace Ioex {
  struct Mesh
  {
    Mesh() :   dimensionality(0), nodeCount(0), elementCount(0),
	       blockCount(0), nodesetCount(0), sidesetCount(0) {}

    Mesh(int dim, int nodes, int elements, int blocks, int nsets, int ssets,
	 char* the_title)
      :  dimensionality(dim), nodeCount(nodes), elementCount(elements),
	 blockCount(blocks), nodesetCount(nsets), sidesetCount(ssets)
    {
      std::strncpy(title, the_title, MAX_LINE_LENGTH+1);
      title[MAX_LINE_LENGTH] = '\0';
    }

    char title[MAX_LINE_LENGTH+1];
    int dimensionality;
    int nodeCount;
    int elementCount;
    int blockCount;
    int nodesetCount;
    int sidesetCount;
  };

  struct Block
  {
    Block() : id(0), elementCount(0), nodesPerElement(0), attributeCount(0),
	      offset_(-1)
    {
      std::strcpy(elType, "");
      std::strcpy(name, "");
    }

    Block(const Block &other) : id(other.id), elementCount(other.elementCount),
				nodesPerElement(other.nodesPerElement),
				attributeCount(other.attributeCount), offset_(other.offset_)
    {
      std::strcpy(elType, other.elType);
      std::strcpy(name,   other.name);
    }

    Block(const Ioss::ElementBlock &other);

    Block& operator=(const Block& other);

    ~Block() {}

    bool operator==(const Block&) const;
    bool operator!=(const Block& other) const {return !(*this == other);}

    char elType[MAX_STR_LENGTH+1];
    char name[MAX_STR_LENGTH+1];
    int id;
    int elementCount;
    int nodesPerElement;
    int attributeCount;
    int offset_;
    private:
  };

  struct NodeSet
  {
    NodeSet() : id(0), nodeCount(0), dfCount(0)
    {
      std::strcpy(name, "");
    }
    NodeSet(const Ioss::NodeSet &other);
    bool operator==(const NodeSet&) const;
    bool operator!=(const NodeSet& other) const {return !(*this == other);}

    char name[MAX_STR_LENGTH+1];
    int id;
    int nodeCount;
    int dfCount;
  };

  struct SideSet
  {
    SideSet() : id(0), sideCount(0), dfCount(0)
    {
      std::strcpy(name, "");
    }
    SideSet(const Ioss::SideBlock &other); // 3D
    SideSet(const Ioss::SideSet   &other); // 3D
    bool operator==(const SideSet&) const;
    bool operator!=(const SideSet& other) const {return !(*this == other);}

    char name[MAX_STR_LENGTH+1];
    int id;
    int sideCount;
    int dfCount;
  };

  struct CommunicationMap
  {
    CommunicationMap() : id(0), entityCount(0), type('U') {}
    CommunicationMap(int the_id, int count, char the_type) :
      id(the_id), entityCount(count), type(the_type) {}
    bool operator==(const CommunicationMap&) const;
    bool operator!=(const CommunicationMap& other) const {return !(*this == other);}
    int id;
    int entityCount;
    char type; // 'n' for node, 'e' for element
  };

  struct CommunicationMetaData
  {
    CommunicationMetaData() : processorId(0), processorCount(0),
			      globalNodes(0), globalElements(0),
			      nodesInternal(0), nodesBorder(0), nodesExternal(0),
			      elementsInternal(0), elementsBorder(0) {}

    std::vector<CommunicationMap> nodeMap;
    std::vector<CommunicationMap> elementMap;
    int processorId;
    int processorCount;
    int globalNodes;
    int globalElements;
    int nodesInternal;
    int nodesBorder;
    int nodesExternal;
    int elementsInternal;
    int elementsBorder;

    private:
    CommunicationMetaData(const CommunicationMetaData &);
  };

  class Redefine
  {
  public:
    explicit Redefine(int exoid);
    ~Redefine();

  private:
    Redefine(const Redefine& from); // do not implement
    Redefine& operator=(const Redefine& from); // do not implement
    int exodusFilePtr;
  };

  class Internals
  {
  public:
    explicit Internals(int exoid);

    int write_meta_data(const Mesh &mesh,
			const std::vector<Block>   &blocks,
			const std::vector<NodeSet> &nodesets,
			const std::vector<SideSet> &sidesets,
			const CommunicationMetaData &comm);

  /*!  A restart file may contain an attribute which contains
   *   information about the processor count and current processor id
   *   * when the file was written.  This code checks whether that
   *   information matches the current processor count and id.  If it
   *   * exists, but doesn't match, a warning message is printed.
   *   Eventually, this will be used to determine whether certain
   *   decomposition-related data in the file is valid or has been
   *   invalidated by a join/re-spread to a different number of
   *   processors.
   */
    bool check_processor_info(int processor_count, int processor_id);

    void update_last_time_attribute(double value);
    bool read_last_time_attribute(double *value);

  private:
    Internals(const Internals& from); // do not implement
    Internals& operator=(const Internals& from); // do not implement

    int put_metadata(const Mesh &mesh,
		     const CommunicationMetaData &comm);
    int put_metadata(const std::vector<Block> &blocks);
    int put_metadata(const std::vector<NodeSet> &nodesets);
    int put_metadata(const std::vector<SideSet> &sidesets);

    int put_non_define_data(const Mesh &mesh,
			    const CommunicationMetaData &comm);
    int put_non_define_data(const std::vector<Block> &blocks);
    int put_non_define_data(const std::vector<NodeSet> &nodesets);
    int put_non_define_data(const std::vector<SideSet> &sidesets);

    int exodusFilePtr;
    int nodeMapVarID[3];
    int elementMapVarID[2];
    int commIndexVar;
    int elemCommIndexVar;
  };
}
#endif /* IOSS_Ioex_Internals_h */
