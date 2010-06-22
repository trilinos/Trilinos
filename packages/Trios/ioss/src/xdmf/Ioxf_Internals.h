/*--------------------------------------------------------------------
 *    Copyright 2001 Sandia Corporation.
 *    Under the terms of Contract DE-AC04-94AL85000, there is a
 *    non-exclusive license for use of this work by or on behalf
 *    of the U.S. Government.  Export of this program may require
 *    a license from the United States Government.
 *--------------------------------------------------------------------
 */
// Author: gdsjaar@sandia.gov

#ifndef IOSS_Ioxf_Internals_h
#define IOSS_Ioxf_Internals_h

#include <vector>
#include <cstring>

#ifndef MAX_LINE_LENGTH
#define MAX_LINE_LENGTH 128
#endif

#ifndef MAX_STR_LENGTH
#define MAX_STR_LENGTH 32
#endif

namespace Ioss {
  class ElementBlock;
  class NodeSet;
  class FaceBlock;
  class EdgeBlock;
}
namespace Ioxf {
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
    SideSet(const Ioss::FaceBlock &other); // 3D
    SideSet(const Ioss::EdgeBlock &other); // 2D
    bool operator==(const SideSet&) const;
    bool operator!=(const SideSet& other) const {return !(*this == other);}

    char name[MAX_STR_LENGTH+1];
    int id;
    int sideCount;
    int dfCount;
    int elemCount;
    int nodesPerSideSet;
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
}
#endif /* IOSS_Ioxf_Internals_h */
