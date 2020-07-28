// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef Iodw_MetaData_h
#define Iodw_MetaData_h

#include <map>
#include <string>
#include <vector>

namespace Iodw {

  namespace meta {

    struct Property
    {
    };

    struct Field
    {
    };

    struct GroupingEntity
    {
      std::vector<Property> properties;
      std::vector<Field>    fields;
    };

    struct EntityBlock : public GroupingEntity
    {
    };

    struct EntitySet : public GroupingEntity
    {
    };

    struct NodeBlock : public EntityBlock
    {
    };
    struct SideBlock : public EntityBlock
    {
    };
    struct ElementBlock : public EntityBlock
    {
    };
    struct EdgeBlock : public EntityBlock
    {
    };
    struct FaceBlock : public EntityBlock
    {
    };

    struct NodeSet : public EntitySet
    {
    };
    struct SideSet : public EntitySet
    {
    };
    struct ElementSet : public EntitySet
    {
    };
    struct EdgeSet : public EntitySet
    {
    };
    struct FaceSet : public EntitySet
    {
    };

    struct Region
    {
      using Name       = std::string;
      using Key        = std::string;
      using NameKeyMap = std::map<Name, Key>;

      NameKeyMap name_key_map;

      std::vector<ElementBlock> element_blocks;
      std::vector<NodeBlock>    node_blocks;

      std::vector<ElementSet> element_sets;
      std::vector<NodeSet>    node_sets;
    };

  } // namespace meta

} // namespace Iodw

#endif // Iodw_MetaData_h
