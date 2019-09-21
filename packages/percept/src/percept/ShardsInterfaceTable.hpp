// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_ShardsInterfaceTable_hpp
#define percept_ShardsInterfaceTable_hpp

#include <string.h>
#include <map>

#include <Shards_Array.hpp>
#include <Shards_ArrayVector.hpp>
#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Stencils.hpp>
#include <stk_mesh/base/TopologyDimensions.hpp>


  namespace percept {
    namespace interface_table {

      typedef void (*stk_set_cell_topology_fptr)(stk::mesh::Part & );

 
      enum  ElemTypes { 
        shards_Node                     ,  
        shards_Particle                 ,  

        shards_Line_2                   ,  
        shards_ShellLine_2              ,  
        shards_ShellLine_3              ,  
        shards_Beam_2                   ,  
        shards_Beam_3                   ,  

        shards_Triangle_3               ,  
        shards_Triangle_6               ,  
        shards_ShellTriangle_3          ,  
        shards_ShellTriangle_6          ,  

        shards_Quadrilateral_4          ,  
        shards_Quadrilateral_8          ,  
        shards_Quadrilateral_9          ,  
        shards_ShellQuadrilateral_4     ,  
        shards_ShellQuadrilateral_8     ,  
        shards_ShellQuadrilateral_9     ,  

        shards_Pentagon_5               ,  
        shards_Hexagon_6                ,  

        shards_Tetrahedron_4            ,  
        shards_Tetrahedron_10           ,  

        shards_Pyramid_5                ,  
        shards_Pyramid_13               ,  
        shards_Pyramid_14               ,  

        shards_Wedge_6                  ,  
        shards_Wedge_15                 ,  
        shards_Wedge_18                 ,  

        shards_Hexahedron_8             ,  
        shards_Hexahedron_20            ,  
        shards_Hexahedron_27            ,
 

        NUM_ELEM_TYPES };

      typedef struct
      {
        //unsigned elemEnumType;
        ElemTypes elemEnumType;
        const CellTopologyData *cellTopoData;
        const char *name;
        unsigned vertex_count;
        unsigned node_count;
        int sweptElemType;
        stk_set_cell_topology_fptr setCellTopoFptr;
      }  elemInfoType;


    } // namespace interface_table

    using namespace interface_table;

    struct lstr {
      bool operator()(const char *a, const char *b) const { return strcmp(a,b)<0; }
    };

    class ShardsInterfaceTable
    {
      ShardsInterfaceTable() {
        checkTable();
        buildMaps();
      }
      void checkTable();
      void buildMaps();
      static  std::map<const char *, int, lstr> s_nameToIdMap;

    public:
      static ShardsInterfaceTable s_singleton;
      static  elemInfoType *s_elemInfo;
      int lookupShardsId(const char *);

    }; //class ShardsInterfaceTable

  } // namespace percept

#endif
