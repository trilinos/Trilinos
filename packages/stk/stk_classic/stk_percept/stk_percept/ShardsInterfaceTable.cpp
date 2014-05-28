#include "ShardsInterfaceTable.hpp"
#include <iostream>
#include <stdexcept>

using namespace shards;

namespace stk {
  namespace percept {


    static  interface_table::elemInfoType elemInfo[NUM_ELEM_TYPES] = {
      { shards_Node                 ,  getCellTopologyData< shards::Node  >(), 0,                 shards::Node::vertex_count,                  shards::Node::node_count                  , shards_Line_2 , 0},
      { shards_Particle             ,  getCellTopologyData< shards::Particle >(), 0,               shards::Particle::vertex_count,              shards::Particle::node_count              , -1 , 0},

      { shards_Line_2               ,  getCellTopologyData< shards::Line<2> >(), 0,                shards::Line<2>::vertex_count,               shards::Line<2>::node_count               , shards_Quadrilateral_4 ,
        stk::mesh::fem::set_cell_topology< shards::Line<2> > },

      { shards_ShellLine_2          ,  getCellTopologyData< shards::ShellLine<2> >(), 0,           shards::ShellLine<2>::vertex_count,          shards::ShellLine<2>::node_count          , -1 , 0},
      { shards_ShellLine_3          ,  getCellTopologyData< shards::ShellLine<3> >(), 0,           shards::ShellLine<3>::vertex_count,          shards::ShellLine<3>::node_count          , -1 , 0},
      { shards_Beam_2               ,  getCellTopologyData< shards::Beam<2> >(), 0,                shards::Beam<2>::vertex_count,               shards::Beam<2>::node_count               , -1 , 0},
      { shards_Beam_3               ,  getCellTopologyData< shards::Beam<3> >(), 0,                shards::Beam<3>::vertex_count,               shards::Beam<3>::node_count               , -1 , 0},

      { shards_Triangle_3           ,  getCellTopologyData< shards::Triangle<3> >(), 0,            shards::Triangle<3>::vertex_count,           shards::Triangle<3>::node_count           , shards_Wedge_6 ,
        stk::mesh::fem::set_cell_topology< shards::Triangle<3> > },
      { shards_Triangle_6           ,  getCellTopologyData< shards::Triangle<6> >(), 0,            shards::Triangle<6>::vertex_count,           shards::Triangle<6>::node_count           , -1 , 0},
      { shards_ShellTriangle_3      ,  getCellTopologyData< shards::ShellTriangle<3> >(), 0,       shards::ShellTriangle<3>::vertex_count,      shards::ShellTriangle<3>::node_count      , -1 , 0},
      { shards_ShellTriangle_6      ,  getCellTopologyData< shards::ShellTriangle<6> >(), 0,       shards::ShellTriangle<6>::vertex_count,      shards::ShellTriangle<6>::node_count      , -1 , 0},

      { shards_Quadrilateral_4      ,  getCellTopologyData< shards::Quadrilateral<4> >(), 0,       shards::Quadrilateral<4>::vertex_count,      shards::Quadrilateral<4>::node_count      , shards_Hexahedron_8 ,
        stk::mesh::fem::set_cell_topology< shards::Quadrilateral<4> > },
      { shards_Quadrilateral_8      ,  getCellTopologyData< shards::Quadrilateral<8> >(), 0,       shards::Quadrilateral<8>::vertex_count,      shards::Quadrilateral<8>::node_count      , -1 , 0},
      { shards_Quadrilateral_9      ,  getCellTopologyData< shards::Quadrilateral<9> >(), 0,       shards::Quadrilateral<9>::vertex_count,      shards::Quadrilateral<9>::node_count      , -1 , 0},
      { shards_ShellQuadrilateral_4 ,  getCellTopologyData< shards::ShellQuadrilateral<4> >(), 0,  shards::ShellQuadrilateral<4>::vertex_count, shards::ShellQuadrilateral<4>::node_count , -1 , 0},
      { shards_ShellQuadrilateral_8 ,  getCellTopologyData< shards::ShellQuadrilateral<8> >(), 0,  shards::ShellQuadrilateral<8>::vertex_count, shards::ShellQuadrilateral<8>::node_count , -1 , 0},
      { shards_ShellQuadrilateral_9 ,  getCellTopologyData< shards::ShellQuadrilateral<9> >(), 0,  shards::ShellQuadrilateral<9>::vertex_count, shards::ShellQuadrilateral<9>::node_count , -1 , 0},

      { shards_Pentagon_5           ,  getCellTopologyData< shards::Pentagon<5> >(), 0,            shards::Pentagon<5>::vertex_count,           shards::Pentagon<5>::node_count           , -1 , 0},
      { shards_Hexagon_6            ,  getCellTopologyData< shards::Hexagon<6> >(), 0,             shards::Hexagon<6>::vertex_count,            shards::Hexagon<6>::node_count            , -1 , 0},

      { shards_Tetrahedron_4        ,  getCellTopologyData< shards::Tetrahedron<4> >(), 0,         shards::Tetrahedron<4>::vertex_count,        shards::Tetrahedron<4>::node_count        , -1 ,
        stk::mesh::fem::set_cell_topology< shards::Tetrahedron<4> >},
      { shards_Tetrahedron_10       ,  getCellTopologyData< shards::Tetrahedron<10> >(), 0,        shards::Tetrahedron<10>::vertex_count,       shards::Tetrahedron<10>::node_count       , -1 , 0},

      { shards_Pyramid_5            ,  getCellTopologyData< shards::Pyramid<5> >(), 0,             shards::Pyramid<5>::vertex_count,            shards::Pyramid<5>::node_count            , -1 , 0},
      { shards_Pyramid_13           ,  getCellTopologyData< shards::Pyramid<13> >(), 0,            shards::Pyramid<13>::vertex_count,           shards::Pyramid<13>::node_count           , -1 , 0},
      { shards_Pyramid_14           ,  getCellTopologyData< shards::Pyramid<14> >(), 0,            shards::Pyramid<14>::vertex_count,           shards::Pyramid<14>::node_count           , -1 , 0},

      { shards_Wedge_6              ,  getCellTopologyData< shards::Wedge<6> >(), 0,               shards::Wedge<6>::vertex_count,              shards::Wedge<6>::node_count              , -1 ,
        stk::mesh::fem::set_cell_topology< shards::Wedge<6> >},
      { shards_Wedge_15             ,  getCellTopologyData< shards::Wedge<15> >(), 0,              shards::Wedge<15>::vertex_count,             shards::Wedge<15>::node_count             , -1 , 0},
      { shards_Wedge_18             ,  getCellTopologyData< shards::Wedge<18> >(), 0,              shards::Wedge<18>::vertex_count,             shards::Wedge<18>::node_count             , -1 , 0},

      { shards_Hexahedron_8         ,  getCellTopologyData< shards::Hexahedron<8> >(), 0,          shards::Hexahedron<8>::vertex_count,         shards::Hexahedron<8>::node_count         , -1 ,
        stk::mesh::fem::set_cell_topology< shards::Hexahedron<8> > },
      { shards_Hexahedron_20        ,  getCellTopologyData< shards::Hexahedron<20> >(), 0,         shards::Hexahedron<20>::vertex_count,        shards::Hexahedron<20>::node_count        , -1 , 0},
      { shards_Hexahedron_27        ,  getCellTopologyData< shards::Hexahedron<27> >(), 0,         shards::Hexahedron<27>::vertex_count,        shards::Hexahedron<27>::node_count        , -1 , 0}

    };

    //    ShardsInterfaceTable::s_elemInfo = elemInfo;
    interface_table::elemInfoType *ShardsInterfaceTable::s_elemInfo = elemInfo;

    std::map<const char *, int, lstr> ShardsInterfaceTable::s_nameToIdMap;

    ShardsInterfaceTable ShardsInterfaceTable::s_singleton;

    void ShardsInterfaceTable::buildMaps()
    {
      //static std::map<const char *, int, lstr> nameToIdMap;
      for (int i = 0; i < NUM_ELEM_TYPES; i++)
        {
          s_elemInfo[i].name = s_elemInfo[i].cellTopoData->name;
        }
      for (int i = 0; i < NUM_ELEM_TYPES; i++)
        {
          s_nameToIdMap[s_elemInfo[i].name] = i;
        }
      //s_nameToIdMap = nameToIdMap;
    }

    int ShardsInterfaceTable::lookupShardsId(const char *name)
    {
      return s_nameToIdMap[name];
    }

    void ShardsInterfaceTable::checkTable()
    {
      for (int i = 0; i < NUM_ELEM_TYPES; i++)
        {
          //s_elemInfo[i]
          if (s_elemInfo[i].elemEnumType != i)
            {
              std::ostringstream msg;
              msg << "ShardsInterfaceTable error - check table's code\n";
              throw std::runtime_error( msg.str() );
            }
        }
    }

  } // namespace percept
} // namespace stk
