#include <gtest/gtest.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>

#include <stk_unit_test_utils/ioUtils.hpp>

#include "TetFixture.hpp"

TEST_F(DGTetFixture, tet)
{
    std::vector<stk::mesh::EntityIdVector> tet_conn = {
            {1, 2, 3, 4}, // id 1
            {2, 3, 4, 5}  // id 2
    };

    std::vector< std::vector<double> > node_coords= {
            {0, 0, 0}, // 1
            {1, 0, 0}, // 2
            {0, 1, 0}, // 3
            {0.5, 0.5, 1.0}, // 6...just kidding, it's 4
            {1.0, 1.0, 1.0}
    };

    setup_mesh(tet_conn, node_coords);

    stk::unit_test_util::write_mesh_using_stk_io("mike.g", get_bulk(), get_bulk().parallel());

    //////////////////////////////////////////////////////////////////////////////////////

    stk::mesh::EntityVector elements;
    stk::mesh::get_selected_entities(get_meta().locally_owned_part(), get_bulk().buckets(stk::topology::ELEM_RANK), elements);

    std::cerr << "num elements: " << elements.size() << std::endl;

    stk::mesh::create_exposed_block_boundary_sides(get_bulk(), get_meta().locally_owned_part(), {get_skin_part()});

    unsigned num_faces = get_bulk().num_faces(elements[0]);
    const stk::mesh::Entity* faces = get_bulk().begin_faces(elements[0]);

    std::cerr << "num faces: " << num_faces << std::endl;

    for(unsigned i=0;i<num_faces;i++)
    {
        stk::mesh::Entity face = faces[i];
        unsigned num_nodes = get_bulk().num_nodes(face);
        const stk::mesh::Entity* nodes = get_bulk().begin_nodes(face);
        for(unsigned j=0;j<num_nodes;++j)
        {
            std::cerr << "Node " << j+1 << " of face " << i+1 << " is " << get_bulk().identifier(nodes[j]) << std::endl;
            double *nodeCoord = static_cast<double*>(stk::mesh::field_data(*get_coord_field(), nodes[j]));
            std::cerr << "Has coordinates: " << nodeCoord[0] << "    " << nodeCoord[1] << "    " << nodeCoord[2] << std::endl;
        }
    }
}

