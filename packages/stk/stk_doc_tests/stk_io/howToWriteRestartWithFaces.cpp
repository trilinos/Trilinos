#include <unistd.h>
#include <gtest/gtest.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/GetEntities.hpp>

extern int gl_argc;
extern char** gl_argv;

namespace
{

TEST(StkIoHowTo, WriteRestartWithFaceBlock)
{
    std::string filename = "output.rst";
    unsigned numStates = 1;
    int outputTimeStep = 1;
    double outputTime = 0.0;
    {
        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
        stk::mesh::Part* part = &meta.declare_part_with_topology("faceBlock", stk::topology::QUAD_4);
        stk::mesh::Field<double>& faceField = meta.declare_field<stk::mesh::Field<double> >(stk::topology::FACE_RANK, "faceField", numStates);
        stk::mesh::put_field_on_mesh(faceField, meta.universal_part(),
                                    (stk::mesh::FieldTraits<stk::mesh::Field<double> >::data_type*) nullptr);
        stk::io::put_face_block_io_part_attribute(*part);
        stk::io::fill_mesh("generated:1x1x1", bulk);
        bool connectFacesToEdges = true;
        stk::mesh::create_all_sides(bulk, meta.universal_part(), {part}, connectFacesToEdges);

        stk::mesh::EntityVector faces;
        stk::mesh::get_selected_entities(meta.universal_part(), bulk.buckets(stk::topology::FACE_RANK), faces);

        for(auto face : faces) {
          double* data = reinterpret_cast<double*>(stk::mesh::field_data(faceField, face));

          *data = bulk.identifier(face);
        }

        stk::io::StkMeshIoBroker ioBroker;
        ioBroker.set_bulk_data(bulk);
        stk::io::write_mesh_with_fields(filename, ioBroker, outputTimeStep, outputTime, stk::io::WRITE_RESTART);
    }

    {
        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
        stk::mesh::Field<double>& faceField = meta.declare_field<stk::mesh::Field<double> >(stk::topology::FACE_RANK, "faceField", numStates);
        stk::mesh::put_field_on_mesh(faceField, meta.universal_part(),
                                    (stk::mesh::FieldTraits<stk::mesh::Field<double> >::data_type*) nullptr);

        stk::io::set_field_role(faceField, Ioss::Field::TRANSIENT);

        stk::io::StkMeshIoBroker stkIo;
        stk::io::fill_mesh_with_fields(filename, stkIo, bulk, stk::io::READ_RESTART);

        int numSteps = stkIo.get_num_time_steps();
        EXPECT_EQ(1, numSteps);

        stk::mesh::EntityVector faces;
        stk::mesh::get_entities(bulk, stk::topology::FACE_RANK, faces);

        for(auto face : faces) {
            double* data = reinterpret_cast<double*>(stk::mesh::field_data(faceField, face));
            double expectedValue = bulk.identifier(face);

            EXPECT_EQ(expectedValue, *data);
        }
    }

    unlink(filename.c_str());
}

}
