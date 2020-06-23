#include "gtest/gtest.h"
#include "mpi.h"
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/base/ForEachEntity.hpp>

namespace {

class EntityLessCoords
{
public:
    EntityLessCoords(const stk::mesh::BulkData& bulk)
    :mesh(bulk)
    ,coordsBase(bulk.mesh_meta_data().coordinate_field())
    { }

    bool operator()(stk::mesh::Entity a, stk::mesh::Entity b)
    {
        if (mesh.entity_rank(a) != stk::topology::NODE_RANK) {
            return stk::mesh::EntityLess(mesh)(a,b);
        }
        double* aCoords = static_cast<double*>(stk::mesh::field_data(*coordsBase, a));
        double* bCoords = static_cast<double*>(stk::mesh::field_data(*coordsBase, b));
        if (less_than(aCoords[0], bCoords[0]))
        {
            return true;
        }
        else if (equal(aCoords[0], bCoords[0]))
        {
            if (less_than(aCoords[1], bCoords[1]))
            {
                return true;
            }
            else if (equal(aCoords[1], bCoords[1]))
            {
                if (less_than(aCoords[2], bCoords[2]))
                {
                    return true;
                }
            }
        }
        return false;
    }

    bool equal(const double& val1, const double& val2) const
    {
        return (abs(val1 - val2) < tolerance);
    }

    bool less_than(const double& val1, const double& val2) const
    {
        if (equal(val1, val2))
            return false;
        return (val1 < val2);
    }

    void set_tolerance(double toleranceIn)
    {
        tolerance = toleranceIn;
    }

private:
    const stk::mesh::BulkData& mesh;
    const stk::mesh::FieldBase* coordsBase;
    double tolerance = 1.0e-6;
};

class EntityCoordSorter : public stk::mesh::EntitySorterBase
{
public:
    virtual void sort(stk::mesh::BulkData &bulk, stk::mesh::EntityVector& entityVector) const
    {
        std::sort(entityVector.begin(), entityVector.end(), EntityLessCoords(bulk));
    }
};


class EntitySortingFixture : public stk::unit_test_util::MeshFixture
{
protected:

    EntitySortingFixture()
    {
        setup_mesh("generated:1x1x4", stk::mesh::BulkData::AUTO_AURA);
    }

    void verify_node_ordering()
    {
        EntityLessCoords entityLessCoords(get_bulk());
        stk::mesh::for_each_entity_run(get_bulk(), stk::topology::NODE_RANK,
            [&entityLessCoords](const stk::mesh::BulkData& bulk, const stk::mesh::MeshIndex& meshIndex)
            {
                 if(meshIndex.bucket_ordinal > 0)
                     EXPECT_TRUE(entityLessCoords((*meshIndex.bucket)[meshIndex.bucket_ordinal-1],
                                                  (*meshIndex.bucket)[meshIndex.bucket_ordinal]));
            }
        );
    }

};

TEST_F(EntitySortingFixture, trivial)
{
    get_bulk().sort_entities(EntityCoordSorter());
    verify_node_ordering();
}

TEST_F(EntitySortingFixture, skin)
{
    get_bulk().sort_entities(EntityCoordSorter());
    verify_node_ordering();
    stk::mesh::Part& surface1 = get_meta().declare_part("Surface 1");
    stk::mesh::create_exposed_block_boundary_sides(get_bulk(),get_meta().universal_part(),{&surface1});
    get_bulk().sort_entities(EntityCoordSorter());
    verify_node_ordering();
}

} // namespace
