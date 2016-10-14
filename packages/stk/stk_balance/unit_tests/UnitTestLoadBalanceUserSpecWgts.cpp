#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>

#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>

#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

namespace
{

class UserSpecficWeights;

////////////////////////////////////////////////////////////////////////////////////////////

void test_weighted_balance(stk::mesh::BulkData& stkMeshBulkData, UserSpecficWeights* fixture);
std::vector<double> get_default_weights(const stk::mesh::BulkData& bulkData, UserSpecficWeights* fixture);
void balance_mesh_with_weights(stk::mesh::BulkData& bulkData, UserSpecficWeights* fixture);
void test_weight_is_balanced(const stk::mesh::BulkData& stkMeshBulkData, UserSpecficWeights* fixture);
void test_weight_is_not_balanced(const stk::mesh::BulkData& stkMeshBulkData, UserSpecficWeights* fixture);
void transform_mesh_around_y_axis(stk::mesh::BulkData& bulkData, double angle_in_degrees);

////////////////////////////////////////////////////////////////////////////////////////////

class UsingSmallEdgeWeights: public stk::balance::UserSpecifiedVertexWeightsSetting
{
    public:
    virtual double getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const { return 1.0e-6; }
};

class UserSpecficWeights : public stk::unit_test_util::MeshFixture
{
public:
    virtual std::vector<double> get_all_element_weights(size_t num_elements) = 0;
    virtual std::string get_decomp_method() const = 0;
    virtual int num_procs_for_test() const { return 2; }
};

class UserSpecficWeightsAnyMethod : public UserSpecficWeights
{
public:
    virtual std::vector<double> get_all_element_weights(size_t num_elements)
    {
        std::vector<double> weights(num_elements/2,1);
        weights[0] = num_elements/2-1;
        weights.resize(num_elements, 0);
        return weights;
    }

    virtual std::string get_decomp_method() const { return method; }
    void set_decomp_method(const std::string& input_method) { method = input_method; }
    std::string method = "rcb";
};

class UserSpecficWeightsParMetisWithManyZeroWeightsParMetis : public UserSpecficWeights
{
public:
    virtual std::vector<double> get_all_element_weights(size_t num_elements)
    {
        std::vector<double> weights(num_elements,0);
        weights[0] = 1;
        weights[num_elements-1] = 1;
        scale_weights(weights);
        return weights;
    }

    virtual std::string get_decomp_method() const { return "parmetis"; }
private:
    void scale_weights(std::vector<double>& weights)
    {
        const size_t scale = 1e5;
        for(size_t i=0;i<weights.size();++i)
        {
            if(weights[i]==0)
                weights[i]=1;
            else
                weights[i]*=scale;
        }
    }
};

/////////////////////////

TEST_F(UserSpecficWeightsAnyMethod, testLoadBalanceSomeHaveZeroWeightsAuraParMetis)
{
    setup_mesh("generated:1x1x20", stk::mesh::BulkData::AUTO_AURA);
    set_decomp_method("parmetis");
    test_weighted_balance(get_bulk(), this);
}

TEST_F(UserSpecficWeightsAnyMethod, testLoadBalanceSomeHaveZeroWeightsWithoutAuraParMetis)
{
    setup_mesh("generated:1x1x20", stk::mesh::BulkData::NO_AUTO_AURA);
    set_decomp_method("parmetis");
    test_weighted_balance(get_bulk(), this);
}

// Weirdness of following:
// 1) RCB - does not balance
// 2) RIB - does not balance with rotated mesh
// 3) RCB - does not balance with rotated mesh

TEST_F(UserSpecficWeightsAnyMethod, DISABLED_testLoadBalanceSomeHaveZeroWeightsAuraRCB)
{
    setup_mesh("generated:1x1x20", stk::mesh::BulkData::NO_AUTO_AURA);
    set_decomp_method("rcb");
    test_weighted_balance(get_bulk(), this);
}

TEST_F(UserSpecficWeightsAnyMethod, DISABLED_testLoadBalanceSomeHaveZeroWeightsRotatedMeshAuraRIB)
{
    setup_mesh("generated:1x1x20", stk::mesh::BulkData::NO_AUTO_AURA);
    double angle_in_degrees = 45;
    transform_mesh_around_y_axis(get_bulk(), angle_in_degrees);
    set_decomp_method("rib");
    test_weighted_balance(get_bulk(), this);
}

TEST_F(UserSpecficWeightsAnyMethod, DISABLED_testLoadBalanceSomeHaveZeroWeightsRotatedMeshAuraRCB)
{
    setup_mesh("generated:1x1x20", stk::mesh::BulkData::NO_AUTO_AURA);
    double angle_in_degrees = 45;
    transform_mesh_around_y_axis(get_bulk(), angle_in_degrees);
    set_decomp_method("rcb");
    test_weighted_balance(get_bulk(), this);
}

/////////////////////////

// weirdness: use integer weights for load balancing

/////////////////////////

// weirdness: vertices that have non-zero weights (element 1 and 20) are connected to other elements
// which have zero weight. Somehow this puts all both elements on same proc. But if you scale the weights, good stuff.

TEST_F(UserSpecficWeightsParMetisWithManyZeroWeightsParMetis, DISABLED_testLoadBalanceWithoutAuraParMetis)
{
    setup_mesh("generated:1x1x20", stk::mesh::BulkData::NO_AUTO_AURA);
    test_weighted_balance(get_bulk(), this);
}

const stk::mesh::FieldBase * get_coordinate_field(const stk::mesh::MetaData& meta_data, const std::string& coordinateFieldName)
{
    const stk::mesh::FieldBase * coord = meta_data.get_field(stk::topology::NODE_RANK, coordinateFieldName);
    ThrowRequireMsg(coord != nullptr, "Null coordinate field for name=" << coordinateFieldName << ". Contact sierra-help@sandia.gov for support.");
    return coord;
}

void make_coordinate_data_parallel_consistent(stk::mesh::BulkData& bulkData, const stk::mesh::FieldBase* coord)
{
    std::vector< const stk::mesh::FieldBase *> fields;
    fields.push_back(coord);
    stk::mesh::communicate_field_data(bulkData, fields);
}

std::string get_coordinates_field_name()
{
    return "coordinates";
}

double convert_angle_from_degrees_to_radians(const double angle_in_degrees)
{
    return 3.14159*angle_in_degrees/180;
}

std::pair<double, double> get_rotated_coords_about_y(const double angle_in_degrees, double* xyz)
{
    double angle_in_radians = convert_angle_from_degrees_to_radians(angle_in_degrees);
    double new_x = cos(angle_in_radians)*xyz[0] + sin(angle_in_radians)*xyz[2];
    double new_z = -sin(angle_in_radians)*xyz[1] + cos(angle_in_radians)*xyz[2];
    return std::make_pair(new_x, new_z);
}

void update_coords_for_node(const double angle_in_degrees, const stk::mesh::FieldBase* coord, stk::mesh::Entity node)
{
    double *xyz = static_cast<double *>(stk::mesh::field_data(*coord, node));
    std::pair<double, double> new_coords = get_rotated_coords_about_y(angle_in_degrees, xyz);
    xyz[0] = new_coords.first;
    xyz[2] = new_coords.second;
}

void transform_coordinates_of_buckets(const double angle_in_degrees, const stk::mesh::FieldBase* coord, const stk::mesh::BucketVector& buckets)
{
    for(size_t i=0;i<buckets.size();++i)
    {
        const stk::mesh::Bucket& bucket = *buckets[i];
        for(size_t j=0;j<bucket.size();++j)
            update_coords_for_node(angle_in_degrees, coord, bucket[j]);
    }
}

void transform_mesh_around_y_axis(stk::mesh::BulkData& bulkData, double angle_in_degrees)
{
    const stk::mesh::FieldBase* coord = get_coordinate_field(bulkData.mesh_meta_data(), get_coordinates_field_name());
    const stk::mesh::BucketVector& buckets = bulkData.get_buckets(stk::topology::NODE_RANK, bulkData.mesh_meta_data().locally_owned_part());
    transform_coordinates_of_buckets(angle_in_degrees, coord, buckets);
    make_coordinate_data_parallel_consistent(bulkData, coord);
}

void test_weighted_balance(stk::mesh::BulkData& stkMeshBulkData, UserSpecficWeights* fixture)
{
    if(stkMeshBulkData.parallel_size() == fixture->num_procs_for_test())
    {
        test_weight_is_not_balanced(stkMeshBulkData, fixture);
        balance_mesh_with_weights(stkMeshBulkData, fixture);
        test_weight_is_balanced(stkMeshBulkData, fixture);
    }
}

void balance_mesh_with_weights(stk::mesh::BulkData& bulkData, UserSpecficWeights* fixture)
{
    stk::balance::UserSpecifiedVertexWeightsSetting graphSettings;
    std::vector<double> weights = get_default_weights(bulkData, fixture);
    graphSettings.setVertexWeights(weights);
    graphSettings.setDecompMethod(fixture->get_decomp_method());
    graphSettings.setCoordinateFieldName(get_coordinates_field_name());
    stk::balance::balanceStkMesh(graphSettings, bulkData);
}

double get_total_weight(size_t num_elements, UserSpecficWeights* fixture)
{
    std::vector<double> weights = fixture->get_all_element_weights(num_elements);
    double total=0;
    for(size_t i=0;i<weights.size();++i)
    {
        total += weights[i];
    }
    return total;
}

std::vector<double> get_default_weights(const stk::mesh::BulkData& bulkData, UserSpecficWeights* fixture)
{
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(bulkData, counts);
    std::vector<double> weights_all_elements = fixture->get_all_element_weights(counts[stk::topology::ELEM_RANK]);
    std::vector<double> weights;
    const stk::mesh::BucketVector& buckets = bulkData.get_buckets(stk::topology::ELEM_RANK, bulkData.mesh_meta_data().locally_owned_part());
    for(size_t i = 0; i < buckets.size(); ++i)
    {
        const stk::mesh::Bucket& bucket = *buckets[i];
        for(size_t j = 0; j < bucket.size(); ++j)
        {
            stk::mesh::EntityId id = bulkData.identifier(bucket[j]);
            unsigned index = id - 1;
            weights.push_back(weights_all_elements[index]);
        }
    }
    return weights;
}

size_t get_global_num_elements(const stk::mesh::BulkData& bulkData)
{
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(bulkData, counts);
    return counts[stk::topology::ELEM_RANK];
}

double get_weight_this_proc(const stk::mesh::BulkData& bulkData, UserSpecficWeights* fixture)
{
    std::vector<double> weights_all_elements = fixture->get_all_element_weights(get_global_num_elements(bulkData));
    double weight_this_proc = 0;
    const stk::mesh::BucketVector& buckets = bulkData.get_buckets(stk::topology::ELEM_RANK, bulkData.mesh_meta_data().locally_owned_part());
    unsigned num_elements = 0;
    for(size_t i = 0; i < buckets.size(); ++i)
    {
        const stk::mesh::Bucket& bucket = *buckets[i];
        for(size_t j = 0; j < bucket.size(); ++j)
        {
            stk::mesh::EntityId id = bulkData.identifier(bucket[j]);
            unsigned index = id - 1;
            weight_this_proc += weights_all_elements[index];
            num_elements++;
        }
    }
    return weight_this_proc;
}

void test_weight_is_not_balanced(const stk::mesh::BulkData& bulkData, UserSpecficWeights* fixture)
{
    double weight_this_proc = get_weight_this_proc(bulkData, fixture);
    double weight_gold = get_total_weight(get_global_num_elements(bulkData), fixture)/bulkData.parallel_size();
    EXPECT_NE(weight_gold, weight_this_proc);
}

void test_weight_is_balanced(const stk::mesh::BulkData& bulkData, UserSpecficWeights* fixture)
{
    double weight_this_proc = get_weight_this_proc(bulkData, fixture);
    double weight_gold = get_total_weight(get_global_num_elements(bulkData), fixture)/bulkData.parallel_size();
    EXPECT_EQ(weight_gold, weight_this_proc);
}

}
