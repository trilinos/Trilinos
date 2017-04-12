#ifndef STK_BALANCE_MTONREBALANCER_HPP
#define STK_BALANCE_MTONREBALANCER_HPP

#include <string>
#include <vector>
#include <stk_mesh/base/Types.hpp>
#include <stk_io/IossBridge.hpp>

namespace stk { namespace mesh { class MetaData; }}
namespace stk { namespace mesh { class BulkData; }}
namespace stk { namespace mesh { class FieldBase; }}
namespace stk { namespace mesh { class Bucket; }}

namespace stk {
namespace balance {
namespace internal {

class MtoNRebalancer
{
public:
    MtoNRebalancer(stk::mesh::BulkData &bulkData, int num_target_procs);
    ~MtoNRebalancer();

    void generate_n_proc_decomp();
    void move_subdomains_such_that_entire_subdomain_doesnt_span_proc_boundaries(const std::vector<unsigned>& target_proc_to_starting_proc);
    stk::io::EntitySharingInfo get_node_sharing_info(unsigned subdomain);
    void create_subdomain_and_write(const std::string &filename, unsigned subdomain, int global_num_nodes, int global_num_elems, const stk::io::EntitySharingInfo &nodeSharingInfo, int numSteps = -1, double timeStep = 0.0);
    bool does_this_proc_own_subdomain(unsigned subdomainOwner);

private:
    void move_entities_into_mapped_subdomain_parts(const std::vector<unsigned>& mappings);
    void change_parts_on_entities_on_all_subdomains(const std::vector<unsigned>& subdomain_proc_mapping);
    void move_entities_into_subdomain_part(size_t i, const stk::mesh::EntityVector &entities);
    stk::mesh::PartVector get_parts_to_add_for_subdomain(size_t subdomain_num);
    stk::mesh::Part* get_subdomain_part(size_t subdomain_num);
    std::string getSubdomainPartName(int subdomainId);
    std::vector<stk::mesh::Entity> get_entities_for_subdomain(size_t subdomain_num);
    stk::mesh::EntityVector get_entitites_for_subdomain_using_field_from_buckets(size_t subdomain_num, const stk::mesh::BucketVector& buckets);
    void add_owned_entities_from_bucket_using_target_decomp_field(const stk::mesh::Bucket& bucket, size_t subdomain_num, stk::mesh::EntityVector& entities);
    void add_entities_from_bucket_using_target_decomp_field(const stk::mesh::Bucket& bucket, size_t subdomain_num, stk::mesh::EntityVector& entities);
    void declare_all_subdomain_parts();
    void store_off_target_proc_on_elements_before_moving_subdomains();
    stk::mesh::EntityVector get_nodes_shared_between_subdomains(int this_subdomain_index, int other_subdomain_index);
    void fill_shared_node_proc_info(stk::mesh::EntityVector& shared_nodes, std::vector<int>& procs_for_shared_nodes, int this_subdomain_num, int other_subdomain_num);
    void fill_shared_node_info_for_this_subdomain(const unsigned this_subdomain_num, stk::mesh::EntityVector& shared_nodes, std::vector<int>& procs_for_shared_nodes);
    stk::mesh::FieldBase* get_target_decomp_field();
    const std::string get_target_decomp_field_name();
    stk::mesh::MetaData& get_meta();
    stk::mesh::BulkData& get_bulk();

private:
    MtoNRebalancer( const MtoNRebalancer& other );
    MtoNRebalancer& operator=( const MtoNRebalancer& other );

    stk::mesh::BulkData &mBulkData;
    int mNumTargetProcs;
    stk::mesh::EntityProcVec decomp;

};

}}}
#endif
