#ifndef MESH_FIXTURE_MXN
#define MESH_FIXTURE_MXN

#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_util/parallel/DebugTool.hpp>

#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
#include "SubMeshCreator.hpp"

#include <stk_unit_test_utils/MeshFixture.hpp>
#include <test_utils/StkMeshInterfaceToNemesis.hpp>

#include <test_utils/BalanceTestUtilities.hpp>
#include <stk_mesh/base/Comm.hpp>

inline void fill_decomp(const int num_partitions, stk::mesh::BulkData& bulk, stk::mesh::EntityProcVec &decomp)
{
    stk::balance::BasicZoltan2Settings graphSettings;
    std::vector<stk::mesh::Selector> selectors = { bulk.mesh_meta_data().locally_owned_part() };
    stk::balance::internal::callZoltan2(graphSettings, num_partitions, decomp, bulk, selectors);
}

inline stk::mesh::EntityProcVec get_element_decomp(const int num_partitions, stk::mesh::BulkData& bulk)
{
    stk::mesh::EntityProcVec decomp;
    fill_decomp(num_partitions, bulk, decomp);
    return decomp;
}

class MeshFixtureMxNRebalance : public stk::unit_test_util::MeshFixture
{
protected:
    MeshFixtureMxNRebalance()
    : MeshFixture(), decomp() {}

    void setup_and_test_balance_of_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        write_rebalanced_mxn();
        test_decomp();
    }

    void rebalance_mxn()
    {
        setup_initial_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        decomp = get_element_decomp(get_num_procs_target_decomp(), get_bulk());
    }

    void write_rebalanced_mxn()
    {
        rebalance_mxn();
        write_subdomain_files();
    }

    void setup_initial_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        create_target_decomp_field_on_entire_mesh();
        setup_mesh(get_input_mesh_file_name(), auraOption);
    }

    void create_field_on_entire_mesh(const std::string& fieldName)
    {
        stk::mesh::Field<double> &field = get_meta().declare_field<stk::mesh::Field<double> >(stk::topology::ELEMENT_RANK, fieldName, 1);
        stk::mesh::put_field(field, get_meta().universal_part());
    }

    void create_target_decomp_field_on_entire_mesh()
    {
        create_field_on_entire_mesh(get_target_decomp_field_name());
    }

    const std::string get_target_decomp_field_name()
    {
        return "TargetDecomp";
    }

    void test_decomp()
    {
        test_subdomain_files();
        test_decomp_is_balanced();
    }

    void write_subdomain_files()
    {
        std::vector<unsigned> target_proc_to_starting_proc = create_mappings(stk::parallel_machine_size(get_comm()));
        move_subdomains_such_that_entire_subdomain_doesnt_span_proc_boundaries(target_proc_to_starting_proc);
        write_out_files(target_proc_to_starting_proc);
        add_nemesis_info_to_files(target_proc_to_starting_proc);
    }

    void move_subdomains_such_that_entire_subdomain_doesnt_span_proc_boundaries(const std::vector<unsigned>& target_proc_to_starting_proc)
    {
        store_off_target_proc_on_elements_before_moving_subdomains();
        stk::balance::internal::rebalance(get_bulk(), target_proc_to_starting_proc, decomp);
        move_entities_into_mapped_subdomain_parts(target_proc_to_starting_proc);
    }

    void store_off_target_proc_on_elements_before_moving_subdomains()
    {
        balance_utils::putVectorDataIntoField(get_bulk(), decomp, get_target_decomp_field());
    }

    void test_subdomain_files()
    {
        test_num_elements_per_file();
        //clean_up_output_files();
    }

    std::vector<unsigned> create_mappings(unsigned num_procs)
    {
        std::vector<unsigned> mappings(get_num_procs_target_decomp(), std::numeric_limits<unsigned>::max());
        for(unsigned i = 0; i < get_num_procs_target_decomp(); i++)
            mappings[i] = i%num_procs;
        return mappings;
    }

    void test_decomp_is_balanced()
    {
        std::vector<unsigned> num_elements_assigned_to_each_proc_local = get_num_elements_assigned_to_each_proc();
        std::vector<unsigned> num_elements_assigned_to_each_proc_global = sum_across_procs(num_elements_assigned_to_each_proc_local);
        test_num_elements_per_proc(num_elements_assigned_to_each_proc_global);
    }

    void test_num_elements_per_proc(const std::vector<unsigned>& num_elements_assigned_to_each_proc_global)
    {
        for(size_t i=0;i<num_elements_assigned_to_each_proc_global.size();++i)
            EXPECT_EQ(get_num_elements_per_proc_for_target_decomp(), num_elements_assigned_to_each_proc_global[i]);
    }

    std::vector<unsigned> sum_across_procs(std::vector<unsigned>& num_elements_assigned_to_each_proc_local)
    {
        std::vector<unsigned> num_elements_assigned_to_each_proc_global(get_num_procs_target_decomp(),0);
        stk::all_reduce_sum(get_comm(), num_elements_assigned_to_each_proc_local.data(), num_elements_assigned_to_each_proc_global.data(), get_num_procs_target_decomp());
        return num_elements_assigned_to_each_proc_global;
    }

    std::vector<unsigned> get_num_elements_assigned_to_each_proc()
    {
        std::vector<unsigned> num_elements_assigned_to_each_proc_local(get_num_procs_target_decomp(),0);
        for(size_t i=0;i<decomp.size();++i)
            num_elements_assigned_to_each_proc_local[decomp[i].second]++;
        return num_elements_assigned_to_each_proc_local;
    }

    std::string getSubdomainPartName(int subdomainId)
    {
        std::ostringstream out;
        out << "subdomain_" << subdomainId;
        return out.str();
    }

    stk::mesh::FieldBase* get_target_decomp_field()
    {
        return get_meta().get_field(stk::topology::ELEMENT_RANK, get_target_decomp_field_name());
    }

    void move_entities_into_mapped_subdomain_parts(const std::vector<unsigned>& mappings)
    {
        declare_all_subdomain_parts();
        get_bulk().modification_begin();
        change_parts_on_entities_on_all_subdomains(mappings);
        get_bulk().modification_end();
    }

    void change_parts_on_entities_on_all_subdomains(const std::vector<unsigned>& subdomain_proc_mapping)
    {
        for(size_t i = 0; i < get_num_procs_target_decomp(); i++)
        {
            if(subdomain_proc_mapping[i] == static_cast<unsigned>(stk::parallel_machine_rank(get_comm())))
            {
                std::vector<stk::mesh::Entity> entities = get_entities_for_subdomain(i);
                move_entities_into_subdomain_part(i, entities);
            }
        }
    }

    void move_entities_into_subdomain_part(size_t i, const stk::mesh::EntityVector &entities)
    {
        stk::mesh::PartVector partVector = get_parts_to_add_for_subdomain(i);
        for(size_t j = 0; j < entities.size(); j++)
            get_bulk().change_entity_parts(entities[j], partVector);
    }

    std::vector<stk::mesh::Entity> get_entities_for_subdomain(size_t subdomain_num)
    {
        const stk::mesh::BucketVector &buckets = get_bulk().buckets(stk::topology::ELEMENT_RANK);
        return get_entitites_for_subdomain_using_field_from_buckets(subdomain_num, buckets);
    }

    stk::mesh::EntityVector get_entitites_for_subdomain_using_field_from_buckets(size_t subdomain_num, const stk::mesh::BucketVector& buckets)
    {
        std::vector<stk::mesh::Entity> entities;
        for(size_t j = 0; j < buckets.size(); j++)
            add_owned_entities_from_bucket_using_target_decomp_field(*buckets[j], subdomain_num, entities);
        return entities;
    }

    void add_owned_entities_from_bucket_using_target_decomp_field(const stk::mesh::Bucket& bucket, size_t subdomain_num, stk::mesh::EntityVector& entities)
    {
        if(bucket.owned())
            add_entities_from_bucket_using_target_decomp_field(bucket, subdomain_num, entities);
    }

    void add_entities_from_bucket_using_target_decomp_field(const stk::mesh::Bucket& bucket, size_t subdomain_num, stk::mesh::EntityVector& entities)
    {
        double *bucketSubdomainData = static_cast<double*>(stk::mesh::field_data(*get_target_decomp_field(), bucket));
        for(size_t k = 0; k < bucket.size(); k++)
        {
            if(bucketSubdomainData[k] == static_cast<double>(subdomain_num))
                entities.push_back(bucket[k]);
        }
    }

    void declare_all_subdomain_parts()
    {
        for(size_t i=0;i<get_num_procs_target_decomp();++i)
        {
            std::string partNameForSubdomain = getSubdomainPartName(i);
            get_meta().declare_part(partNameForSubdomain, stk::topology::ELEMENT_RANK);
        }
    }

    stk::mesh::Part* get_subdomain_part(size_t subdomain_num)
    {
        std::string partNameForSubdomain = getSubdomainPartName(subdomain_num);
        return get_meta().get_part(partNameForSubdomain);
    }

    stk::mesh::PartVector get_parts_to_add_for_subdomain(size_t subdomain_num)
    {
        stk::mesh::Part* part = get_subdomain_part(subdomain_num);
        return {part};
    }

    std::string getFilename(const std::string& filename, int numProcsDecomp, int subdomainId)
    {
        int width = balance_utils::getWidth(numProcsDecomp);
        std::ostringstream os;
        os << filename << "." << numProcsDecomp << "." << std::setfill('0') << std::setw(width) << subdomainId;
        return os.str();
    }

    virtual const std::string get_output_filename() const
    {
        return "subdomain.exo";
    }

    void write_out_files(const std::vector<unsigned>& mappings)
    {
        for(size_t i = 0; i < get_num_procs_target_decomp(); i++)
        {
            if(mappings[i] == static_cast<unsigned>(stk::parallel_machine_rank(get_comm())))
                write_subdomain_file(i);
        }
    }

    void write_subdomain_file(size_t subdomain_num)
    {
        stk::mesh::MetaData newMeta(get_meta().spatial_dimension());
        stk::mesh::BulkData newBulkData(newMeta, MPI_COMM_SELF);
        stk::balance::internal::createNewSubMesh(get_meta(), get_bulk(), *get_meta().get_part(getSubdomainPartName(subdomain_num)), newMeta, newBulkData);
        write_submesh_to_file(newBulkData, subdomain_num);
    }

    void write_submesh_to_file(stk::mesh::BulkData& newBulkData, size_t subdomain_num)
    {
        std::string localFilename = getFilename(get_output_filename(), get_num_procs_target_decomp(), subdomain_num);
        stk::io::write_mesh(localFilename, newBulkData);
    }

    void test_num_elements_per_file()
    {
        if(stk::parallel_machine_rank(get_comm())==0)
        {
            for(size_t i=0;i<get_num_procs_target_decomp();++i)
                test_num_elements_this_subdomain(getFilename(get_output_filename(), get_num_procs_target_decomp(), i));
        }
    }

    void test_num_elements_this_subdomain(const std::string& subdomain_filename)
    {
        stk::mesh::MetaData feta;
        stk::mesh::BulkData bulkData(feta, MPI_COMM_SELF);
        stk::io::fill_mesh(subdomain_filename, bulkData);
        EXPECT_EQ(get_num_elements_per_proc_for_target_decomp(), stk::mesh::count_selected_entities(feta.universal_part(), bulkData.buckets(stk::topology::ELEM_RANK)));
    }

    void clean_up_output_files()
    {
        if(stk::parallel_machine_rank(get_comm()) == 0)
            balance_utils::clearFiles(get_output_filename(), get_num_procs_target_decomp());
    }

    size_t get_num_elements_per_proc_for_target_decomp() const
    {
        return get_x()*get_y()*get_z()/get_num_procs_target_decomp();
    }

    virtual const std::string get_input_mesh_file_name() const
    {
        std::ostringstream os;
        os << "generated:" << get_x() << "x" << get_y() << "x" << get_z() << std::endl;
        return os.str();
    }

    virtual const unsigned get_x() const { return 0; }
    virtual const unsigned get_y() const { return 0; }
    virtual const unsigned get_z() const { return 0; }

    virtual const unsigned get_num_procs_initial_decomp() const = 0;
    virtual const unsigned get_num_procs_target_decomp() const = 0;

    void fill_subdomain_independent_nemesis_info()
    {
        StkMeshInterfaceToNemesis meshInterfaceToNem(get_bulk(), get_comm());
        nemesis_info = meshInterfaceToNem.getNemesisInfo();
    }

    stk::mesh::EntityVector get_nodes_shared_between_subdomains(int this_subdomain_index, int other_subdomain_index)
    {
        stk::mesh::Selector selected_nodes = *get_subdomain_part(this_subdomain_index) & *get_subdomain_part(other_subdomain_index);
        stk::mesh::EntityVector nodes;
        stk::mesh::get_selected_entities(selected_nodes, get_bulk().buckets(stk::topology::NODE_RANK), nodes);
        return nodes;
    }

    void fill_shared_node_proc_info(stk::mesh::EntityVector& shared_nodes, std::vector<int>& procs_for_shared_nodes, int this_subdomain_num, int other_subdomain_num)
    {
        stk::mesh::EntityVector nodes = get_nodes_shared_between_subdomains(this_subdomain_num, other_subdomain_num);
        shared_nodes.insert(shared_nodes.end(), nodes.begin(), nodes.end());
        procs_for_shared_nodes.resize(shared_nodes.size(), other_subdomain_num);
    }

    void fill_shared_node_info_for_this_subdomain(const unsigned this_subdomain_num, stk::mesh::EntityVector& shared_nodes, std::vector<int>& procs_for_shared_nodes)
    {
        for(unsigned other_subdomain_num=0;other_subdomain_num<get_num_procs_target_decomp();++other_subdomain_num)
        {
            if(this_subdomain_num!=other_subdomain_num)
                fill_shared_node_proc_info(shared_nodes, procs_for_shared_nodes, this_subdomain_num, other_subdomain_num);
        }
    }

    int open_exodus_file(const std::string& filename)
    {
        // EX_WRITE + EXLARG
        int cpu_word_size = sizeof(double);
        int io_word_size = sizeof(double);
        float version = 0.0;
        return ex_open(filename.c_str(), EX_WRITE, &cpu_word_size, &io_word_size, &version);
    }

    stk::mesh::EntityVector get_unique_shared_node_list(const stk::mesh::EntityVector& shared_nodes)
    {
        stk::mesh::EntityVector unique_shared_nodes = shared_nodes;
        stk::util::sort_and_unique(unique_shared_nodes);
        return unique_shared_nodes;
    }

    stk::mesh::EntityVector get_internal_nodes(int subdomain_num, stk::mesh::EntityVector& unique_shared_nodes)
    {
        stk::mesh::EntityVector internal_nodes;
        nemesis_info.fill_internal_nodes(get_subdomain_part(subdomain_num), unique_shared_nodes, internal_nodes, get_bulk());
        return internal_nodes;
    }

    void set_nemesis_node_sharing_info(int exoid, stk::mesh::EntityVector &shared_nodes, stk::mesh::EntityVector &unique_shared_nodes,
            stk::mesh::EntityVector &internal_nodes)
    {
        nemesis_info.num_border_nodes_per_subdomain = unique_shared_nodes.size();
        nemesis_info.num_internal_nodes_per_subdomain = internal_nodes.size();
        ThrowRequireWithSierraHelpMsg(nemesis_info.num_border_nodes_per_subdomain+nemesis_info.num_internal_nodes_per_subdomain>0);
        nemesis_info.fill_node_num_map(exoid, nemesis_info.get_number_subdomain_nodes());
        nemesis_info.set_nemesis_node_ids(shared_nodes, internal_nodes, unique_shared_nodes, get_bulk());
    }

    void set_nemesis_info_on_nodes_using_stk_mesh_entities(int exoid, int subdomain_num)
    {
        stk::mesh::EntityVector shared_nodes = get_shared_node_list_from_this_subdomain_to_other_subdomains(subdomain_num);
        stk::mesh::EntityVector unique_shared_nodes = get_unique_shared_node_list(shared_nodes);
        stk::mesh::EntityVector internal_nodes = get_internal_nodes(subdomain_num, unique_shared_nodes);
        set_nemesis_node_sharing_info(exoid, shared_nodes, unique_shared_nodes, internal_nodes);
    }

    void set_nemesis_info_on_elements(int exoid, int subdomain_num)
    {
        nemesis_info.fill_elem_num_map(exoid, stk::mesh::count_selected_entities(*get_subdomain_part(subdomain_num), get_bulk().buckets(stk::topology::ELEM_RANK)));
        nemesis_info.set_nemesis_elem_ids(get_subdomain_part(subdomain_num), get_bulk());
    }

    void fill_nemesis_info_for_subdomain_using_exodus_file_for_maps(int exoid, int subdomain_num)
    {
        set_nemesis_info_on_nodes_using_stk_mesh_entities(exoid, subdomain_num);
        set_nemesis_info_on_elements(exoid, subdomain_num);
    }

    void close_exodus_file(int exoid)
    {
        ex_close(exoid);
    }

    void write_nemesis_info_to_file(int subdomain_num)
    {
        std::string filename = getFilename(get_output_filename(), get_num_procs_target_decomp(), subdomain_num);
        int exoid = open_exodus_file(filename);
        nemesis_info.add_nemesis_info_to_exodus_file(exoid, get_subdomain_part(subdomain_num), subdomain_num, get_num_procs_target_decomp());
        close_exodus_file(exoid);
    }

    void set_subdomain_nemesis_info(int subdomain_num)
    {
        std::string filename = getFilename(get_output_filename(), get_num_procs_target_decomp(), subdomain_num);
        int exoid = open_exodus_file(filename);
        fill_nemesis_info_for_subdomain_using_exodus_file_for_maps(exoid, subdomain_num);
        close_exodus_file(exoid);
    }

    stk::mesh::EntityVector get_shared_node_list_from_this_subdomain_to_other_subdomains(int subdomain_num)
    {
        stk::mesh::EntityVector shared_nodes;
        fill_shared_node_info_for_this_subdomain(subdomain_num, shared_nodes, nemesis_info.procs_for_shared_nodes);
        nemesis_info.set_node_communication_map_counts(shared_nodes.size());
        return shared_nodes;
    }

    void fill_and_write_nemesis_data(int subdomain_num)
    {
        set_subdomain_nemesis_info(subdomain_num);
        write_nemesis_info_to_file(subdomain_num);
        nemesis_info.clear_subdomain_data();
    }

    void add_nemesis_info_to_files(std::vector<unsigned> target_proc_to_starting_proc)
    {
        fill_subdomain_independent_nemesis_info();
        for(unsigned subdomain_num=0;subdomain_num<get_num_procs_target_decomp();++subdomain_num)
            if(target_proc_to_starting_proc[subdomain_num]==static_cast<unsigned>(get_bulk().parallel_rank()))
                fill_and_write_nemesis_data(subdomain_num);
    }

private:
    stk::mesh::EntityProcVec decomp;
    NemesisInfo nemesis_info;
};

#endif
