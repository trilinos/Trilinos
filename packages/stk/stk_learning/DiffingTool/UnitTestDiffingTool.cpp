#include <gtest/gtest.h>

#include <vector>
#include <algorithm>
#include <stdlib.h>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldTraits.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/CreateFaces.hpp>

#include <DiffingToolLib/DiffingTool.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/environment/EnvData.hpp>

#include <stk_util/util/SortAndUnique.hpp>
#include <stk_util/diag/StringUtil.hpp>

#include <stk_util/parallel/DebugTool.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker

#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/getOption.h>

#include "stk_unit_test_utils/unittestMeshUtils.hpp"

using namespace stk::diff;

namespace 
{

template <typename T>
void compare_vectors(std::vector<T> items1, std::vector<T>& items2, int line_number)
{
    ASSERT_EQ(items1.size(), items2.size()) << "Vector sizes diff. Called from line: " << line_number << std::endl;;
    for(size_t i=0;i<items1.size();++i)
    {
        EXPECT_EQ(items1[i], items2[i]) << "Vectors diff at index = " << i << ", called from line: " << line_number << std::endl;
    }
}

std::vector<std::string> filter_out_empty_strings(const std::vector<std::string>& strings)
{
    std::vector<std::string> filtered_strings = strings;
    filtered_strings.erase(std::remove(filtered_strings.begin(), filtered_strings.end(), ""), filtered_strings.end());
    return filtered_strings;
}

std::string join_string_vector_ignoring_empty_strings(const std::vector<std::string>& strings, char sep)
{
    std::vector<std::string> filtered_strings = filter_out_empty_strings(strings);
    return stk::util::join(filtered_strings, sep);
}

std::vector<std::string> split_strings_from_stream(std::stringstream& ss, std::string& scratch_space, char delim)
{
    std::vector<std::string> split_strings;
    while (std::getline(ss, scratch_space, delim))
    {
        if(scratch_space != "")
        {
            split_strings.push_back(scratch_space);
        }
    }
    return split_strings;
}

std::vector<std::string> split(const std::string &s, char delim)
{
    std::stringstream ss(s);
    std::string scratch_space;
    return split_strings_from_stream(ss, scratch_space, delim);
}

void find_and_replace_within_string_vector(std::vector<std::string>& split_strings, const std::string& string_to_find, const std::string& new_string)
{
    size_t size_of_string_to_find = string_to_find.size();
    for(size_t i=0;i<split_strings.size();++i)
    {
        if(split_strings[i].substr(0, size_of_string_to_find) == string_to_find)
        {
            split_strings[i] = new_string;
        }
    }
}

void replace_all_custom_ghostings_with_fmwk_recv_ghost_part(std::vector<std::string> &part_names)
{
    const std::string string_to_find = "{custom_ghosting_";
    const std::string replace_found_string_with_this = "{FMWK_RECV_GHOST_PART}";
    find_and_replace_within_string_vector(part_names, string_to_find, replace_found_string_with_this);
    stk::util::sort_and_unique(part_names);
}

void replace_fmwk_stk_shared_part_with_shares(std::vector<std::string> &part_names)
{
    const std::string string_to_find = "{FMWK_STK_SHARED_PART}";
    const std::string replace_found_string_with_this = "{SHARES}";
    find_and_replace_within_string_vector(part_names, string_to_find, replace_found_string_with_this);
}

void replace_fmwk_send_ghost_part_with_send_ghost(std::vector<std::string> &part_names)
{
    const std::string string_to_find = "{FMWK_SEND_GHOST_PART}";
    const std::string replace_found_string_with_this = "{SEND_GHOST}";
    find_and_replace_within_string_vector(part_names, string_to_find, replace_found_string_with_this);
}

void replace_fmwk_locally_owned_context_bit_with_owns(std::vector<std::string> &part_names)
{
    const std::string string_to_find = "LOCALLY_OWNED_CONTEXT_BIT";
    const std::string replace_found_string_with_this = "{OWNS}";
    find_and_replace_within_string_vector(part_names, string_to_find, replace_found_string_with_this);
}

void yank_out_part(std::vector<std::string> &part_names, const std::string& part_to_yank_out)
{
    find_and_replace_within_string_vector(part_names, part_to_yank_out, "");
    part_names.erase(std::remove(part_names.begin(), part_names.end(), ""), part_names.end());
}

void adjust_stk_part_names(std::vector<std::string>& stk_part_names)
{
    replace_all_custom_ghostings_with_fmwk_recv_ghost_part(stk_part_names);
    yank_out_part(stk_part_names, "LOCALLY_OWNED_CONTEXT_BIT");
    yank_out_part(stk_part_names, "sides_created_during_death");
}

void adjust_fmwk_part_names(std::vector<std::string>& fmwk_part_names)
{
    yank_out_part(fmwk_part_names, "{SHARES}");
    yank_out_part(fmwk_part_names, "{OWNS}");
    replace_fmwk_send_ghost_part_with_send_ghost(fmwk_part_names);
    replace_fmwk_locally_owned_context_bit_with_owns(fmwk_part_names);
    replace_fmwk_stk_shared_part_with_shares(fmwk_part_names);
    stk::util::sort_and_unique(fmwk_part_names);
}

void adjust_stk_part_names_per_bucket(std::vector<std::string>& stk_part_names)
{
    for(size_t i=0;i<stk_part_names.size();++i)
    {
        std::vector<std::string> part_names_for_bucket = split(stk_part_names[i], ' ');
        adjust_stk_part_names(part_names_for_bucket);
        stk_part_names[i] = join_string_vector_ignoring_empty_strings(part_names_for_bucket, ' ');
    }
}

void adjust_fmwk_part_names_per_bucket(std::vector<std::string>& fmwk_part_names)
{
    for(size_t i=0;i<fmwk_part_names.size();++i)
    {
        std::vector<std::string> part_names_for_bucket = split(fmwk_part_names[i], ' ');
        adjust_fmwk_part_names(part_names_for_bucket);
        fmwk_part_names[i] = join_string_vector_ignoring_empty_strings(part_names_for_bucket, ' ');
    }
}

TEST(ParallelDebugTool, testStringManips)
{
    {
        std::string part_name_before = "a b c {custom_ghosting_10} d e";
        std::string gold_result = "a b c d e {FMWK_RECV_GHOST_PART}";

        std::vector<std::string> part_names = split(part_name_before, ' ');
        adjust_stk_part_names(part_names);
        std::vector<std::string> gold_part_names = split(gold_result, ' ');
        compare_vectors(part_names, gold_part_names, __LINE__);
    }

    {
        std::string part_name_before = "a b c {custom_ghosting_2} {custom_ghosting_10} d e";
        std::string gold_result = "a b c d e {FMWK_RECV_GHOST_PART}";

        std::vector<std::string> part_names = split(part_name_before, ' ');
        adjust_stk_part_names(part_names);
        std::vector<std::string> gold_part_names = split(gold_result, ' ');
        compare_vectors(part_names, gold_part_names, __LINE__);
    }

    {
        std::string part_name_before = "a b c {SHARES} d e";
        std::string gold_result = "a b c d e";

        std::vector<std::string> part_names = split(part_name_before, ' ');
        adjust_fmwk_part_names(part_names);
        std::vector<std::string> gold_part_names = split(gold_result, ' ');
        compare_vectors(part_names, gold_part_names, __LINE__);
    }

    {
        std::string part_name_before = "a b c {FMWK_SEND_GHOST_PART} d e";
        std::string gold_result = "a b c d e {SEND_GHOST}";

        std::vector<std::string> part_names = split(part_name_before, ' ');
        adjust_fmwk_part_names(part_names);
        std::vector<std::string> gold_part_names = split(gold_result, ' ');
        compare_vectors(part_names, gold_part_names, __LINE__);
    }

    {
        std::string part_name_before = "a b c {FMWK_STK_SHARED_PART} d e";
        std::string gold_result = "a b c d e {SHARES}";

        std::vector<std::string> part_names = split(part_name_before, ' ');
        adjust_fmwk_part_names(part_names);
        std::vector<std::string> gold_part_names = split(gold_result, ' ');
        compare_vectors(part_names, gold_part_names, __LINE__);
    }

    {
        std::string stk_part_namesA  = "a b c {custom_ghosting_2} {custom_ghosting_3} {SEND_GHOST} {SHARES}";
        std::string fmwk_part_namesB = "a b c {SHARES} {FMWK_RECV_GHOST_PART} {FMWK_SEND_GHOST_PART} {FMWK_STK_SHARED_PART}";

        std::vector<std::string> stk_part_names = split(stk_part_namesA, ' ');
        adjust_stk_part_names(stk_part_names);
        std::vector<std::string> fmwk_part_names = split(fmwk_part_namesB, ' ');
        adjust_fmwk_part_names(fmwk_part_names);

        compare_vectors(stk_part_names, fmwk_part_names, __LINE__);
    }
}

void setup_env_data(stk::EnvData& env_data, stk::ParallelMachine global_comm, int color)
{
    env_data.m_parallelComm = splitComm(color, global_comm);
    env_data.m_worldComm = global_comm;
}

void communicate_part_differences_to_apps(stk::ParallelMachine comm, int numLocalDiffs)
{
    get_global_part_differences(comm, numLocalDiffs);
}

void wait_for_apps_to_pack_and_send_data(stk::CommSparse& comm)
{
    comm.allocate_buffers();
    comm.communicate();
}

std::vector<std::string> unpack_buffer_to_string_vector(stk::CommBuffer& buff)
{
    std::vector<std::string> partNames;
    while(buff.remaining())
    {
        partNames.push_back(unpack_string(buff));
    }
    return partNames;
}

std::vector<std::string> get_sorted_part_names(stk::CommBuffer& buff)
{
    std::vector<std::string> partNames = unpack_buffer_to_string_vector(buff);
    std::sort(partNames.begin(), partNames.end());
    return partNames;
}

void create_output_message(std::ostringstream &os, std::vector<std::string>& part_diffs)
{
    if(!part_diffs.empty())
    {
        os << join_string_vector_ignoring_empty_strings(part_diffs, ' ') << std::endl;
    }
    else
    {
        os << std::endl;
    }
}

std::vector<std::string> part_name_set_difference(const std::vector<std::string> & partNames1, const std::vector<std::string> & partNames2)
{
    std::vector<std::string> part_diffs(partNames1.size());
    auto iter = std::set_difference(partNames1.begin(), partNames1.end(), partNames2.begin(), partNames2.end(), part_diffs.begin());
    part_diffs.resize(iter - part_diffs.begin());
    return part_diffs;
}

int part_name_set_difference_message(const std::vector<std::string> & partNames1, const std::vector<std::string> & partNames2, std::ostringstream& os)
{
    std::vector<std::string> part_diffs = part_name_set_difference(partNames1, partNames2);
    create_output_message(os, part_diffs);
    int num_diffs = part_diffs.size();
    return num_diffs;
}

int get_parts_in_first_not_in_second(const std::vector<std::string> & partNames1, const std::vector<std::string> & partNames2, std::ostringstream & globalss)
{
    std::ostringstream os("Parts in first app, not in second app: ", std::ios_base::ate);
    int num_diffs = part_name_set_difference_message(partNames1, partNames2, os);
    if (num_diffs > 0)
    {
        globalss << os.str();
    }
    return num_diffs;
}

int get_parts_in_second_not_in_first(const std::vector<std::string> & partNames1, const std::vector<std::string> & partNames2, std::ostringstream & globalss)
{
    std::ostringstream os("Parts in second app, not in first app: ", std::ios_base::ate);
    int num_diffs = part_name_set_difference_message(partNames2, partNames1, os);
    if (num_diffs > 0)
    {
        globalss << os.str();
    }
    return num_diffs;
}

int compare_part_names(const std::vector<std::string> & partNames1, const std::vector<std::string> & partNames2, std::ostringstream & globalss)
{
    return get_parts_in_first_not_in_second(partNames1, partNames2, globalss) +
           get_parts_in_second_not_in_first(partNames1, partNames2, globalss);
}

void unpack_part_names(stk::EnvData & env_data, std::vector<std::string> & partNames1, std::vector<std::string> & partNames2)
{
    stk::CommSparse comm(env_data.m_worldComm);
    wait_for_apps_to_pack_and_send_data(comm);
    std::pair<int,int> sourceProcs = getSourceProcs(env_data.m_worldComm);
    partNames1 = get_sorted_part_names(comm.recv_buffer(sourceProcs.first));
    partNames2 = get_sorted_part_names(comm.recv_buffer(sourceProcs.second));
}

void output_message(const std::string &message, int localProcId)
{
    if (!message.empty())
    {
        std::ostringstream massaged_output;
        std::vector<std::string> tokens = split(message,'\n');
        for(std::string& s: tokens)
        {
            massaged_output << "[" << localProcId << "] " << s << std::endl;
        }
        std::cerr << massaged_output.str();
        {
            std::ostringstream oss;
            oss << "diffingOutput." << localProcId;
            std::ofstream out(oss.str(), std::ios_base::app);
            out << massaged_output.str();
        }
    }
}

enum DiffingOption { STK_VS_STK, FMWK_VS_STK, FMWK_DEATH_VS_STK_DEATH};

void compare_parts(stk::EnvData & env_data,  DiffingOption diffing_option)
{
    std::ostringstream globalss;

    std::vector<std::string> stkPartNames, fmwkPartNames;
    unpack_part_names(env_data, stkPartNames, fmwkPartNames);
    if(diffing_option == FMWK_VS_STK)
    {
        adjust_stk_part_names(stkPartNames);
        adjust_fmwk_part_names(fmwkPartNames);
    }
    else if(diffing_option == FMWK_DEATH_VS_STK_DEATH)
    {
        adjust_stk_part_names(stkPartNames);
        adjust_stk_part_names(fmwkPartNames);
    }

    int num_diffs = compare_part_names(stkPartNames, fmwkPartNames, globalss);
    communicate_part_differences_to_apps(env_data.m_worldComm, num_diffs);
    if (num_diffs > 0)
    {
        output_message(globalss.str(), stk::parallel_machine_rank(env_data.m_parallelComm));
    }
}

std::pair<int,int> getXYdim()
{
    int xdim = 2;
    int ydim = 1;
    return std::make_pair(xdim, ydim);
}

std::string get_generated_mesh_spec(int numProcs)
{
    std::ostringstream os;
    std::pair<int,int> xydim = getXYdim();
    os << "generated:" << xydim.first << "x" << xydim.second << "x" << numProcs;
    return os.str();
}

class TestMesh
{
public:
    TestMesh(stk::ParallelMachine comm)
    :m_meta_data(), m_bulk_data(m_meta_data, comm)
    {
        stk::unit_test_util::fill_mesh_using_stk_io(get_generated_mesh_spec(m_bulk_data.parallel_size()), m_bulk_data);
        stk::mesh::create_edges(m_bulk_data);
        stk::mesh::create_faces(m_bulk_data);
    }

    stk::mesh::BulkData & get_bulk_data()
    {
        return m_bulk_data;
    }

private:
    stk::mesh::MetaData m_meta_data;
    stk::mesh::BulkData m_bulk_data;
};

void expect_mpmd(stk::EnvData & env_data)
{
    int num_global_procs = stk::parallel_machine_size(env_data.m_worldComm);
    int num_local_procs = stk::parallel_machine_size(env_data.m_parallelComm);
    EXPECT_TRUE(num_global_procs%3 == 0 && num_local_procs==num_global_procs/3);
}

void expect_num_elements(const stk::mesh::BulkData & bulk)
{
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(bulk, counts);
    unsigned zdim = bulk.parallel_size();
    std::pair<int,int> xydim = getXYdim();
    size_t expected_num_elems = xydim.first*xydim.second*zdim;
    EXPECT_EQ(expected_num_elems, counts[stk::topology::ELEM_RANK]);
}

std::vector<size_t> unpack_bucket_counts(stk::CommBuffer& buf)
{
    std::vector<size_t> bucket_counts;
    while(buf.remaining())
    {
        size_t count = 0;
        buf.unpack<size_t>(count);
        bucket_counts.push_back(count);
    }
    return bucket_counts;
}

int compare_bucket_counts_message(int procId, const std::vector<size_t>& bucket_counts1, const std::vector<size_t>& bucket_counts2, std::ostringstream& globalss)
{
    std::ostringstream localss;
    int num_diffs = 0;
    if (bucket_counts1.size() != bucket_counts2.size())
    {
        num_diffs = 1;
        localss << "Number of entity ranks doesn't match: " << bucket_counts1.size() << " != " << bucket_counts2.size() << std::endl;
    }
    else
    {
        for(size_t i=0; i<bucket_counts1.size(); ++i)
        {
            if (bucket_counts1[i] != bucket_counts2[i])
            {
                num_diffs++;
                localss << "Number of buckets for entity rank "<< i << " doesn't match: " << bucket_counts1[i] << " != " << bucket_counts2[i] << std::endl;
            }
        }
    }
    if (num_diffs > 0)
    {
        globalss << localss.str();
    }
    return num_diffs;
}

bool compare_bucket_part_names_message(size_t num_entity_ranks1, size_t num_entity_ranks2, stk::EnvData& env_data, std::ostringstream& globalss, DiffingOption diffing_option)
{
    size_t num_equal_counts = std::min(num_entity_ranks1, num_entity_ranks2);
    std::pair<int,int> sourceProcs = getSourceProcs(env_data.m_worldComm);
    bool diffs_occurred = false;

    for(size_t i=0;i<num_equal_counts;++i)
    {
        stk::CommSparse comm(env_data.m_worldComm);
        wait_for_apps_to_pack_and_send_data(comm);

        std::vector<std::string> stk_part_names_for_buckets_of_rank = unpack_buffer_to_string_vector(comm.recv_buffer(sourceProcs.first));
        std::vector<std::string> fmwk_part_names_for_buckets_of_rank = unpack_buffer_to_string_vector(comm.recv_buffer(sourceProcs.second));

        if(diffing_option == FMWK_VS_STK)
        {
            adjust_stk_part_names_per_bucket(stk_part_names_for_buckets_of_rank);
            adjust_fmwk_part_names_per_bucket(fmwk_part_names_for_buckets_of_rank);
        }
        else if(diffing_option == FMWK_DEATH_VS_STK_DEATH)
        {
            adjust_stk_part_names_per_bucket(stk_part_names_for_buckets_of_rank);
            adjust_stk_part_names_per_bucket(fmwk_part_names_for_buckets_of_rank);
        }

        size_t min_bucket_compare = std::min(stk_part_names_for_buckets_of_rank.size(), fmwk_part_names_for_buckets_of_rank.size());
        for(size_t j=0;j<min_bucket_compare;++j)
        {
            std::vector<std::string> names1 = split(stk_part_names_for_buckets_of_rank[j], ' ');
            std::vector<std::string> names2 = split(fmwk_part_names_for_buckets_of_rank[j], ' ');

            {
                std::ostringstream os;
                os << "Buckets in first app, not in second app, with parts for rank " << i << " for bucket index " << j << "\n";
                int num_local_diffs = part_name_set_difference_message(names1, names2, os);
                if (num_local_diffs > 0)
                {
                    globalss << os.str();
                    diffs_occurred = true;
                }
            }

            {
                std::ostringstream os;
                os << "Buckets in second app, not in first app, with parts for rank " << i << " for bucket index " << j << " \n";
                int num_local_diffs = part_name_set_difference_message(names2, names1, os);
                if (num_local_diffs > 0)
                {
                    globalss << os.str();
                    diffs_occurred = true;
                }
            }
        }
    }

    return diffs_occurred;
}

int compare_bucket_counts(stk::EnvData& env_data, size_t& num_ranks1, size_t& num_ranks2, std::ostringstream& globalss)
{
    stk::CommSparse comm(env_data.m_worldComm);
    wait_for_apps_to_pack_and_send_data(comm);
    std::pair<int,int> sourceProcs = getSourceProcs(env_data.m_worldComm);
    std::vector<size_t> bucket_counts1 = unpack_bucket_counts(comm.recv_buffer(sourceProcs.first));
    std::vector<size_t> bucket_counts2 = unpack_bucket_counts(comm.recv_buffer(sourceProcs.second));
    num_ranks1 = bucket_counts1.size();
    num_ranks2 = bucket_counts2.size();
    int localProcId = stk::parallel_machine_rank(env_data.m_parallelComm);
    int num_diffs = compare_bucket_counts_message(localProcId, bucket_counts1, bucket_counts2, globalss);
    return get_global_bucket_count_differences(env_data.m_worldComm, num_diffs);
}

void build_message_for_missing_ranks(stk::EnvData& env_data, size_t num_ranks1, size_t num_ranks2, std::ostringstream& globalss)
{
    std::pair<int,int> sourceProcs = getSourceProcs(env_data.m_worldComm);
    int which_app = -1;
    size_t max = 0;
    int which_proc = -1;
    size_t num_equal_counts = std::min(num_ranks1, num_ranks2);
    if(num_ranks1 != num_equal_counts)
    {
        which_app = 0;
        max = num_ranks1;
        which_proc = sourceProcs.first;
    }
    else if (num_ranks2 != num_equal_counts)
    {
        which_app = 1;
        max = num_ranks2;
        which_proc = sourceProcs.second;
    }

    if(which_app != -1)
    {
        for(size_t i=num_equal_counts;i<max;++i)
        {
            stk::CommSparse comm(env_data.m_worldComm);
            wait_for_apps_to_pack_and_send_data(comm);

            stk::CommBuffer & buf = comm.recv_buffer(which_proc);
            std::vector<std::string> part_names_for_buckets_of_rank = get_sorted_part_names(buf);

            globalss << "Proc " << which_proc << " of app # " << which_app << " has extra entity ranks\n";
        }
    }
}

int bool_to_int(bool flag)
{
    return flag ? 1 : 0;
}

void compare_bucket_part_memberships(stk::EnvData& env_data, DiffingOption diffing_option)
{
    std::ostringstream globalss;
    size_t num_ranks1 = 0, num_ranks2 = 0;
    int num_count_issues = compare_bucket_counts(env_data, num_ranks1, num_ranks2, globalss);

    if(num_count_issues>0)
    {
        std::pair<int,int> sourceProcs = getSourceProcs(env_data.m_worldComm);
        for(size_t i=0;i<num_ranks1;++i)
        {
            stk::CommSparse comm(env_data.m_worldComm);
            wait_for_apps_to_pack_and_send_data(comm);
            std::vector<std::string> stk_part_names_for_entities = unpack_buffer_to_string_vector(comm.recv_buffer(sourceProcs.first));
            std::vector<std::string> fmwk_part_names_for_entities = unpack_buffer_to_string_vector(comm.recv_buffer(sourceProcs.second));

            if(diffing_option == FMWK_VS_STK)
            {
                adjust_stk_part_names_per_bucket(stk_part_names_for_entities);
                adjust_fmwk_part_names_per_bucket(fmwk_part_names_for_entities);
            }
            else if(diffing_option == FMWK_DEATH_VS_STK_DEATH)
            {
                adjust_stk_part_names_per_bucket(stk_part_names_for_entities);
                adjust_stk_part_names_per_bucket(fmwk_part_names_for_entities);
            }

            if(stk_part_names_for_entities.size() != fmwk_part_names_for_entities.size())
            {
                globalss << "Num entities of rank " << i << " does not match! First app: " << stk_part_names_for_entities.size()
                        << " and second app: " << fmwk_part_names_for_entities.size() << std::endl;
            }
            else
            {
                size_t num_entities = stk_part_names_for_entities.size();
                for(size_t j=0;j<num_entities;++j)
                {
                    std::vector<std::string> names1 = split(stk_part_names_for_entities[j], ' ');
                    std::vector<std::string> names2 = split(fmwk_part_names_for_entities[j], ' ');

                    {
                        std::ostringstream os;
                        os << "Parts on entity in first app, not in second app, for rank " << i << " for entity index " << j << "\n";
                        os << "App 1: " << stk_part_names_for_entities[j] << std::endl;
                        os << "App 2: " << fmwk_part_names_for_entities[j] << std::endl;
                        int num_local_diffs = part_name_set_difference_message(names1, names2, os);
                        if (num_local_diffs > 0)
                        {
                            globalss << os.str();
                        }
                    }

                    {
                        std::ostringstream os;
                        os << "Parts on entity in second app, not in first app, for rank " << i << " for entity index " << j << "\n";
                        os << "App 1: " << stk_part_names_for_entities[j] << std::endl;
                        os << "App 2: " << fmwk_part_names_for_entities[j] << std::endl;
                        int num_local_diffs = part_name_set_difference_message(names2, names1, os);
                        if (num_local_diffs > 0)
                        {
                            globalss << os.str();
                        }
                    }
                }
            }
        }
    }

    bool diffs_exist = compare_bucket_part_names_message(num_ranks1, num_ranks2, env_data, globalss, diffing_option);
    if (diffs_exist)
    {
        build_message_for_missing_ranks(env_data, num_ranks1, num_ranks2, globalss);
    }
    int diffs = get_global_bucket_part_membership_differences(env_data.m_worldComm, bool_to_int(diffs_exist));
    if (diffs > 0)
    {
        output_message(globalss.str(), stk::parallel_machine_rank(env_data.m_parallelComm));
    }
}

stk::mesh::Part& put_elem1_in_new_part_for_app1(stk::mesh::BulkData& bulk)
{
    stk::mesh::Part& newPart = bulk.mesh_meta_data().declare_part("new_part");
    stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, 1);
    bulk.modification_begin();
    if (bulk.is_valid(elem1) && bulk.bucket(elem1).owned())
    {
        bulk.change_entity_parts(elem1, {&newPart}, {});
    }
    bulk.modification_end();
    return newPart;
}

enum DifferentPartOption {
    KeepPartsTheSame,
    AddDifferentPartDoNotIgnore,
    AddDifferentPartIgnore
};

void call_parts_match(DifferentPartOption differentPartOption, stk::mesh::BulkData& bulk, stk::EnvData& env_data, stk::mesh::Part* newPart)
{
    switch(differentPartOption) {
        case AddDifferentPartDoNotIgnore:
            EXPECT_FALSE(stk::diff::parts_match(bulk, env_data));
            break;
        case AddDifferentPartIgnore:
            EXPECT_TRUE(stk::diff::parts_match_except(bulk, env_data, newPart));
            break;
        default:
            EXPECT_TRUE(stk::diff::parts_match(bulk, env_data));
            break;
    }
}

void run_mock_app(int color, DifferentPartOption differentPartOption)
{
    stk::EnvData & env_data = stk::EnvData::instance();
    setup_env_data(env_data, MPI_COMM_WORLD, color);
    expect_mpmd(env_data);

    TestMesh mesh(env_data.m_parallelComm);
    stk::mesh::BulkData & bulk = mesh.get_bulk_data();
    expect_num_elements(bulk);

    stk::mesh::Part* newPart = nullptr;

    if (differentPartOption != KeepPartsTheSame && color == 0)
    {
        newPart = &put_elem1_in_new_part_for_app1(bulk);
    }

    bool continue_runs = true;
    size_t num_mock_steps = 4;
    for(size_t i=0;i<num_mock_steps;++i)
    {
        communicate_run_state(env_data, continue_runs);
        call_parts_match(differentPartOption, bulk, env_data, newPart);

        if (differentPartOption == KeepPartsTheSame)
        {
            EXPECT_TRUE(bucket_part_memberships_match(bulk, env_data));
        }
        else
        {
            EXPECT_TRUE(!bucket_part_memberships_match(bulk, env_data));
        }
    }

    continue_runs = false;
    communicate_run_state(env_data, continue_runs);
}

TEST(ParallelDebugTool, mock_app_1_with_different_part)
{
    int color = 0;
    DifferentPartOption partOption = AddDifferentPartDoNotIgnore;
    run_mock_app(color, partOption);
}

TEST(ParallelDebugTool, mock_app_2_with_different_part)
{
    int color = 1;
    DifferentPartOption partOption = AddDifferentPartDoNotIgnore;
    run_mock_app(color, partOption);
}

TEST(ParallelDebugTool, mock_app_1_ignore_different_part)
{
    int color = 0;
    DifferentPartOption partOption = AddDifferentPartIgnore;
    run_mock_app(color, partOption);
}

TEST(ParallelDebugTool, mock_app_2_ignore_different_part)
{
    int color = 1;
    DifferentPartOption partOption = AddDifferentPartIgnore;
    run_mock_app(color, partOption);
}

TEST(ParallelDebugTool, mock_app_1)
{
    int color = 0;
    DifferentPartOption partOption = KeepPartsTheSame;
    run_mock_app(color, partOption);
}

TEST(ParallelDebugTool, mock_app_2)
{
    int color = 1;
    DifferentPartOption partOption = KeepPartsTheSame;
    run_mock_app(color, partOption);
}

bool check_if_continuing(stk::EnvData &env_data)
{
    std::pair<int,int> sourceProcs = getSourceProcs(env_data.m_worldComm);
    stk::CommSparse comm(env_data.m_worldComm);
    wait_for_apps_to_pack_and_send_data(comm);

    int receive_int1 = 0;
    comm.recv_buffer(sourceProcs.first).unpack<int>(receive_int1);
    int global_rec_int1 = parallel_sum(env_data.m_parallelComm, receive_int1);

    if(receive_int1 !=0)
    {
        ThrowRequireMsg(global_rec_int1 != 0, "Receiving mixed messages in diffing tool.\n");
    }
    else
    {
        ThrowRequireMsg(global_rec_int1 == 0, "Receiving mixed messages in diffing tool.\n");
    }

    int receive_int2 = 0;
    comm.recv_buffer(sourceProcs.second).unpack<int>(receive_int2);
    int global_rec_int2 = parallel_sum(env_data.m_parallelComm, receive_int2);

    if(receive_int2 != 0)
    {
        ThrowRequireMsg(global_rec_int2 != 0, "Receiving mixed messages in diffing tool.\n");
    }
    else
    {
        ThrowRequireMsg(global_rec_int2 == 0, "Receiving mixed messages in diffing tool.\n");
    }

    ThrowRequireMsg(global_rec_int1 == global_rec_int2, "Receiving mixed messages in diffing tool.\n");

    bool continue_runs = true;
    if(global_rec_int1 == 0)
    {
        continue_runs = false;
    }
    return continue_runs;
}

void run_diffing_tool(stk::EnvData &env_data)
{
    bool continue_work = true;
    while(continue_work)
    {
        continue_work = check_if_continuing(env_data);
        if(continue_work)
        {
            compare_parts(env_data, DiffingOption::STK_VS_STK);
            compare_bucket_part_memberships(env_data, DiffingOption::STK_VS_STK);
        }
    }
}

TEST(ParallelDebugTool, mockDiffingToolUnit)
{
    int color = 2;
    stk::EnvData &env_data = stk::EnvData::instance();
    setup_env_data(env_data, MPI_COMM_WORLD, color);
    expect_mpmd(env_data);
    run_diffing_tool(env_data);
}

TEST(ParallelDebugTool, mockDiffingToolApp)
{
    int color = 2;
    stk::EnvData &env_data = stk::EnvData::instance();
    setup_env_data(env_data, MPI_COMM_WORLD, color);
    expect_mpmd(env_data);

    bool continue_work = true;
    while(continue_work)
    {
        continue_work = check_if_continuing(env_data);
        if(continue_work)
        {
            compare_parts(env_data, DiffingOption::FMWK_VS_STK);
            compare_bucket_part_memberships(env_data, DiffingOption::FMWK_VS_STK);
        }
    }
}

TEST(ParallelDebugTool, mockDiffingToolAppComparingDeathAlgorithms)
{
    int color = 2;
    stk::EnvData &env_data = stk::EnvData::instance();
    setup_env_data(env_data, MPI_COMM_WORLD, color);
    expect_mpmd(env_data);

    bool continue_work = true;
    while(continue_work)
    {
        continue_work = check_if_continuing(env_data);
        if(continue_work)
        {
            compare_parts(env_data, DiffingOption::FMWK_DEATH_VS_STK_DEATH);
            compare_bucket_part_memberships(env_data, DiffingOption::FMWK_DEATH_VS_STK_DEATH);
        }

    }
}
}
