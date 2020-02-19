// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef MESH_FIXTURE_MXN
#define MESH_FIXTURE_MXN

#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>

#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>

#include <stk_unit_test_utils/MeshFixture.hpp>
#include <test_utils/StkMeshInterfaceToNemesis.hpp>

#include <test_utils/BalanceTestUtilities.hpp>
#include <stk_mesh/base/Comm.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_balance/internal/balanceMtoN.hpp>
#include <stk_balance/internal/MxNutils.hpp>
#include <stk_balance/internal/entityDataToField.hpp>

class MeshFixtureMxNRebalance : public stk::unit_test_util::MeshFixture
{
protected:
    MeshFixtureMxNRebalance() : MeshFixture() {}

    void setup_and_test_balance_of_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        write_rebalanced_mxn();
        test_decomp();
    }

    void write_rebalanced_mxn()
    {
        setup_initial_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        EXPECT_TRUE(stk::balance::internal::rebalanceMtoN(get_bulk(), *targetDecompField, get_num_procs_target_decomp(), get_output_filename()));
    }

    // FIXME
    void setup_initial_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        create_target_decomp_field_on_entire_mesh();
        setup_mesh(get_input_mesh_file_name(), auraOption);
    }

    void create_field_on_entire_mesh(const std::string& fieldName)
    {
        targetDecompField = &get_meta().declare_field<stk::mesh::Field<double> >(stk::topology::ELEMENT_RANK, fieldName, 1);
        stk::mesh::put_field_on_mesh(*targetDecompField, get_meta().universal_part(), static_cast<double*>(nullptr));
    }

    void create_target_decomp_field_on_entire_mesh()
    {
        create_field_on_entire_mesh(get_target_decomp_field_name());
    }

    const std::string get_target_decomp_field_name()
    {
        return "Target_Decomp";
    }

    void test_decomp()
    {
        // Make sure all procs have written their files before p0 tries to read them
        MPI_Barrier(get_comm());

        test_subdomain_files();
        test_decomp_is_balanced();
    }

    void test_subdomain_files()
    {
        test_num_elements_per_file();
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
        stk::balance::GraphCreationSettings graphSettings;
        stk::mesh::EntityProcVec localDecomp = stk::balance::internal::get_element_decomp(get_num_procs_target_decomp(), get_bulk(), graphSettings);

        std::vector<unsigned> num_elements_assigned_to_each_proc_local(get_num_procs_target_decomp(),0);
        for(size_t i=0;i<localDecomp.size();++i)
            num_elements_assigned_to_each_proc_local[localDecomp[i].second]++;
        return num_elements_assigned_to_each_proc_local;
    }

    // FIXME
    std::string getSubdomainPartName(int subdomainId)
    {
        std::ostringstream out;
        out << "subdomain_" << subdomainId;
        return out.str();
    }

    // FIXME
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

    virtual std::string get_output_filename() const
    {
        return "subdomain.exo";
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

    size_t get_num_elements_per_proc_for_target_decomp() const
    {
        return get_x()*get_y()*get_z()/get_num_procs_target_decomp();
    }

    virtual std::string get_input_mesh_file_name() const
    {
        std::ostringstream os;
        os << "generated:" << get_x() << "x" << get_y() << "x" << get_z() << std::endl;
        return os.str();
    }

    virtual unsigned get_x() const { return 0; }
    virtual unsigned get_y() const { return 0; }
    virtual unsigned get_z() const { return 0; }
    virtual unsigned get_num_procs_initial_decomp() const = 0;
    virtual unsigned get_num_procs_target_decomp() const = 0;

    stk::mesh::Field<double> *targetDecompField;
};

#endif
