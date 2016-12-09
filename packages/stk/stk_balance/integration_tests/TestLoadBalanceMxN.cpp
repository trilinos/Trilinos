#include <test_utils/OptionsForTesting.hpp>
#include <test_utils/MeshFixtureMxNRebalance.hpp>

namespace
{

class TestBalanceMxNRebalanceUsingInputFiles : public MeshFixtureMxNRebalance
{
protected:
    TestBalanceMxNRebalanceUsingInputFiles()
    : MeshFixtureMxNRebalance(), num_procs_initial_decomp(0),
    num_procs_target_decomp(0), mOutputFileName("subomain.exo"), mInputFileName("")
    {
    }

    virtual const unsigned get_num_procs_initial_decomp() const { return num_procs_initial_decomp; }
    virtual const unsigned get_num_procs_target_decomp()  const { return num_procs_target_decomp; }
    virtual const std::string get_output_filename() const { return mOutputFileName; }
    virtual const std::string get_input_mesh_file_name() const { return mInputFileName; }

    void set_options()
    {
        Options options = getOptionsForTest("none");
        test_options_are_set(options);
        init_member_data(options);
    }

    void init_member_data(const Options &options)
    {
        num_procs_initial_decomp = stk::parallel_machine_size(get_comm());
        num_procs_target_decomp = options.getNumTargetProcs();
        mInputFileName = options.getMeshFileName();
        mOutputFileName = options.getOutputFilename();
    }

private:
    void test_options_are_set(const Options &options)
    {
        ASSERT_TRUE(options.getMeshFileName() != "none");
        ASSERT_TRUE(options.getNumTargetProcs() != 0);
    }

private:
    int num_procs_initial_decomp;
    int num_procs_target_decomp;
    std::string mOutputFileName;
    std::string mInputFileName;
};


TEST_F(TestBalanceMxNRebalanceUsingInputFiles, MxN_decompositionWithoutAura)
{
    set_options();
    write_rebalanced_mxn();
}

}

