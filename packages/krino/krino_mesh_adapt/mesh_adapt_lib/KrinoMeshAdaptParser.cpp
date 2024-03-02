// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <KrinoMeshAdaptInputData.hpp>
#include <KrinoMeshAdaptParser.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_util/environment/Env.hpp>
#include <exception>                         // for exception
#include <string>                            // for string
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace krino
{

static std::string make_flag(const stk::CommandLineOption &option)
{
    return option.name+","+option.abbreviation;
}

bool MeshAdaptParser::read_command_line( int argc, char *argv[], stk::ParallelMachine comm)
{
    MeshAdaptInputData defaultValues;

    const stk::CommandLineOption inmeshOption{"input_mesh", "i", "Filename of input genesis mesh."};
    const stk::CommandLineOption outmeshOption{"output_mesh", "o", "Filename of output genesis refined mesh to write."};

    const stk::CommandLineOption numUniformRefinementLevelsOption{"--number_refines", "n", "Number of uniform refinement passes.  Default is 1."};
    const stk::CommandLineOption force64BitOption{"force_64_bit", "", "Force all new nodes and elements past 32bit limit."};
    const stk::CommandLineOption assert32BitOption{"assert_32_bit", "", "Assert that all nodes and elements are within the 32bit limit."};
    const stk::CommandLineOption noAutoComposeOption{"no_auto_compose", "", "Do not automatically recompose the parallel files."};

    stk::CommandLineParserParallel commandLine(comm);
    commandLine.disallow_unrecognized();

    commandLine.add_required<std::string>(inmeshOption);
    commandLine.add_required<std::string>(outmeshOption);
    commandLine.add_optional<int>(numUniformRefinementLevelsOption, mInputData.algorithmParams.numUniformRefinementLevels);
    commandLine.add_flag(make_flag(force64BitOption), force64BitOption.description);
    commandLine.add_flag(make_flag(assert32BitOption), assert32BitOption.description);
    commandLine.add_flag(make_flag(noAutoComposeOption), noAutoComposeOption.description);

    stk::CommandLineParser::ParseState state = commandLine.parse(argc, const_cast<const char**>(argv));
    if(state == stk::CommandLineParser::ParseComplete)
    {
        mInputData.meshIn = commandLine.get_option_value<std::string>(inmeshOption.name);
        mInputData.meshOut = commandLine.get_option_value<std::string>(outmeshOption.name);

        mInputData.algorithmParams.force64Bit = commandLine.is_option_provided(force64BitOption.name);
        mInputData.algorithmParams.assert32Bit = commandLine.is_option_provided(assert32BitOption.name);
        mInputData.algorithmParams.autoCompose = !commandLine.is_option_provided(noAutoComposeOption.name);

        if(commandLine.is_option_parsed(numUniformRefinementLevelsOption.name))
        {
            mInputData.algorithmParams.numUniformRefinementLevels = commandLine.get_option_value<unsigned>(numUniformRefinementLevelsOption.name);
            if (mInputData.algorithmParams.numUniformRefinementLevels < 1)
            {
                sierra::Env::outputP0() << "ERROR: Number of uniform refinement passes specified via --"+numUniformRefinementLevelsOption.name+" must be >= 1." << std::endl;
                state = stk::CommandLineParser::ParseError;
            }
        }
    }

    if(state != stk::CommandLineParser::ParseComplete)
        sierra::Env::outputP0() << commandLine.get_usage() << std::endl;

    return state == stk::CommandLineParser::ParseComplete;
}

}
