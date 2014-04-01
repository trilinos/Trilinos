#include <gtest/gtest.h>
#include <mpi.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <Ioss_SubSystem.h>
#include <stk_util/util/ParameterList.hpp>

inline bool fieldWithNameChangedIsOutput(stk::io::StkMeshIoBroker &stkIo, MPI_Comm communicator, const size_t resultsOutputIndex, const std::string &goldFieldName)
{
    double dummyTime = 0;
    stkIo.process_output_request(resultsOutputIndex, dummyTime);
    Ioss::Region *outputRegion = stkIo.get_output_io_region(resultsOutputIndex).get();
    Ioss::NodeBlock *nodeBlockAssociatedWithDisplacementField = outputRegion->get_node_blocks()[0];
    Ioss::NameList fieldNames;
    nodeBlockAssociatedWithDisplacementField->field_describe(Ioss::Field::TRANSIENT, &fieldNames);

    return (goldFieldName == fieldNames[0]);
}

inline void validate_parameters_equal_value(const stk::util::Parameter &parameter,
					    const stk::util::Parameter &gold_parameter)
{
  ASSERT_EQ(parameter.type, gold_parameter.type);
  switch(parameter.type)
    {
    case stk::util::ParameterType::INTEGER:
      {
	ASSERT_EQ(boost::any_cast<int>(parameter.value),
		  boost::any_cast<int>(gold_parameter.value));
	break;
      }
    case stk::util::ParameterType::INT64:
      {
	ASSERT_EQ(boost::any_cast<int64_t>(parameter.value),
		  boost::any_cast<int64_t>(gold_parameter.value));
	break;
      }
    case stk::util::ParameterType::DOUBLE:
      {
	ASSERT_EQ(boost::any_cast<double>(parameter.value),
		  boost::any_cast<double>(gold_parameter.value));
	break;
      }
    case stk::util::ParameterType::FLOAT:
      {
	ASSERT_EQ(boost::any_cast<float>(parameter.value),
		  boost::any_cast<float>(gold_parameter.value));
	break;
      }
    case stk::util::ParameterType::DOUBLEVECTOR:
      {
	std::vector<double> vec = boost::any_cast<std::vector<double> >(parameter.value);
	std::vector<double> gvec = boost::any_cast<std::vector<double> >(gold_parameter.value);
	ASSERT_EQ(vec.size(), gvec.size());
	for (size_t j = 0; j < vec.size(); ++j) {
	  ASSERT_EQ(vec[j], gvec[j]);
	}
	break;
      }
    case stk::util::ParameterType::FLOATVECTOR:
      {
	std::vector<float> vec = boost::any_cast<std::vector<float> >(parameter.value);
	std::vector<float> gvec = boost::any_cast<std::vector<float> >(gold_parameter.value);
	ASSERT_EQ(vec.size(), gvec.size());
	for (size_t j = 0; j < vec.size(); ++j) {
	  ASSERT_EQ(vec[j], gvec[j]);
	}
	break;
      }
    case stk::util::ParameterType::INTEGERVECTOR:
      {
	std::vector<int> vec = boost::any_cast<std::vector<int> >(parameter.value);
	std::vector<int> gvec = boost::any_cast<std::vector<int> >(gold_parameter.value);
	ASSERT_EQ(vec.size(), gvec.size());
	for (size_t j = 0; j < vec.size(); ++j) {
	  ASSERT_EQ(vec[j], gvec[j]);
	}
	break;
      }
    case stk::util::ParameterType::INT64VECTOR:
      {
	std::vector<int64_t> vec = boost::any_cast<std::vector<int64_t> >(parameter.value);
	std::vector<int64_t> gvec = boost::any_cast<std::vector<int64_t> >(gold_parameter.value);
	ASSERT_EQ(vec.size(), gvec.size());
	for (size_t j = 0; j < vec.size(); ++j) {
	  ASSERT_EQ(vec[j], gvec[j]);
	}
	break;
      }
    default:
      ASSERT_EQ(1,0) << "Invalid type found in validate_parameters_equal_value";
    }
}
