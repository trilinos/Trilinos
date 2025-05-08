#include "stk_middle_mesh/application_interface.hpp"
#include "stk_transfer/ConservativeTransfer.hpp"
#include "stk_unit_test_utils/ConservativeTransferUserExample.hpp"

#include "stk_middle_mesh_util/create_stk_mesh.hpp"
#include "stk_middle_mesh_util/exodus_writer.hpp"
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include "stk_util/parallel/Parallel.hpp"



void check_conservation(mesh::FieldPtr<double> functionValsSend, mesh::FieldPtr<double> functionValsRecv,
                        std::shared_ptr<ConservativeTransferUserForTest> sendCallback,
                        std::shared_ptr<ConservativeTransferUserForTest> recvCallback)
{
  assert(functionValsSend);
  assert(functionValsRecv);

  double integralValSend = sendCallback->integrate_function(functionValsSend)[0];
  double integralValRecv = recvCallback->integrate_function(functionValsRecv)[0];
        
  if (std::abs(integralValSend - integralValRecv) > 1e-12)
    throw std::runtime_error("transfer was not conservative");
}

template <typename Tfunc>
void set_field(mesh::FieldPtr<double> fieldPtr, Tfunc func)
{
  auto& field = *fieldPtr;
  for (auto& vert : fieldPtr->get_mesh()->get_vertices())
    if (vert)
    {
      utils::Point pt = vert->get_point_orig(0);
      field(vert, 0, 0) = func(pt);
    }
}


template <typename Tfunc>
void check_field(mesh::FieldPtr<double> fieldPtr, Tfunc func)
{
  auto& field = *fieldPtr;
  for (auto& vert : field.get_mesh()->get_vertices())
    if (vert)
    {
      utils::Point pt = vert->get_point_orig(0);
      double valExpected = func(pt);
      std::cout << "field = " << field(vert, 0, 0) << ", val expected = " << valExpected << ", error = " << std::abs(field(vert, 0, 0) - valExpected) << std::endl;
      if (std::abs(field(vert, 0, 0) - valExpected) > 1e-12)
        throw std::runtime_error("field transfer was not exact");
    }
}


std::function<double(const utils::Point&)> function_factory(const std::string& functionName)
{
  if (functionName == "constant")
  {
    return [](const utils::Point& /*pt*/) { return 1; };  
  } else if (functionName == "linear")
  {
    return [](const utils::Point& pt) { return pt.x + 2*pt.y + 3*pt.z; };
  } else if (functionName == "quadratic")
  {
    return [](const utils::Point& pt) { return pt.x*pt.x + 2*pt.y*pt.y + 3*pt.z; };
  } else if (functionName == "exponential")
  {
    return [](const utils::Point& pt) { return std::exp(pt.x + pt.y + pt.z ); };
  } else
    throw std::runtime_error("unrecognized function name: " + functionName);
}

void write_output(mesh::FieldPtr<double> field1, mesh::FieldPtr<double> field2, const std::string& functionName)
{
  std::shared_ptr<mesh::Mesh> inputMesh1 = field1->get_mesh();
  std::shared_ptr<mesh::Mesh> inputMesh2 = field2->get_mesh();

  auto func = function_factory(functionName);
  auto field1Exact = mesh::create_field<double>(inputMesh1, mesh::FieldShape(1, 0, 0), 1);
  auto field2Exact = mesh::create_field<double>(inputMesh2, mesh::FieldShape(1, 0, 0), 1);
  set_field(field1Exact, func);
  set_field(field2Exact, func);
  
  auto field1Adaptor = std::make_shared<stk_interface::impl::FieldOutputAdaptorDouble>(field1, "field");
  auto field1ExactAdaptor = std::make_shared<stk_interface::impl::FieldOutputAdaptorDouble>(field1Exact, "fieldExact");
  stk_interface::impl::ExodusWriter writer1(inputMesh1, {field1Adaptor, field1ExactAdaptor});
  writer1.write("mesh1_conservative.exo");

  auto field2Adaptor = std::make_shared<stk_interface::impl::FieldOutputAdaptorDouble>(field2, "field");
  auto field2ExactAdaptor = std::make_shared<stk_interface::impl::FieldOutputAdaptorDouble>(field2Exact, "fieldExact");
  stk_interface::impl::ExodusWriter writer2(inputMesh2, {field2Adaptor, field2ExactAdaptor});
  writer2.write("mesh2_conservative.exo");   
}


int main(int argc, char* argv[])
{
  stk::initialize(&argc, &argv);

  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    throw std::runtime_error("mock app only works on 1 process");

  int sendOneToTwoDefault = true;
  bool sendOneToTwo = stk::get_command_line_option(argc, argv, "send-one-to-two", sendOneToTwoDefault);

  std::string defaultFileName1 = "generated:3x3x1|sideset:Z|bbox:0,0,0,1,1,1";
  std::string defaultFileName2 = "generated:4x4x1|sideset:z|bbox:0,0,1,1,1,2";

  std::string defaultFunctionName = "linear";
  std::string functionName = stk::get_command_line_option(argc, argv, "function-name", defaultFunctionName);

  std::string meshFileName1 = stk::get_command_line_option(argc, argv, "mesh1", defaultFileName1);
  std::string meshFileName2 = stk::get_command_line_option(argc, argv, "mesh2", defaultFileName2);

  std::string defaultPartName1 = "surface_1";
  std::string defaultPartName2 = "surface_1";
  std::string partName1 = stk::get_command_line_option(argc, argv, "part-name1", defaultPartName1);
  std::string partName2 = stk::get_command_line_option(argc, argv, "part-name2", defaultPartName2);

  int defaultNumIters = 64;
  int numIters = stk::get_command_line_option(argc, argv, "num-iters", defaultNumIters);

  {
    stk_interface::StkMeshCreator creator1(meshFileName1, "NONE", MPI_COMM_WORLD);
    std::shared_ptr<mesh::Mesh> inputMesh1 = creator1.create_mesh_from_part(partName1).mesh;

    stk_interface::StkMeshCreator creator2(meshFileName2, "NONE", MPI_COMM_WORLD);
    std::shared_ptr<mesh::Mesh> inputMesh2 = creator2.create_mesh_from_part(partName2).mesh;

    auto transferCallback1 = std::make_shared<ConservativeTransferUserForTest>(inputMesh1);
    auto transferCallback2 = std::make_shared<ConservativeTransferUserForTest>(inputMesh2);
    MPI_Comm unionComm = MPI_COMM_WORLD;
    stk::transfer::ConservativeTransfer transfer(unionComm, inputMesh1, inputMesh2, transferCallback1, transferCallback2);

    auto func = function_factory(functionName);
    auto functionVals1 = mesh::create_field<double>(inputMesh1, mesh::FieldShape(1, 0, 0), 1);
    auto functionVals2 = mesh::create_field<double>(inputMesh2, mesh::FieldShape(1, 0, 0), 1);

    auto functionValsSend = sendOneToTwo ? functionVals1 : functionVals2;
    auto functionValsRecv = sendOneToTwo ? functionVals2 : functionVals1;
    auto transferCallbackSend = sendOneToTwo ? transferCallback1 : transferCallback2;
    auto transferCallbackRecv = sendOneToTwo ? transferCallback2 : transferCallback1;
    set_field(functionValsSend, func);
    
    for (int i=0; i < numIters; ++i)
    {
      transfer.start_transfer(functionValsSend, functionValsRecv);
      transfer.finish_transfer();

      transfer.start_transfer(functionValsRecv, functionValsSend);
      transfer.finish_transfer();

      check_conservation(functionValsSend, functionValsRecv, transferCallbackSend, transferCallbackRecv);
    }

    if (functionName == "linear")
      check_field(functionValsRecv, func);

    write_output(functionVals1, functionVals2, functionName);
  }

  stk::finalize();
}
