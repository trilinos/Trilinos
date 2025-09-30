#include "stk_middle_mesh/communication_api.hpp"
#include "stk_middle_mesh_util/create_stk_mesh.hpp"
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include "stk_util/parallel/Parallel.hpp"
#include "stk_unit_test_utils/ConservativeTransferUserExample.hpp"


void check_conservation(MPI_Comm unionComm, 
                        mesh::FieldPtr<double> functionValsSend, mesh::FieldPtr<double> functionValsRecv,
                        std::shared_ptr<ConservativeTransferUserForTest> sendCallback,
                        std::shared_ptr<ConservativeTransferUserForTest> recvCallback)
{
  int myRank = utils::impl::comm_rank(unionComm);
  if (functionValsSend)
  {
    double integralVal = sendCallback->integrate_function(functionValsSend)[0];

    MPI_Send(&integralVal, 1, MPI_DOUBLE, 1 - myRank, 666, unionComm);
  } else
  {
    double integralVal = recvCallback->integrate_function(functionValsRecv)[0];   

    double integralValSender;
    MPI_Recv(&integralValSender, 1, MPI_DOUBLE, 1 - myRank, 666, unionComm, MPI_STATUS_IGNORE);

    std::cout << "integralValSender = " << integralValSender << ", integralValRecver = " << integralVal << ", diff = " << std::abs(integralValSender - integralVal) << std::endl;
    if (std::abs(integralValSender - integralVal) > 1e-12)
      throw std::runtime_error("transfer was not conservative");
  }
}

template <typename Tfunc>
mesh::FieldPtr<double> create_field(std::shared_ptr<mesh::Mesh> mesh, Tfunc func)
{
  auto fieldPtr = mesh::create_field<double>(mesh, mesh::FieldShape(1, 0, 0), 1);

  auto& field = *fieldPtr;
  for (auto& vert : mesh->get_vertices())
    if (vert)
    {
      utils::Point pt = vert->get_point_orig(0);
      field(vert, 0, 0) = func(pt);
    }

  return fieldPtr;
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
      //std::cout << "field = " << field(vert, 0, 0) << ", val expected = " << valExpected << ", error = " << std::abs(field(vert, 0, 0) - valExpected) << std::endl;
      if (std::abs(field(vert, 0, 0) - valExpected) > 1e-12)
        throw std::runtime_error("field transfer was not exact");
    }
}


std::function<double(const utils::Point&)> function_factory(const std::string& functionName)
{
  if (functionName == "linear")
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

int main(int argc, char* argv[])
{
  stk::initialize(&argc, &argv);

  int defaultColor = -1;
  int color = stk::get_command_line_option(argc, argv, "app-color", defaultColor);
  STK_ThrowRequireMsg(color == 0 || color == 1, "app-color command line argument must be provided, and must either 0 or 1");

  int sendColor = stk::get_command_line_option(argc, argv, "send-color", defaultColor);
  STK_ThrowRequireMsg(sendColor == 0 || sendColor == 1, "send-color command line argument must be provided, and must either 0 or 1");

  std::string defaultFileName;
  if (color == 0)
    defaultFileName = "generated:3x3x1|sideset:Z|bbox:0,0,0,1,1,1";
  else
    defaultFileName = "generated:4x4x1|sideset:z|bbox:0,0,1,1,1,2";

  std::string meshFileName = stk::get_command_line_option(argc, argv, "mesh", defaultFileName);

  std::string defaultPartName = "surface_1";
  std::string partName = stk::get_command_line_option(argc, argv, "part-name", defaultPartName);

  std::string defaultFunctionName = "linear";
  std::string functionName = stk::get_command_line_option(argc, argv, "function-name", defaultFunctionName);

  {
    MPI_Comm meshComm;
    MPI_Comm_split(MPI_COMM_WORLD, color, 0, &meshComm);

    stk_interface::StkMeshCreator creator(meshFileName, "NONE", meshComm);
    std::shared_ptr<mesh::Mesh> inputMesh = creator.create_mesh_from_part(partName).mesh;

    std::shared_ptr<mesh::Mesh> inputMesh1 = color == 0 ? inputMesh : nullptr;
    std::shared_ptr<mesh::Mesh> inputMesh2 = color == 0 ? nullptr   : inputMesh;

    auto transferCallback = std::make_shared<ConservativeTransferUserForTest>(inputMesh);
    std::shared_ptr<ConservativeTransferUserForTest> transferCallback1 = color == 0 ? transferCallback : nullptr;
    std::shared_ptr<ConservativeTransferUserForTest> transferCallback2 = color == 0 ? nullptr : transferCallback;
    if (color == 0)
      std::cout << "mesh1 numel = " << inputMesh1->get_elements().size() << std::endl;

    if (color == 1)
      std::cout << "mesh2 numel = " << inputMesh2->get_elements().size() << std::endl;    
      
    MPI_Comm unionComm = MPI_COMM_WORLD;
    stk::transfer::ConservativeTransfer transfer(unionComm, inputMesh1, inputMesh2, transferCallback1, transferCallback2);

    bool amISender = color == sendColor;
    auto func = function_factory(functionName);
    mesh::FieldPtr<double> functionValsSend, functionValsRecv;
    if (amISender) {
      functionValsSend = create_field(inputMesh, func);
      functionValsRecv = nullptr;
    } else {
      functionValsSend = nullptr;
      functionValsRecv = mesh::create_field<double>(inputMesh, mesh::FieldShape(1, 0, 0), 1);
    }

    
    transfer.start_transfer(functionValsSend, functionValsRecv);
    transfer.finish_transfer();

    auto transferCallbackSend = amISender ? transferCallback : nullptr;
    auto transferCallbackRecv = amISender ? nullptr : transferCallback;
    check_conservation(unionComm, functionValsSend, functionValsRecv, transferCallbackSend, transferCallbackRecv);

    if (!amISender && functionName == "linear")
      check_field(functionValsRecv, func);

    MPI_Comm_free(&meshComm);
  }

  stk::finalize();
}
