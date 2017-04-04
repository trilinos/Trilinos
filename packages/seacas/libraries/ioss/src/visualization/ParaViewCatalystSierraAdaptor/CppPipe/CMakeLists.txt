cmake_minimum_required(VERSION 2.8.8)
project(CppPipe)

find_package(ParaView 4.1 REQUIRED COMPONENTS vtkPVCatalyst)
include("${PARAVIEW_USE_FILE}")
ADD_DEFINITIONS("-DUSE_VTK_JSON")

set(CppPipe_SRCS
    vtkCppPipe.cxx
    CppPipeRepresentation.cxx
    CppPipeImageSet.cxx
    CppPipeProxies.cxx
    CppPipeMultiCamera8.cxx
    CppPipeCamera.cxx
    CppPipeSlice.cxx
    CppPipeClip.cxx
    CppPipeBoxClip.cxx
    CppPipeThreshold.cxx
    CppPipeContour.cxx
)

add_library(CppPipe SHARED ${CppPipe_SRCS})

target_link_libraries(CppPipe vtkPVCatalyst vtkjsoncpp)
