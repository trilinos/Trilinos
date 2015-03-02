// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
//

#ifndef stk_util_DebugTool_hpp
#define stk_util_DebugTool_hpp

#include <execinfo.h>
#include <cxxabi.h>
#include <stk_util/stk_config.h>
#include <string>
#include <mpi.h>


inline int getProcId(MPI_Comm comm)
{
    int procId = -1;
    MPI_Comm_rank(comm, &procId);
    return procId;
}
inline int getProcSize(MPI_Comm comm)
{
    int procSize = -1;
    MPI_Comm_size(comm, &procSize);
    return procSize;
}

inline MPI_Comm splitComm(int color, MPI_Comm globalComm)
{
    int myGlobalId = getProcId(globalComm);
    MPI_Comm localComm;
    MPI_Comm_split(MPI_COMM_WORLD, color, myGlobalId, &localComm);
    return localComm;
}

inline int getDestinationProc(MPI_Comm globalComm)
{
    int numLocalProcs = getProcSize(globalComm) / 3;
    int myProc = getProcId(globalComm) % numLocalProcs;
    int destinationProc = myProc + 2*numLocalProcs;
    return destinationProc;
}

const int STRING_SIZE_TAG = 111;
const int STRING_TAG = 112;
const int MESHDATA_SIZE_TAG = 113;
const int MESHDATA_TAG = 114;
const int FIELDDATA_SIZE_TAG = 115;
const int FIELDDATA_TAG = 116;
//const char *testStringValue = "test string";
const int oneForNullTerminator = 1;

inline std::string demangleFunctionNames(char** symbollist, int addrlen)
{
    std::string mangledNamesString("");
#if defined(__GNUC__) && !defined(__ICC)
    size_t funcnamesize = 256;
    char* funcname = (char*) malloc(funcnamesize);

    for(int i = 1; i < addrlen; i++)
    {
        char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

        for(char *p = symbollist[i]; *p; ++p)
        {
            if(*p == '(')
                begin_name = p;
            else if(*p == '+')
                begin_offset = p;
            else if(*p == ')' && begin_offset)
            {
                end_offset = p;
                break;
            }
        }

        if(begin_name && begin_offset && end_offset && begin_name < begin_offset)
        {
            *begin_name++ = '\0';
            *begin_offset++ = '\0';
            *end_offset = '\0';

            int status;
            char* ret = abi::__cxa_demangle(begin_name, funcname, &funcnamesize, &status);
            if(status == 0)
            {
                mangledNamesString += std::string(ret) + " " + begin_offset + "\n";
            }
            else
            {
                mangledNamesString += std::string(begin_name) + " " + begin_offset + "\n";
            }
        }
        else
        {
            mangledNamesString += symbollist[i] + std::string("\n");
        }
    }

    free(funcname);
#endif
    return mangledNamesString;
}

inline std::string getStackTrace()
{
#if defined(__GNUC__) && !defined(__ICC)
    void *trace[16];
    char **mangledFunctionNames = (char **) NULL;
    int trace_size = 0;

    trace_size = backtrace(trace, 16);
    mangledFunctionNames = backtrace_symbols(trace, trace_size);

    std::string demangledNames = demangleFunctionNames(mangledFunctionNames, trace_size);
    free(mangledFunctionNames);

    return demangledNames;
#else
    return "none";
#endif
}

inline void sendStackTrace(MPI_Comm globalComm)
{
//    std::string stringToSend(testStringValue);
    std::string stringToSend(getStackTrace());
    int lengthOfString = stringToSend.length();
    lengthOfString += oneForNullTerminator;

    MPI_Send(&lengthOfString, 1, MPI_INT, getDestinationProc(globalComm), STRING_SIZE_TAG, globalComm);
    MPI_Send(const_cast<char*>(stringToSend.c_str()), lengthOfString, MPI_CHAR, getDestinationProc(globalComm), STRING_TAG, globalComm);
}

inline void receiveString(int source, MPI_Comm localComm, MPI_Comm globalComm)
{
    MPI_Status status;
    int lengthOfString = -1;
    MPI_Recv(&lengthOfString, 1, MPI_INT, source, STRING_SIZE_TAG, globalComm, &status);

//    EXPECT_EQ(std::string(testStringValue).length(), static_cast<size_t>(lengthOfString-oneForNullTerminator));

    char *receiveString = new char[lengthOfString];
    MPI_Recv(receiveString, lengthOfString, MPI_CHAR, source, STRING_TAG, globalComm, &status);

//    EXPECT_EQ(std::string(testStringValue), std::string(receiveString));

    if(getProcId(localComm) == 0)
    {
        std::cerr << "Stack Trace from proc " << source << ":\n" << receiveString << std::endl;
    }

    delete[] receiveString;
}

inline void receiveStackTrace(MPI_Comm localComm, MPI_Comm globalComm)
{
    int source0 = getProcId(localComm);
    receiveString(source0, localComm, globalComm);
    int source1 = getProcId(localComm) + getProcSize(localComm);
    receiveString(source1, localComm, globalComm);
}

//inline void sendMeshAndFieldData(std::vector<uint64_t> &meshData, std::vector<double> &meshFieldData, MPI_Comm globalComm)
//{
//    size_t sizeMeshData = meshData.size();
//    MPI_Send(&sizeMeshData, 1, MPI_UINT64_T, getDestinationProc(globalComm), MESHDATA_SIZE_TAG, globalComm);
//    MPI_Send(&meshData[0], sizeMeshData, MPI_UINT64_T, getDestinationProc(globalComm), MESHDATA_TAG, globalComm);
//    size_t sizeFieldData = meshFieldData.size();
//    MPI_Send(&sizeFieldData, 1, MPI_UINT64_T, getDestinationProc(globalComm), FIELDDATA_SIZE_TAG, globalComm);
//    MPI_Send(&meshFieldData[0], sizeFieldData, MPI_DOUBLE, getDestinationProc(globalComm), FIELDDATA_TAG, globalComm);
//}
//
//inline void receiveMeshAndFieldData(std::vector<uint64_t> &meshData, std::vector<double> &meshFieldData,
//                                    int source, MPI_Comm globalComm)
//{
//    MPI_Status status;
//    int sizeMeshData = -1;
//    MPI_Recv(&sizeMeshData, 1, MPI_UINT64_T, source, MESHDATA_SIZE_TAG, globalComm, &status);
//
//    meshData.clear();
//    meshData.resize(sizeMeshData);
//    MPI_Recv(&meshData[0], sizeMeshData, MPI_UINT64_T, source, MESHDATA_TAG, globalComm, &status);
//
//    int sizeMeshFieldData = -1;
//    MPI_Recv(&sizeMeshFieldData, 1, MPI_UINT64_T, source, FIELDDATA_SIZE_TAG, globalComm, &status);
//
//    meshFieldData.clear();
//    meshFieldData.resize(sizeMeshFieldData);
//    MPI_Recv(&meshFieldData[0], sizeMeshFieldData, MPI_DOUBLE, source, FIELDDATA_TAG, globalComm, &status);
//}

//inline void receiveMeshData(std::vector<uint64_t> &meshData, std::vector<double> &meshFieldData,
//                                    int source, MPI_Comm globalComm)
//{
//
//}
//
//inline void sendMeshData(sierra::Fmwk::MeshBulkData &fmwk_mesh, stk::mesh::Selector selector, const stk::mesh::FieldVector &fields, MPI_Comm globalComm)
//{
//    stk::mesh::BulkData &mesh = fmwk_mesh.get_stk_mesh();
//    std::vector<uint64_t> meshData;
//    std::vector<double> meshFieldData;
//    for(stk::mesh::EntityRank entity_rank = stk::topology::NODE_RANK; entity_rank <= stk::topology::CONSTRAINT_RANK; ++entity_rank)
//    {
//        stk::mesh::EntityVector entities;
//        entities.reserve(1000);
//        const stk::mesh::BucketVector &buckets = mesh.buckets(entity_rank);
//
//        meshData.push_back(buckets.size());
//        for(size_t bucketIndex=0; bucketIndex<buckets.size(); bucketIndex++)
//        {
//            const stk::mesh::Bucket &buck = *buckets[bucketIndex];
//            if (selector(buck))
//            {
//                meshData.push_back(buck.size());
//                const stk::mesh::PartVector& parts = buck.supersets();
//                meshData.push_back(parts.size());
//                for(size_t partIndex=0; partIndex<parts.size(); ++partIndex)
//                {
//                    if (!fmwk_mesh.mesh_meta_data().is_special_framework_part(*parts[partIndex]))
//                    {
//                        meshData.push_back(parts[partIndex]->mesh_meta_data_ordinal());
//                    }
//                }
//                for(size_t entityIndex=0; entityIndex<buck.size(); entityIndex++)
//                {
//                    entities.push_back(buck[entityIndex]);
//                    meshData.push_back(mesh.entity_key(buck[entityIndex]));
//                    meshData.push_back(mesh.parallel_owner_rank(entities[entityIndex]));
//                    unsigned totalNumConnected = 0;
//                    for(stk::mesh::EntityRank connected_rank = stk::topology::NODE_RANK; connected_rank <= stk::topology::CONSTRAINT_RANK; ++connected_rank)
//                    {
//                        if(connected_rank != entity_rank)
//                        {
//                            totalNumConnected += mesh.num_connectivity(buck[entityIndex], connected_rank);
//                        }
//                    }
//                    meshData.push_back(totalNumConnected);
//                    for(stk::mesh::EntityRank connected_rank = stk::topology::NODE_RANK; connected_rank <= stk::topology::CONSTRAINT_RANK; ++connected_rank)
//                    {
//                        if(connected_rank != entity_rank)
//                        {
//                            unsigned numConnected = mesh.num_connectivity(buck[entityIndex], connected_rank);
//                            const stk::mesh::Entity *connectedEntities = mesh.begin(buck[entityIndex], connected_rank);
//                            for(unsigned connIndex=0; connIndex<numConnected; connIndex++)
//                            {
//                                meshData.push_back(mesh.entity_key(connectedEntities[connIndex]));
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        // stk::mesh::get_selected_entities(selector, mesh.buckets(entity_rank), entities);
//        meshData.push_back(entities.size());
//        for(size_t entityIndex=0; entityIndex<entities.size(); entityIndex++)
//        {
//            size_t numFieldsToSend = 0;
//            for(size_t fieldIndex=0; fieldIndex<fields.size(); fieldIndex++)
//            {
//                if(fields[fieldIndex]->entity_rank() == entity_rank)
//                {
//                    int numScalars = stk::mesh::field_scalars_per_entity(*fields[fieldIndex], entities[entityIndex]);
//                    if(numScalars > 0)
//                    {
//                        numFieldsToSend++;
//                    }
//                }
//            }
//            meshData.push_back(numFieldsToSend);
//            for(size_t fieldIndex=0; fieldIndex<fields.size(); fieldIndex++)
//            {
//                if(fields[fieldIndex]->entity_rank() == entity_rank)
//                {
//                    int numScalars = stk::mesh::field_scalars_per_entity(*fields[fieldIndex], entities[entityIndex]);
//                    if(numScalars > 0)
//                    {
//                        meshData.push_back(fields[fieldIndex]->mesh_meta_data_ordinal());
//                        meshData.push_back(numScalars);
//
//                        void *fieldData = stk::mesh::field_data(*fields[fieldIndex], entities[entityIndex]);
//                        if(fields[fieldIndex]->data_traits().is_floating_point && fields[fieldIndex]->data_traits().size_of == 8)
//                        {
//                            double *doubleValues = static_cast<double*>(fieldData);
//                            for(int i=0; i<numScalars; i++)
//                            {
//                                meshFieldData.push_back(doubleValues[i]);
//                            }
//                        }
//                        else if(fields[fieldIndex]->data_traits().is_integral && fields[fieldIndex]->data_traits().size_of == 4)
//                        {
//                            int *intValues = static_cast<int *>(fieldData);
//                            for(int i=0; i<numScalars; i++)
//                            {
//                                meshFieldData.push_back(static_cast<double>(intValues[i]));
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//    sendMeshAndFieldData(meshData, meshFieldData, globalComm);
//}

#endif

