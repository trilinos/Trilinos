
#ifndef BALANCETESTUTILITIES_HPP_
#define BALANCETESTUTILITIES_HPP_

#include <vector>
#include <string>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <test_utils/NemesisInfo.hpp>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class FieldBase; } }

namespace balance_utils
{
void putFieldDataOnMesh(stk::mesh::BulkData &bulkData, const std::vector<int>& coloring);
int getWidth(int numProcsDecomp);
void clearFiles(const std::string &baseFilename, int numProcs);
std::string getFilename(const std::string& filename, int numProcsDecomp, int subdomainId);
void putVectorDataIntoField(stk::mesh::BulkData &bulkData, const stk::mesh::EntityProcVec& vector_of_data, stk::mesh::FieldBase *field);
void putDecompFieldDataOnMesh(stk::mesh::BulkData &bulkData, int proc_rank);
void putEntityProcOnMeshField(stk::mesh::BulkData &bulkData, const stk::mesh::EntityProcVec& coloring);
}


#endif /* BALANCETESTUTILITIES_HPP_ */
