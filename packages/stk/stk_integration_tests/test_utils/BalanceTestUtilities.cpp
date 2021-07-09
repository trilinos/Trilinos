#include "BalanceTestUtilities.hpp"
#include <iomanip>
#include <unistd.h>
#include <stk_mesh/base/GetEntities.hpp>
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_mesh/base/Field.hpp>      // for Field
#include "stk_balance/internal/entityDataToField.hpp"

namespace balance_utils
{

stk::mesh::EntityProcVec convert_vector_data_into_entity_proc_vec(stk::mesh::BulkData &bulkData, const std::vector<int>& coloring)
{
    stk::mesh::EntityProcVec entity_coloring(coloring.size());
    const stk::mesh::BucketVector &buckets = bulkData.buckets(stk::topology::ELEMENT_RANK);
    unsigned elem_counter = 0;
    for(size_t i = 0; i < buckets.size(); i++)
    {
        const stk::mesh::Bucket &bucket = *buckets[i];
        if(bucket.owned())
        {
            for(size_t j = 0; j < bucket.size(); j++)
            {
                entity_coloring[elem_counter] = std::make_pair(bucket[j], coloring[elem_counter]);
                ++elem_counter;
            }
        }
    }
    return entity_coloring;
}

void putEntityProcOnMeshField(stk::mesh::BulkData &bulkData, const stk::mesh::EntityProcVec& coloring)
{
    const std::string fieldName2 = "Coloring";
    stk::mesh::Field<double> &field2 = *bulkData.mesh_meta_data().get_field<stk::mesh::Field<double> >(stk::topology::ELEMENT_RANK, fieldName2);
    stk::balance::internal::put_entity_data_to_field(coloring, &field2);
    putDecompFieldDataOnMesh(bulkData, bulkData.parallel_rank());
}

void putFieldDataOnMesh(stk::mesh::BulkData &bulkData, const std::vector<int>& coloring)
{
    stk::mesh::EntityProcVec entity_coloring = convert_vector_data_into_entity_proc_vec(bulkData, coloring);
    putEntityProcOnMeshField(bulkData, entity_coloring);
}

void putDecompFieldDataOnMesh(stk::mesh::BulkData &bulkData, int proc_rank)
{
    const stk::mesh::BucketVector &buckets = bulkData.buckets(stk::topology::ELEMENT_RANK);

    const std::string fieldName1 = "DomainProc";
    stk::mesh::Field<double> &field1 = *bulkData.mesh_meta_data().get_field<stk::mesh::Field<double> >(stk::topology::ELEMENT_RANK, fieldName1);

    for(size_t i = 0; i < buckets.size(); i++)
    {
        const stk::mesh::Bucket &bucket = *buckets[i];
        if(bucket.owned())
        {
            for(size_t j = 0; j < bucket.size(); j++)
            {
                double *decomp_data = stk::mesh::field_data(field1, bucket[j]);
                *decomp_data = static_cast<double>(proc_rank);
            }
        }
    }
}

int getWidth(int numProcsDecomp)
{
    return std::log10(static_cast<double>(numProcsDecomp - 1))+1;
}

std::string getFilename(const std::string& filename, int numProcsDecomp, int subdomainId)
{
    int width = balance_utils::getWidth(numProcsDecomp);
    std::ostringstream os;
    os << filename << "." << numProcsDecomp << "." << std::setfill('0') << std::setw(width) << subdomainId;
    return os.str();
}

void clearFiles(const std::string &baseFilename, int numProcs)
{
    for(int i = 0; i < numProcs; i++)
    {
        unlink(getFilename(baseFilename, numProcs, i).c_str());
    }
}

}

