#ifndef STK_BALANCE_LAST_STEP_FIELD_WRITER_H
#define STK_BALANCE_LAST_STEP_FIELD_WRITER_H

#include <string>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include "stk_mesh/base/Field.hpp"

#include <stk_util/parallel/Parallel.hpp>

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_io/FillMesh.hpp>

namespace stk { namespace balance { namespace internal {

class LastStepFieldWriter
{
public:

    LastStepFieldWriter(stk::mesh::BulkData& bulk, const std::string& inputFilename)
    : LastStepFieldWriter(bulk)
    {
        stk::io::StkMeshIoBroker stkIo;
        fill_mesh_and_time_data(stkIo, inputFilename);
    }

    ~LastStepFieldWriter() {}

    void write_mesh(const std::string& filename,
                    stk::io::DatabasePurpose databasePurpose = stk::io::WRITE_RESULTS)
    {
        if(has_field_data())
        {
            stk::io::StkMeshIoBroker broker(bulkData.parallel());
            broker.set_bulk_data(bulkData);
            size_t fh = broker.create_output_mesh(filename, databasePurpose);
            add_transient_fields(broker, fh);

            broker.begin_output_step(fh, maxTime);
            broker.write_defined_output_fields(fh);
            broker.end_output_step(fh);
        }
        else
        {
            stk::io::write_mesh(filename, bulkData, databasePurpose);
        }
    }

    void set_output_time(double outputTime)
    {
        maxTime = outputTime;
        numSteps = 1;
    }

    double get_max_time() const { return maxTime; }

protected:
    LastStepFieldWriter(stk::mesh::BulkData& bulk) : bulkData(bulk), numSteps(-1), maxTime(0.0) {}

    bool has_field_data() const { return numSteps>0; }

    void fill_time_data(stk::io::StkMeshIoBroker& stkIo)
    {
        numSteps = stkIo.get_num_time_steps();
        if(numSteps>0)
        {
            stkIo.read_defined_input_fields(numSteps);
            maxTime = stkIo.get_max_time();
        }
    }

    void add_transient_fields(stk::io::StkMeshIoBroker& stkIo, size_t fileHandle) const
    {
        const stk::mesh::FieldVector& out_fields = stkIo.meta_data().get_fields();

        for(size_t i=0;i<out_fields.size();++i)
        {
            const Ioss::Field::RoleType* fieldRole = stk::io::get_field_role(*out_fields[i]);
            if(fieldRole != nullptr && *fieldRole == Ioss::Field::TRANSIENT)
            {
                stkIo.add_field(fileHandle, *out_fields[i]);
            }
        }
    }

    void fill_mesh_and_time_data(stk::io::StkMeshIoBroker& stkIo, const std::string& inputFilename)
    {
        stk::io::fill_mesh_preexisting(stkIo, inputFilename, bulkData);
        fill_time_data(stkIo);
    }

    stk::mesh::BulkData& bulkData;
    int numSteps;
    double maxTime;
};

class LastStepFieldWriterAutoDecomp : public LastStepFieldWriter
{
public:
    LastStepFieldWriterAutoDecomp(stk::mesh::BulkData& bulk, const std::string& inputFilename)
    : LastStepFieldWriter(bulk)
    {
        stk::io::StkMeshIoBroker stkIo;
        stkIo.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));
        fill_mesh_and_time_data(stkIo, inputFilename);
    }
};

}}}

#endif
