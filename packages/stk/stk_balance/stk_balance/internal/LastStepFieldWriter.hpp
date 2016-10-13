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
        fill_mesh_and_time_data(inputFilename);
    }

    ~LastStepFieldWriter() {}

    void write_mesh(const std::string& filename,
                    stk::io::DatabasePurpose databasePurpose = stk::io::WRITE_RESULTS)
    {
        size_t fh = stkIo.create_output_mesh(filename, databasePurpose);

        if(has_field_data())
        {
            add_transient_fields(fh);

            stkIo.begin_output_step(fh, maxTime);
            stkIo.write_defined_output_fields(fh);
            stkIo.end_output_step(fh);
        }
        else
        {
            stkIo.write_output_mesh(fh);        }
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

    void fill_time_data()
    {
        numSteps = stkIo.get_num_time_steps();
        if(numSteps>0)
        {
            stkIo.read_defined_input_fields(numSteps);
            maxTime = stkIo.get_max_time();
        }
    }

    void add_transient_fields(size_t fileHandle)
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

    void fill_mesh_and_time_data(const std::string& inputFilename)
    {
        stk::io::fill_mesh_preexisting(stkIo, inputFilename, bulkData);
        fill_time_data();
    }

    stk::mesh::BulkData& bulkData;
    int numSteps;
    double maxTime;
    stk::io::StkMeshIoBroker stkIo;
};

class LastStepFieldWriterAutoDecomp : public LastStepFieldWriter
{
public:
    LastStepFieldWriterAutoDecomp(stk::mesh::BulkData& bulk, const std::string& inputFilename)
    : LastStepFieldWriter(bulk)
    {
        stkIo.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));
        fill_mesh_and_time_data(inputFilename);
    }
};

}}}

#endif
