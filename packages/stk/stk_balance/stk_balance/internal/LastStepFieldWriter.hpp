// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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
    write_mesh_with_fields(fh);
  }

  void write_mesh_with_additional_fields(const std::string& filename,
                                         const stk::mesh::FieldVector& fieldsCreatedDuringRuntime,
                                         stk::io::DatabasePurpose databasePurpose = stk::io::WRITE_RESULTS)
  {
    if(fieldsCreatedDuringRuntime.size()>0)
    {
      mHaveFields = true;
    }

    size_t fh = stkIo.create_output_mesh(filename, databasePurpose);

    for(stk::mesh::FieldBase * out_field : fieldsCreatedDuringRuntime)
      stkIo.add_field(fh, *out_field);

    write_mesh_with_fields(fh);
  }


  void set_output_time(double outputTime)
  {
    maxTime = outputTime;
    numSteps = 1;
  }

  double get_max_time() const { return maxTime; }

protected:
  LastStepFieldWriter(stk::mesh::BulkData& bulk) : bulkData(bulk), numSteps(-1), maxTime(0.0) {}

  void write_mesh_with_fields(size_t fileHandle)
  {
    if(has_field_data())
    {
      const stk::mesh::FieldVector& fieldsFromInputFile = stkIo.meta_data().get_fields();
      add_transient_fields(fileHandle, fieldsFromInputFile);

      stkIo.begin_output_step(fileHandle, maxTime);
      stkIo.write_defined_output_fields(fileHandle);
      stkIo.end_output_step(fileHandle);
    }
    else
    {
      stkIo.write_output_mesh(fileHandle);
    }
  }

  bool has_field_data() const
  {
    return numSteps>0 || mHaveFields;
  }

  void fill_time_data()
  {
    numSteps = stkIo.get_num_time_steps();
    if(numSteps>0)
    {
      stkIo.read_defined_input_fields(numSteps);
      maxTime = stkIo.get_max_time();
    }
  }

  void add_transient_fields(size_t fileHandle, const stk::mesh::FieldVector& out_fields)
  {
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
  bool mHaveFields = false;
};

class AllStepFieldWriterAutoDecomp : public LastStepFieldWriter
{
public:
  AllStepFieldWriterAutoDecomp(stk::mesh::BulkData& bulk, const std::string& inputFilename)
    : LastStepFieldWriter(bulk)
  {
    stkIo.property_add(Ioss::Property("DECOMPOSITION_METHOD", "LINEAR"));
    fill_mesh_and_time_data(inputFilename);
  }
};

}}}

#endif
