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

#ifndef DEBUGWRITER_HPP
#define DEBUGWRITER_HPP

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>                     // for BulkData
#include <stk_mesh/base/DumpMeshInfo.hpp>
#include <stddef.h>                                  // for size_t
#include <fstream>                                   // for operator<<, etc
#include <stk_mesh/base/FieldBase.hpp>               // for FieldBase, etc
#include <stk_mesh/base/Types.hpp>                   // for EntityVector, etc
#include <vector>                                    // for vector
#include "stk_mesh/base/Bucket.hpp"                  // for Bucket
#include "stk_mesh/base/Entity.hpp"                  // for Entity
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/NgpFieldParallel.hpp"
#include "stk_mesh/base/Part.hpp"                    // for Part
#include "stk_mesh/base/Selector.hpp"                // for Selector
#include "stk_topology/topology.hpp"                 // for operator<<, etc
#include "stk_util/ngp/NgpSpaces.hpp"

namespace stk {
namespace debug {

std::string getFileName(const std::string &baseFileName, int numProcs, int localProc);

class DebugWriter
{
public:
    DebugWriter(const std::string& callingFile, int lineNumber, int numProcs, int localProc, const std::string basename = "selected")
    : DebugWriter(basename, callingFile, lineNumber, numProcs, localProc)
    {
    }

    virtual ~DebugWriter()
    {
        m_out.close();
    }

    template <typename Value>
    void write_to_file(const std::string &tagString, Value value)
    {
        m_out  << tagString << ": " << value << std::endl;
    }

    void write_selected_entities(const stk::mesh::BulkData& meshBulk, stk::mesh::Selector selector,
                                 const stk::mesh::FieldVector &fields)
    {
        if(meshBulk.in_modifiable_state())
        {
            m_out << "Mesh is in a modifiable state." << std::endl;
        }

        stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(meshBulk.mesh_meta_data().entity_rank_count());
        for(stk::mesh::EntityRank entityRank = stk::topology::NODE_RANK; entityRank < endRank; ++entityRank)
        {
            write_buckets(meshBulk, selector, entityRank);

            write_connectivity_and_field_data(meshBulk, selector, fields, entityRank);
        }
    }

    void write_meta_data(const stk::mesh::BulkData& meshBulk)
    {
        stk::mesh::impl::dump_all_meta_info(meshBulk.mesh_meta_data(), m_out);
    }

protected:
    DebugWriter(const std::string &baseFileName, const std::string& callingFile, int lineNumber, int numProcs, int localProc)
    : m_out(getFileName(baseFileName, numProcs, localProc).c_str(), std::ios_base::app)
    {
        m_out.precision(16);
        writeLocationInformation(m_out, callingFile, lineNumber);
    }

    virtual bool do_printing_for_part(stk::mesh::Part &part) const
    {
        return !stk::mesh::is_auto_declared_part(part);
    }

private:
    using DoubleField = stk::mesh::ConstFieldData<double, stk::ngp::HostSpace, stk::mesh::Layout::Auto>;
    using IntField    = stk::mesh::ConstFieldData<int, stk::ngp::HostSpace, stk::mesh::Layout::Auto>;

    void writeLocationInformation(std::ofstream &out, const std::string& callingFile, int lineNumber)
    {
        static int counter = 0;
        unsigned slashLocation = callingFile.find_last_of("/");
        std::string fileBaseName = callingFile.substr(slashLocation);
        out << "Called from: " << fileBaseName << ": " << lineNumber << ", counter: " << counter << std::endl;
        counter++;
    }

    void write_buckets(const stk::mesh::BulkData& meshBulk,
                       stk::mesh::Selector selector,
                       stk::mesh::EntityRank entityRank)
    {
        const stk::mesh::BucketVector &buckets = meshBulk.buckets(entityRank);

        m_out << "num " << entityRank << " buckets is " << buckets.size() << std::endl;
        for(size_t bucketIndex = 0; bucketIndex < buckets.size(); bucketIndex++)
        {
            const stk::mesh::Bucket &buck = *buckets[bucketIndex];
            if(selector(buck))
            {
                m_out << "bucket size = " << buck.size() << ", rank = " << buck.entity_rank() << ", parts: ";

                const stk::mesh::PartVector& parts = buck.supersets();
                for(size_t partIndex = 0; partIndex < parts.size(); ++partIndex)
                {
                    if(do_printing_for_part(*parts[partIndex]))
                    {
                        m_out << " " << parts[partIndex]->name();
                    }
                }
                m_out << std::endl;
                for(size_t entityIndex = 0; entityIndex < buck.size(); entityIndex++)
                {
                    m_out << meshBulk.entity_key(buck[entityIndex]) << std::endl;
                }
            }
        }
    }

    void write_connectivity_and_field_data(const stk::mesh::BulkData& meshBulk,
                                           stk::mesh::Selector selector,
                                           const stk::mesh::FieldVector &fields,
                                           stk::mesh::EntityRank entityRank)
    {

        stk::mesh::EntityVector entities;
        stk::mesh::get_selected_entities(selector, meshBulk.buckets(entityRank), entities);

        std::vector<DoubleField> doubleFields;
        std::vector<IntField> intFields;
        std::vector<std::string> doubleFieldNames;
        std::vector<std::string> intFieldNames;
        for (stk::mesh::FieldBase* field : fields)
        {
            if (is_double_field(field))
            {
                doubleFields.push_back(field->data<double, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Auto>());
                doubleFieldNames.push_back(field->name());
            } else if (is_int_field(field))
            {
                intFields.push_back(field->data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Auto>());
                intFieldNames.push_back(field->name());
            }
        }

        for(size_t entityIndex = 0; entityIndex < entities.size(); entityIndex++)
        {
            stk::mesh::Entity entity = entities[entityIndex];
            write_entity(meshBulk, entity);

            write_connectivity(meshBulk, entity);

            write_fields(doubleFields, doubleFieldNames, entity, entityRank);
            write_fields(intFields,    intFieldNames,    entity, entityRank);
        }
    }

    void write_entity(const stk::mesh::BulkData& stkMesh, stk::mesh::Entity entity)
    {
        m_out << stkMesh.entity_key(entity);
        m_out << " owned by proc " << stkMesh.parallel_owner_rank(entity);
    }

    void write_connectivity(const stk::mesh::BulkData& stkMesh, stk::mesh::Entity entity)
    {
        m_out << " downward connected to";
        for(stk::mesh::EntityRank connectedRank = stkMesh.entity_rank(entity);
                connectedRank >= stk::topology::NODE_RANK; --connectedRank)
        {
            const unsigned numConnected = stkMesh.num_connectivity(entity, connectedRank);
            const stk::mesh::Entity *connectedEntities = stkMesh.begin(entity, connectedRank);
            for(unsigned connI = 0; connI < numConnected; ++connI)
            {
                m_out << " " << stkMesh.entity_key(connectedEntities[connI]);
            }
        }
        m_out << std::endl;
    }

    template <typename FieldDataType>
    void write_fields(const std::vector<FieldDataType> &fields,
                      const std::vector<std::string>& fieldNames,
                      stk::mesh::Entity entity,
                      stk::mesh::EntityRank entityRank)
    {
        static_assert(std::is_same_v<FieldDataType, DoubleField> || std::is_same_v<FieldDataType, IntField>);

        for (size_t i=0; i < fields.size(); ++i)
        {
            if(fieldNames[i] != "nodeClosureCount" && fieldNames[i] != "faceClosureCount")
            {
                if(fields[i].entity_rank() == entityRank)
                {
                    auto entityValues = fields[i].entity_values(entity);
                    if(entityValues.num_scalars() > 0)
                    {
                        m_out << "    " << fieldNames[i];
                        write_field_values(entityValues);
                    }
                }
            }
        }
        m_out << std::endl;
    }


    bool is_double_field(stk::mesh::FieldBase *field)
    {
        return (field->data_traits().is_floating_point && field->data_traits().size_of == 8);
    }

    bool is_int_field(stk::mesh::FieldBase *field)
    {
        return (field->data_traits().is_integral && field->data_traits().size_of == 4);
    }

    template <typename EntityValues>
    void write_field_values(EntityValues values)
    {
        //for(int i = 0; i < numScalars; i++)
        for (stk::mesh::CopyIdx copy : values.copies())
        {
            for (stk::mesh::ComponentIdx comp : values.components())
            {
                m_out << " " << values(copy, comp);
            }
        }
    }

protected:
    std::ofstream m_out;
};

class DebugWriterWithInternalParts : public DebugWriter
{
public:
    DebugWriterWithInternalParts(const std::string& callingFile, int lineNumber, int numProcs, int localProc, const std::string basename = "internals")
    : DebugWriter(basename, callingFile, lineNumber, numProcs, localProc)
    {
    }

protected:
    virtual bool do_printing_for_part(stk::mesh::Part & /*part*/) const
    {
        return true;
    }
};

inline
std::string getFileName(const std::string &baseFileName, int numProcs, int localProc)
{
    std::ostringstream os;
    os << baseFileName << "." <<numProcs<<"."<<localProc;
    std::string filename = os.str();
    static bool first_time = true;
    if (first_time)
    {
        std::ofstream out(filename.c_str(), std::ios_base::out);
        out.close();
        first_time = false;
    }
    return filename;
}

}//namespace debug
}//namespace stk

#endif

