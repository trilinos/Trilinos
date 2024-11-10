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

#ifndef SKIN_BOUNDARY_ERROR_REPORTER_HPP
#define SKIN_BOUNDARY_ERROR_REPORTER_HPP

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/SideSetEntry.hpp>
#include <stk_mesh/base/SideSetEntryCompare.hpp>

namespace stk
{
namespace mesh
{
namespace impl
{

class SkinBoundaryErrorReporter
{
public:
    SkinBoundaryErrorReporter(std::ostream &stream, const stk::mesh::BulkData &bulkData)
    : m_stream(stream), m_bulkData(bulkData), m_cmp(bulkData) {}

    void add_entry(stk::mesh::Entity face, const SideSetEntry &sideSetEntry)
    {
        FaceToSideSetMap::iterator iter = m_faceToSideSetMap.find(face);
        if(iter == m_faceToSideSetMap.end())
        {
            add_new_entry(face, sideSetEntry);
        }
        else
        {
            iter->second.insert(sideSetEntry);
        }
    }

    void report(const stk::mesh::EntityVector & skinnedSides, const stk::mesh::EntityVector &sidesetSides, const Part& skinnedPart)
    {
        std::ostringstream os;

        os << "" << std::endl;
        os << "Skin report for skinned part named: " << skinnedPart.name() << " on processor " << m_bulkData.parallel_rank() << std::endl;
        os << "-----------------------------------" << std::endl;

        report_difference_from_skinned_part_to_sideset(skinnedSides, sidesetSides, os);
        report_difference_from_sideset_to_skinned_part(skinnedSides, sidesetSides, os);

        m_stream << os.str();
    }

private:
    SkinBoundaryErrorReporter();

    void add_new_entry(stk::mesh::Entity face, const SideSetEntry &sideSetEntry)
    {
        std::set<SideSetEntry, SideSetEntryLess> setEntry(m_cmp);
        setEntry.insert(sideSetEntry);
        m_faceToSideSetMap.insert( std::pair< stk::mesh::Entity, std::set<SideSetEntry, SideSetEntryLess> > (face, setEntry));
    }

    stk::mesh::EntityVector get_entity_vector_difference(const stk::mesh::EntityVector & A, const stk::mesh::EntityVector &B)
    {
        stk::mesh::EntityVector difference(A.size() + B.size());
        stk::mesh::EntityVector::iterator it;

        it = std::set_difference(A.begin(), A.end(), B.begin(), B.end(), difference.begin());
        difference.resize(it - difference.begin());

        return difference;
    }

    void report_difference_from_skinned_part_to_sideset(const stk::mesh::EntityVector & skinnedSides, const stk::mesh::EntityVector &sidesetSides, std::ostringstream &os)
    {
        stk::mesh::EntityVector diffFromPartToSideSet = get_entity_vector_difference(skinnedSides, sidesetSides);

        os << "    List of sides in skinned part but not in skinned sideset" << std::endl;
        for(size_t i = 0; i < diffFromPartToSideSet.size(); ++i)
        {
            os << "        (" << i << ") " << m_bulkData.identifier(diffFromPartToSideSet[i]) << std::endl;
        }
        os << "    -----------------------------------" << std::endl;
    }

    void report_difference_from_sideset_to_skinned_part(const stk::mesh::EntityVector & skinnedSides, const stk::mesh::EntityVector &sidesetSides, std::ostringstream &os)
    {
        stk::mesh::EntityVector diffFromSideSetToPart = get_entity_vector_difference(sidesetSides, skinnedSides);

        os << "    List of sides in skinned sideset but not in skinned part" << std::endl;
        for(size_t i = 0; i < diffFromSideSetToPart.size(); ++i)
        {
            const stk::mesh::Bucket* bptr = m_bulkData.bucket_ptr(diffFromSideSetToPart[i]);
            os << "        (" << i << ") side-entity " << m_bulkData.identifier(diffFromSideSetToPart[i]) << " " << bptr->topology() << " parts: { ";
            for(const stk::mesh::Part* part : bptr->supersets()) {
              os << part->name() << " ";
            }
            os << "}" << std::endl;
            report_sideset_info(diffFromSideSetToPart[i], os);
        }
        os << "    -----------------------------------" << std::endl;
    }

    void report_sideset_entry_info(const std::set<SideSetEntry, SideSetEntryLess> &setList, std::ostringstream &os)
    {
        int count = 0;
        for(const SideSetEntry & entry : setList)
        {
            os << "            Sideset Info[" << count << "] = (" << m_bulkData.identifier(entry.element) << "," << m_bulkData.bucket(entry.element).topology() << "," << (uint32_t)entry.side << ")"<< std::endl;
            count++;
        }
    }

    void report_sideset_info(stk::mesh::Entity face, std::ostringstream &os)
    {
        FaceToSideSetMap::iterator iter = m_faceToSideSetMap.find(face);
        if(iter != m_faceToSideSetMap.end())
            report_sideset_entry_info(iter->second, os);
    }

    typedef std::map< stk::mesh::Entity, std::set<SideSetEntry, SideSetEntryLess> > FaceToSideSetMap;

    std::ostream &m_stream;
    const stk::mesh::BulkData &m_bulkData;
    FaceToSideSetMap m_faceToSideSetMap;
    SideSetEntryLess m_cmp;
};

} } }

#endif

