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

#ifndef STK_MESH_DIAGNOSTIC_OBSERVER_HPP
#define STK_MESH_DIAGNOSTIC_OBSERVER_HPP

#include <stk_mesh/base/ModificationObserver.hpp>
#include <unordered_set>

namespace stk { namespace mesh {

class BulkData;

enum MeshDiagnosticFlag
{
    NO_RULE   = 0     ,
    RULE_1    = 1L<<1 ,
    RULE_2    = 1L<<2 ,
    RULE_3    = 1L<<3 ,
    ALL_RULES = RULE_1 | RULE_2 | RULE_3,
    SOLO_SIDES = 1L<<4
};

class MeshDiagnosticObserver : public stk::mesh::ModificationObserver
{
public:
    MeshDiagnosticObserver(stk::mesh::BulkData& bulkData)
      : ModificationObserver(ModificationObserverPriority::STK_INTERNAL_LOW_PRIORITY),
        m_initialized(false),
        m_throwOnError(false),
        m_bulkData(bulkData),
        m_diagnosticFlag(NO_RULE) {}

    MeshDiagnosticObserver(stk::mesh::BulkData& bulkData, stk::mesh::MeshDiagnosticFlag diagnosticFlag)
      : ModificationObserver(ModificationObserverPriority::STK_INTERNAL_LOW_PRIORITY),
        m_initialized(false),
        m_throwOnError(false),
        m_bulkData(bulkData),
        m_diagnosticFlag(diagnosticFlag) {}

    virtual void finished_modification_end_notification() override;
    void enable_rule(stk::mesh::MeshDiagnosticFlag flag);
    unsigned get_number_of_errors() const {return m_errorDatabase.size();}
    void throw_if_errors_exist() {m_throwOnError = true;}
private:
    bool m_initialized;
    bool m_throwOnError;
    stk::mesh::BulkData &m_bulkData;
    std::unordered_set<std::string> m_errorDatabase;
    stk::mesh::MeshDiagnosticFlag m_diagnosticFlag;

    void gather_new_errors(std::ofstream & out, const std::vector<std::string> & errorList);
};

} }

#endif
