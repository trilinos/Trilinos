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
//

#ifndef STK_STK_IO_STK_IO_DYNAMICTOPOLOGY_HPP_
#define STK_STK_IO_STK_IO_DYNAMICTOPOLOGY_HPP_

#include <Ioss_DynamicTopology.h>
#include <Ioss_DynamicTopologyObserver.h>
#include <Ioss_DynamicTopologyBroker.h>
#include <stk_mesh/base/Types.hpp>

namespace stk {
  namespace io {

    enum class FileOption {
      NO_DYNAMIC_TOPOLOGY_FILE_CONTROL = 0,
      USE_DYNAMIC_TOPOLOGY_MULTI_FILE  = 1,
      USE_DYNAMIC_TOPOLOGY_GROUP_FILE  = 2,
    };

    Ioss::FileControlOption get_ioss_file_control_option(FileOption fileOption);

    class StkMeshIoBroker;

    class DynamicTopologyObserver : public Ioss::DynamicTopologyObserver
    {
    public:
      DynamicTopologyObserver(StkMeshIoBroker& stkIoBroker, size_t index,
                              const Ioss::FileControlOption fileControlOption_);

      virtual ~DynamicTopologyObserver() {}

      void define_model() override;

      void write_model() override;

      void define_transient() override;

      Ioss::FileControlOption get_control_option() const override;

      bool needs_new_output_file() const override;

      bool get_output_refined_mesh() const { return m_outputRefinedMesh; }

      void set_output_refined_mesh(bool flag) { m_outputRefinedMesh = flag; }

      void initialize_region() override;

    private:
      DynamicTopologyObserver();

      StkMeshIoBroker& m_stkIoBroker;
      size_t m_outputFileIndex;
      Ioss::FileControlOption m_fileControlOption;

      bool m_outputRefinedMesh{true};

      unsigned count_output_field_variables(stk::mesh::EntityRank rank) const;
    };

  }
}


#endif /* STK_STK_IO_STK_IO_DYNAMICTOPOLOGY_HPP_ */
