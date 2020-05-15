// Copyright(C) 1999-2010 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
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
//     * Neither the name of NTESS nor the names of its
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

#ifndef IOSS_Ioad_AdiosWrapper_h
#define IOSS_Ioad_AdiosWrapper_h

#include <Ioss_PropertyManager.h>
#include <adios2.h>
#include <string>

namespace Ioad {

  class AdiosWrapper : private adios2::ADIOS, private adios2::IO, private adios2::Engine
  {
  public:
    AdiosWrapper(MPI_Comm communicator, const std::string &filename, bool is_input,
                 unsigned long rank, const Ioss::PropertyManager &properties);
    AdiosWrapper(AdiosWrapper &&wrapper);
    ~AdiosWrapper();
    adios2::StepStatus BeginStep();
    void               EndStep();
    template <typename T>
    void DefineMetaVariable(const std::string &meta_name, const std::string &variable_name = "");

    template <typename T>
    void DefineVariable(const std::string &name, const adios2::Dims &shape = adios2::Dims(),
                        const adios2::Dims &start        = adios2::Dims(),
                        const adios2::Dims &count        = adios2::Dims(),
                        const bool          constantDims = false);
    template <typename T> void Put(const std::string &name, const T *value);

    template <typename T> void DefineAttribute(const std::string &name, const T &value);

    template <typename T> void InquireAndPut(const std::string &name, const T *value);
    template <typename T>
    T GetAttribute(const std::string &attribute_name, bool ignore_missing = false,
                   T default_value = T());

    template <typename T> void GetSync(adios2::Variable<T> var, T *data);

    template <typename T> void GetSync(std::string var_name, T *data);

    template <typename T> void GetSync(adios2::Variable<T> var, T &data);
    template <typename T> void GetSync(std::string var_name, T &data);

    template <typename T>
    void PutMetaVariable(const std::string &meta_name, T value,
                         const std::string &variable_name = "");
    template <typename T>
    T GetMetaVariable(const std::string &meta_name, const std::string &variable_name = "");
    std::pair<std::string, std::string> DecodeMetaName(std::string name) const;
    std::string                         EncodeMetaVariable(const std::string &meta_name,
                                                           const std::string &variable_name = "") const;

    bool IsStreaming() const { return m_IsStreaming; };

    using adios2::Engine::AllStepsBlocksInfo;
    using adios2::IO::AvailableVariables;
    using adios2::IO::InquireAttribute;
    using adios2::IO::InquireVariable;

  private:
    adios2::IO     IOInit(const Ioss::PropertyManager &properties, bool is_input);
    adios2::Engine EngineInit(const std::string &filename, bool is_input);

    const std::string m_MetaSeparator{"::"};

    const int      m_Rank;
    const MPI_Comm m_Communicator;

    bool m_OpenStep;
    bool m_IsStreaming;
    int  count_real_begin = 0;
  };

} // end of namespace Ioad

#include "adios/AdiosWrapper.hpp"

#endif
