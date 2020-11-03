// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "adios/AdiosWrapper.h"

#include <algorithm>

namespace Ioad {

  AdiosWrapper::AdiosWrapper(MPI_Comm comm, const std::string &filename, bool is_input,
                             unsigned long rank, const Ioss::PropertyManager &properties)
      : adios2::ADIOS(comm), adios2::IO(IOInit(properties, is_input)), adios2::Engine(EngineInit(
                                                                           filename, is_input)),
        m_Rank(rank), m_Communicator(comm), m_OpenStep(false)
  {
  }

  adios2::IO AdiosWrapper::IOInit(const Ioss::PropertyManager &properties, bool is_input)
  {
    adios2::IO     bpio                 = this->ADIOS::DeclareIO("io");
    std::string    engine               = "BPFile";
    std::string    transport            = "File";
    std::string    library              = "POSIX";
    adios2::Params transport_parameters = {{"Library", "POSIX"}};
    // Set engine based on properties
    if (is_input && properties.exists("Engine_in")) {
      engine = properties.get("Engine_in").get_string();
    }
    if (!is_input && properties.exists("Engine_out")) {
      engine = properties.get("Engine_out").get_string();
    }
    if (engine == "") {
      // By default, use BPFile. This is the default in ADIOS, but this enforces that
      // if the default changes in ADIOS, it doesn't change here.
      engine = "bpfile";
    }
    std::transform(engine.begin(), engine.end(), engine.begin(), ::tolower);
    bpio.SetEngine(engine);
    // TODO: See if there is a way to avoid hardcoding this list
    // of engines (and possibly engines have multiple corresponding
    // strings).
    if (engine == "bpfile" || engine == "hdf5") {
      m_IsStreaming = false;
    }
    else {
      m_IsStreaming = true;
    }
    // Set transport based on properties
    if (properties.exists("Transport")) {
      transport = properties.get("Transport").get_string();
    }
    // if(transport == "File") // Leave transport parameters to default value.
    if (transport == "WAN") {
      transport_parameters = {{"Library", "ZMQ"}, {"IPAddress", "127.0.0.1"}};
    }
    else if (transport == "InSituMPI") {
      transport_parameters = {};
    }
    else if (transport == "SST") {
      transport_parameters = {{"MarshalMethod", "BP"}};
    }
    // Set transport parameters based on properties
    if (properties.exists("Library")) {
      transport_parameters["Library"] = properties.get("Library").get_string();
    }
    if (properties.exists("IPAddress")) {
      transport_parameters["IPAddress"] = properties.get("IPAddress").get_string();
    }
    if (properties.exists("MarshalMethod")) {
      transport_parameters["MarshalMethod"] = properties.get("MarshalMethod").get_string();
    }
    if (properties.exists("verbose")) {
      transport_parameters["verbose"] = properties.get("verbose").get_string();
    }
    bpio.AddTransport(transport, transport_parameters);
    // TODO: Add support for passing parameters, such as "num_threads".
    return bpio;
  }

  adios2::Engine AdiosWrapper::EngineInit(const std::string &filename, bool is_input)
  {
    adios2::Mode mode = adios2::Mode::Read;
    if (!is_input) {
      mode = adios2::Mode::Write;
    }
    return this->IO::Open(filename, mode);
  }

  AdiosWrapper::~AdiosWrapper()
  {
    EndStep();
    this->Engine::Close();
  }

  adios2::StepStatus AdiosWrapper::BeginStep()
  {
    if (!m_OpenStep) {
      adios2::StepStatus status = this->Engine::BeginStep();
      if (status == adios2::StepStatus::OK) {
        m_OpenStep = true;
      }
      return status;
    }
    // Either we  called `BeginStep()` successfully or we didn't need to. Either
    // way, there was no issue.
    return adios2::StepStatus::OK;
  }

  void AdiosWrapper::EndStep()
  {
    if (m_OpenStep) {
      this->Engine::EndStep();
      m_OpenStep = false;
    }
  }

  std::string AdiosWrapper::EncodeMetaVariable(const std::string &meta_name,
                                               const std::string &variable_name) const
  {
    if (!variable_name.empty()) {
      return variable_name + m_MetaSeparator + meta_name;
    }
    else {
      return meta_name;
    }
  }

  std::pair<std::string, std::string> AdiosWrapper::DecodeMetaName(std::string name) const
  {
    std::size_t pos = 0;
    std::string meta;
    pos = name.rfind(m_MetaSeparator);
    if (pos != std::string::npos && pos != name.size() - 1) {
      meta = name.substr(pos + m_MetaSeparator.size());
      name = name.substr(0, pos);
    }
    return std::make_pair(name, meta);
  }

} // namespace Ioad
