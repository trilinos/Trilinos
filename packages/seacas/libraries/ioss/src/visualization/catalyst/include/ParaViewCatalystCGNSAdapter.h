// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef __PARAVIEW_CATALYST_CGNS_ADAPTER_H
#define __PARAVIEW_CATALYST_CGNS_ADAPTER_H

#include "CatalystParserInterface.h"
#include <stdint.h>
#include <string>
#include <vector>

class ParaViewCatalystCGNSAdapterBase
{
public:
  ParaViewCatalystCGNSAdapterBase(){};
  virtual ~ParaViewCatalystCGNSAdapterBase(){};
  virtual void CreateNewPipeline(const char *catalyst_python_filename,
                                 const char *catalyst_sierra_block_json) = 0;
  virtual void CleanupCatalyst() = 0;
  virtual void PerformCoProcessing() = 0;
  virtual void SetTimeData(double currentTime, int timeStep) = 0;
  virtual void CreateBase(int base_id,
                          const std::string& base_name) = 0;
  virtual void AddStructuredZoneData(int base_id,
                                     int zone_id,
                                     const std::string& zone_name,
                                     const std::string& data_name,
                                     int ni,
                                     int nj,
                                     int nk,
                                     int comp_count,
                                     bool is_cell_field,
                                     char field_suffix_separator,
                                     double* data,
                                     int size) = 0;
  virtual int parseFile(const std::string &filepath,
                        CatalystParserInterface::parse_info &pinfo) = 0;
  virtual int parseString(const std::string &s, CatalystParserInterface::parse_info &pinfo) = 0;
};

typedef ParaViewCatalystCGNSAdapterBase *(*ParaViewCatalystCGNSAdapterBaseSignature)();

extern "C" {
ParaViewCatalystCGNSAdapterBase *ParaViewCatalystCGNSAdapterCreateInstance();
}

class ParaViewCatalystCGNSAdapterImplementation;

class ParaViewCatalystCGNSAdapter : public ParaViewCatalystCGNSAdapterBase
{
public:
  ParaViewCatalystCGNSAdapter() {}
  virtual ~ParaViewCatalystCGNSAdapter() {}
  virtual void CreateNewPipeline(const char *catalyst_python_filename,
                                 const char *catalyst_sierra_block_json);
  virtual void CleanupCatalyst();
  virtual void PerformCoProcessing();
  virtual void SetTimeData(double currentTime, int timeStep);
  virtual void CreateBase(int base_id,
                          const std::string& base_name);
  virtual void AddStructuredZoneData(int base_id,
                                     int zone_id,
                                     const std::string& zone_name,
                                     const std::string& data_name,
                                     int ni,
                                     int nj,
                                     int nk,
                                     int comp_count,
                                     bool is_cell_field,
                                     char field_suffix_separator,
                                     double* data,
                                     int size);
  virtual int parseFile(const std::string &filepath,
                        CatalystParserInterface::parse_info &pinfo);
  virtual int parseString(const std::string &s, CatalystParserInterface::parse_info &pinfo);
};

#endif /* __PARAVIEW_CATALYST_CGNS_ADAPTER_H */
