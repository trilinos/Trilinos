/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef __PARAVIEW_CATALYST_IOSS_ADAPTER_H
#define __PARAVIEW_CATALYST_IOSS_ADAPTER_H

#include "CatalystParserInterface.h"
#include <stdint.h>
#include <string>
#include <vector>

// Base class needed for Ioss's dynamic library
// registration.

class ParaViewCatalystIossAdapterBase
{
public:
  ParaViewCatalystIossAdapterBase(){};
  virtual ~ParaViewCatalystIossAdapterBase(){};
  virtual std::string getName() const { return "ParaViewCatalystIossAdapterBase"; }
  virtual void        DeletePipeline(const char *results_output_filename) = 0;
  virtual void        CleanupCatalyst()                                   = 0;
  virtual void        CreateNewPipeline(
             const char *catalyst_python_filename, const char *catalyst_sierra_block_json,
             const char *catalyst_sierra_separator_character, const char *catalyst_sierra_input_deck_name,
             int UnderscoreVectors, int ApplyDisplacements, const char *restart_tag, int enable_logging,
             int debug_level, const char *results_output_filename, const char *catalyst_output_directory,
             std::vector<std::string> &catalyst_sierra_data)                                   = 0;
  virtual void PerformCoProcessing(const char *              results_output_filename,
                                   std::vector<int> &        error_and_warning_codes,
                                   std::vector<std::string> &error_and_warning_messages)       = 0;
  virtual void SetTimeData(double currentTime, int timeStep,
                           const char *results_output_filename)                                = 0;
  virtual void CreateGlobalVariable(std::vector<std::string> &component_names, const double *data,
                                    const char *results_output_filename)                       = 0;
  virtual void CreateGlobalVariable(std::vector<std::string> &component_names, const int *data,
                                    const char *results_output_filename)                       = 0;
  virtual void InitializeGlobalPoints(int num_points, int dimension, const double *data,
                                      const char *results_output_filename)                     = 0;
  virtual void InitializeElementBlocks(const std::vector<int> &element_block_id_list,
                                       const char *            results_output_filename)                    = 0;
  virtual void CreateElementBlock(const char *elem_block_name, int elem_block_id,
                                  const std::string &elem_type, int nodes_per_elem, int num_elem,
                                  const int64_t *global_elem_ids, int *connectivity,
                                  const char *results_output_filename)                         = 0;
  virtual void CreateElementBlock(const char *elem_block_name, int elem_block_id,
                                  const std::string &elem_type, int nodes_per_elem, int num_elem,
                                  const int64_t *global_elem_ids, int64_t *connectivity,
                                  const char *results_output_filename)                         = 0;
  virtual void CreateNodeSet(const char *node_set_name, int node_set_id, int num_ids,
                             const int *data, const char *results_output_filename)             = 0;
  virtual void CreateNodeSet(const char *node_set_name, int node_set_id, int num_ids,
                             const int64_t *data, const char *results_output_filename)         = 0;
  virtual void CreateSideSet(/*const char* side_set_name,*/
                             const char *ss_owner_name, int side_set_id, int num_ids,
                             const int *element_ids, const int *face_ids,
                             const char *results_output_filename) = 0;
  virtual void CreateSideSet(/*const char* side_set_name,*/
                             const char *ss_owner_name, int side_set_id, int num_ids,
                             const int64_t *element_ids, const int64_t *face_ids,
                             const char *results_output_filename) = 0;
  virtual void CreateElementVariable(std::vector<std::string> &component_names, int elem_block_id,
                                     const double *data, const char *results_output_filename)  = 0;
  virtual void CreateElementVariable(std::vector<std::string> &component_names, int elem_block_id,
                                     const int *data, const char *results_output_filename)     = 0;
  virtual void CreateElementVariable(std::vector<std::string> &component_names, int elem_block_id,
                                     const int64_t *data, const char *results_output_filename) = 0;
  virtual void CreateNodalVariable(std::vector<std::string> &component_names, const double *data,
                                   const char *results_output_filename)                        = 0;
  virtual void CreateNodalVariable(std::vector<std::string> &component_names, const int *data,
                                   const char *results_output_filename)                        = 0;
  virtual void CreateNodalVariable(std::vector<std::string> &component_names, const int64_t *data,
                                   const char *results_output_filename)                        = 0;
  virtual void ReleaseMemory(const char *results_output_filename)                              = 0;
  virtual void logMemoryUsageAndTakeTimerReading(const char *results_output_filename)          = 0;
  virtual int  parseFile(const std::string &                  filepath,
                         CatalystParserInterface::parse_info &pinfo)                           = 0;

  virtual int parseString(const std::string &s, CatalystParserInterface::parse_info &pinfo) = 0;
};

typedef ParaViewCatalystIossAdapterBase *(*ParaViewCatalystIossAdapterBaseSignature)();

extern "C" {
ParaViewCatalystIossAdapterBase *ParaViewCatalystIossAdapterCreateInstance();
}

// ParaViewCatalystIossAdapter is a link between Ioss's mesh
// data structure and a ParaView vtkMultiBlockDataSet data structure.

class ParaViewCatalystIossAdapterImplementation;

class ParaViewCatalystIossAdapter : public ParaViewCatalystIossAdapterBase
{
public:
  ParaViewCatalystIossAdapter() {}
  virtual ~ParaViewCatalystIossAdapter() {}

  // Description:
  // Deletes pipeline with name results_output_filename and any associated
  // logging data.
  virtual void DeletePipeline(const char *results_output_filename);

  // Description:
  // Cleanup ParaView Catalyst and free resources.
  virtual void CleanupCatalyst();

  // Description:
  // Initializes ParaView Catalyst to perform in-situ co-processing
  // with the Python file catalyst_python_filename.  This method can
  // be called multiple times with different co-processing Python scripts.
  // If initialization fails, co-processing will not occur in any other
  // methods on this class.
  // Additional arguments:
  //   UnderscoreVectors - joined vector variable names end in an underscore.
  //   ApplyDisplacements - a nodal variable named DISPL or displ is applied to
  //                        the mesh node coordinates each time-step.
  //   restart_tag - if not empty, contains the current restart iteration string, ie s0001
  //   enable_logging - turn on logging in the adapter. Default is off.
  //   debug_level - enable catalyst debug output 0, 1, 2. Default is 0.
  //   results_output_filename - filename associated with the Ioss results output block.
  //   catalyst_output_directory - name of the output directory for storing Catalyst output.
  //                               Default is CatalystOutput.
  //   catalyst_sierra_data - string data vector for development and debugging.
  virtual void CreateNewPipeline(
      const char *catalyst_python_filename, const char *catalyst_sierra_block_json,
      const char *catalyst_sierra_separator_character, const char *catalyst_sierra_input_deck_name,
      int UnderscoreVectors, int ApplyDisplacements, const char *restart_tag, int enable_logging,
      int debug_level, const char *results_output_filename, const char *catalyst_output_directory,
      std::vector<std::string> &catalyst_sierra_data);

  // Description:
  // Calls the ParaView Catalyst pipeline to run co-processing for this time iteration.
  virtual void PerformCoProcessing(const char *              results_output_filename,
                                   std::vector<int> &        error_and_warning_codes,
                                   std::vector<std::string> &error_and_warning_messages);

  // Description:
  // Sets time data for this ParaView Catalyst co-processing iteration.
  // currentTime is the current Ioss simulation time and timeStep is
  // the current time iteration count.
  virtual void SetTimeData(double currentTime, int timeStep, const char *results_output_filename);

  // Description:
  // Creates global vtkPoints array.
  virtual void InitializeGlobalPoints(int num_points, int dimension, const double *data,
                                      const char *results_output_filename);

  // Description:
  // Creates empty element blocks for the ids contained in element_block_id_list
  virtual void InitializeElementBlocks(const std::vector<int> &element_block_id_list,
                                       const char *            results_output_filename);

  // Description:
  // Clears all nodal and element variables from the vtkMultiBlockDataSet.
  // Clears the global vtkPoints.
  virtual void ReleaseMemory(const char *results_output_filename);

  // Description:
  // Creates a global variable on the vtkMultiBlockDataSet.
  virtual void CreateGlobalVariable(std::vector<std::string> &component_names, const double *data,
                                    const char *results_output_filename);
  virtual void CreateGlobalVariable(std::vector<std::string> &component_names, const int *data,
                                    const char *results_output_filename);

  // Description:
  // Creates an element block on the vtkMultiBlockDataSet.
  virtual void CreateElementBlock(const char *elem_block_name, int elem_block_id,
                                  const std::string &elem_type, int nodes_per_elem, int num_elem,
                                  const int64_t *global_elem_ids, int *connectivity,
                                  const char *results_output_filename);
  virtual void CreateElementBlock(const char *elem_block_name, int elem_block_id,
                                  const std::string &elem_type, int nodes_per_elem, int num_elem,
                                  const int64_t *global_elem_ids, int64_t *connectivity,
                                  const char *results_output_filename);

  // Description:
  // Creates a node set on the vtkMultiBlockDataSet.
  virtual void CreateNodeSet(const char *node_set_name, int node_set_id, int num_ids,
                             const int *data, const char *results_output_filename);
  virtual void CreateNodeSet(const char *node_set_name, int node_set_id, int num_ids,
                             const int64_t *data, const char *results_output_filename);

  // Description:
  // Creates a side set (side block) on the vtkMultiBlockDataSet.
  void CreateSideSet(/*const char* side_set_name,*/
                     const char *ss_owner_name, int side_set_id, int num_ids,
                     const int *element_ids, const int *face_ids,
                     const char *results_output_filename);
  void CreateSideSet(/*const char* side_set_name,*/
                     const char *ss_owner_name, int side_set_id, int num_ids,
                     const int64_t *element_ids, const int64_t *face_ids,
                     const char *results_output_filename);

  // Description:
  // Creates an element variable on the vtkMultiBlockDataSet.
  virtual void CreateElementVariable(std::vector<std::string> &component_names, int elem_block_id,
                                     const double *data, const char *results_output_filename);
  virtual void CreateElementVariable(std::vector<std::string> &component_names, int elem_block_id,
                                     const int *data, const char *results_output_filename);
  virtual void CreateElementVariable(std::vector<std::string> &component_names, int elem_block_id,
                                     const int64_t *data, const char *results_output_filename);

  // Description:
  // Creates a nodal variable on the vtkMultiBlockDataSet.
  virtual void CreateNodalVariable(std::vector<std::string> &component_names, const double *data,
                                   const char *results_output_filename);
  virtual void CreateNodalVariable(std::vector<std::string> &component_names, const int *data,
                                   const char *results_output_filename);
  virtual void CreateNodalVariable(std::vector<std::string> &component_names, const int64_t *data,
                                   const char *results_output_filename);

  // Description:
  // Collects memory usage information from all processors and
  // writes the min, max, and mean to the log file.  Also writes the
  // min, max, and mean of the elapsed time since this method was
  // last called.
  virtual void logMemoryUsageAndTakeTimerReading(const char *results_output_filename);

  int parseFile(const std::string &filepath, CatalystParserInterface::parse_info &pinfo);

  int parseString(const std::string &s, CatalystParserInterface::parse_info &pinfo);

  virtual std::string getName() const { return "ParaViewCatalystIossAdapter"; }
};

#endif /* __PARAVIEW_CATALYST_IOSS_ADAPTER_H */
