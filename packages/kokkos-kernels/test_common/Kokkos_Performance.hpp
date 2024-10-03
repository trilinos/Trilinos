//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_PERFORMANCE_HPP
#define KOKKOS_PERFORMANCE_HPP

#include "Kokkos_Core.hpp"
#include <yaml-cpp/yaml.h>
#include <iostream>
#include <fstream>
#include <unistd.h>

namespace KokkosKernels {

// Note that this class was directly derived from code in the
// Teuchos_XMLPerfTestArchive.hpp class so there are many
// places which follow those methods exactly.

class Performance {
 public:
  /**
   * \brief Performance class manages the archive.
   *
   * This will generate a starting point for a machine configuration. Users
   * should add new entries via set_machine_config if necessary.
   * For example Kokkos users might want to provide the name of the user
   * Kokkos NodeType or Kokkos DeviceType. Contains information mostly
   * extracted from /proc/cpuinfo if possible.
   * Entries are:
   * - Compiler: The compiler name.
   * - Compiler_Version: A compiler version number.
   * - CPU_Name: The CPUs model name.
   * - CPU_Sockets: Number of CPU sockets in the system.
   * - CPU_Cores_Per_Socket: Number of CPU cores per socket.
   * - CPU_Total_HyperThreads: Total number of threads in a node.
   */
  Performance() : machine_configuration_node(get_machine_configuration()) {}

  /**
   * \brief set_config adds/changes a test config parameter
   *
   * \param name [in] The name used to identify the parameter in the archive.
   * \param val [in] The value assinged to the parameter.
   */
  template <typename T>
  void set_config(const std::string& name, const T& val) {
    test_configuration_node[name] = val;
  }

  /**
   * \brief set_result adds/changes a test result with a tolerance
   *
   * \param name [in] The name used to identify the result.
   * \param val [in] The recorded result which came from running the test.
   * \param tolerance [in] abs((new-old)/old)>tolerance triggers test failure.
   */
  void set_result(const std::string& name, double val, double tolerance) {
    validate_input_result_name(name);
    results_node[name] = Performance::Tolerance(val, tolerance).as_string();
  }

  /**
   * \brief set_result adds/changes a test result
   *
   * \param name [in] The name used to identify the result.
   * \param time [in] The recorded time which came from running the test.
   * \param tolerance_low [in] Lower bound of tolerance.
   * \param tolerance_high [in] Upper bound of tolerance.
   */
  void set_result(const std::string& name, double val, double tolerance_low, double tolerance_high) {
    validate_input_result_name(name);
    results_node[name] = Performance::Tolerance(val, tolerance_low, tolerance_high).as_string();
  }

  /**
   * \brief set_result adds/changes a test result for exact comparison
   *
   * \param name [in] The name used to identify the result.
   * \param val [in] The recorded result which came from running the test.
   */
  template <typename T>
  void set_result(const std::string& name, const T& val) {
    validate_input_result_name(name);
    results_node[mark_name_with_exact_code(name)] = val;
  }

  /**
   * \brief set_machine_config adds/changes a machine configuration
   *
   * \param name [in] The name used to identify the machine config parameter.
   * \param val [in] The setting for the machine config parameter.
   */
  template <typename T>
  void set_machine_config(const std::string& name, const T& val) {
    machine_configuration_node[name] = val;
  }

  /**
   * \brief Result codes after creating/comparing a test entry
   */
  enum Result { Failed, Passed, NewMachine, NewConfiguration, NewTest, NewTestConfiguration, UpdatedTest, Unknown };

  /**
   * \brief Processes the test and update the yaml archive
   *
   * This should be called after all necessary inserts, such as set_config,
   * set_time, and set_result.
   *
   * \param test_name [in] Named used to match test entries.
   * \param archive_name [in] The local yaml path to generate the archive.
   * \param host_name [in] An optional hostname to be used instead of
   *   the one provided by the OS.
   *
   * \return Whether a matching test is found, or if it was added to an
   *   archive.
   *
   * Will search for a matching machine name with matching machine
   * configuration and matching test configuration. If one is found the
   * result values will be compared, if not a new test entry is
   * generated and the result written back to the file.
   *
   * Here is the list of valid return values:
   *
   * - Failed: Matching configuration found, but results are
   *     deviating more than the allowed tolerance.
   * - Passed: Matching configuration found, and results are
   *     within tolerances.
   * - NewMachine: The test archive didn't contain an entry with
   *     the same machine name. A new entry was generated.
   * - NewConfiguration: No matching machine configuration was
   *     found. A new entry was generated.
   * - NewTest: No matching testname was found. A new entry was
   *     generated and added to the archive.
   * - NewTestConfiguration: A matching testname was found, but
   *     different parameters were used. A new entry was generated.
   * - UpdatedTest: A matching test was found but more result
   *     values were given then previously found. The entry is updated.
   *     This will only happen if all the old result values are present in
   *     the new ones, and are within their respective tolerances.
   */
  Result run(const std::string& archive_name, const std::string& test_name, const std::string& host_name = "") const;

  /**
   * \brief print_archive will std::cout the yaml archive for inspection.
   *
   * \param archive_name [in] The local yaml path to print.
   */
  static void print_archive(const std::string& archive_name);

  /**
   * \brief erase_archive will delete the archive.
   *
   * \param archive_name [in] The local yaml path to print.
   */
  static void erase_archive(const std::string& archive_name);

 private:
  typedef YAML::Node node_t;

  /* Tolerance is an internal helper struct.
   * The tolerance can be either expressed as relative or through an upper and
   * lower bound. This is now private to the Performance class and follows the
   * original Teuchos implementation.
   */
  struct Tolerance {
    double value;
    double lower;
    double upper;
    double tolerance;
    bool use_tolerance;
    Tolerance();
    Tolerance(double val, double tol);
    Tolerance(double val, double low, double up);
    Tolerance(std::string str);
    bool operator==(const Tolerance& rhs);
    std::string as_string();
    void from_string(const std::string& valtol_str);
  };

  // Set up the machine configuration - users can modify the default setup
  node_t get_machine_configuration() const;

  // Compares two nodes and determines if they are have the same
  // members but does not compare the values of those members.
  bool hasSameElements(const node_t& a, const node_t& b) const;

  // does string include the exact key character '*'
  bool string_includes_exact_code(const std::string& name) const;

  // append the '*' key character to the name
  std::string mark_name_with_exact_code(const std::string& name) const;

  // make sure the input name doesn't have the character '*'
  void validate_input_result_name(const std::string& name) const;

  // private data members
  node_t machine_configuration_node;  // stores machine config settings
  node_t test_configuration_node;     // stores test config settings
  node_t results_node;                // stores result settings
};

Performance::node_t Performance::get_machine_configuration() const {
  // Get CPUName, Number of Sockets, Number of Cores, Number of Hyperthreads
  std::string cpuname("Undefined");
  unsigned int threads          = 0;
  unsigned int cores_per_socket = 0;
  unsigned int highest_socketid = 0;

  std::ifstream cpuinfo("/proc/cpuinfo");
  std::string line;
  if ((cpuinfo.rdstate() & cpuinfo.failbit)) {
#ifndef __clang__  // TODO decide how to best handle this generically
    std::cout << "Failed to open /proc/cpuinfo\n";
#endif
  }
  while (!cpuinfo.eof() && !(cpuinfo.rdstate() & cpuinfo.failbit)) {
    getline(cpuinfo, line);
    if (line.find("model name") < line.size()) {
      cpuname = line.substr(line.find(":") + 2);
      threads++;
    }
    if (line.find("physical id") < line.size()) {
      unsigned int socketid = atoi(line.substr(line.find(":") + 2).c_str());
      highest_socketid      = highest_socketid > socketid ? highest_socketid : socketid;
    }
    if (line.find("cpu cores") < line.size()) {
      cores_per_socket = atoi(line.substr(line.find(":") + 2).c_str());
    }
  }

  std::string compiler_name = "Unknown";
  int compiler_version      = 0;

#if defined __clang__
  compiler_name    = "Clang";
  compiler_version = __clang_major__ * 100 + __clang_minor__ * 10 + __clang_patchlevel__;
#endif

#if defined __GNUC__ && !defined KOKKOS_COMPILER_NAME && !defined __clang__
  compiler_name    = "Gnu GCC";
  compiler_version = __GNUC__ * 100 + __GNUC_MINOR__ * 10 + __GNUC_PATCHLEVEL__;
#endif

#if defined __PGIC__ && !defined KOKKOS_COMPILER_NAME
  compiler_name    = "PGI C++";
  compiler_version = __PGIC__ * 100 + __PGIC_MINOR__ * 10 + __PGIC_PATCHLEVEL__;
#endif

  node_t machine_config;
  machine_config["Compiler"]               = compiler_name;
  machine_config["Compiler_Version"]       = compiler_version;
  machine_config["CPU_Name"]               = cpuname;
  machine_config["CPU_Sockets"]            = highest_socketid + 1;
  machine_config["CPU_Cores_Per_Socket"]   = cores_per_socket;
  machine_config["CPU_Total_HyperThreads"] = threads;
  return machine_config;
}

Performance::Result Performance::run(const std::string& archive_name, const std::string& test_name,
                                     const std::string& host_name) const {
  // These are abitrary category names used in the yaml
  const std::string test_configuration_string    = "TestConfiguration";
  const std::string test_results_string          = "TestResults";
  const std::string machine_configuration_string = "MachineConfiguration";
  const std::string tests_string                 = "Tests";

  // Now create the test entry - combincation of configuration and times/results
  Performance::node_t new_test_entry;  // the entry will have two bits added below
  new_test_entry[test_configuration_string] = test_configuration_node;
  new_test_entry[test_results_string]       = results_node;

  // Run the archiver which will either add the results, or compare them to
  // prior results if they already exist. This method will open the yaml,
  // import everything, do appropriate comparisons, then write out a new yaml.
  Performance::node_t database;

  Result return_value = Performance::Passed;
  bool is_new_config  = true;

  // Open YAML File whhich stores test database
  if (std::ifstream(archive_name)) {
    database = YAML::LoadFile(archive_name);
  }

  // Get host_setting which is read by default or set to optional host_name
  char host_setting[256];
  memset(host_setting, 0, 256);
  if (host_name.empty()) {
    gethostname(host_setting, 255);
  } else {
    strncat(host_setting, host_name.c_str(), 255);
  }

  // Does host_setting exist?
  if (database[host_setting]) {
    Performance::node_t machine = database[host_setting];

    // Find matching machine configuration
    for (size_t machine_index = 0; machine_index < machine.size(); ++machine_index) {
      Performance::node_t configuration = machine[machine_index];
      if (!configuration[machine_configuration_string] || !configuration[tests_string]) {
        throw std::logic_error(
            "Configuration must has child MachineConfiguration and a child "
            "\"Tests\".");
      }

      Performance::node_t machine_configuration = configuration[machine_configuration_string];
      Performance::node_t old_tests             = configuration[tests_string];
      if (hasSameElements(machine_configuration, machine_configuration_node)) {
        is_new_config = false;

        // Find existing test with same name as the new test
        if (old_tests[test_name]) {
          Performance::node_t old_test_array = old_tests[test_name];
          int match_test_index               = -1;
          for (size_t entry_index = 0; entry_index < old_test_array.size(); ++entry_index) {
            Performance::node_t old_test_entry = old_test_array[entry_index];
            if (hasSameElements(old_test_entry[test_configuration_string], new_test_entry[test_configuration_string])) {
              match_test_index = static_cast<int>(entry_index);
            }
          }
          if (match_test_index == -1) {
            database[host_setting][machine_index][tests_string][test_name].push_back(new_test_entry);
            return_value = Performance::NewTestConfiguration;
          } else {
            bool deviation                     = false;
            Performance::node_t old_test_entry = old_test_array[match_test_index];
            Performance::node_t old_results    = old_test_entry[test_results_string];
            Performance::node_t new_results    = new_test_entry[test_results_string];
            // Compare all entries
            for (YAML::const_iterator old_r = old_results.begin(); old_r != old_results.end(); ++old_r) {
              Performance::node_t result_entry = old_r->second;
              // Finding entry with same name
              std::string result_name = old_r->first.Scalar();
              bool exists             = new_results[result_name];
              if (exists) {
                std::string oldv_str      = old_r->second.Scalar();
                std::string old_test_name = test_name;
                std::ostringstream new_result_entry_name_stream;
                new_result_entry_name_stream << new_results[result_name];
                std::string new_result_data = new_result_entry_name_stream.str();

                // based on name does result use tolerance?
                // if it has the '*' key character appended it means it's an
                // exact
                if (!string_includes_exact_code(result_name)) {
                  Performance::Tolerance old_valtol(oldv_str);
                  Performance::Tolerance new_valtol(new_results[result_name].Scalar());
                  if (old_valtol.use_tolerance) {
                    double diff = old_valtol.value - new_valtol.value;
                    diff *= diff;

                    double normalization = old_valtol.value;
                    normalization *= normalization;
                    if (normalization == 0 ? diff > 0
                                           : diff / normalization > old_valtol.tolerance * old_valtol.tolerance) {
                      deviation = true;
                      std::cout << std::endl
                                << "  DeviationA in Test: \"" << old_test_name << "\" for entry \"" << result_name
                                << "\"" << std::endl;
                      std::cout << "    Existing Value: \"" << oldv_str << "\"" << std::endl;
                      std::cout << "    New Value:      \"" << new_result_data << "\"" << std::endl << std::endl;
                    }
                  } else {
                    if ((old_valtol.lower > new_valtol.value) || (old_valtol.upper < new_valtol.value)) {
                      deviation = true;
                      std::cout << std::endl
                                << "  DeviationB in Test: \"" << old_test_name << "\" for entry \"" << result_name
                                << "\"" << std::endl;
                      std::cout << "    Existing Value: \"" << oldv_str << "\"" << std::endl;
                      std::cout << "    New Value:      \"" << new_result_data << "\"" << std::endl << std::endl;
                    }
                  }
                } else {
                  // Compare exact match for every other type of entry
                  if (oldv_str.compare(new_result_data) != 0) {
                    deviation = true;
                    std::cout << std::endl
                              << "  DeviationC in Test: \"" << old_test_name << "\" for entry \"" << result_name << "\""
                              << std::endl;
                    std::cout << "    Existing Value: \"" << oldv_str << "\"" << std::endl;
                    std::cout << "    New Value:      \"" << new_result_data << "\"" << std::endl << std::endl;
                  }
                }
              }
              // An old value was not given in the new test: this is an error;
              if (!exists) {
                std::cout << "  Error New test has same name as an existing "
                             "one, but one of the old entries is missing."
                          << std::endl;
                deviation = true;
              }
            }
            if (deviation) {
              return_value = Performance::Failed;
            } else {
              // Did someone add new values to the test?
              if (new_results.size() != old_results.size()) {
                for (YAML::const_iterator new_r = new_results.begin(); new_r != new_results.end(); ++new_r) {
                  if (!old_results[new_r->first.Scalar()]) {
                    old_results[new_r->first.Scalar()] = (new_r->second);
                  }
                }
                return_value = Performance::UpdatedTest;
              }
            }
          }
        } else {  // End Test Exists
          // Add new test if no match was found
          database[host_setting][machine_index][tests_string][test_name].push_back(new_test_entry);
          return_value = Performance::NewTest;
        }
      }  // End MachineConfiguration Exists
    }    // End loop over MachineConfigurations

    // Did not find matching MachineConfiguration
    if (is_new_config) {
      Performance::node_t machine_entry;
      machine_entry[machine_configuration_string] = machine_configuration_node;
      machine_entry[tests_string][test_name].push_back(new_test_entry);
      database[host_setting].push_back(machine_entry);
      return_value = Performance::NewConfiguration;
    }
  } else {  // Machine Entry does not exist
    Performance::node_t machine_entry;
    machine_entry[machine_configuration_string] = machine_configuration_node;
    machine_entry[tests_string][test_name].push_back(new_test_entry);
    database[host_setting].push_back(machine_entry);
    return_value = Performance::NewMachine;
  }

  if (return_value > Performance::Passed) {
    // write the actual database out
    std::ofstream fout(archive_name.c_str());
    fout << database << std::endl;
  }
  return return_value;
}

bool Performance::hasSameElements(const Performance::node_t& a, const Performance::node_t& b) const {
  if (a.size() != b.size()) {
    return false;
  }

  for (YAML::const_iterator i = a.begin(); i != a.end(); ++i) {
    std::string cat_name = i->first.Scalar();
    // validate we can find this cat in b
    if (!b[cat_name]) {
      return false;
    }
    Performance::node_t sub_a = i->second;
    Performance::node_t sub_b = b[cat_name];
    if (sub_a.Scalar() != sub_b.Scalar()) {
      return false;
    }

    if (!hasSameElements(sub_a, sub_b)) {
      return false;
    }
  }

  return true;
}

bool Performance::string_includes_exact_code(const std::string& name) const {
  return (name.length() != 0 && name[name.length() - 1] == '*');
}

std::string Performance::mark_name_with_exact_code(const std::string& name) const { return name + "*"; }

void Performance::validate_input_result_name(const std::string& name) const {
  if (string_includes_exact_code(name)) {
    throw std::logic_error("The name " + name +
                           " ends with the * key character"
                           " which is reserved for internal use to signify an "
                           "exact value not using"
                           " the tolerance settings.");
  }
}

void Performance::print_archive(const std::string& archiveName) {
  std::cout << YAML::LoadFile(archiveName) << std::endl;
}

void Performance::erase_archive(const std::string& yamlArchive) { std::ofstream(yamlArchive) << std::endl; }

Performance::Tolerance::Tolerance() {
  value         = 0;
  lower         = 0;
  upper         = 0;
  tolerance     = 0;
  use_tolerance = true;
}

Performance::Tolerance::Tolerance(double val, double tol) {
  value         = val;
  lower         = 0;
  upper         = 0;
  tolerance     = tol;
  use_tolerance = true;
}

Performance::Tolerance::Tolerance(double val, double low, double up) {
  value         = val;
  upper         = up;
  lower         = low;
  tolerance     = 0;
  use_tolerance = false;
}

Performance::Tolerance::Tolerance(std::string str) { from_string(str); }

bool Performance::Tolerance::operator==(const Tolerance& rhs) {
  return (value == rhs.value) && (tolerance == rhs.tolerance) && (lower == rhs.lower) && (upper == rhs.upper) &&
         (use_tolerance == rhs.use_tolerance);
}

std::string Performance::Tolerance::as_string() {
  std::ostringstream strs;
  if (use_tolerance)
    strs << value << " , " << tolerance;
  else
    strs << value << " , " << lower << " , " << upper;
  return strs.str();
}

void Performance::Tolerance::from_string(const std::string& valtol_str) {
  std::string value_str = valtol_str.substr(0, valtol_str.find(","));
  value                 = atof(value_str.c_str());
  std::string tol_str   = valtol_str.substr(valtol_str.find(",") + 1);
  if (tol_str.find(",") <= tol_str.length()) {
    use_tolerance         = false;
    std::string lower_str = tol_str.substr(0, tol_str.find(","));
    lower                 = atof(lower_str.c_str());
    std::string upper_str = tol_str.substr(tol_str.find(",") + 1);
    upper                 = atof(upper_str.c_str());
  } else {
    use_tolerance = true;
    tolerance     = atof(tol_str.c_str());
  }
}

}  // namespace KokkosKernels

#endif  // KOKKOS_PERFORMANCE_HPP
