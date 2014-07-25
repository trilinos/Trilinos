// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_XMLPERFTESTARCHIVE_HPP
#define TEUCHOS_XMLPERFTESTARCHIVE_HPP

/// \file Teuchos_XMLPerfTestArchive.hpp
/// \brief Tools for an XML-based performance test archive

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_FileInputSource.hpp>
#include <Teuchos_XMLObject.hpp>
#include <sstream>

//----------------------------------------------------------------------------
//-------- Identify Compiler Version -----------------------------------------
//----------------------------------------------------------------------------

#if defined __ECC || defined __ICC || defined __INTEL_COMPILER
  #define TEUCHOS_COMPILER_NAME "Intel C++"
  #if defined __ICC
    #define TEUCHOS_COMPILER_VERSION __ICC
  #else
    #if defined __INTEL_COMPILER
      #define TEUCHOS_COMPILER_VERSION __INTEL_COMPILER
    #else
      #define TEUCHOS_COMPILER_VERSION __ECC
    #endif
  #endif
  #define TEUCHOS_COMPILER_INTEL 1
#endif

#if defined __IBMC__ || defined __IBMCPP__
  #define TEUCHOS_COMPILER_NAME "IBM C++"
  #if defined __IBMC__
    #define TEUCHOS_COMPILER_VERSION __IBMC__
  #else
    #define TEUCHOS_COMPILER_VERSION __IBMCPP__
  #endif
  #define TEUCHOS_COMPILER_IBM 1
#endif

#if defined __APPLE_CC__
   /* Apple uses GNU as compiler */
  #define TEUCHOS_COMPILER_APPLECC 1
#endif

#if defined __clang__
  #define TEUCHOS_COMPILER_NAME "Clang"
  #define TEUCHOS_COMPILER_VERSION __clang_major__*100+__clang_minor__*10+__clang_patchlevel__
  #define TEUCHOS_COMPILER_CLANG 1
#endif

#if defined __GNUC__ && !defined TEUCHOS_COMPILER_NAME && !defined __clang__
  #define TEUCHOS_COMPILER_NAME "Gnu GCC"
  #define TEUCHOS_COMPILER_VERSION __GNUC__*100+__GNUC_MINOR__*10+__GNUC_PATCHLEVEL__
  #define TEUCHOS_COMPILER_GCC 1
#endif

#if defined __PGIC__ && !defined TEUCHOS_COMPILER_NAME
  #define TEUCHOS_COMPILER_NAME "PGI C++"
  #define TEUCHOS_COMPILER_VERSION __PGIC__*100+__PGIC_MINOR__*10+__PGIC_PATCHLEVEL__
  #define TEUCHOS_COMPILER_PGI 1
#endif

#if defined __NVCC__
  #define TEUCHOS_DEVICE_COMPILER_NAME "NVIDIA NVCC"
  #define TEUCHOS_DEVICE_COMPILER_VERSION __NVCC__
#endif

#if !defined TEUCHOS_COMPILER_NAME
  #define TEUCHOS_COMPILER_NAME "Unknown compiler"
#endif

#if !defined TEUCHOS_COMPILER_VERSION
  #define TEUCHOS_COMPILER_VERSION 0
#endif

#if !defined TEUCHOS_DEVICE_COMPILER_NAME
  #define TEUCHOS_DEVICE_COMPILER_NAME TEUCHOS_COMPILER_NAME
#endif

#if !defined TEUCHOS_DEVICE_COMPILER_VERSION
  #define TEUCHOS_DEVICE_COMPILER_VERSION TEUCHOS_COMPILER_VERSION
#endif

namespace Teuchos {
  /**
   * \brief ValueTolerance is a struct to keep a tuple of value and a tolerance.
   * The tolerance can be either expressed as a relative or through an upper and
   * lower bound.
   */
struct ValueTolerance {
  double value;
  double lower;
  double upper;
  double tolerance;
  bool use_tolerance;

  ValueTolerance();
  ValueTolerance(double val, double tol);
  ValueTolerance(double val, double low, double up);

  ValueTolerance(std::string str);

  bool operator ==(ValueTolerance& rhs);

  std::string as_string();
  void from_string(const std::string& valtol_str);
};


/**
 * \class XMLTestNode
 * \brief Subclass of XMLObject used by the performance archive.
 *
 * This subclass of XMLObject generates an XML list in a style more
 * suitable for a performance test archive. It also provides a number
 * of convenience functions helpful for working with a test archive.
 */
class XMLTestNode : public XMLObject {
public:
  XMLTestNode();
  XMLTestNode(const std::string &tag);
  XMLTestNode(XMLObjectImplem *ptr);
  XMLTestNode(XMLObject obj);
  void addDouble (const std::string& name, double val);
  void addInt (const std::string& name, int val);
  void addBool (const std::string& name, bool val);
  void addValueTolerance(const std::string& name, ValueTolerance val);
  void addString (const std::string& name, std::string val);

  template<class T>
  void addAttribute (const std::string& name, T val) {
    for (size_t i = 0; i < name.length (); i++) {
      if (name[i] == ' ') {
        return;
      }
    }
    std::ostringstream strs;
    strs << val;
    XMLTestNode entry (name);
    entry.addContent (strs.str ());
    XMLObject::addChild (entry);
  }

  bool hasChild(const std::string &name) const;

  void appendContentLine(const size_t& i, const std::string &str);

  XMLTestNode getChild(const std::string &name) const;

  XMLTestNode getChild(const int &i) const;

  const XMLObject* xml_object() const;

  bool hasSameElements(XMLTestNode const & lhs) const;
};

/**
 * \brief PerfTest_MachineConfig generates a basic machine configuration XMLTestNode.
 *
 * \details The function provides a starting point for a machine configuration. Users
 * should add new entries to the returned XMLTestNode to provide test relevant machine
 * configuration entries. For example Kokkos users might want to provide the name of the
 * user Kokkos NodeType or Kokkos DeviceType. The returned config contains information
 * mostly extracted from /proc/cpuinfo if possible. On non unix systems most values
 * will be unknown. Entries are:
 * - Compiler: The compiler name.
 * - Compiler_Version: A compiler version number.
 * - CPU_Name: The CPUs model name.
 * - CPU_Sockets: Number of CPU sockets in the system.
 * - CPU_Cores_Per_Socket: Number of CPU cores per socket.
 * - CPU_Total_HyperThreads: Total number of threads in a node.
 */
XMLTestNode PerfTest_MachineConfig();

/**
 * \brief ReturnValues for PerfTest_CheckOrAdd_Test
 */
enum PerfTestResult {PerfTestFailed, PerfTestPassed,
                     PerfTestNewMachine, PerfTestNewConfiguration,
                     PerfTestNewTest, PerfTestNewTestConfiguration,
                     PerfTestUpdatedTest};

/**
 * \brief Check whether a test is present and match an existing test
 *   in an archive.
 *
 * This function consumes a machine configuration XMLTestNode and a
 * test entry XMLTestNode.  It will attempt to read from an existing
 * file containing a test archive, or generate a new one.  Optionally
 * a hostname override can be provided, which is for example useful
 * when running on clusters, where the cluster name should be used for
 * the test entries instead of the compute node name.
 * PerfTest_CheckOrAdd_Test will go through the test archive and
 * search for a matching machine name with matching machine
 * configuration and matching test configuration. If one is found the
 * result values will be compared, if not a new test entry is
 * generated and the result written back to the file.
 *
 * \param machine_config [in] An XMLTestNode describing the machine
 *   configuration.
 * \param new_test [in] An XMLTestNode describing the test.
 * \param filename [in] The name of a file containing a performance
 *   test archive.
 * \param ext_hostname [in] An optional hostname to be used instead of
 *   the one provided by the OS.
 *
 * \return Whether a matching test is found, or if it was added to an
 *   archive.
 *
 * Here is the list of valid return values:
 *
 * - PerfTestFailed: Matching configuration found, but results are
 *   deviating more than the allowed tolerance.
 * - PerfTestPassed: Matching configuration found, and results are
 *   within tolerances.
 * - PerfTestNewMachine: The test archive didn't contain an entry with
 *   the same machine name. A new entry was generated.
 * - PerfTestNewConfiguration: No matching machine configuration was
 *   found. A new entry was generated.
 * - PerfTestNewTest: No matching testname was found. A new entry was
 *   generated.
 * - PerfTestNewTestConfiguration: A matching testname was found, but
 *   different parameters were used. A new entry was generated.
 * - PerfTestUpdatedTest: A matching test was found but more result
 *   values were given then previously found. The entry is updated.
 *   This will only happen if all the old result values are present in
 *   the new ones, and are within their respective tolerances.
 */
PerfTestResult
PerfTest_CheckOrAdd_Test (XMLTestNode machine_config,
                          XMLTestNode new_test,
                          const std::string filename,
                          const std::string ext_hostname = std::string ());

} // namespace Teuchos

#endif
