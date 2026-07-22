// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <Teuchos_XMLObject.hpp>
#include <Teuchos_FileInputSource.hpp>
#include <Teuchos_XMLPerfTestArchive.hpp>
#ifdef _WIN32
#include <winsock2.h>
#else
#include <unistd.h>
#endif

namespace Teuchos {

ValueTolerance::ValueTolerance() {
  value = 0;
  lower = 0;
  upper = 0;
  tolerance = 0;
  use_tolerance = true;
}

ValueTolerance::ValueTolerance(double val, double tol) {
  value = val;
  lower = 0;
  upper = 0;
  tolerance = tol;
  use_tolerance = true;
}

ValueTolerance::ValueTolerance(double val, double low, double up) {
  value = val;
  upper = up;
  lower = low;
  tolerance = 0;
  use_tolerance = false;
}

ValueTolerance::ValueTolerance(std::string str) {
  from_string(str);
}

bool ValueTolerance::operator ==(ValueTolerance& rhs) {
  return (value == rhs.value) &&
         (tolerance == rhs.tolerance) &&
         (lower == rhs.lower) &&
         (upper == rhs.upper) &&
         (use_tolerance == rhs.use_tolerance);
}

std::string ValueTolerance::as_string(){
  std::ostringstream strs;
  if(use_tolerance)
    strs << value << " , " << tolerance;
  else
    strs << value << " , " << lower << " , " << upper;
  return  strs.str();
}

void ValueTolerance::from_string(const std::string& valtol_str) {
  std::string value_str = valtol_str.substr(0,valtol_str.find(","));
  value = atof(value_str.c_str());
  std::string tol_str = valtol_str.substr(valtol_str.find(",")+1);
  if(tol_str.find(",")<=tol_str.length()) {
    use_tolerance = false;
    std::string lower_str = tol_str.substr(0,tol_str.find(","));
    lower = atof(lower_str.c_str());
    std::string upper_str = tol_str.substr(tol_str.find(",")+1);
    upper = atof(upper_str.c_str());
  } else {
    use_tolerance = true;
    tolerance = atof(tol_str.c_str());
  }
}

  XMLTestNode::XMLTestNode():XMLObject() {}

  XMLTestNode::XMLTestNode(const std::string &tag):XMLObject(tag) {}

  XMLTestNode::XMLTestNode(XMLObjectImplem *ptr):XMLObject(ptr) {}

  XMLTestNode::XMLTestNode(XMLObject obj):XMLObject(obj) {}

  void  XMLTestNode::addDouble (const std::string &name, double val) {
    addAttribute<double>(name,val);
  }

  void  XMLTestNode::addInt (const std::string &name, int val) {
    addAttribute<int>(name,val);
  }

  void  XMLTestNode::addBool (const std::string &name, bool val) {
    addAttribute<bool>(name,val);
  }

  void XMLTestNode::addValueTolerance(const std::string &name, ValueTolerance val){
    addAttribute<std::string>(name,val.as_string());
  }

  void  XMLTestNode::addString (const std::string &name, std::string val) {
    addAttribute<std::string>(name,val);
  }

  bool XMLTestNode::hasChild(const std::string &name) const {
    bool found = false;
    for(int i = 0; i < numChildren(); i++) {
      if(name.compare(XMLObject::getChild(i).getTag()) == 0) {
        found = true;
        i = numChildren();
      }
    }
    return found;
  }

  void XMLTestNode::appendContentLine(const size_t& i, const std::string &str) {
    ptr_->appendContentLine(i,str);
  }

  XMLTestNode XMLTestNode::getChild(const std::string &name) const {
    XMLTestNode child;
    for(int i = 0; i < numChildren(); i++) {
      if(name.compare(XMLObject::getChild(i).getTag()) == 0)
        child = XMLObject::getChild(i);
    }
    return child;
  }

  XMLTestNode XMLTestNode::getChild(const int &i) const {
    return XMLObject::getChild(i);
  }

  const XMLObject* XMLTestNode::xml_object() const {
    return (XMLObject*) this;
  }

  bool XMLTestNode::hasSameElements(XMLTestNode const & lhs) const {

    if((numChildren()!=lhs.numChildren()) ||
       (numContentLines()!= lhs.numContentLines()) ||
       (getTag().compare(lhs.getTag())!=0)) return false;

    for(int i = 0; i<numChildren(); i++) {
      const XMLTestNode child = XMLObject::getChild(i);
      if( (!lhs.hasChild(child.getTag())) ||
          (!child.hasSameElements(lhs.getChild(child.getTag()))) ) return false;
    }

    for(int i = 0; i<numContentLines(); i++)
      if(getContentLine(i).compare(lhs.getContentLine(i))!=0) return false;

    return true;
  }

XMLTestNode PerfTest_MachineConfig() {

  // Get CPUName, Number of Sockets, Number of Cores, Number of Hyperthreads
  std::string cpuname("Undefined");
  unsigned int threads = 0;
  unsigned int cores_per_socket = 0;
  unsigned int highest_socketid = 0;

  {
    std::ifstream cpuinfo("/proc/cpuinfo");
    std::string line;
    if((cpuinfo.rdstate()&cpuinfo.failbit)) std::cout<<"Failed to open filen\n";
    while (!cpuinfo.eof() && !(cpuinfo.rdstate()&cpuinfo.failbit)) {
      getline (cpuinfo,line);
      if (line.find("model name") < line.size()) {
        cpuname = line.substr(line.find(":")+2);
        threads++;
      }
      if (line.find("physical id") < line.size()) {
        unsigned int socketid = atoi(line.substr(line.find(":")+2).c_str());
        highest_socketid = highest_socketid>socketid?highest_socketid:socketid;
      }
      if (line.find("cpu cores") < line.size()) {
        cores_per_socket = atoi(line.substr(line.find(":")+2).c_str());
      }
    }
  }


  XMLTestNode machine_config("MachineConfiguration");

  machine_config.addString("Compiler", TEUCHOS_COMPILER_NAME);
  machine_config.addInt("Compiler_Version",  TEUCHOS_COMPILER_VERSION);
  machine_config.addString("CPU_Name", cpuname);
  machine_config.addInt("CPU_Sockets", highest_socketid+1);
  machine_config.addInt("CPU_Cores_Per_Socket", cores_per_socket);
  machine_config.addInt("CPU_Total_HyperThreads", threads);
  return machine_config;
}

PerfTestResult
PerfTest_CheckOrAdd_Test (XMLTestNode machine_config,
                          XMLTestNode new_test,
                          const std::string filename,
                          const std::string ext_hostname)
{
  XMLTestNode database;
  PerfTestResult return_value = PerfTestPassed;
  bool is_new_config = true;

  // Open Database File
  //
  // FIXME (mfh 09 Apr 2014) This actually opens the file twice.
  if (std::ifstream (filename.c_str ())) {
    database = FileInputSource (filename).getObject ();
  }

  // Get Current Hostname
  char hostname[256];
  memset (hostname, 0, 256);
  if (ext_hostname.empty ()) {
    gethostname (hostname, 255);
  } else {
    strncat (hostname, ext_hostname.c_str (), 255);
  }

  XMLTestNode new_test_entry = new_test.getChild ("TestEntry");

  if (database.isEmpty ()) {
    database = XMLTestNode ("PerfTests");
  }
  // Does hostname exist?
  if (database.hasChild (hostname)) {
    XMLTestNode machine = database.getChild (hostname);

    // Find matching machine configuration
    for (int i = 0; i < machine.numChildren (); ++i) {
      XMLTestNode configuration = machine.getChild (i);
      TEUCHOS_TEST_FOR_EXCEPTION(
        configuration.getTag ().compare ("Configuration") != 0,
        std::runtime_error, "Unexpected Tag \"" << configuration.getTag ()
        << "\"; only children with Tag = \"Configuration\" are allowed in a "
        "MachineEntry.");

      TEUCHOS_TEST_FOR_EXCEPTION(
        ! configuration.hasChild ("MachineConfiguration") ||
        ! configuration.hasChild ("Tests"),
        std::runtime_error,
        "A Configuration needs to have a child \"MachineConfiguration\" and a "
        "child \"Tests\".");

      XMLTestNode machine_configuration = configuration.getChild ("MachineConfiguration");
      XMLTestNode old_tests = configuration.getChild ("Tests");

      if (machine_configuration.hasSameElements (machine_config)) {
        is_new_config = false;

        // Find existing test with same tag as the new test
        if (old_tests.hasChild (new_test.getTag ())) {

          XMLTestNode old_test = old_tests.getChild (new_test.getTag ());

          int new_test_config = -1;
          for (int k = 0; k < old_test.numChildren (); ++k) {
            XMLTestNode old_test_entry = old_test.getChild (k);

            TEUCHOS_TEST_FOR_EXCEPTION(
              ! old_test_entry.hasChild ("TestConfiguration") ||
              ! new_test_entry.hasChild ("TestResults"),
              std::runtime_error, "A TestEntry needs to have a child "
              "\"TestConfiguration\" and a child \"TestResults\".");

            if (old_test_entry.getChild ("TestConfiguration").hasSameElements (new_test_entry.getChild ("TestConfiguration"))) {
              new_test_config = k;
            }
          }

          if (new_test_config < 0) {
            old_test.addChild (new_test_entry);
            return_value = PerfTestNewTestConfiguration;
          } else {
            bool deviation = false;
            XMLTestNode old_test_entry = old_test.getChild (new_test_config);
            XMLTestNode old_results = old_test_entry.getChild ("TestResults");
            XMLTestNode new_results = new_test_entry.getChild ("TestResults");

            // Compare all entries
            for (int old_r = 0; old_r < old_results.numChildren (); ++old_r) {
              XMLTestNode result_entry = old_results.getChild (old_r);

              // Finding entry with same name
              bool exists = new_results.hasChild (result_entry.getTag ());

              if (exists) {
                std::string oldv_str = result_entry.getContentLine (0);

                // If it is a time or result compare numeric values with tolerance
                if((result_entry.getTag().find("Time")==0) || (result_entry.getTag().find("Result")==0)) {
                  ValueTolerance old_valtol(oldv_str);
                  ValueTolerance new_valtol(new_results.getChild(result_entry.getTag()).getContentLine(0));

                  if(old_valtol.use_tolerance) {
                    double diff = old_valtol.value - new_valtol.value;
                    diff*=diff;

                    double normalization = old_valtol.value;
                    normalization*=normalization;

                    if(normalization==0?diff>0:diff/normalization>old_valtol.tolerance*old_valtol.tolerance) {
                      deviation = true;
                      std::cout << std::endl
                          << "DeviationA in Test: \"" << old_test.getTag()
                          << "\" for entry \"" <<  result_entry.getTag() << "\"" << std::endl;
                      std::cout << "  Existing Value: \"" << oldv_str << "\"" << std::endl;
                      std::cout << "  New Value:      \"" << new_results.getChild(result_entry.getTag()).getContentLine(0) << "\"" << std::endl << std::endl;
                    }
                  } else {
                    if( (old_valtol.lower>new_valtol.value) || (old_valtol.upper<new_valtol.value)) {
                      deviation = true;
                      std::cout << std::endl
                          << "DeviationB in Test: \"" << old_test.getTag()
                          << "\" for entry \"" <<  result_entry.getTag() << "\"" << std::endl;
                      std::cout << "  Existing Value: \"" << oldv_str << "\"" << std::endl;
                      std::cout << "  New Value:      \"" << new_results.getChild(result_entry.getTag()).getContentLine(0) << "\"" << std::endl << std::endl;
                    }
                  }
                } else {
                  // Compare exact match for every other type of entry
                  if(oldv_str.compare(new_results.getChild(result_entry.getTag()).getContentLine(0))!=0) {
                    deviation = true;
                    std::cout << std::endl
                        << "DeviationC in Test: \"" << old_test.getTag()
                        << "\" for entry \"" <<  result_entry.getTag() << "\"" << std::endl;
                    std::cout << "  Existing Value: \"" << oldv_str << "\"" << std::endl;
                    std::cout << "  New Value:      \"" << new_results.getChild(result_entry.getTag()).getContentLine(0) << "\"" << std::endl << std::endl;
                  }
                }
              }
              // An old value was not given in the new test: this is an error;
              if(!exists) {
                std::cout << "Error New test has same name as an existing one, but one of the old entries is missing." << std::endl;
                deviation = true;
              }
            }

            if(deviation) { return_value = PerfTestFailed; }
            else {
              // Did someone add new values to the test?
              if(new_results.numChildren()!=old_results.numChildren()) {
                for(int new_r = 0; new_r < new_results.numChildren() ; new_r++) {
                  if(!old_results.hasChild(new_results.getChild(new_r).getTag())) {
                    old_results.addChild(new_results.getChild(new_r));
                  }
                }

                return_value = PerfTestUpdatedTest;
              }
            }
          }
        } else { // End Test Exists
          // Add new test if no match was found
          old_tests.addChild(new_test);
          return_value = PerfTestNewTest;
        }
      } // End MachineConfiguration Exists
    } // End loop over MachineConfigurations

    // Did not find matching MachineConfiguration
    if(is_new_config) {
      XMLTestNode config("Configuration");
      config.addChild(machine_config);
      XMLTestNode tests("Tests");
      tests.addChild(new_test);

      config.addChild(tests);
      machine.addChild(config);

      return_value = PerfTestNewConfiguration;
    }
  } else { // Machine Entry does not exist
    XMLTestNode machine(hostname);

    XMLTestNode config("Configuration");
    config.addChild(machine_config);
    XMLTestNode tests("Tests");
    tests.addChild(new_test);
    config.addChild(tests);

    machine.addChild(config);

    database.addChild(machine);

    return_value = PerfTestNewMachine;
  }


  if(return_value>PerfTestPassed) {
    std::ofstream fout(filename.c_str());
    fout << database << std::endl;
  }

  return return_value;
}
}
