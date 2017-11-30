// Copyright(C) 2015-2017 National Technology & Engineering Solutions of
// Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above
//   copyright notice, this list of conditions and the following
//   disclaimer in the documentation and/or other materials provided
//   with the distribution.
//
// * Neither the name of NTESS nor the names of its
//   contributors may be used to endorse or promote products derived
//   from this software without specific prior written permission.
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

#include <fstream>
#include <iostream>

#include "aprepro.h"

int main(int argc, char *argv[])
{
  SEAMS::Aprepro           aprepro;
  std::vector<std::string> input_files;

  bool quiet = false;

  // Parse all options...
  for (int ai = 1; ai < argc; ++ai) {
    std::string arg = argv[ai];
    if ((arg[0] == '-' && arg[1] == 'q') || (arg[0] == '-' && arg[1] == '-' && arg[2] == 'q')) {
      quiet = true;
    }

    if (arg[0] == '-') { // Parse "--arg [val]" or "--arg=val" or "--arg"
      std::string val = ai + 1 < argc ? argv[ai + 1] : "";
      ai += aprepro.set_option(arg, val);
    }
    else if (arg.find('=') != std::string::npos) { // Parse var=value option.
      size_t      index = arg.find_first_of('=');
      std::string var   = arg.substr(0, index);
      std::string value = arg.substr(index + 1);
      if (value[0] == '\"' || value[0] == '\'') {
        value = value.substr(1, value.length() - 2);
        aprepro.add_variable(var, value, true); // Make it immutable
      }
      else {
        double dval = strtod(value.c_str(), nullptr);
        aprepro.add_variable(var, dval, true);
      }
    }
    else {
      input_files.emplace_back(argv[ai]);
    }
  }

  // Size of input_files should be either 0, 1, or 2:
  // 0 -- interactive, output to std::cout
  // 1 -- read from input_files[0], output to std::cout
  // 2 -- read from  input_files[0], output to input_files[1]

  if (input_files.empty()) {
    if (!quiet) {
      auto comment = aprepro.getsym("_C_")->value.svar;
      std::cout << comment << " Algebraic Preprocessor -- Aprepro, version " << aprepro.version()
                << "\n";
    }
    aprepro.ap_options.interactive = true;
    try {
      aprepro.parse_stream(std::cin, "standard input");
    }
    catch (std::exception &e) {
      std::cerr << "Aprepro terminated due to exception: " << e.what() << '\n';
    }
  }
  else {
    std::fstream infile(input_files[0].c_str());
    if (!infile.good()) {
      std::cerr << "APREPRO: ERROR: Could not open file: " << input_files[0] << '\n';
      return 0;
    }

    // Read and parse a file.  The entire file will be parsed and
    // then the output can be obtained in an std::ostringstream via
    // Aprepro::parsing_results()
    try {
      bool result = aprepro.parse_stream(infile, input_files[0]);

      if (result) {
        if (input_files.size() > 1) {
          std::ofstream ofile(input_files[1].c_str());
          if (!quiet) {
            auto comment = aprepro.getsym("_C_")->value.svar;
            ofile << comment << " Algebraic Preprocessor (Aprepro) version " << aprepro.version()
                  << "\n";
          }
          ofile << aprepro.parsing_results().str();
        }
        else {
          if (!quiet) {
            auto comment = aprepro.getsym("_C_")->value.svar;
            std::cout << comment << " Algebraic Preprocessor (Aprepro) version "
                      << aprepro.version() << "\n";
          }
          std::cout << aprepro.parsing_results().str();
        }
      }
    }
    catch (std::exception &e) {
      std::cerr << "Aprepro terminated due to exception: " << e.what() << '\n';
    }
  }
  if (aprepro.ap_options.debugging || aprepro.ap_options.dumpvars) {
    aprepro.dumpsym("variable", false);
  }
}
