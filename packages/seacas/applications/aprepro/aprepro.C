// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <cstring>
#include <fstream>
#include <iostream>

#include "aprepro.h"

int main(int argc, char *argv[])
{
  SEAMS::Aprepro           aprepro;
  std::vector<std::string> input_files;

  bool quiet = false;

  int exit_status = EXIT_SUCCESS;

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
        try {
          double dval = std::stod(value);
          aprepro.add_variable(var, dval, true);
        }
        catch (std::exception &e) {
          // If cannot convert to double; make it a string variable...
          aprepro.add_variable(var, value, true); // Make it immutable
        }
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

      if (aprepro.ap_options.errors_fatal && aprepro.get_error_count() > 0) {
        exit_status = EXIT_FAILURE;
      }
      if ((aprepro.ap_options.errors_and_warnings_fatal) &&
          (aprepro.get_error_count() + aprepro.get_warning_count() > 0)) {
        exit_status = EXIT_FAILURE;
      }
    }
    catch (std::exception &e) {
      std::cerr << "Aprepro terminated due to exception: " << e.what() << '\n';
    }
  }
  else {
    std::fstream infile(input_files[0], std::fstream::in);
    if (!infile.good()) {
      std::cerr << "APREPRO: ERROR: Could not open file: " << input_files[0] << '\n'
                << "                Error Code: " << strerror(errno) << '\n';
      return 0;
    }

    // Read and parse a file.  The entire file will be parsed and
    // then the output can be obtained in an std::ostringstream via
    // Aprepro::parsing_results()
    bool writeResults = true;
    try {
      bool result = aprepro.parse_stream(infile, input_files[0]);

      if (aprepro.ap_options.errors_fatal && aprepro.get_error_count() > 0) {
        writeResults = false;
      }
      if ((aprepro.ap_options.errors_and_warnings_fatal) &&
          (aprepro.get_error_count() + aprepro.get_warning_count() > 0)) {
        writeResults = false;
      }

      if (writeResults) {
        if (result) {
          if (input_files.size() > 1) {
            std::ofstream ofile(input_files[1]);
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
      else {
        exit_status = EXIT_FAILURE;
        std::cerr << "There were " << aprepro.get_error_count() << " errors and "
                  << aprepro.get_warning_count() << " warnings."
                  << "\n";
        if (aprepro.ap_options.errors_and_warnings_fatal) {
          std::cerr << "Errors and warnings are fatal. No output has been written"
                    << "\n";
        }
        else if (aprepro.ap_options.errors_fatal) {
          std::cerr << "Errors are fatal. No output has been written."
                    << "\n";
        }
        else {
          std::cerr << "Neither errors nor warnings are fatal. "
                    << "If you see this message, then there is a bug in Aprepro. "
                    << "No output has been written."
                    << "\n";
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
  return exit_status;
}
