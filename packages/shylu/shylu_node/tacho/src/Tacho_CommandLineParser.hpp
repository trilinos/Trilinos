// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_COMMAND_LINE_PARSER_HPP__
#define __TACHO_COMMAND_LINE_PARSER_HPP__

/// \file Tacho_CommandLineParser.hpp
/// \brief Command line parser
/// \author Kyungjoo Kim (kyukim@sandia.gov)

// "std" includes
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <tuple>

namespace Tacho {

template <typename T> struct Option {
  std::string _desc;
  std::string _type;
  std::string _val;
  Option(std::string desc, T val);
  Option(std::string val) : _val(val){};
  T value() const;
};
template <> Option<bool>::Option(std::string desc, bool val) : _desc(desc), _type("bool"), _val(std::to_string(val)) {}
template <> Option<int>::Option(std::string desc, int val) : _desc(desc), _type("int"), _val(std::to_string(val)) {}
template <> Option<std::string>::Option(std::string desc, std::string val) : _desc(desc), _type("string"), _val(val) {}

struct CommandLineParser {
private:
  typedef std::tuple<std::string, std::string, std::string, void *> option_type;

  std::string _desc;
  std::map<std::string, option_type> _map;

public:
  CommandLineParser(){};
  CommandLineParser(const std::string desc) : _desc(desc) {}

  template <typename T> void set_option(const std::string opt, const std::string desc, T *val) {
    const auto in = Option<T>(desc, *val);
    _map[opt] = std::make_tuple(in._desc, in._type, in._val, val);
  }
  void print_option(const std::string argv0) {
    std::cout << "Usage: " << argv0 << " [options]\n";
    std::cout << "  options:\n";

    const std::string prefix = "--";
    for (auto it = _map.begin(); it != _map.end(); ++it) {
      auto key = prefix + it->first;
      auto val = it->second;
      std::cout << std::left << "  " << std::setw(30) << key << std::setw(10) << std::get<1>(val) << std::get<0>(val)
                << "\n";
      if (std::get<1>(val) == "bool") {
        if (it->first != "help" && it->first != "echo-command-line")
          std::cout << std::setw(42) << " "
                    << "(default: " << key << "=" << (std::get<2>(val) == "0" ? "false" : "true") << ")\n";
      } else {
        std::cout << std::setw(42) << " "
                  << "(default: " << key << "=" << std::get<2>(val) << ")\n";
      }
    }
    std::cout << "Description:\n";
    std::cout << "  " << _desc << "\n\n";
  }
  bool parse(int argc, char **argv) {
    bool help = false, echo = false;
    this->set_option<bool>("help", "Print this help message", &help);
    this->set_option<bool>("echo-command-line", "Echo the command-line but continue as normal", &echo);

    // check help
    for (int i = 1; i < argc; ++i) {
      std::string s = argv[i];
      if (s == "--help")
        help = true;
      if (s == "--echo-command-line")
        echo = true;
    }

    if (help) {
      print_option(argv[0]);
    } else {
      // parse
      for (int i = 1; i < argc; ++i) {
        std::string s = argv[i];
        for (auto it = _map.begin(); it != _map.end(); ++it) {
          const auto t = it->second;
          // find option starting with --
          if (s.find("--") != std::string::npos && s.length() != 2) {
            std::string desc = std::get<0>(t);
            std::string type = std::get<1>(t);

            size_t pos = s.find("=");
            if (pos != std::string::npos) {
              // --opt=val
              std::string key = s.substr(2, pos - 2);
              if (key == it->first) {
                std::string sval = s.substr(pos + 1, s.length());
                void *tval = std::get<3>(t);
                if (tval != NULL) {
                  if (type == "int") {
                    *((int *)tval) = atoi(sval.c_str());
                  } else if (type == "string") {
                    *((std::string *)tval) = sval;
                  } else if (type == "bool") {
                    *((bool *)tval) = (sval == "true");
                  } else {
                    std::cout << " int somethng wrong\n";
                  }
                  _map[it->first] = std::make_tuple(desc, type, sval, (void *)NULL);
                }
              }
            } else {
              // --opt (bool)
              std::string key = s.substr(2, s.length());
              if (key == it->first) {
                void *tval = std::get<3>(t);
                if (tval != NULL) {
                  if (type == "bool") {
                    *((bool *)tval) = true;
                  } else {
                    std::cout << " bool somethng wrong\n";
                  }
                  _map[it->first] = std::make_tuple(desc, type, "1", (void *)NULL);
                }
              }
            }
          }
        }
      }

      // print out unused options
      if (echo) {
        // echo command
        std::cout << "Echoing the command-line:\n\n";
        for (int i = 0; i < argc; ++i)
          std::cout << argv[i] << " ";
        std::cout << "\n\n";

        // used options
        std::cout << "Used options:\n\n";
        for (auto it = _map.begin(); it != _map.end(); ++it) {
          std::string key = it->first;
          std::string type = std::get<1>(it->second);
          std::string val = std::get<2>(it->second);
          bool used = std::get<3>(it->second) == NULL;
          if (used) {
            std::cout << "  ";
            if (type == "bool")
              std::cout << "--" << key << "\n";
            else
              std::cout << "--" << key << "=" << val << "\n";
          }
        }
        std::cout << "\n";

        // not used options
        std::cout << "Not used options:\n\n";
        for (auto it = _map.begin(); it != _map.end(); ++it) {
          std::string key = it->first;
          std::string type = std::get<1>(it->second);
          std::string val = std::get<2>(it->second);
          bool used = std::get<3>(it->second) == NULL;
          if (!used) {
            std::cout << "  ";
            if (type == "bool")
              std::cout << "--" << key << "\n";
            else
              std::cout << "--" << key << "=" << val << " (default)\n";
          }
        }
        std::cout << "\n";
      }
    }
    return help;
  }
};
} // namespace Tacho

#endif
