#ifndef PHALANX_EXCEPTIONS_HPP
#define PHALANX_EXCEPTIONS_HPP

#include <string>
#include <stdexcept>

namespace PHX {

  struct circular_dag_exception : public std::runtime_error {
    circular_dag_exception(const std::string message);
  };

  struct missing_evaluator_exception : public std::runtime_error {
    missing_evaluator_exception(const std::string message);
  };

  struct multiple_evaluator_for_field_exception : public std::runtime_error {
    multiple_evaluator_for_field_exception(const std::string message);
  };

}

#endif
