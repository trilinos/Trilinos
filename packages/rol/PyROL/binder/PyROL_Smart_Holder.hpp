#pragma once
#include <pybind11/smart_holder.h>
#include <ROL_Vector.hpp>
#include <ROL_Vector_SimOpt.hpp>
#include <ROL_Objective.hpp>
#include <ROL_Objective_SimOpt.hpp>
#include <ROL_Reduced_Objective_SimOpt.hpp>
#include <ROL_Constraint.hpp>
#include <ROL_Constraint_SimOpt.hpp>

PYBIND11_SMART_HOLDER_TYPE_CASTERS(ROL::Vector<double>)
PYBIND11_SMART_HOLDER_TYPE_CASTERS(ROL::Vector_SimOpt<double>)
PYBIND11_SMART_HOLDER_TYPE_CASTERS(ROL::Objective<double>)
PYBIND11_SMART_HOLDER_TYPE_CASTERS(ROL::Objective_SimOpt<double>)
PYBIND11_SMART_HOLDER_TYPE_CASTERS(ROL::Reduced_Objective_SimOpt<double>)
PYBIND11_SMART_HOLDER_TYPE_CASTERS(ROL::Constraint<double>)
PYBIND11_SMART_HOLDER_TYPE_CASTERS(ROL::Constraint_SimOpt<double>)
