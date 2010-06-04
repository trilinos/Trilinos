/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <heartbeat/Iohb_Layout.h>
#include <string>
#include <iomanip>

namespace Iohb {
  Layout::Layout(bool show_labels, int precision)
    : showLabels(show_labels), precision_(precision), count_(0), layout_()
  {}

  Layout::~Layout()
  {}

  std::ostream& operator<<(std::ostream& o, Layout& lo)
  {
    o << lo.layout_.str();
    return o;
  }

  void Layout::add_literal(const std::string& label)
  {
    layout_ << label;
  }

  void Layout::add(const std::string& name, double value)
  {
    if (count_++ > 0) {
      layout_ << ", ";
    }

    if (showLabels && name != "") {
      layout_ << name;
      layout_ << "=";
    }
    layout_.setf(std::ios::scientific);
    layout_.setf(std::ios::showpoint);
    layout_ << std::setprecision(precision_) << value;
  }

  void Layout::add(const std::string& name, int value)
  {
    if (count_++ > 0) {
      layout_ << ", ";
    }

    if (showLabels && name != "") {
      layout_ << name;
      layout_ << "=";
    }
    layout_ << value;
  }

  void Layout::add(const std::string& name, long value)
  {
    if (count_++ > 0) {
      layout_ << ", ";
    }

    if (showLabels && name != "") {
      layout_ << name;
      layout_ << "=";
    }
  layout_ << value;
  }

  void Layout::add(const std::string& name, const std::string &value)
  {
    if (count_++ > 0) {
      layout_ << ", ";
    }

    if (showLabels && name != "") {
      layout_ << name;
      layout_ << "=";
    }
    layout_ << value;
  }

  void Layout::add(const std::string& name, std::vector<double> &value)
  {
    if (value.size() == 1) {
      add(name, value[0]);
    } else {
      if (count_++ > 0) {
	layout_ << ", ";
      }

      if (showLabels && name != "") {
	layout_ << name;
	layout_ << "=";
      }
      layout_.setf(std::ios::scientific);
      layout_.setf(std::ios::showpoint);
      for (std::vector<double>::size_type i=0; i < value.size(); i++) {
	layout_ << std::setprecision(precision_) << value[i];
	if (i < value.size()-1)
	  layout_ << ", ";
      }

    }
  }

  void Layout::add(const std::string& name, std::vector<int> &value)
  {
    if (value.size() == 1) {
      add(name, value[0]);
    } else {
      if (count_++ > 0) {
	layout_ << ", ";
      }

      if (showLabels && name != "") {
	layout_ << name;
	layout_ << "=";
      }
      for (std::vector<int>::size_type i=0; i < value.size(); i++) {
	layout_ << value[i];
	if (i < value.size()-1)
	  layout_ << ", ";
      }
    }
  }

  void Layout::add(const std::string& name, std::vector<long> &value)
  {
    if (value.size() == 1) {
      add(name, value[0]);
    } else {
      if (count_++ > 0) {
	layout_ << ", ";
      }

      if (showLabels && name != "") {
	layout_ << name;
	layout_ << "=";
      }
      for (std::vector<long>::size_type i=0; i < value.size(); i++) {
	layout_ << value[i];
	if (i < value.size()-1)
	  layout_ << ", ";
      }
    }
  }

  void Layout::add(const std::string& name, std::vector<std::string> &value)
  {
    if (value.size() == 1)
      add(name, value[0]);
  }
}
