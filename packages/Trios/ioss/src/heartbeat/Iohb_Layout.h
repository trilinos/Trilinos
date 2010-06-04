/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef SIERRA_Iohb_Layout_h
#define SIERRA_Iohb_Layout_h

#include <iostream>
#include <sstream>
#include <vector>
#include <string>

namespace Iohb {
  class Layout
    {
    public:
      explicit Layout(bool show_labels=true, int precision=5);
      ~Layout();

      friend std::ostream& operator<<(std::ostream&, Layout&);

      void add_literal(const std::string& label);
      void add(const std::string& name, double value);
      void add(const std::string& name, int value);
      void add(const std::string& name, long value);
      void add(const std::string& name, const std::string &value);

      void add(const std::string& name, std::vector<double> &value);
      void add(const std::string& name, std::vector<int> &value);
      void add(const std::string& name, std::vector<long> &value);
      void add(const std::string& name, std::vector<std::string> &value);

    private:
      Layout(const Layout&); // do not implement
      Layout& operator=(const Layout&); // do not implement

      bool showLabels;
      int precision_;

      int count_; // Number of fields on current line...
      std::ostringstream layout_;
    };
}

#endif // SIERRA_Iohb_Layout_h
