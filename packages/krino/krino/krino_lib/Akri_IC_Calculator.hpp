// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_IC_Calculator_h
#define Akri_IC_Calculator_h

namespace krino { class LevelSet; }

namespace krino {

// abstract class used to define an IC Calculator.

class IC_Calculator {
public:
  IC_Calculator() {}
  virtual ~IC_Calculator() {}

  // compute signed distance for all nodes with LevelSet defined
  virtual void compute_signed_distance(const LevelSet &ls) const = 0;
};

class IC_Binder : public IC_Calculator {
public:
  IC_Binder(const double interface_size, const double smooth_bridge_size,
            const double smooth_bridge_offset, const double other_ls_scale_factor,
            const int binder_type, const bool root_smooth_bridge) :
    my_interface_size(interface_size),
    my_smooth_bridge_size(smooth_bridge_size),
    my_smooth_bridge_offset(smooth_bridge_offset),
    my_other_ls_scale_factor(other_ls_scale_factor),
    my_binder_type(binder_type),
    my_root_smooth_bridge(root_smooth_bridge) {}
  virtual ~IC_Binder() {}

  // compute signed distance for all nodes with LevelSet defined
  void compute_signed_distance(const LevelSet &ls) const override;
private:
  double my_interface_size;
  double my_smooth_bridge_size;
  double my_smooth_bridge_offset;
  double my_other_ls_scale_factor;
  int my_binder_type;
  bool my_root_smooth_bridge;
};

} // namespace krino

#endif // Akri_IC_Calculator_h
