template< typename Scalar , class DeviceType >
struct compute_acceleration_velocity_displacement;

template<typename Scalar>
struct compute_acceleration_velocity_displacement<Scalar, KOKKOS_MACRO_DEVICE>
{
	typedef KOKKOS_MACRO_DEVICE     device_type ;
	typedef typename Kokkos::MDArrayView<Scalar,device_type> array_type ;

  compute_acceleration_velocity_displacement(
      const array_type     & arg_internal_force,
      const array_type     & arg_nodal_mass,
      const array_type     & arg_acceleration,
      const array_type     & arg_velocity,
      const array_type     & arg_displacement,
      const Scalar           arg_dt,
      const int              arg_current_state,
      const int              arg_previous_state
      )
    : internal_force(arg_internal_force)
    , nodal_mass(arg_nodal_mass)
    , acceleration(arg_acceleration)
    , velocity(arg_velocity)
    , displacement(arg_displacement)
    , dt(arg_dt)
    , current_state(arg_current_state)
    , previous_state(arg_previous_state)
  {}


  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( int inode ) const
  {
    Scalar v_new[3];

    acceleration(inode,0,current_state) = internal_force(inode,0) / nodal_mass(inode);
    acceleration(inode,1,current_state) = internal_force(inode,1) / nodal_mass(inode);
    acceleration(inode,2,current_state) = internal_force(inode,2) / nodal_mass(inode);

    velocity(inode,0,current_state) = v_new[0] = velocity(inode,0,previous_state) + dt*acceleration(inode,0,previous_state);
    velocity(inode,1,current_state) = v_new[1] = velocity(inode,1,previous_state) + dt*acceleration(inode,1,previous_state);
    velocity(inode,2,current_state) = v_new[2] = velocity(inode,2,previous_state) + dt*acceleration(inode,2,previous_state);

    displacement(inode,0,current_state) = displacement(inode,0,previous_state) + dt*v_new[0];
    displacement(inode,1,current_state) = displacement(inode,1,previous_state) + dt*v_new[1];
    displacement(inode,2,current_state) = displacement(inode,2,previous_state) + dt*v_new[2];

#if 0
    std::cout << "acceleration (previous) = (";
    std::cout << acceleration(inode,0,previous_state) << ",";
    std::cout << acceleration(inode,1,previous_state) << ",";
    std::cout << acceleration(inode,2,previous_state) << ")" << std::endl;

    std::cout << "acceleration (current) = (";
    std::cout << acceleration(inode,0,current_state) << ",";
    std::cout << acceleration(inode,1,current_state) << ",";
    std::cout << acceleration(inode,2,current_state) << ")" << std::endl;

    std::cout << "velocity (previous) = (";
    std::cout << velocity(inode,0,previous_state) << ",";
    std::cout << velocity(inode,1,previous_state) << ",";
    std::cout << velocity(inode,2,previous_state) << ")" << std::endl;

    std::cout << "velocity (current) = (";
    std::cout << velocity(inode,0,current_state) << ",";
    std::cout << velocity(inode,1,current_state) << ",";
    std::cout << velocity(inode,2,current_state) << ")" << std::endl;

    std::cout << "displacement (previous) = (";
    std::cout << displacement(inode,0,previous_state) << ",";
    std::cout << displacement(inode,1,previous_state) << ",";
    std::cout << displacement(inode,2,previous_state) << ")" << std::endl;

    std::cout << "displacement (current) = (";
    std::cout << displacement(inode,0,current_state) << ",";
    std::cout << displacement(inode,1,current_state) << ",";
    std::cout << displacement(inode,2,current_state) << ")" << std::endl;
#endif
  }

      const array_type      internal_force;
      const array_type      nodal_mass;
      const array_type      acceleration;
      const array_type      velocity;
      const array_type      displacement;
      const Scalar          dt;
      const int             current_state;
      const int             previous_state;
};
