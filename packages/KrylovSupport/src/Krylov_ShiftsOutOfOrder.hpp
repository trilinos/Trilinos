#ifndef __Krylov_ShiftsOutOfOrder_hpp
#define __Krylov_ShiftsOutOfOrder_hpp


///
/// \note Shifts can only be out of order for modified Leja ordering.
struct ShiftsOutOfOrder : public std::exception {
  ShiftsOutOfOrder (const int where, const int when)
    : where_ (where), when_ (when)
  {}
  const char* what() const { 
    if (errMsg_.size() == 0)
      {
	// This method really shouldn't throw, but ostringstream
	// operations might throw due to out-of-memory.  However, if we
	// run out of memory with this little work, the program is
	// completely trashed, so I'm not going to worry about it.
	std::ostringstream os;
	os << "Shifts out of order at index " << where_ << " of inputShifts, at "
	  "step " << when_ << " of computing the modified Leja ordering";
	errMsg_ = os.str();
      }
    return errMsg_.c_str();
  }
  std::string errMsg_;
  int where_, when_;
};

} // namespace Krylov
#endif // __Krylov_ShiftsOutOfOrder_hpp
