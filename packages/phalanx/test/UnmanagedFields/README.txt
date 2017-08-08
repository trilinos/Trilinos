Each subdirectory tests the unmanaged field support for one of the
three Field types (MDField, Field and Kokkos::View). The tests mimic
eachother for consistency. MDField checks both dynamic and static
fields. We have copied this pattern over into the other tests even
though they do not need the dynamic testing.

This test also covers the FieldManager.getFieldData() methods (both
const and non-const versions).
