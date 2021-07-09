Due to licensing restrictions, we've added a simple testing framework instead of
leveraging https://github.com/google/benchmark.
# To add a blas3 routine to this directory
1. Create `KokkosBlas3_ROUTINE_perf_test.hpp`
2. Populate `KokkosBlas3_ROUTINE_perf_test.hpp` with:
```bash
// Forward declarations
void do_ROUTINE_serial_blas(options_t options);
void do_ROUTINE_serial_batched(options_t options);
void do_ROUTINE_parallel_blas(options_t options);
void do_ROUTINE_parallel_batched(options_t options);

// ROUTINE invoke table
void (*do_ROUTINE_invoke[LOOP_N][TEST_N])(options_t) = {
  do_ROUTINE_serial_blas,
  do_ROUTINE_serial_batched,
  do_ROUTINE_parallel_blas,
  do_ROUTINE_parallel_batched
};
```
3. Update the definitions in `KokkosBlas3_common.hpp`, where the comment `//ADD MORE BLAS3 ROUTINES HERE` is.
4. Add the `do_ROUTINE_invoke` table to the `routine_table` in
   `KokkosBlas3_perf_test.cpp`, where the comment `//ADD MORE BLAS3 ROUTINES HERE` is.
5. Update the commandline argument processing in
   `KokkosBlas3_perf_test.cpp` to specify how to run ROUTINE.
6. To run the new routine by default, append `ROUTINE,` to `#define DEFAULT_BLAS_ROUTINES` in `KokkosBlas3_common.hpp`.
