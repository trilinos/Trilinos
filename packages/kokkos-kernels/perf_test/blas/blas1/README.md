Due to licensing restrictions, we've added a simple testing framework instead of
leveraging https://github.com/google/benchmark.
# To add a blas routine to this directory
1. Create `KokkosBlas_ROUTINE_perf_test.hpp`
2. Populate `KokkosBlas_ROUTINE_perf_test.hpp` with:
```bash
// Forward declarations
void do_ROUTINE_serial_blas(options_t options);
void do_ROUTINE_serial_batched(options_t options);
void do_ROUTINE_parallel_blas(options_t options);
void do_ROUTINE_parallel_batched(options_t options);

// ROUTINE invoke table
void (*do_ROUTINE_invoke[LOOP_N][TEST_N])(options_t) = {
  {do_ROUTINE_serial_blas, do_ROUTINE_serial_batched},
  {do_ROUTINE_parallel_blas, do_ROUTINE_parallel_batched}
};
```
. Update the definitions in `KokkosBlas_common.hpp`, where the comment `//ADD MORE BLAS ROUTINES HERE` is.
4. Add a conditional to invoke the new routine via `do_ROUTINE_invoke` in
   `KokkosBlas_trmm_perf_test.hpp`, where the comment `//ADD MORE BLAS ROUTINES HERE` is.
5. Update the commandline argument processing in
   `KokkosBlas_trmm_perf_test.hpp` to specify how to run ROUTINE.
6. Append `ROUTINE,` to `#define DEFAULT_BLAS_ROUTINES` in `KokkosBlas_common.hpp`.
