
For the test "periodic_32bit_int_limit":

The exodus mesh must be generated using cubit with the journal file
peridoic_32bit_int_limit.jou. Then the gids must be shifted above the
32bit int limit using the SEACAS "grepos" tool/executable. Run grepos
with:

grepos -64 periodic_32bit_int_limit_UNSHIFTED.exo periodic_32bit_int_limit.exo

Once in grepos do the following to shift the
indices:

increment nodemap 7000000000
increment elemmap 7000000000
exit

On exit, the new file will be written with the shifted indices.
