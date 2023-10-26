exe_path=$(head -n 1 "binary_dir.txt")

${exe_path}/SPMV/KokkosBatched_Test_SPMV -A ../data/A.mm -B ../data/B.mm -X ../output/X_SPMV -timers ../output/timers_SPMV -n1 10 -n2 100 -team_size -1 -implementation 3 -l -vector_length 8 -N_team 8