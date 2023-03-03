exe_path=$(head -n 1 "binary_dir.txt")

${exe_path}/CG/KokkosBatched_Test_CG -A ../data/A.mm -B ../data/B.mm -X ../output/X_CG -timers ../output/timers_CG -n1 10 -n2 100 -team_size -1 -implementation 0 -l -n_iterations 20 -tol 1e-8 -vector_length 8 -N_team 8