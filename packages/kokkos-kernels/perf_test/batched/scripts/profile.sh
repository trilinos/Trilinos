exec=$1
b=$2
testfile=$exec.$b.txt
#rm -f $testfile

echo $exec blk $b # > $testfile
nvprof --metrics \
    achieved_occupancy,\
dram_read_throughput,dram_write_throughput,gst_throughput,gld_throughput,l2_l1_read_throughput,l2_tex_read_throughput,l2_tex_write_throughput,local_load_throughput,local_store_throughput,shared_load_throughput,shared_store_throughput,\
gld_efficiency,shared_efficiency,\
l1_cache_hit_rate,tex_cache_hit_rate,\
sm_efficiency,warp_execution_efficiency,branch_efficiency,\
stall_inst_fetch,stall_memory_dependency,stall_exec_dependency,stall_memory_throttle,stall_pipe_busy,stall_not_selected \
    ./$exec -B $b |& tee $testfile



