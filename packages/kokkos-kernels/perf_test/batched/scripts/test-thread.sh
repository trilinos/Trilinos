exec=$1
numacmd="numactl --membind 1"
testfile=$exec.thread.txt
rm -f $testfile

echo $exec > $testfile

for th in 1 2 4 8 16 34 68; do
    echo "$numacmd ./$exec --kokkos-threads=$th  >> $testfile"
    $numacmd ./$exec --kokkos-threads=$th  >> $testfile
done
