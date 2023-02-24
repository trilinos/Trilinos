exec=$1
numacmd="numactl --membind 1"
testfile=$exec.workset.hot.txt
rm -f $testfile

echo $exec > $testfile

for w in 1024 2048 4096 8192 16384 32768 65536 131072; do
    echo "$numacmd ./$exec --kokkos-threads=68 >> $testfile"
    $numacmd ./$exec --kokkos-threads=68 -N $w -hot-cache  >> $testfile
done
