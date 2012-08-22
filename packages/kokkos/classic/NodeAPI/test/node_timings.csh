#!/bin/tcsh 

set SerialNodeNT    = "1"
set TBBNodeNT       = "1 2"
set TPINodeNT       = "1 2"
set ThrustGPUNodeNT = "1"
set sizes           = "10000 100000 1000000 10000000"

echo "Legend: "
echo "node problem_size num_cpu_threads"
echo "float Kokkos init time"
echo "float native init time"
echo "float Kokkos sum  time"
echo "float native sum  time"
echo "int Kokkos init time"
echo "int native init time"
echo "int Kokkos sum  time"
echo "int native sum  time"
echo ""

foreach size ($sizes)
  foreach node (SerialNode TBBNode TPINode ThrustGPUNode)
    foreach nt (`eval echo \$${node}NT`)
      echo $node $size $nt
      ./Kokkos_${node}TestAndTiming.exe --test-size=${size} --num-iters=20 --show-test-details=ALL --num-threads=${nt} --not-unit-test="*MemoryInitTest*" | grep "Kokkos<\|Native" | sed "s/.*time: //"
    end
  end
end
