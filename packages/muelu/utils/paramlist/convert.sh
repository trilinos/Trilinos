#!/bin/bash
if (( $# != 1 )); then
  echo "Usage: ./convert.sh <in> <out>"
  exit
fi
in=$1
if [ ! -f $in ]; then
  echo "File \"$in\" does not exist."
  exit
fi

sed -i 's/"PDE equations"/"number of equations"/'                               $in
sed -i 's/"numDesiredLevel"/"max levels"/'                                      $in
sed -i 's/"implicit transpose"/"transpose: use implicit"/'                      $in
sed -i '/rcb/! s/"algorithm"/"aggregation: drop scheme"/'                       $in
sed -i '/multijagged/! s/"algorithm"/"aggregation: drop scheme"/'               $in
sed -i 's/"algorithm"/"aggregation: drop scheme"/'                              $in
sed -i 's/"aggregation threshold"/"aggregation: drop tol"/'                     $in
sed -i 's/"MinNodesPerAggregate"/"aggregation: min agg size"/'                  $in
sed -i 's/"MaxNodesPerAggregate"/"aggregation: max agg size"/'                  $in
sed -i 's/"MaxNeighAlreadySelected"/"aggregation: max selected neighbors"/'     $in
sed -i 's/"Dirichlet detection threshold"/"aggregation: Dirichlet threshold"/'  $in
sed -i 's/"maxCoarseSize"/"coarse: max size"/'                                  $in
sed -i '/="0"/! s/"startLevel"/"repartition: start level"/'                     $in
sed -i 's/"minRowsPerProcessor"/"repartition: min rows per proc"/'              $in
sed -i 's/"nonzeroImbalance"/"repartition: max imbalance"/'                     $in
sed -i 's/"remapPartitions"/"repartition: remap parts"/'                        $in
sed -i 's/"numRemapValues"/"repartition: remap num values"/'                    $in
sed -i 's/"alwaysKeepProc0"/"repartition: keep proc 0"/'                        $in
sed -i 's/"Damping factor"/"sa: damping factor"/'                               $in
sed -i 's/"lumping"/"filtered matrix: use lumping"/'                            $in
sed -i 's/"Ordering"/"aggregation: ordering"/'                                  $in
sed -i 's/value="Natural"/value="natural"/'                                     $in
sed -i 's/value="Graph"/value="graph"/'                                         $in
sed -i 's/value="Random"/value="random"/'                                       $in
sed -i 's/"laplacian"/"distance laplacian"/'                                    $in
sed -i 's/"aggregation: visualize"/"aggregation: export visualization data"/'   $in
sed -i 's/"print"/"export data"/'                                               $in
