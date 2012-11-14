#!/bin/sh

if [ $# -lt 1 ]; then
  name=graph.dot
else
  name=$1
fi

sed -i 's/label=Graph\]/label=\"Graph\"\]/' $name
sed -i 's/\\"/"/g' $name
sed -i 's/"</</' $name
sed -i 's/>"/>/' $name
dot -Tps $name -o `echo $name | sed 's/.dot/.ps/'`
