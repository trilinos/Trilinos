#!/bin/sh
#MULTX="1 2 3 4 5 6 7 8 9 10"
MULTX="1 2 5 10 20 40 80 100"
TEMPLATE=muelu.pbs.template
SPREF=muelu
DPREF=run
CDIR=`pwd`
BINARY_EXE="MueLu_ScalingTest.exe"
CPN=24
#NXBASE=60
#MULTIPLIES=100
MATRIX=Laplace1D

# Glory / Hopper
if [ $MATRIX == "Laplace1D" ]; then
  NXBASE=22000000
elif [ $MATRIX == "Laplace2D" ]; then
  NXBASE=2000
else
  NXBASE=150
fi

# Redsky
#if [ $MATRIX == "Laplace1D" ]; then
#  NXBASE=4000000
#elif [ $MATRIX == "Laplace2D" ]; then
#  NXBASE=2000
#else
#  NXBASE=150
#fi
      


# What timing info am I extracting?
declare -a LABELS=(build setup solve, its)
declare -a TIMELINES=("Matrix Build" "MueLu Setup" "Belos Solve" "Number of iterations")

###############################
calc(){
  awk "BEGIN{ printf $* }" | awk "{printf \"%5.1f\", \$1 }"
}
###############################
template_nodes(){
  MULT=$1

  if [ $MATRIX == "Laplace1D"   ]; then 
    expr $MULT
  elif [ $MATRIX == "Laplace2D" ]; then 
    expr $MULT \* $MULT 
  else
    expr $MULT \* $MULT \* $MULT
  fi
}

###############################
ocoee_build(){
  NODES=$1; DIR=$2; BINARY=$3; CPN=$4 NX=$5
  echo "Building $DIR..."
  CORES=`expr $CPN \* $NODES`   

  # Make directory, if needed
  if [ -d $DIR ]; then :
  else mkdir $DIR; fi
 
  # Link Binary file, if needed
  if [ -e $DIR/$BINARY ]; then :
  else cd $DIR; ln -s $CDIR/$BINARY .; cd $CDIR; fi

  # Build SLURM script
  cat $TEMPLATE | sed "s/_NODES_/$NODES/" | sed "s/_NX_/$NX/g" \
      | sed "s/_CPN_/$CPN/" \
      | sed "s/_CORES_/$CORES/" \
      | sed "s/_MTYPE_/$MATRIX/" \
      | sed "s#_BINARY_#./$BINARY#" \
      > $DIR/$SPREF-$NODES.pbs

}

###############################
ocoee_clean(){
  NODES=$1; DIR=$2;
  rm -rf $DIR
}

###############################
ocoee_run(){
  NODES=$1; DIR=$2; BATCH=$3
  echo "Running $DIR..." 

  cd $DIR
  qsub $SPREF-$NODES.pbs
  cd $CDIR
}

###############################
ocoee_analyze(){
  NODES=$1; DIR=`echo $2 | awk '{ printf "%20s", $1 }'`


  cd $DIR
  EXEC=`ls screen.out* 2> /dev/null | wc -l`
  if [ $EXEC -eq 0 ]; then
    echo "$DIR: not run"
  else
    FNAME=`ls -rt screen.out* |tail -n 2`   

   for ((M=0;M<${#TIMELINES[@]};M++)); do
     T=${TIMELINES[${M}]}
     #TIME[$M]=`grep $T screen.out.* | cut -f3 -d')' | cut -f1 -d'(' | awk '{printf "%11.1f", $1 }'`
     TIME[$M]=`grep "$T" screen.out.* | cut -f2 -d':' | awk '{printf "%11.1f", $1 }'`
   done

   NORM=`grep Problem screen.out.* | cut -f2 -d':' | awk '{printf "%6.4e", $1}'`

   if [ $NODES == 1 ] ; then 
     for ((M=0;M<${#TIMELINES[@]};M++)); do
       BASE_TIME[$M]=${TIME[$M]}
     done
   fi 

   for ((M=0;M<${#TIMELINES[@]};M++)); do
     if [ "${TIME[$M]}" == "" ]; then EFF[$M]=0.0;
     else EFF[$M]=`calc 100 \* ${BASE_TIME[$M]} / ${TIME[$M]}`; fi
   done

   OUTSTR=""
   for ((M=0;M<${#TIMELINES[@]};M++)); do
     OUTSTR="$OUTSTR ${TIME[$M]} ${EFF[$M]}% "
   done
   echo "$DIR: $NORM $OUTSTR"
  fi
  cd ..
}


###############################
###############################
###############################
OPT=$1
# Sanity Check
if [ "$OPT" != "b" ] && [ "$OPT" != "r" ] && [ "$OPT" != "c" ] && [ "$OPT" != "a" ]; then   echo "Syntax: $0 [b|r|c|a]"; exit 1; fi

# Analyze header
if [ "$OPT" == "a" ]; then 
  OUTSTR=""
  for ((M=0;M<${#LABELS[@]};M++)); do
    TL=`echo ${LABELS[$M]} | awk '{ printf "%7s-time    eff", $1 }'`  
    OUTSTR="$OUTSTR $TL"
  done
    echo "% file              :       norm$OUTSTR"
fi


# Main Loop
for I in $MULTX; do
    NODES=`template_nodes $I`
    BIN="$BINARY_EXE"
    DIR=${DPREF}_$NODES
    NX=`expr $I \* $NXBASE`
 
      # Build Mode
    if [ "$OPT" == "b" ]; then
	ocoee_build $NODES $DIR $BIN $CPN $NX
    elif [ "$OPT" == "r" ]; then
	ocoee_run $NODES $DIR $BATCH
    elif [ "$OPT" == "c" ]; then
	ocoee_clean $NODES $DIR
    elif [ "$OPT" == "a" ]; then
	ocoee_analyze $NODES $DIR
    fi
done


