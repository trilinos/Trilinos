#!/bin/sh
#MULTX="1 2 3 4 5 6 7 8 9 10"
#MULTX="1 2 3 4 6 8 10 12 14 16"
MULTX="1 2"
#MULTX="3"
TEMPLATE=petra.pbs.template
SPREF=epetra
DPREF=run
CDIR=`pwd`
BINARY_EXE="MueLu_RAPScalingTest_UsingCMSMMM.exe"
BINARY_EXE2="MueLu_RAPScalingTest_UsingMLMMM.exe"
CPN=24
#NXBASE=60
#MULTIPLIES=100

NXBASE=350
MULTIPLIES=50

MATRIX=Laplace2D

# What timing info am I extracting?
#declare -a LABELS=(epetra reuse ix ml-2mm)
#declare -a TIMELINES=("Epetra " "Reuse" "All I&X  " "ML:2matmult")
declare -a LABELS=(RAP-kernel)
declare -a TIMELINES=("RAP kernel")

###############################
calc(){
  awk "BEGIN{ printf $* }" | awk "{printf q\"%5.1f\", \$1 }"
}
###############################
template_nodes(){
  MULT=$1

  expr $MULT \* $MULT
}

###############################
arcadia_build(){
  NODES=$1; DIR=$2; BINARY=$3; CPN=$4; NX=$5 ; BINARY2=$6
  echo "Building $DIR..."

  # Make directory, if needed
  if [ -d $DIR ]; then :
  else mkdir $DIR; fi

  # Link Binary file, if needed
  if [ -e $DIR/$BINARY ]; then :
  else cd $DIR; ln -s $CDIR/$BINARY .; cd $CDIR; fi

  # Link Binary file, if needed
  if [ -e $DIR/$BINARY2 ]; then :
  else cd $DIR; ln -s $CDIR/$BINARY2 .; cd $CDIR; fi

  # Build SLURM script
  #cat $TEMPLATE | sed "s/_NODES_/$NODES/" | sed "s/_NX_/$NX/g" \
  #    | sed "s/_NMULT_/$MULTIPLIES/" | sed "s/_CPN_/$CPN/" \
  #    | sed "s/_MTYPE_/$MATRIX/" \
  #    > $DIR/$SPREF-$NODES.slurm

  # Build PBS Script
  CORES=`expr $NODES \* $CPN`
  cat $TEMPLATE | sed "s/_NODES_/$NODES/" | sed "s/_NX_/$NX/g" \
      | sed "s/_NMULT_/$MULTIPLIES/" | sed "s/_CPN_/$CPN/" \
      | sed "s/_CORES_/$CORES/" \
      | sed "s/_MTYPE_/$MATRIX/" \
      > $DIR/$SPREF-$NODES.pbs


}

###############################
arcadia_clean(){
  NODES=$1; DIR=$2;
  rm -rf $DIR
}

###############################
arcadia_run(){
  NODES=$1; DIR=$2; BATCH=$3
  echo "Running $DIR..."

  cd $DIR
  #sbatch $SPREF-$NODES.slurm
  qsub $SPREF-$NODES.pbs
  cd $CDIR
}

###############################
arcadia_analyze(){
  NODES=$1; DIR=`echo $2 | awk '{ printf "%20s", $1 }'`


  cd $DIR
  EXEC=`ls screen_1.cms 2> /dev/null | wc -l`
  if [ $EXEC -eq 0 ]; then
    echo "$DIR: not run"
  else
    FNAME=`ls -rt screen_1.cms |tail -n 2`

   for ((M=0;M<${#TIMELINES[@]};M++)); do
     M2=`expr $M + ${#TIMELINES[@]}`
     T=${TIMELINES[${M}]}
     TIME[$M]=`grep "$T" screen_1.cms | cut -f3 -d')' | cut -f1 -d'(' | awk '{printf "%11.2f", $1 }'`
     TIME[$M2]=`grep "$T" screen_1.ml | cut -f3 -d')' | cut -f1 -d'(' | awk '{printf "%11.1f", $1 }'`
   done

   if [ $NODES == 1 ] ; then
     for ((M=0;M<${#TIMELINES[@]};M++)); do
       M2=`expr $M + ${#TIMELINES[@]}`
       BASE_TIME[$M]=${TIME[$M]}
       BASE_TIME[$M2]=${TIME[$M2]}
     done
   fi

   for ((M=0;M<${#TIMELINES[@]};M++)); do
     M2=`expr $M + ${#TIMELINES[@]}`
     if [ "${TIME[$M]}" == "" ]; then EFF[$M]=0.0;
     else EFF[$M]=`calc 100 \* ${BASE_TIME[$M]} / ${TIME[$M]}`; fi

     if [ "${TIME[$M2]}" == "" ]; then EFF[$M2]=0.0;
     else EFF[$M2]=`calc 100 \* ${BASE_TIME[$M2]} / ${TIME[$M2]}`; fi

   done

   OUTSTR=""
   for ((M=0;M<${#TIMELINES[@]};M++)); do
     OUTSTR="$OUTSTR ${TIME[$M]} ${EFF[$M]}% "
     OUTSTR="$OUTSTR ${TIME[$M2]} ${EFF[$M2]}% "
   done
   echo "$DIR: $OUTSTR"
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
    TL=`echo ${LABELS[$M]} | awk '{ printf "cs%5s-time    eff", $1 }'`
    OUTSTR="$OUTSTR $TL"
    TL=`echo ${LABELS[$M]} | awk '{ printf "ml%5s-time    eff", $1 }'`
    OUTSTR="$OUTSTR $TL"
  done
    echo "% file              :       $OUTSTR"
fi



# Main Loop
for I in $MULTX; do
    NODES=`template_nodes $I`
    BIN=$BINARY_EXE
    DIR=${DPREF}_$NODES
    NX=`expr $I \* $NXBASE`
    BIN2=$BINARY_EXE2

      # Build Mode
    if [ "$OPT" == "b" ]; then
	arcadia_build $NODES $DIR $BIN $CPN $NX $BIN2
    elif [ "$OPT" == "r" ]; then
	arcadia_run $NODES $DIR $BATCH
    elif [ "$OPT" == "c" ]; then
	arcadia_clean $NODES $DIR
    elif [ "$OPT" == "a" ]; then
	arcadia_analyze $NODES $DIR
    fi
done


