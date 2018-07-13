#!/bin/bash -le

TRILINOS_UPDATE_BRANCH=$1
TRILINOS_PRISTINE_BRANCH=$2
KOKKOS_UPDATE_BRANCH=$3

if [ -z $TRILINOS_UPDATE_BRANCH ]
then
  TRILINOS_UPDATE_BRANCH=develop
fi

if [ -z $TRILINOS_PRISTINE_BRANCH ]
then
  TRILINOS_PRISTINE_BRANCH=develop
fi

export TRILINOS_UPDATED_PATH=${PWD}/trilinos-update
export TRILINOS_PRISTINE_PATH=${PWD}/trilinos-pristine

#rm -rf ${KOKKOSKERNELS_PATH}
#rm -rf ${TRILINOS_UPDATED_PATH}
#rm -rf ${TRILINOS_PRISTINE_PATH}

#Already done:
if [ ! -d "${TRILINOS_UPDATED_PATH}" ]; then
  git clone https://github.com/trilinos/trilinos ${TRILINOS_UPDATED_PATH}
fi
if [ ! -d "${TRILINOS_PRISTINE_PATH}" ]; then
  git clone https://github.com/trilinos/trilinos ${TRILINOS_PRISTINE_PATH}
fi

cd ${TRILINOS_UPDATED_PATH}
git checkout $TRILINOS_UPDATE_BRANCH
git reset --hard origin/$TRILINOS_UPDATE_BRANCH
git pull
cd ..

if [ -n $KOKKOS_UPDATE_BRANCH ]
then
  python ${TRILINOS_PRISTINE_PATH}/packages/kokkos/scripts/snapshot.py ${KOKKOS_PATH} ${TRILINOS_UPDATED_PATH}/packages
fi
python ${TRILINOS_PRISTINE_PATH}/packages/kokkos/scripts/snapshot.py ${KOKKOSKERNELS_PATH} ${TRILINOS_UPDATED_PATH}/packages


cd ${TRILINOS_UPDATED_PATH}
echo ""
echo ""
echo "Trilinos State:"
git log --pretty=oneline --since=7.days
cd ..

cd ${TRILINOS_PRISTINE_PATH}
git status
echo "Checkout $TRILINOS_PRISTINE_BRANCH"
git checkout $TRILINOS_PRISTINE_BRANCH
echo "Pull"
git pull
cd ..

cd ${TRILINOS_PRISTINE_PATH}
echo ""
echo ""
echo "Trilinos Pristine State:"
git log --pretty=oneline --since=7.days
cd ..
