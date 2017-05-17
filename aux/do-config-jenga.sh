#!/bin/bash -ex

# Modify these paths for your system
TRILINOS_DIR=${TRILINOS_INSTALL_DIR}
GOAL_SRC_DIR=/lore/granzb/goal
GOAL_INSTALL_DIR=/lore/granzb/goal/install
GOAL_DATA_DIR=/lore/granzb/goal-data

cmake $GOAL_SRC_DIR \
-DCMAKE_CXX_COMPILER="mpicxx" \
-DCMAKE_INSTALL_PREFIX=$GOAL_INSTALL_DIR \
-DBUILD_TESTING=ON \
-DTrilinos_PREFIX=$TRILINOS_DIR \
-DGoal_DATA=$GOAL_DATA_DIR \
-DGoal_FAD_SIZE=60 \
-DGoal_ENABLE_ALL_APPS=ON \
-DGoal_CXX_FLAGS="-g -O2 -std=c++11 -Wall -Wno-deprecated-declarations -Wno-unused-parameter -Wno-unused-local-typedefs" \
2>&1 | tee config_log
