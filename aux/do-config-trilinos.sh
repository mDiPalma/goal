#!/bin/bash -ex

# Modify these paths for your system.
TPLS=
TRILINOS_SRC_DIR=
TRILINOS_INSTALL_DIR=
MPI_DIR=
PARMETIS_DIR=
BOOST_DIR=
YAML_DIR=
HAVE_LL=

cmake $TRILINOS_SRC_DIR \
\
-DCMAKE_INSTALL_PREFIX:PATH=$TRILINOS_INSTALL_DIR \
-DCMAKE_BUILD_TYPE:STRING=NONE \
-DCMAKE_C_FLAGS:STRING="-g -O2" \
-DCMAKE_CXX_FLAGS:STRING="-g -O2 -Wno-deprecated-declarations -Wno-sign-compare" \
-DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
\
-DBUILD_SHARED_LIBS:BOOL=ON \
\
-DTrilinos_ENABLE_SECONDARY_TESTED_CODE:BOOL=OFF \
-DTrilinos_ENABLE_EXPORT_MAKEFILES:BOOL=OFF \
-DTrilinos_ASSERT_MISSING_PACKAGES:BOOL=OFF \
-DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
-DTrilinos_WARNINGS_AS_ERRORS_FLAGS:STRING="" \
\
-DTrilinos_ENABLE_Teuchos:BOOL=ON \
-DTrilinos_ENABLE_Shards:BOOL=ON \
-DTrilinos_ENABLE_Sacado:BOOL=ON \
-DTrilinos_ENABLE_Belos:BOOL=ON \
-DTrilinos_ENABLE_Tpetra:BOOL=ON \
-DTrilinos_ENABLE_Ifpack2:BOOL=ON \
-DTrilinos_ENABLE_MiniTensor:BOOL=ON \
\
-DTrilinos_ENABLE_Kokkos:BOOL=ON \
-DTrilinos_ENABLE_KokkosCore:BOOL=ON \
-DKokkos_ENABLE_Serial:BOOL=ON \
-DKokkos_ENABLE_OpenMP:BOOL=OFF \
-DKokkos_ENABLE_Pthread:BOOL=OFF \
\
-DTrilinos_ENABLE_Phalanx:BOOL=ON \
-DPhalanx_INDEX_SIZE_TYPE:STRING="INT" \
-DPhalanx_SHOW_DEPRECATED_WARNINGS:BOOL=OFF \
-DPhalanx_KOKKOS_DEVICE_TYPE:STRING="SERIAL" \
\
-DTrilinos_ENABLE_Zoltan:BOOL=ON \
-DZoltan_ENABLE_ULLONG_IDS:BOOL=ON \
\
-DTrilinos_ENABLE_Xpetra:BOOL=OFF \
-DTrilinos_ENABLE_NOX:BOOL=OFF \
-DTrilinos_ENABLE_Thyra:BOOL=OFF \
-DTrilinos_ENABLE_Zoltan2:BOOL=OFF \
-DTrilinos_ENABLE_MueLu:BOOL=OFF \
-DTrilinos_ENABLE_Amesos2:BOOL=OFF \
-DTrilinos_ENABLE_Stratimikos:BOOL=OFF \
-DTrilinos_ENABLE_Rythmos:BOOL=OFF \
-DTrilinos_ENABLE_Stokhos:BOOL=OFF \
-DTrilinos_ENABLE_Piro:BOOL=OFF \
-DTrilinos_ENABLE_Teko:BOOL=OFF \
-DTrilinos_ENABLE_STKIO:BOOL=OFF \
-DTrilinos_ENABLE_STKMesh:BOOL=OFF \
-DTrilinos_ENABLE_SEACAS:BOOL=OFF \
-DTrilinos_ENABLE_Epetra:BOOL=OFF \
-DTrilinos_ENABLE_EpetraExt:BOOL=OFF \
-DTrilinos_ENABLE_Ifpack:BOOL=OFF \
-DTrilinos_ENABLE_AztecOO:BOOL=OFF \
-DTrilinos_ENABLE_Amesos:BOOL=OFF \
-DTrilinos_ENABLE_Anasazi:BOOL=OFF \
-DTrilinos_ENABLE_ML:BOOL=OFF \
-DTrilinos_ENABLE_Intrepid:BOOL=OFF \
\
-DTPL_ENABLE_MPI:BOOL=ON \
-DMPI_BASE_DIR:PATH=$MPI_DIR \
\
-DTPL_ENABLE_Boost:BOOL=ON \
-DBoost_INCLUDE_DIRS:FILEPATH="$BOOST_DIR/include" \
-DBoost_LIBRARY_DIRS:FILEPATH="$BOOST_DIR/lib" \
\
-DTPL_ENABLE_ParMETIS:STRING=ON \
-DParMETIS_INCLUDE_DIRS:PATH="$PARMETIS_DIR/include" \
-DParMETIS_LIBRARY_DIRS:PATH="$PARMETIS_DIR/lib" \
\
-DTPL_ENABLE_METIS:STRING=ON \
-DMETIS_INCLUDE_DIRS:PATH="$PARMETIS_DIR/include" \
-DMETIS_LIBRARY_DIRS:PATH="$PARMETIS_DIR/lib" \
\
-DTPL_ENABLE_yaml-cpp:BOOL=ON \
-Dyaml-cpp_INCLUDE_DIRS:PATH="$YAML_DIR/include" \
-Dyaml-cpp_LIBRARY_DIRS:PATH="$YAML_DIR/lib" \
\
-DTPL_ENABLE_Netcdf:BOOL=OFF \
-DTPL_ENABLE_HDF5:BOOL=OFF \
-DTPL_ENABLE_SuperLU:BOOL=OFF \
-DTPL_ENABLE_X11:BOOL=OFF \
-DTPL_ENABLE_Matio:BOOL=OFF \
\
-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
-DTpetra_INST_FLOAT:BOOL=OFF \
-DTpetra_INST_INT_INT:BOOL=ON \
-DTpetra_INST_DOUBLE:BOOL=ON \
-DTpetra_INST_COMPLEX_FLOAT:BOOL=OFF \
-DTpetra_INST_COMPLEX_DOUBLE:BOOL=OFF \
-DTpetra_INST_INT_LONG:BOOL=OFF \
-DTpetra_INST_INT_UNSIGNED:BOOL=OFF \
-DTpetra_INST_INT_LONG_LONG:BOOL=$HAVE_LL \
\
2>&1 | tee config_log
