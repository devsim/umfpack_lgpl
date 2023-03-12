#!/bin/bash
set -e
set -u
mkdir -p build
(cd build && ${CMAKE} ${ARCH_ARG} -DCMAKE_C_COMPILER=${CC} -DCMAKE_BUILD_TYPE=RELEASE ../umfpack)
