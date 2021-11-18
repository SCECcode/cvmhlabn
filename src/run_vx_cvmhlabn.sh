#!/bin/bash

if [ $# -lt 2 ]; then
	printf "Usage: %s: <infile> <outfile>\n" $(basename $0) >&2    
        exit 1
fi

SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"
IN_FILE=$1
OUT_FILE=$2

if [[ -z "${UCVM_INSTALL_PATH}" ]]; then
  if [[ -f "${UCVM_INSTALL_PATH}/model/cvmhlabn/lib" ]]; then
    env DYLD_LIBRARY_PATH=${UCVM_INSTALL_PATH}/model/cvmhlabn/lib ${SCRIPT_DIR}/vx_cvmhlabn < ${IN_FILE} > ${OUT_FILE}
    if [ $? -ne 0 ]; then
        exit 1
    fi
    exit 0
  fi
fi

env DYLD_LIBRARY_PATH=../src ${SCRIPT_DIR}/vx_cvmhlabn < ${IN_FILE} > ${OUT_FILE}
if [ $? -ne 0 ]; then
    exit 1
fi

exit 0
