#!/bin/bash

cd test_cvmh

make run_unit | test result_unit.txt

p=`grep -c FAIL result_unit.txt` 
if [ $p != 0 ]; then
   echo "something wrong.."
   exit 1 
else
   echo "okay"
   exit 0
fi

make run_accept | test result_accept.txt

p=`grep -c FAIL result_accept.txt` 
if [ $p != 0 ]; then
   echo "something wrong.."
   exit 1 
else
   echo "okay"
   exit 0
fi
