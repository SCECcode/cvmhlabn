#!/bin/bash

cd test_cvmh

make run_unit | test result_unit.txt

make run_accept | test result_accept.txt

p=`grep -c FAIL result.txt` 

if [ $p != 0 ]; then
   echo "something wrong.."
   exit 1 
else
   echo "okay"
   exit 0
fi


