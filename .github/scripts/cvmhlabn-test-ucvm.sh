#!/bin/bash

## testing cvmhlabn validation with ucvm

cd test_validation

make run_validate | tee result_validate.txt

p=`grep -c FAIL result_validate.txt` 
if [ $p != 0 ]; then
   echo "something wrong.."
   exit 1 
fi

make run_validate_with_ucvm | tee result_validate_ucvm.txt
p=`grep -c FAIL result_validate_ucvm.txt` 
if [ $p != 0 ]; then
   echo "something wrong.."
   exit 1 
fi

exit 0 

