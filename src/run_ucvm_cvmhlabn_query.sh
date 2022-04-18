##
##

IN_FILE=$1

if [ "x${UCVM_INSTALL_PATH}" != "x" ] ; then
  SCRIPT_DIR=${UCVM_INSTALL_PATH}/bin
  source $SCRIPT_DIR/../conf/ucvm_env.sh
  ucvm_query -f ${SCRIPT_DIR}/../conf/ucvm.conf -m cvmhlabn  < ${IN_FILE}  
  else
    echo "need to have UCVM_INSTALL_PATH set!!!"
fi


