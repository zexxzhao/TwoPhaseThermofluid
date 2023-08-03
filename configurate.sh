
CONFIGFILE=make.config

if [ -f ${CONFIGFILE} ]; then
    echo "rm ${CONFIGFILE}"
    bash -c "rm ${CONFIGFILE}"
fi

echo "##########################" > ${CONFIGFILE}
echo "# System definitions" >> ${CONFIGFILE}
echo "##########################" >> ${CONFIGFILE}

echo >> ${CONFIGFILE}
echo "FC=mpif90">> ${CONFIGFILE}
echo "FFLAGS:=-O2 -g -fPIC -fbounds-check -fcheck=mem -fallow-argument-mismatch" >> ${CONFIGFILE}

echo >> ${CONFIGFILE}
echo "CC=mpicc">> ${CONFIGFILE}
echo "CFLAGS:=-O2 -g -fPIC " >> ${CONFIGFILE}

echo >> ${CONFIGFILE}
echo  "LIBS=" >> ${CONFIGFILE}
echo  "INC=" >> ${CONFIGFILE}

echo >> ${CONFIGFILE}
echo >> ${CONFIGFILE}
echo "##########################" >> ${CONFIGFILE}
echo "# PATH" >> ${CONFIGFILE}
echo "##########################" >> ${CONFIGFILE}


echo "DIRECTORY:=${PWD}">> ${CONFIGFILE}
echo "DEST:=${PWD}/bin">> ${CONFIGFILE}
echo >> ${CONFIGFILE}
