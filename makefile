

include make.config
SUBDIR:=src genmesh


all:
	for DIR in ${SUBDIR}; do\
		cd $${DIR} && ${MAKE} && cd ..; \
	done


info:
	for DIR in ${SUBDIR}; do\
		cd $${DIR} && ${MAKE} info && cd ..; \
	done

install: ${BIN_INSTALL} ${LIB_INSTALL} ${INC_INSTALL}
	for DIR in ${SUBDIR}; do\
		cd $${DIR} && ${MAKE} install && cd ..; \
	done

${BIN_INSTALL}:
	mkdir ${BIN_INSTALL}

${LIB_INSTALL}:
	mkdir ${LIB_INSTALL}

${INC_INSTALL}:
	mkdir ${INC_INSTALL}

clean:
	for DIR in ${SUBDIR}; do\
		cd $${DIR} && ${MAKE} clean && cd ..; \
	done
	
distclean: clean
	for DIR in ${SUBDIR}; do\
		cd $${DIR} && ${MAKE} distclean && cd ..; \
	done