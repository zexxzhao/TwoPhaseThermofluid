

include make.config
SUBDIR:=src genmesh


all: ${DEST}
	for DIR in ${SUBDIR}; do\
		cd $${DIR} && ${MAKE} && cd ..; \
	done

${DEST}:
	mkdir ${DEST}

clean:
	for DIR in ${SUBDIR}; do\
		cd $${DIR} && ${MAKE} clean && cd ..; \
	done
	rm -r ${DEST}
	
