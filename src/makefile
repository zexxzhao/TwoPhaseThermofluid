
include ../make.config

NAME:=ns
TARGET:=${DEST}/${NAME}.exe 

SUBDIR:=legacy util

all: ${TARGET}

${TARGET}: ${SUBDIR}
	${CC} ${CFLAGS} obj/*.o ${LIBS} ${MPI_COMPILE_FLAGS} ${MPI_LINK_FLAGS} -lgfortran -lm -ldl -o $@ 

${SUBDIR}:
	for DIR in $@; do\
		cd $${DIR} && ${MAKE} && cd ..; \
	done



clean:
	for DIR in ${SUBDIR}; do\
		cd $${DIR} && ${MAKE} clean && cd ..; \
	done
	rm -f ${TARGET}

.PHONY: all ${SUBDIR} clean

