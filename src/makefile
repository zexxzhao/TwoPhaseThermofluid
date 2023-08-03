
include ../make.config

NAME:=ns
TARGET:=${DEST}/${NAME}.exe 

SUBDIR:=legacy util

all: ${TARGET}

${TARGET}: ${SUBDIR}
	${CC} ${CFLAGS} -o $@ obj/*.o ${LIBS} -lgfortran -lm -ldl -lmpichfort

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

