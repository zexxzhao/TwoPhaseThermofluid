
include ../make.config

NAME:=ns
TARGET:=${DEST}/${NAME}.exe 
OBJDIR:=${DIRECTORY}/src/obj

SUBDIR:=legacy util

all: ${TARGET}

info:
	@echo "FC    :" ${FC}  
	@echo "FFLAGS:" ${FFLAGS}  
	@echo "CC    :" ${CC}  
	@echo "CFLAGS:" ${CFLAGS}  
	@echo "NAME  :" ${NAME}
	@echo "DEST  :" ${DEST} 
	@echo "LIBS  :" ${LIBS}	
	@echo "OBJDIR  :" ${OBJDIR}


${TARGET}: ${SUBDIR} | ${DEST}
	${CC} ${CFLAGS} ${OBJDIR}/*.o ${LIBS} ${MPI_FLAGS} -lm -ldl -o $@ 

${SUBDIR}: ${OBJDIR}
	for DIR in $@; do\
		${MAKE} -C $${DIR} OBJDIR=${OBJDIR}; \
	done

${OBJDIR}:
	mkdir -p ${OBJDIR}

clean:
	for DIR in ${SUBDIR}; do\
		cd $${DIR} && ${MAKE} clean && cd ..; \
	done
	rm -f ${TARGET}
	rm -rf ${OBJDIR}

.PHONY: all ${SUBDIR} info clean

