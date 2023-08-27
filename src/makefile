
include ../make.config

#NAME:=ns
BUILD_DIR:=build

SUBDIR:=legacy util sys# program

LIBNAME:=FEM

STATIC_LIBS:=${BUILD_DIR}/lib${LIBNAME}_legacy${STATIC_LIB_SUFFIX}\
	${BUILD_DIR}/lib${LIBNAME}_util${STATIC_LIB_SUFFIX}

SONAME:=${BUILD_DIR}/lib${LIBNAME}${SHARED_LIB_SUFFIX}

# static lib from each subdir

F90EXEFILE:=$(wildcard program_dir/*.f90)
CEXEFILE:=$(wildcard program_dir/*.c)
EXE_NAME=${F90EXEFILE:%.f90=%${EXE_SUFFIX}}\
	${CEXEFILE:%.c=%${EXE_SUFFIX}}
#remove directory
EXE_NAME:=$(notdir ${EXE_NAME})


all: ${SONAME}  

info:
	@echo "FC    :" ${FC}  
	@echo "FFLAGS:" ${FFLAGS}  
	@echo "CC    :" ${CC}  
	@echo "CFLAGS:" ${CFLAGS}  
	@echo "NAME  :" ${NAME}
	@echo "DEST  :" ${DEST} 
	@echo "LIBS  :" ${LIBS}	
	@echo "STATIC:" ${STATIC_LIBS}
	@echo "SHARED:" ${SONAME}
	@echo "EXE   :" ${EXE_NAME}


install_subdir: all
	for DIR in ${SUBDIR}; do\
		cd $${DIR}_dir && ${MAKE} install BUILD_DIR=../${BUILD_DIR} && cd ..; \
	done
	@echo "Installing shared library $< ..."
	@mkdir -p ${LIB_INSTALL}
	@cp -f ${SONAME} ${LIB_INSTALL}

install: install_subdir ${EXE_NAME:%=${BUILD_DIR}/%}
	@echo "Installing executables $< ..."
	@mkdir -p ${BIN_INSTALL}
	@cp -f ${EXE_NAME:%${EXE_SUFFIX}=${BUILD_DIR}/%${EXE_SUFFIX}} ${BIN_INSTALL}

clean:
	for DIR in ${SUBDIR}; do\
		cd $${DIR}_dir && ${MAKE} clean BUILD_DIR=../${BUILD_DIR} && cd ..; \
	done
	cd program_dir && ${MAKE} clean BUILD_DIR=../${BUILD_DIR} && cd ..
	rm -f ${SONAME} ${EXE_NAME}
	rm -rf ${BUILD_DIR}

distclean: clean
	rm -r ${LIB_INSTALL} ${BIN_INSTALL} ${INC_INSTALL}

.PHONY: all info clean distclean install ${SUBDIR:%=%_dir}

##################################
# Rules
##################################

${SONAME}: ${STATIC_LIBS} | ${BUILD_DIR}
	@echo "Creating shared library $@ ..."
	mkdir -p tmp_objs
	cd tmp_objs && \
	for lib in ${STATIC_LIBS}; do\
		${AR} -x ../$${lib} || exit 1; \
	done 
	ld -shared -o $@ tmp_objs/*.o
	rm -rf tmp_objs

${BUILD_DIR}/lib${LIBNAME}_%${STATIC_LIB_SUFFIX}: ${BUILD_DIR}
	${MAKE} -C $*_dir BUILD_DIR=../${BUILD_DIR}

# the target build/driver.exe can't find this rule. What's wrong?
#${BUILD_DIR}/%${EXE_SUFFIX}: ${SONAME} program_dir/%.o
${BUILD_DIR}/%.exe: ${SONAME} exe_prerequisite
	@echo "Creating executable $@ ..."
	${CC} ${CFLAGS} ${MPI_FLAGS} -I${INC_INSTALL} \
		-L${LIB_INSTALL} -lFEM -lgfortran -lm -ldl ${LIBS} -o $@ program_dir/$*.o

exe_prerequisite:
	cd program_dir && ${MAKE} BUILD_DIR=../${BUILD_DIR} && cd ..


${BUILD_DIR}:
	mkdir -p $@