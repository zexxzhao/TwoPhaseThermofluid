##################################
# Configuration
##################################
include ../../make.config

##################################
#  Files
##################################
include file.list


MODOBJFILES = $(MODSRCFILES:%.f90=%.o)

OBJFILES = $(F90SRCFILES:%.f90=%.o)\
           $(F77SRCFILES:%.f=%.o)

TARGET:=${BUILD_DIR}/libFEM_legacy${STATIC_LIB_SUFFIX}

##################################
# Targets
##################################


all: ${TARGET} 

info: 

install: all 
	@echo "Installing F90 MOD $< ..."
	@mkdir -p ${INC_INSTALL}
	@cp -f *.mod ${INC_INSTALL}


clean:
	rm -f ${OBJFILES} 
	rm -f ${MODOBJFILES}
	rm -f *.mod 
	rm -f ${TARGET}

distclean: clean

${OBJFILES}: ${MODOBJFILES}
##################################
# Rules
##################################

include ../make.rules

.PHONY.: all clean distclean install

${TARGET}: ${OBJFILES} ${MODOBJFILES}
	@echo "Creating static library $@ ..."
	${AR} ${ARFLAGS} -o $@ $^
	@ranlib $@

