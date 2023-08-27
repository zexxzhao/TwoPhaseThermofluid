include ../../make.config


CFILE:=main.c
F90FILE:=driver.f90
OBJFILE:=${CFILE:%.c=%.o}\
	${F90FILE:%.f90=%.o}

INC+=-I${INC_INSTALL}

all: ${OBJFILE}

info:

install:

clean:
	rm -f ${OBJFILE} ${CFILE:%.c=%.d} ${F90FILE:%.f90=%.d}

distclean: clean

.PHONY: all clean distclean install


##################################
# Rules
##################################

include ../make.rules
