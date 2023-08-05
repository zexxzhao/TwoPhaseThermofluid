##################################
# Configuration
##################################
include ../../make.config

##################################
#  Files
##################################
include file.list


MODOBJFILES = $(MODSRCFILES:%.f90=${OBJDIR}/%.o)
MODFILES = $(MODSRCFILES:%.f90=${OBJDIR}/%.mod)

OBJFILES = $(F90SRCFILES:%.f90=${OBJDIR}/%.o)\
           $(F77SRCFILES:%.f=${OBJDIR}/%.o)

##################################
# Targets
##################################

all: ${OBJFILES}

${OBJFILES}: ${MODOBJFILES}

clean:
	rm -f ${OBJFILES} 


include ../make.rules

.PHONY.: all clean

