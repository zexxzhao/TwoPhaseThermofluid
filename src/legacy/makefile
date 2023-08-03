##################################
# Configuration
##################################
include ../../make.config

##################################
#  Files
##################################
include file.list


MODFILES = $(MODSRCFILES:%.f90=${OBJDIR}/%.mod)

OBJFILES = $(MODSRCFILES:%.f90=${OBJDIR}/%.o) \
		   $(F90SRCFILES:%.f90=${OBJDIR}/%.o)\
           $(F77SRCFILES:%.f=${OBJDIR}/%.o)

##################################
# Targets
##################################

all: ${OBJFILES}

clean:
	rm -f ${OBJFILES} 


include ../make.rules

.PHONY.: all clean

