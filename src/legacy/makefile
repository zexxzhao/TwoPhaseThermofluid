##################################
# Configuration
##################################
include ../../make.config
OBJ=${DIRECTORY}/src/obj


##################################
#  Files
##################################
include file.list


MODFILES = $(MODSRCFILES:%.f90=${OBJ}/%.mod)

OBJFILES = $(MODSRCFILES:%.f90=${OBJ}/%.o) \
		   $(F90SRCFILES:%.f90=${OBJ}/%.o)\
           $(F77SRCFILES:%.f=${OBJ}/%.o)

##################################
# Targets
##################################

#default: ${DEST}/${NAME}.exe 

all: ${OBJFILES}

info:
	@echo "FC    :" ${FC}  
	@echo "FFLAGS:" ${FFLAGS}  
	@echo "CC    :" ${CC}  
	@echo "CFLAGS:" ${CFLAGS}  
	@echo "NAME  :" ${NAME}
	@echo "DEST  :" ${DEST} 
	@echo "LIBS  :" ${LIBS}	
	@echo "OBJ   :" ${OBJ} 

#${DEST}/${NAME}.exe:${OBJFILES} ${MODFILES}| ${DEST}
#	rm -rf ${DEST}/${NAME}.exe
#	${FC}   ${OBJFILES} ${LIBS} -o  ${DEST}/${NAME}.exe
#	@ls -lh ${DEST}/${NAME}.exe

${OBJ}:
	mkdir   ${OBJ}


${DEST}:	
	mkdir   ${DEST}
	
clean:
	rm -f ${OBJFILES} 
	
include ../make.rules

.PHONY.: all clean

