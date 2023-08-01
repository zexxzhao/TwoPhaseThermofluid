##################################
# Configuration
##################################
include make.config

##################################
# Targets
##################################

default: ${DEST}/${NAME}.exe 

info:
	@echo "FC    :" ${FC}  
	@echo "FFLAGS:" ${FFLAGS}  
	@echo "CC    :" ${CC}  
	@echo "CFLAGS:" ${CFLAGS}  
	@echo "NAME  :" ${NAME}
	@echo "DEST  :" ${DEST} 
	@echo "LIBS  :" ${LIBS}	
	@echo "OBJ   :" ${OBJ} 

${DEST}/${NAME}.exe:${OBJFILES} | ${DEST}
	rm -rf ${DEST}/${NAME}.exe
	${FC}   ${OBJFILES} ${LIBS} -o  ${DEST}/${NAME}.exe
	@ls -lh ${DEST}/${NAME}.exe

${OBJ}:
	mkdir   ${OBJ}

${DEST}:	
	mkdir   ${DEST}
	
clean:
	rm -rf obj *.mod ${DEST}/${NAME}.exe
	
##################################
# Rules
##################################


${OBJ}/%.o: %.f90 | ${OBJ}
	$(FC) $(FFLAGS) -c $< -o $@

${OBJ}/%.o: %.f | ${OBJ}
	$(FC) $(FFLAGS) -c $< -o $@



.PHONY.: clean

