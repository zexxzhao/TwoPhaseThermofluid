
include ../../make.config

#OBJDIR=${DIRECTORY}/src/obj


HEADERS:=cfile_state.h\
	config.h\
	sparsity_pattern.h\
	utility.h

FILE:=cfile_state.c\
	config.c\
	sparsity_pattern.c\
	utility.c

OBJFILES:=${FILE:%.c=%.o}

TARGET:=${BUILD_DIR}/libFEM_util${STATIC_LIB_SUFFIX}

INC+=-I..

##################################
# Targets
##################################

all: ${TARGET}

info:

install: all
	@echo "Installing header files $< ..."
	@mkdir -p ${INC_INSTALL}/util_dir
	cp -f ${HEADERS} ${INC_INSTALL}/util_dir

clean:
	rm -f ${OBJFILES} ${FILE:%.c=%.d}
	rm -f ${TARGET}

-include ${OBJFILES:.o=.d}

##################################
# Rules
##################################

.PHONY: clean

include ../make.rules

${TARGET}: ${OBJFILES}
	@echo "Creating static library $@ ..."
	${AR} ${ARFLAGS} -o $@ $^
	@ranlib $@	