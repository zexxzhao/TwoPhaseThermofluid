
include ../../make.config
include ../make.rules

OBJ=${DIRECTORY}/src/obj


FILENAME:=utility.c
OBJFILES:=${FILENAME:%.c=${OBJ}/%.o}

all: ${OBJFILES}

clean:
	rm -f ${OBJFILES}

.PHONY: clean

include ../make.rules
-include ${OBJFILES:.o=.d}