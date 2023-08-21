
include ../../make.config

#OBJDIR=${DIRECTORY}/src/obj


FILENAME:=config.c utility.c #main.c
OBJFILES:=${FILENAME:%.c=${OBJDIR}/%.o}

all: ${OBJFILES}

clean:
	rm -f ${OBJFILES}

.PHONY: clean

include ../make.rules
-include ${OBJFILES:.o=.d}
