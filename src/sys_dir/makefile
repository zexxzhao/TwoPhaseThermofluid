include ../../make.config

HEADERS=defines.h

INC+=-I..
all:

install:
	@echo "Installing header files $< ..."
	@mkdir -p ${INC_INSTALL}/sys_dir
	@cp -f ${HEADERS} ${INC_INSTALL}/sys_dir

clean:

distclean: clean

.PHONY: all clean distclean install