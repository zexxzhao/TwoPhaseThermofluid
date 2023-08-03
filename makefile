

SUBDIR:=src genmesh

all:
	for DIR in ${SUBDIR}; do\
		cd $${DIR} && ${MAKE} && cd ..; \
	done

clean:
	for DIR in ${SUBDIR}; do\
		cd $${DIR} && ${MAKE} clean && cd ..; \
	done
	
