
FILEIST = applyBC.c

OBJIST = $(FILEIST:.c=.o)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf *.o *.mod *.pc *.il *.d *.obj *.pdb *.lib *.exp *.a

