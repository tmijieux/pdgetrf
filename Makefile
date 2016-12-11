TARGET=pdgetrf perf_driver test_driver
CFLAGS=-std=gnu99 -g -Wall -Wextra
LDFLAGS=
GENGETOPT=gengetopt
CC=mpiicc
#CC=mpicc

CFLAGS+=#-fdiagnostics-color=auto

ifdef DEBUG
CFLAGS+=-ggdb -O0 -DDEBUG=1
else
CFLAGS+=-O3
endif

ifeq ($(strip $(BLASLIB)),)
LDFLAGS+=-lopenblas -lm
else
LDFLAGS+= -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
CFLAGS+=-DTDP_USE_MKL  -DMKL_ILP64 -I${MKLROOT}/include
endif


SRC = 	perf/perf.c \
	error.c \
	util.c \
	getrf.c \
	gesv.c \
	trsv.c

test_SRC = test_driver.c $(SRC)
perf_SRC = perf_driver.c $(SRC)
pdgetrf_SRC = main.c $(SRC)

test_OBJ=$(test_SRC:.c=.o)
perf_OBJ=$(perf_SRC:.c=.o)
pdgetrf_OBJ=$(pdgetrf_SRC:.c=.o)

DEP=$(SRC:.c=.d) main.d perf_driver.d

all: $(TARGET)

-include $(DEP)

pdgetrf: $(pdgetrf_OBJ)
	$(CC) $^ -o $@ $(LDFLAGS)

perf_driver: $(perf_OBJ)
	$(CC) $^ -o $@ $(LDFLAGS)

test_driver: $(test_OBJ)
	$(CC) $^ -o $@ $(LDFLAGS)

%.o: %.c
	@$(CC) -MM $(CFLAGS) $*.c > $*.d
	$(CC) -c $(CFLAGS) $*.c -o $*.o

clean:
	$(RM) $(OBJ) $(DEP) *.d *.o

mrproper: clean
	$(RM) $(TARGET)

genopt: matprod.ggo
	$(GENGETOPT) -u"MatrixA MatrixB" < $^

