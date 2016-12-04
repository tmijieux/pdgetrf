TARGET=dgetrf perf_driver
CFLAGS=-std=gnu99 -g -Wall -Wextra $(shell pkg-config --cflags glib-2.0) -fopenmp
LDFLAGS=-lm $(shell pkg-config --libs glib-2.0) 
GENGETOPT=gengetopt
CC=mpicc

CFLAGS+=-fdiagnostics-color=auto

ifdef DEBUG
CFLAGS+=-ggdb -O0 -DDEBUG=1
else
CFLAGS+=-O3
endif

ifeq ($(strip $(BLASLIB)),)
LDFLAGS+=-lopenblas
else
LDFLAGS+= -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
	${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a \
	-Wl,--end-group ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a \
	-ldl -lpthread -lm -fopenmp
CFLAGS+=-DMKL
endif


SRC = 	perf/perf.c \
	error.c \
	util.c \
	getrf.c

test_SRC = main.c $(SRC)
perf_SRC = perf_driver.c $(SRC)

test_OBJ=$(test_SRC:.c=.o)
perf_OBJ=$(perf_SRC:.c=.o)

DEP=$(SRC:.c=.d) main.d perf_driver.d

all: $(TARGET)

-include $(DEP)

dgetrf: $(test_OBJ)
	$(CC) $^ -o $@ $(LDFLAGS)

perf_driver: $(perf_OBJ)
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

