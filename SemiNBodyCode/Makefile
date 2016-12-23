CC=gcc
CFLAGS= -std=c99 -O3 -fPIC

SOURCES=two_planet_perturbation.c kepler_solve.c
OBJECTS=$(SOURCES:.c=.o)

all: $(OBJECTS)
	$(CC) $(CFLAGS) -lm  -shared -o libSemiNbody.so $^
testme: testme.c $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
clean:
	rm -rf *.o *.so
run: testme
	./testme > data.txt
