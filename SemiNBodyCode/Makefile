CC=gcc
CFLAGS= -std=c99 -O3 -fPIC -Wall  

SOURCES=two_planet_perturbation.c kepler_solve.c
OBJECTS=$(SOURCES:.c=.o)

all: $(OBJECTS)
	$(CC) $(CFLAGS) -lm  -shared -o libSemiNbody.so $^
%.o: %.c
	$(CC)  -c -o $@ $< $(CFLAGS)
testme: testme.c $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^
clean:
	rm -rf *.o *.so
run: testme
	./testme > data.txt
