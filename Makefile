CC=gcc
CFLAGS= -std=c99 -O3 -fPIC

SOURCES=two_planet_perturbation.c kepler_solve.c LaplaceCoefficients.c ActionAngleIntegrate.c
OBJECTS=$(SOURCES:.c=.o)

all: $(OBJECTS)
	$(CC) $(CFLAGS) -lm  -shared -o libsemianalyticMEGNO.so $^
testme: testme.c $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
ActionAngleIntegrate.o: ActionAngleIntegrate.c LaplaceCoefficients.o
	$(CC) $(CFLAGS) -c $< -o $@
clean:
	rm -rf *.o *.so
run: testme
	./testme > data.txt
resonances_test: resonances_test.c LaplaceCoefficients.o
	$(CC) -std=c99 -o $@ $^
action_angle_test: action_angle_test.c ActionAngleIntegrate.o LaplaceCoefficients.o 
	$(CC) -std=c99 -lm -o  $@ $^
laplace_test: laplace_test.c LaplaceCoefficients.o
	$(CC) -std=c99 -lm -o  $@ $^
secular_test: secular_test.c secular_evolution.o LaplaceCoefficients.o ActionAngleIntegrate.o
	$(CC) -std=c99 -lm  -o  $@ $^
