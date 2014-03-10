OBJECTS = main.o integrate.o params.o quintessence.o intparams.o basicdump.o
INC = -I/usr/local/include
LIBS = -L/usr/local/lib -lm -lgsl -lgslcblas
CC = g++

all: main

%.o: ../%.cpp
	$(CC) -c $(INC) $< -o $@

main: $(OBJECTS)
	$(CC) $(OBJECTS) $(LIBS) -o $@

clean:
	-rm -f $(OBJECTS) main

