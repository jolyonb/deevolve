OBJECTS = main.o model.o checkpathology.o write.o evolve.o initialise.o
INC = -Wall -I/usr/local/include
LIBS = -L/usr/local/lib -lm -lgsl -lgslcblas -L/usr/local/boost -lboost_filesystem -lboost_system -lboost_timer
CC = g++

all: sfevol

%.o: source/%.cpp
	$(CC) -c $(INC) $< -o $@

sfevol: $(OBJECTS)
	$(CC) $(OBJECTS) $(LIBS) -o $@

clean:
	-rm -f $(OBJECTS) main

