OBJECTS = main.o  
INC = -Wall -I/usr/local/include
LIBS = -L/usr/local/lib -lm -lgsl -lgslcblas -L/usr/local/boost -lboost_filesystem -lboost_system -lboost_timer
CC = g++

all: mcmc

%.o: src/%.cpp
	$(CC) -c $(INC) $< -o $@

mcmc: $(OBJECTS)
	$(CC) $(OBJECTS) $(LIBS) -o $@

clean:
	-rm -f $(OBJECTS) mcmc

