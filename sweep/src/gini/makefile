OBJECTS = main.o  
INC = -Wall -I/usr/local/include
LIBS = -L/usr/local/lib -lm -lgsl -lgslcblas -L/usr/local/boost -lboost_filesystem -lboost_system -lboost_timer
CC = g++

all: gini

%.o: %.cpp
	$(CC) -c $(INC) $< -o $@

gini: $(OBJECTS)
	$(CC) $(OBJECTS) $(LIBS) -o $@

clean:
	-rm -f $(OBJECTS) gini

