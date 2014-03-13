OBJECTS = main.o integrate.o quintessence.o basicdump.o lambdaCDM.o simplecheck.o model.o
INC = -I/usr/local/include
LIBS = -L/usr/local/lib -lm -lgsl -lgslcblas -L/usr/local/boost -lboost_filesystem -lboost_system -lboost_timer
CC = g++

all: main

%.o: ../%.cpp
	$(CC) -c $(INC) $< -o $@

main: $(OBJECTS)
	$(CC) $(OBJECTS) $(LIBS) -o ../$@

clean:
	-rm -f $(OBJECTS) ../main

