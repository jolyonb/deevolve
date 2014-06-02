# Compilation flags

OPTIMIZE = -O3 -march=native -flto
WARNINGS = -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable

CC = g++

CFLAGS = $(OPTIMIZE) $(WARNINGS) -std=gnu++03
INC = -I/usr/local/include
LIBS = -L/usr/local/lib -lm -lgsl -lgslcblas -L/usr/local/boost -lboost_filesystem -lboost_system -lboost_timer

# Make rules

OBJECTS = main.o integrate.o quintessence.o basicdump.o lambdaCDM.o simplecheck.o model.o linearw.o kessence.o process.o KGB.o params.o evolve.o fr.o inireader.o

.PHONY: all
all: main

main: $(OBJECTS)
	$(CC) $(OBJECTS) $(LIBS) -o ../$@

clean:
	-rm -rf $(OBJECTS) ../main .deps

ifeq ($(USER),jolyon)
# Fancy dependency stuff: stores the dependencies of each source file, so that if headers are updated, it knows which to recompile

DEPDIR = .deps

%.o : ../source/%.cpp
	@mkdir -p .deps && $(CC) -MM $(CFLAGS) $< -MP -MF $(DEPDIR)/$*.P
	$(COMPILE.c) -o $@ -c $(INC) $<

-include $(OBJECTS:%.o=$(DEPDIR)/%.P)
else
%.o: ../source/%.cpp
	$(CC) -c $(INC) $< -o $@
endif
