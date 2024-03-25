CC = g++

SRC = $(sort $(wildcard *.cc))
OFS = $(SRC:%.cc=%.o)
DEP = $(SRC:%.cc=%.d)
OBJ = $(OFS)

BFS = main_having

AFS = $(OFS) $(BFS)

ALL = $(AFS)

all : $(ALL)

-include $(DEP)

main_having : $(OBJ)
	$(CC) -o $@ $^

%.o : %.cc
	$(CC) -O3 -c -MMD -std=c++20 -o $@ $*.cc

.PHONY : clean
clean :
	rm -f *.o *.d *.gch a.out $(ALL)

