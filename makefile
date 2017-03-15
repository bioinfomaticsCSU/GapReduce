CC=g++

CPPFLAGS = -g -Wall -O3

splitScaffoldSet:	SplitScaffoldToContig.o ScaffoldSet.o
	$(CC) -o $@ $^
	
fillGapRead:	FillGapRead.o GapRead.o ScaffoldSet.o
	$(CC) -o $@ $^ -lm -ldl -I include/ -L lib/ -lbamtools

SplitScaffoldToContig.o: SplitScaffoldToContig.cpp ScaffoldSet.h
	$(CC) -c SplitScaffoldToContig.cpp

FillGapRead.o: FillGapRead.cpp GapRead.h ScaffoldSet.h
	$(CC) -c FillGapRead.cpp -lm -ldl -I include/ -L lib/ -lbamtools

GapRead.o: GapRead.cpp ScaffoldSet.h
	$(CC) -c GapRead.cpp -lm -ldl -I include/ -L lib/ -lbamtools
	
ScaffoldSet.o: ScaffoldSet.cpp ScaffoldSet.h
	$(CC) -c ScaffoldSet.cpp

GetFinalScaffoldFill:	GetFinalScaffoldSet.o ScaffoldSet.o
	$(CC) -o $@ $^
	
GetFinalScaffoldSet.o:	GetFinalScaffoldSet.cpp ScaffoldSet.h
	$(CC) -c GetFinalScaffoldSet.cpp
	
all: splitScaffoldSet fillGapRead GetFinalScaffoldFill
	make all -C graph/
	cp graph/fillGapRegion ./
	rm -f *.o
clean:
	rm -f *.o
	rm splitScaffoldSet fillGapRead fillGapRegion GetFinalScaffoldFill

