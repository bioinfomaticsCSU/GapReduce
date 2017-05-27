CC=g++

CPPFLAGS = -g -Wall -O3

splitScaffoldSet:	SplitScaffoldToContig.o ScaffoldSet.o 
	$(CC) -o $@ $^ -std=gnu++98
	
fillGapRead:	FillGapRead.o GapRead.o ScaffoldSet.o 
	$(CC) -o $@ $^ -lm -ldl -I include/ -L. -lbamtools -std=gnu++98

SplitScaffoldToContig.o: SplitScaffoldToContig.cpp ScaffoldSet.h 
	$(CC) -c SplitScaffoldToContig.cpp -std=gnu++98

FillGapRead.o: FillGapRead.cpp GapRead.h ScaffoldSet.h 
	$(CC) -c FillGapRead.cpp -lm -ldl -I include/ -L. -lbamtools -std=gnu++98

GapRead.o: GapRead.cpp ScaffoldSet.h
	$(CC) -c GapRead.cpp -lm -ldl -I include/ -L. -lbamtools -std=gnu++98
	
ScaffoldSet.o: ScaffoldSet.cpp ScaffoldSet.h
	$(CC) -c ScaffoldSet.cpp -std=gnu++98

GetFinalScaffoldFill:	GetFinalScaffoldSet.o ScaffoldSet.o
	$(CC) -o $@ $^ -std=gnu++98
	
GetFinalScaffoldSet.o:	GetFinalScaffoldSet.cpp ScaffoldSet.h
	$(CC) -c GetFinalScaffoldSet.cpp -std=gnu++98
	
all: splitScaffoldSet fillGapRead GetFinalScaffoldFill
	make all -C graph/
	cp graph/fillGapRegion ./
	rm -f *.o
clean:
	rm -f *.o
	rm splitScaffoldSet fillGapRead fillGapRegion GetFinalScaffoldFill

