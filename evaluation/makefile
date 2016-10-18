CC=g++

CPPFLAGS = -g -Wall -O3

splitScaffoldSet:	SplitScaffoldToContig.o ScaffoldSet.o
	$(CC) -o $@ $^

SplitScaffoldToContig.o: SplitScaffoldToContig.cpp ScaffoldSet.h
	$(CC) -c SplitScaffoldToContig.cpp
	
ScaffoldSet.o: ScaffoldSet.cpp ScaffoldSet.h
	$(CC) -c ScaffoldSet.cpp

getGapRegionInScaffold:	GetGapRegionInScaffold.o getcontigcoordinate.o
	$(CC) -o $@ $^
	
GetGapRegionInScaffold.o:	GetGapRegionInScaffold.cpp getcontigcoordinate.h
	$(CC) -c GetGapRegionInScaffold.cpp
	
getcontigcoordinate.o:	getcontigcoordinate.cpp getcontigcoordinate.h
	$(CC) -c getcontigcoordinate.cpp
	
getPerfectScaffold:	get_perfect_scaffold.o
	$(CC) -o $@ $^
	
get_perfect_scaffold.o: get_perfect_scaffold.cpp
	$(CC) -c get_perfect_scaffold.cpp
	
all: splitScaffoldSet getGapRegionInScaffold getPerfectScaffold 
	rm -f *.o
clean:
	rm -f *.o
	rm splitScaffoldSet getGapRegionInScaffold getPerfectScaffold

