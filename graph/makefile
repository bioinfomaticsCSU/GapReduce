CC=g++

CPPFLAGS = -g -Wall -O3

fillGapRegion:	main.o dataSet.o graph.o fillGap.o
	$(CC) -o $@ $^

main.o: main.cpp dataSet.h graph.h fillGap.h 
	$(CC) -c main.cpp
	
dataSet.o: dataSet.cpp dataSet.h
	$(CC) -c dataSet.cpp
	
graph.o: graph.cpp dataSet.h graph.h
	$(CC) -c graph.cpp
	
fillGap.o: fillGap.h dataSet.h graph.h fillGap.h
	$(CC) -c fillGap.cpp
	
all: fillGapRegion 
	rm -f *.o
clean:
	rm -f *.o
	rm fillGapRegion

