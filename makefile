CC=g++

CPPFLAGS = -g -Wall -O3

splitScaffoldSet:	SplitScaffoldToContig.o ScaffoldSet.o 
	$(CC) -o $@ $^ -I $(BAMTOOLS_HOME_INCLUDE)/ $(BAMTOOLS_HOME_LIB)/libbamtools.a -lz
	
fillGapRead:	FillGapRead.o GapRead.o ScaffoldSet.o 
	$(CC) -o $@ $^ -lm -ldl -I $(BAMTOOLS_HOME_INCLUDE)/ $(BAMTOOLS_HOME_LIB)/libbamtools.a -lz

SplitScaffoldToContig.o: SplitScaffoldToContig.cpp ScaffoldSet.h 
	$(CC) -c SplitScaffoldToContig.cpp

FillGapRead.o: FillGapRead.cpp GapRead.h ScaffoldSet.h 
	$(CC) -c FillGapRead.cpp -lm -ldl

GapRead.o: GapRead.cpp ScaffoldSet.h
	$(CC) -c GapRead.cpp -lm -ldl
	
ScaffoldSet.o: ScaffoldSet.cpp ScaffoldSet.h
	$(CC) -c ScaffoldSet.cpp

GetFinalScaffoldFill:	GetFinalScaffoldSet.o ScaffoldSet.o
	$(CC) -o $@ $^ -I $(BAMTOOLS_HOME_INCLUDE)/ $(BAMTOOLS_HOME_LIB)/libbamtools.a -lz
	
GetFinalScaffoldSet.o:	GetFinalScaffoldSet.cpp ScaffoldSet.h
	$(CC) -c GetFinalScaffoldSet.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ $(BAMTOOLS_HOME_LIB)/libbamtools.a -lz
	
all: splitScaffoldSet fillGapRead GetFinalScaffoldFill
	make all -C graph/
	cp graph/fillGapRegion ./
	rm -f *.o
clean:
	rm -f *.o
	rm splitScaffoldSet fillGapRead fillGapRegion GetFinalScaffoldFill

