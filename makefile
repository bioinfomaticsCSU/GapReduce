CC=g++

CPPFLAGS = -g -Wall -O3

splitScaffoldSet:	SplitScaffoldToContig.o ScaffoldSet.o 
	$(CC) -o $@ $^ -I $(BAMTOOLS_HOME_INCLUDE)/ $(BAMTOOLS_HOME_LIB)/libbamtools.a
	
fillGapRead:	FillGapRead.o GapRead.o ScaffoldSet.o 
	$(CC) -o $@ $^ -lm -ldl -I $(BAMTOOLS_HOME_INCLUDE)/ $(BAMTOOLS_HOME_LIB)/libbamtools.a

SplitScaffoldToContig.o: SplitScaffoldToContig.cpp ScaffoldSet.h 
	$(CC) -c SplitScaffoldToContig.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ $(BAMTOOLS_HOME_LIB)/libbamtools.a

FillGapRead.o: FillGapRead.cpp GapRead.h ScaffoldSet.h 
	$(CC) -c FillGapRead.cpp -lm -ldl -I $(BAMTOOLS_HOME_INCLUDE)/ $(BAMTOOLS_HOME_LIB)/libbamtools.a

GapRead.o: GapRead.cpp ScaffoldSet.h
	$(CC) -c GapRead.cpp -lm -ldl -I $(BAMTOOLS_HOME_INCLUDE)/ $(BAMTOOLS_HOME_LIB)/libbamtools.a
	
ScaffoldSet.o: ScaffoldSet.cpp ScaffoldSet.h
	$(CC) -c ScaffoldSet.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ $(BAMTOOLS_HOME_LIB)/libbamtools.a

GetFinalScaffoldFill:	GetFinalScaffoldSet.o ScaffoldSet.o
	$(CC) -o $@ $^ -I $(BAMTOOLS_HOME_INCLUDE)/ $(BAMTOOLS_HOME_LIB)/libbamtools.a
	
GetFinalScaffoldSet.o:	GetFinalScaffoldSet.cpp ScaffoldSet.h
	$(CC) -c GetFinalScaffoldSet.cpp -I $(BAMTOOLS_HOME_INCLUDE)/ $(BAMTOOLS_HOME_LIB)/libbamtools.a
	
all: splitScaffoldSet fillGapRead GetFinalScaffoldFill
	make all -C graph/
	cp graph/fillGapRegion ./
	rm -f *.o
clean:
	rm -f *.o
	rm splitScaffoldSet fillGapRead fillGapRegion GetFinalScaffoldFill

