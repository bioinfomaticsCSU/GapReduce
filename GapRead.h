#ifndef FILLGAPREAD_H_INCLUDED 
#define FILLGAPREAD_H_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

//#include "filler.h"

#include "ScaffoldSet.h"

using namespace std;

typedef struct SingleRead{
    long int readIndex;
    char * read;
    bool orientation;
    long int contigPosition;
    long int contigIndex;
}SingleRead;

typedef struct PairedRead{
    long int readIndex;
    char * leftRead;
    bool leftOrientation;
    long int leftContigPosition;
    long int leftContigIndex;
    char * rightRead;
    bool rightOrientation;
    long int rightContigPosition;
    long int rightContigIndex;
}PairedRead;

typedef struct GapRegionSingleMapReadSet{
    SingleRead * leftContigSingleReadSet;
    SingleRead * rightContigSingleReadSet;
    long int leftContigSingleReadCount;
    long int rightContigSingleReadCount;
}GapRegionSingleMapReadSet;

typedef struct GapRegionPairedMapReadSet{
    PairedRead * pairedReadSet;
    long int pairedReadCount;
}GapRegionPairedMapReadSet;

typedef struct GapRegionReadSet{
    GapRegionSingleMapReadSet * singleMapReadSet;
    GapRegionPairedMapReadSet * pairedMapReadSet;
    long int * gapToContigIndex;
    long int gapCount;
}GapRegionReadSet;


GapRegionReadSet * GetReadSetInGapRegion(ScaffoldSetHead * scaffoldSetHead, char * endContigFile, char * endContigBamLeft, char * endContigBamRight, long int maxInsertSize, long int minInsertSize, long int intervalLength, bool isPairedRead);
void WriteGapRegionReadSet(GapRegionReadSet * gapRegionReadSet, char * fileName);
void FillSingleGapByRead(ScaffoldSetHead * scaffoldSetHead, GapRegionReadSet * gapRegionReadSet, long int libraryIndex);

int KMPIndexOfContigOfMisMatch(char * contig, char * pattern);
void SearchSharedKmerRead(ScaffoldSetHead * scaffoldSetHead, GapRegionReadSet * gapRegionReadSet);
bool ReverseComplement(char * temp1, char * temp2);

#endif
