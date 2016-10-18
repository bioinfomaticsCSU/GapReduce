#ifndef SCAFFOLDSET_H_INCLUDED 
#define SCAFFOLDSET_H_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
using namespace std;

typedef struct MapRegionGapIndex{
    
    long int gapIndex;
    long int distance;
    MapRegionGapIndex * next;
    
    
}MapRegionGapIndex;

typedef struct ContigToGapIndex{
    long int nextGapContigIndex;
    long int leftGapIndex;
    long int rightGapIndex;
}ContigToGapIndex;

typedef struct GapToContigIndex{
    long int leftContigIndex;
    long int rightContigIndex;
    long int gapDistance;
}GapToContigIndex;

typedef struct ScaffoldSet{
    char * scaffold;
    char * scaffoldName;
    long int scaffoldLength;
    long int * contigStartCoordinate;
    long int * contigEndCoordinate;
    long int * contigIndex;
    long int * gapDistance;
    long int contigCount;
    long int gapCount;
    ScaffoldSet(){
        scaffold = NULL;
        scaffoldName = NULL;
        scaffoldLength = 0;
        contigStartCoordinate = NULL;
        contigEndCoordinate = NULL;
        contigIndex = NULL;
        gapDistance = NULL;
        contigCount = 0;
        gapCount = 0;
    }

}ScaffoldSet;

typedef struct ScaffoldSetHead{
    ScaffoldSet * scaffoldSet;
    ContigToGapIndex * contigToGapIndex;
    GapToContigIndex * gapToContigIndex;
    long int scaffoldCount; 
    long int contigCount;
    long int gapCount;
    ScaffoldSetHead(){
        scaffoldSet = NULL;
        contigToGapIndex = NULL;
        scaffoldCount = 0;
        contigCount = 0;
        gapCount = 0;
    }

}ScaffoldSetHead;

typedef struct ContigSet{
    char * contig;
    long int contigLength;
    ContigSet(){
        contig = NULL;
        contigLength = 0;
    }
}ContigSet;

typedef struct ContigSetHead{
    ContigSet * contigSet;
    long int contigCount;
    ContigSetHead(){
        contigSet = NULL;
        contigCount = 0;
    }

}ContigSetHead;

MapRegionGapIndex * GetGapReadFromReadPositionOfContig(ScaffoldSetHead * scaffoldSetHead, long int contigIndex, long int readPosition, long int maxInsetSize, long int minInsertSize, bool isPairedRead, bool isReverse);

ScaffoldSetHead * GetScaffoldSetFromScaffoldFile(char * scaffoldFileName);
void SplitScaffoldToContig(ScaffoldSetHead * scaffoldSetHead, long int minGapDistance, char * gapInformationFile);
void WriteContigSetFromScffoldSet(ScaffoldSetHead * scaffoldSetHead, char * contigSetFileName);
void WriteEndContigFromScffoldSet(ScaffoldSetHead * scaffoldSetHead, long int subContigLength, char * endContigFileName);
ContigSetHead * GetContigSet(char * contigSetFile);
#endif 
