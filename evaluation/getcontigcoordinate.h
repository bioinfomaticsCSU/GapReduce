#ifndef GETCONTIGCOORDINATE_H_INCLUDED 
#define GETCONTIGCOORDINATE_H_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

typedef struct GapContigIndex{
    long int gapDistance;
    long int leftContigIndex;
    long int rightContigIndex;
}GapContigIndex;

typedef struct GapContigIndexHead{
    GapContigIndex * gapContigIndex;
    long int gapCount;
}GapContigIndexHead;



typedef struct SubContigCoordinate{
    char * contigName;
    char * referenceName;
    long int referenceIndex;
    long int contigIndex;
    long int referenceStart;
    long int referenceEnd;
    long int contigStart;
    long int contigEnd;
    long int referenceOverlapLength;
    long int contigOverlapLength;
    double identityPercent;
    bool isReverse;
    SubContigCoordinate(){
        referenceName = NULL;
        referenceIndex = -1;
        contigName = NULL;
        contigIndex = -1;
        referenceStart = 0;
        referenceEnd = 0;
        contigStart = 0;
        contigEnd = 0;
        referenceOverlapLength = 0;
        contigOverlapLength = 0;
        identityPercent = 0;
        isReverse = false;
    }
}SubContigCoordinate;

typedef struct SubContigCoordinateSetHead{
    SubContigCoordinate * subContigCoordinateSet;
    long int subContigNumber;
    SubContigCoordinateSetHead(){
        subContigCoordinateSet= NULL;
        subContigNumber = 0;
    }
}SubContigCoordinateSetHead;

typedef struct ReferenceCoordinate{
    char * referenceName;
    char * reference;
    long int referenceLength;
    long int * subContigCoordinateIndexSet;
    long int subContigNumber;
    ReferenceCoordinate(){
        referenceName = NULL;
        reference = NULL;
        referenceLength = 0;
        subContigNumber = 0;
        subContigCoordinateIndexSet = NULL;
    }
}ReferenceCoordinate;

typedef struct ReferenceCoordinateSetHead{
    ReferenceCoordinate * referenceCoordinateSet;
    long int referenceNumber;
    ReferenceCoordinateSetHead(){
        referenceNumber = 0;
        referenceCoordinateSet = NULL;
    }
}ReferenceCoordinateSetHead;

typedef struct ContigSet{
    char * contig;
    long int contigLength;
    char * contigName;
    long int * subContigCoordinateIndexSet;
    long int subContigNumber;
    ContigSet(){
        contig = NULL;
        contigLength = 0;
        contigName = NULL;
        subContigCoordinateIndexSet = NULL;
        subContigNumber = 0;
    }

}ContigSet;


GapContigIndexHead * GetGapContigIndexHead(char * gapContigIndexFile);
ContigSet * GetContigSet(char * contigFileName, long int & contigNumber, SubContigCoordinateSetHead * subContigCoordinateSetHead);
SubContigCoordinateSetHead * GetSubContigCoordinateSetHeadFromCoordinateFile(char * coordinateFileName);
long int SearchContigIndexFromContigName(ContigSet * contigSet, long int contigNumber, char * contigName);
long int SearchReferenceIndexFromContigName(ReferenceCoordinateSetHead * referenceCoordinateSetHead, char * referenceName);

void OutputSubContigSet(ContigSet * contigSet, long int contigNumber, SubContigCoordinateSetHead * subContigCoordinateSetHead);

SubContigCoordinateSetHead * GetSubContigCoordinateSetHeadFromCoordinateFile(char * coordinateFileName);
ReferenceCoordinateSetHead * ConstructReferenceCoordinate(char * referenceFileName, SubContigCoordinateSetHead * subContigCoordinateSetHead);
void OutputContigSetCoordinate(ContigSet * contigSet, SubContigCoordinateSetHead * subContigCoordinateSetHead, long int contigNumber);
void OutputReferenceCoordinateSet(ReferenceCoordinateSetHead * referenceCoordinateSetHead, SubContigCoordinateSetHead * subContigCoordinateSetHead);

void GetMinGapDistanceIndex(SubContigCoordinate * subContigCoordinateSet, long int * leftIndexSet, long int leftCount, long int * rightIndexSet, long int rightCount, long int & leftIndex, long int & rightIndex, long int gapDistance);
void OutPutGapInScaffold(GapContigIndexHead * gapContigIndexHead, ReferenceCoordinateSetHead * referenceCoordinateSetHead, SubContigCoordinateSetHead * subContigCoordinateSetHead, ContigSet * contigSet, long int contigNumber, char * gapFileName);
bool ReverseComplement(char * temp1, char * temp2);

#endif
