#ifndef FILLGAP_H_INCLUDED
#define FILLGAP_H_INCLUDED

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dataSet.h"
#include "graph.h"

using namespace std;

typedef struct Contig{
    char * contig;
    Contig(){
        contig = NULL;
    }
}Contig;

typedef struct ContigSet{
    char * contig;
    struct ContigSet * next;
    ContigSet(){
        contig = NULL;
        next = NULL;
    }
}ContigSet;


char * GetContigFromContigSet(char * contigSetFile, long int contigIndex);
double ScoreTtest(KmerReadIndex * readIndex, ReadSet * readSet, long int mean, long int std, long int gap, bool isLeft);
double ScoreDistribution(long int * distance, long int num, long int mean, long int std);
int ScoreKmerCount(long int kmerCount[][2], double score[][2], long int rowCount);
char * TranverseGraphFromNode(DBGraphHead * deBruijnGraphHead, ReadSetHead * readSetHead, KmerSetHead * largeKmerSetHead, KmerSetHead * kmerSetHead, char * startContig, long int startNodeIndex, long int gapDistance, long int & extendLength, bool isOut);
char * GetGapRegion(char * leftGap, long int leftExtendLength, char * rightGap, long int rightExtendLength, long int gapDistance, long int cutLength);
char * ExtractPathFromGraph(char * leftContig, char * rightContig, long int gapDistance, DBGraphHead * deBruijnGraphHead, ReadSetHead * readSetHead, KmerSetHead * kmerSetHead, KmerSetHead * largeKmerSetHead);

bool SingleKmerReadConsencusContig(char * read, char * kmer, char * contig, bool isOut);
bool KmerConsencus(ReadSet * readSet, long int readLength, char * kmer, KmerReadIndex * temp, char * contig, bool isOut);

#endif
