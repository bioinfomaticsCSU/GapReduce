#ifndef DATASET_H_INCLUDED 
#define DATASET_H_INCLUDED 
#include "fstream"
#include <cstring>
#include <iostream>
#include "stdio.h"
#include "stdlib.h"
#include <pthread.h>

using namespace std;

typedef struct ReadSet{
    char * read;
    long int insertSize;
}ReadSet;

typedef struct ReadSetHead{
    long int readLength;
    long int insertSize;
    long int std;
    ReadSet * leftReadSet;
    ReadSet * rightReadSet;
    long int leftReadCount;
    long int rightReadCount;
    bool isPaired;
}ReadSetHead; 

typedef struct KmerReadIndex{
    long int index;
    KmerReadIndex * next;
}KmerReadIndex;

typedef struct KmerSet{
    char * kmer;
    long int kmerCount;
    KmerReadIndex * readIndex;
}KmerSet;

typedef struct KmerSetHead{
    long int kmerLength;
    long int minKmerFrequency;
    KmerSet * leftKmerSet;
    KmerSet * rightKmerSet;
    long int leftKmerCount;
    long int rightKmerCount;
    
}KmerSetHead;


double RatioOfN(char * contig);
long int file_size(char* filename);
void CopyReadofATGC(char * destination, char * source, long int length);
void AppendRight( char * temp, char * temp1, char * temp2,long int kmerLength);
int KMPIndexOfContig(char * contig, char * pattern);
bool ReverseComplement(char * temp1, char * temp2);
bool ReverseComplement(char * temp);
bool ReverseComplementReadSet(ReadSet * readSet, long int readCount);

ReadSet * GetReadSet(char * readAddress, long int & readCount);
ReadSetHead * GetReadSetHead(char * leftReadAddress, char * rightReadAddress, bool isPaired);


unsigned long int Hash(char * str, unsigned int len, unsigned long int max);
long int SearchKmerOfKmerSet(char * kmer, KmerSet * kmerSet, long int kmerLength, long int kmerCount);
void InsertKmerToKmerSet(char * kmer, long int readIndex, KmerSet * kmerSet, long int kmerLength, long int kmerCount);
KmerSet * GetKmerSet(ReadSet * readSet, long int readCount, long int readLength, long int kmerLength, long int & kmerCount);
KmerSetHead * GetKmerSetHead(ReadSetHead * readSetHead, long int kmerLength,long int minKmerFrequency);


#endif
