#ifndef GRAPH_H_INCLUDED 
#define GRAPH_H_INCLUDED 
#include "fstream"
#include <cstring>
#include <iostream>
#include "stdio.h"
#include "stdlib.h"
#include <pthread.h>

#include "dataSet.h"

using namespace std;

typedef struct GraphNode{
    long int index;
    struct GraphNode * next;
    GraphNode(){
        index = -1;
        next = NULL;
    }
}GraphNode;

typedef struct DBGraph{
    char * contig;
    GraphNode * outNode;
    GraphNode * inNode;
    DBGraph(){
        contig = NULL;
        outNode = NULL;
        inNode = NULL;
    }
}DBGraph;

typedef struct DBGraphHead{
    DBGraph * deBruijnGraph;
    unsigned long int nodeNumber;
    DBGraphHead(){
        deBruijnGraph = NULL;
        nodeNumber = 0;
    }
}DBGraphHead;


int NodeCount(GraphNode * temp);
int DBGSearchNode(GraphNode * tempGraphNode, long int num, long int newNum);
int DBGSearchNode(GraphNode * tempGraphNode, long int num);
int DBGRemoveInNode(DBGraph * graph, long int num);
int DBGRemoveOutNode(DBGraph * graph, long int num);
void AdjustDBG(DBGraph *deBruijnGraph, int kmerLength,long int max);
void SimplePathMerge(DBGraph *deBruijnGraph, int kmerLength, long int max);
void RemoveTip(DBGraph *deBruijnGraph, int kmerLength, long int graphNodeNumber);
void MergeBubble(DBGraph *deBruijnGraph, int kmerLength, long int graphNodeNumber);
void RemoveCycle(DBGraph *deBruijnGraph, int kmerLength, long int graphNodeNumber, KmerSetHead * kmerSetHead);
void RemoveCycle1(DBGraph *deBruijnGraph, int kmerLength, long int graphNodeNumber);
DBGraph * OptimizeDBGraphSpace(DBGraph * deBruijnGraph, unsigned long int & graphNodeCount);
void GetDBGraphFromAddress(DBGraphHead * deBruijnGraphHead, char * str);
long int FindKmerInDeBruijnGraph(char * kmer, long int kmerLength, DBGraph *deBruijnGraph, long int graphNodeCount, long int & mapPosition);



void InsertLeftMateNode(DBGraph * deBruijnGraph, long int right, long int left);
void InsertRightMateNode(DBGraph * deBruijnGraph, long int left, long int right);
long int DBGraphInsertNode(long int previousPos, char * kmer, int kmerLength, DBGraph * deBruijnGraph, unsigned long int graphHashCount);
long int CreateDBGraphFromReadSet(DBGraphHead * deBruijnGraphHead, ReadSet * readSetHead, KmerSetHead * kmerSetHead, bool isLeftRead);
DBGraphHead * CreateDBGraphHead(ReadSetHead * readSetHead, KmerSetHead * kmerSetHead);

void WriteDBGraph(DBGraph * deBruijnGraph, long int count, char * str);






#endif
