#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include <pthread.h>
#include <time.h>
#include "dataSet.h"
#include "graph.h"
#include "fillGap.h"

using namespace std;
 
int main(int argc, char *argv[])
{  
    
    char * gapRegionFileName = new char[30];
    strcpy(gapRegionFileName, "fill_gap_region.fa");
    ofstream ocout;
    ocout.open(gapRegionFileName, ios::app);
    
    long int size = file_size(argv[1]) + file_size(argv[2]);
    if(size <= 0){
        ocout<<">gaptt--"<<argv[1]<<endl;
        char * gapRegionNULL = (char *)malloc(sizeof(char)*(atoi(argv[10]) + 1));
        for(long int i = 0; i < atoi(argv[10]); i++){
            gapRegionNULL[i] = 'N';
        }
        gapRegionNULL[atoi(argv[10])] = '\0';
        ocout<<gapRegionNULL<<endl;
        free(gapRegionNULL);
        gapRegionNULL = NULL;
        return -1;
    }
    
    char * contigSetFile = new char[20];
    strcpy(contigSetFile, "contig_set.fa");
    cout<<"aa"<<endl;
    char * leftContig = GetContigFromContigSet(contigSetFile, atoi(argv[8]));
    cout<<"bb"<<endl;
    char * rightContig = GetContigFromContigSet(contigSetFile, atoi(argv[9]));
    
    ReadSetHead * readSetHead = GetReadSetHead(argv[1], argv[2], atoi(argv[5]));
    readSetHead->insertSize = atoi(argv[3]);
    readSetHead->std = atoi(argv[4]);
    
    cout<<"sdsd"<<endl;
    //cout<<leftContig<<endl;
    //cout<<rightContig<<endl;
    cout<<"sdsd--"<<endl;
    KmerSetHead * kmerSetHead = GetKmerSetHead(readSetHead, atoi(argv[7]), 1);
    KmerSetHead * largeKmerSetHead = GetKmerSetHead(readSetHead, atoi(argv[6]), 1);
    DBGraphHead * deBruijnGraphHead = CreateDBGraphHead(readSetHead, kmerSetHead);
    char * gapRegion = NULL;
    if(deBruijnGraphHead->deBruijnGraph == NULL){
        gapRegion = (char *)malloc(sizeof(char)*(atoi(argv[10]) + 1));
        for(long int i = 0; i < atoi(argv[10]); i++){
            gapRegion[i] = 'N';
        }
    }else{
        gapRegion = ExtractPathFromGraph(leftContig, rightContig, atoi(argv[10]), deBruijnGraphHead, readSetHead, kmerSetHead, largeKmerSetHead);
    }
    cout<<"large-------------kmer!"<<endl;
    double ratio = RatioOfN(gapRegion);
    char * gapRegionOfSmallFrequency = NULL;
    char * finalGapRegion = NULL;
    cout<<"large-------------kmer!--"<<ratio<<endl;
    if(ratio > 0){
        KmerSetHead * kmerSetHeadOfSmallFrequency = GetKmerSetHead(readSetHead, atoi(argv[7]), 0);
        KmerSetHead * largeKmerSetHeadOfSmallFrequency = GetKmerSetHead(readSetHead, atoi(argv[6]), 0);
        DBGraphHead * deBruijnGraphHeadOfSmallFrequency = CreateDBGraphHead(readSetHead, kmerSetHeadOfSmallFrequency);
        if(deBruijnGraphHeadOfSmallFrequency->deBruijnGraph == NULL){
            gapRegionOfSmallFrequency = (char *)malloc(sizeof(char)*(atoi(argv[10]) + 1));
            for(long int i = 0; i < atoi(argv[10]); i++){
                gapRegionOfSmallFrequency[i] = 'N';
            }
        }else{
            gapRegionOfSmallFrequency = ExtractPathFromGraph(leftContig, rightContig, atoi(argv[10]), deBruijnGraphHeadOfSmallFrequency, readSetHead, kmerSetHeadOfSmallFrequency, largeKmerSetHeadOfSmallFrequency);
        }
        
        if(RatioOfN(gapRegionOfSmallFrequency) <= ratio){
            finalGapRegion = gapRegionOfSmallFrequency;
        }else{
            finalGapRegion = gapRegion;
        }
    }else{
        finalGapRegion = gapRegion;
    }
    
    
    ocout<<">gap--"<<argv[1]<<endl;
    ocout<<finalGapRegion<<endl;
    
    /*
    char * file = new char[100];
    long int pid = getpid();
    sprintf(file,"/proc/%ld/status",pid);
    ifstream icin;
    icin.open(file);
    char * line = new char[1000];
    while(icin.getline(line,1000)){
        cout<<line<<endl;
    }
    */
}
