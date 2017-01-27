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
    
    //char * gapRegionFileName = new char[30];
    //strcpy(gapRegionFileName, "fill_gap_region.fa");
    ofstream ocout;
    ocout.open(argv[11], ios::app);
    
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
    
    //char * contigSetFile = new char[20];
    //strcpy(contigSetFile, "contig_set.fa");
    //cout<<"aa"<<endl;
    char * leftContig = GetContigFromContigSet(argv[12], atoi(argv[8]));
    //cout<<"bb"<<endl;
    char * rightContig = GetContigFromContigSet(argv[12], atoi(argv[9]));
    
    ReadSetHead * readSetHead = GetReadSetHead(argv[1], argv[2], atoi(argv[5]));
    readSetHead->insertSize = atoi(argv[3]);
    readSetHead->std = atoi(argv[4]);
    
    //cout<<"cc--"<<atoi(argv[5])<<endl;
    
    long int step = atoi(argv[13]);
    long int gapDistance = atoi(argv[10]);
    char * finalGapRegion = NULL;
    //cout<<"leftContig--:"<<leftContig<<endl;
    //cout<<"rightContig--:"<<rightContig<<endl;
    //cout<<"startfilling a gap:"<<endl;
    for(long int fequency = atoi(argv[14]); fequency >= 0; fequency--){
        double ratio = 0;
        for(long int kmerLength = atoi(argv[6]); kmerLength >= atoi(argv[7]); kmerLength = kmerLength - step){
            if(kmerLength%2 == 0){
                kmerLength = kmerLength + 1;
            }
            ratio = 0;
            KmerSetHead * kmerSetHead = GetKmerSetHead(readSetHead, kmerLength, fequency);
            DBGraphHead * deBruijnGraphHead = CreateDBGraphHead(readSetHead, kmerSetHead);
            char * gapRegion = NULL;
            if(deBruijnGraphHead->deBruijnGraph == NULL){
                gapRegion = (char *)malloc(sizeof(char)*(atoi(argv[10]) + 1));
                for(long int i = 0; i < atoi(argv[10]); i++){
                    gapRegion[i] = 'N';
                }
                gapRegion[atoi(argv[10])] = '\0';
            }else{
                gapRegion = ExtractPathFromGraph(leftContig, rightContig, gapDistance, deBruijnGraphHead, readSetHead, kmerSetHead, NULL);
            }
            
            ratio = RatioOfN(gapRegion);
            if(ratio == 0){
                if(finalGapRegion == NULL){
                    finalGapRegion = gapRegion;
                }
                break;
            }else if(ratio != 1){
                leftContig = MergeGapToContigLeft(leftContig, gapRegion);
                rightContig = MergeGapToContigRight(rightContig, gapRegion);
                finalGapRegion = gapRegionToFinal(finalGapRegion, gapRegion, gapDistance);
                //cout<<"newGapDistance:"<<gapDistance<<endl;
            }else{
                if(finalGapRegion == NULL){
                    finalGapRegion = gapRegion;
                }
            }
            //cout<<"tempFinal:"<<finalGapRegion<<endl;
        }
        if(ratio == 0){
            break;
        }
        
    }
    //cout<<"endingfilling a gap:"<<endl;
    ocout<<">gap--"<<argv[1]<<"--"<<atoi(argv[10])<<endl;
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
