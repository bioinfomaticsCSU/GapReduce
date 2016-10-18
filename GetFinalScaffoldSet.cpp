#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "ScaffoldSet.h"

using namespace std;

int main(int argc, char *argv[]){
    
    ContigSetHead * contigSetHead = GetContigSet(argv[1]);
    ContigSetHead * gapSetHead = GetContigSet(argv[2]);
    cout<<contigSetHead->contigCount<<"--"<<gapSetHead->contigCount<<endl;
                  
    FILE * fp;
    if((fp = fopen(argv[3], "r")) == NULL){
        printf("%s, does not exist!", argv[3]);
        exit(0);
    }
    
    FILE * fpScaffold;
    fpScaffold = fopen(argv[4], "w");
    
    
    long int maxSize = 100;
    char * line = (char*)malloc(sizeof(char)*maxSize);
    long int startIndex = 0;
    long int gapIndex = 0;
    long int scaffoldIndex = 0;
    while((fgets(line, maxSize, fp)) != NULL){ 
        char * p = NULL;
        p = strtok(line, "\t");
        long int endIndex = atoi(p);
        p = strtok(NULL, "\t");
        long int startN = atoi(p);
        p = strtok(NULL, "\t");
        long int endN = atoi(p);
        cout<<endIndex<<":"<<endl;
        fprintf(fpScaffold, ">scaffold_%ld\n", scaffoldIndex);
        if(startN == 1){
            fprintf(fpScaffold, "%s", gapSetHead->contigSet[gapIndex].contig);
            gapIndex++;
        }
        for(long int i = startIndex; i <= endIndex; i++){
            cout<<i<<"--";
            fprintf(fpScaffold, "%s", contigSetHead->contigSet[i].contig);
            if(i != endIndex){
                fprintf(fpScaffold, "%s", gapSetHead->contigSet[gapIndex].contig);
                gapIndex++;
            }
        }
        if(endN == 1){
            fprintf(fpScaffold, "%s", gapSetHead->contigSet[gapIndex].contig);
            gapIndex++;
        }
        cout<<endl;
        fprintf(fpScaffold, "\n");
        startIndex = endIndex + 1;
        scaffoldIndex++;
    }
    
    
    
    
}

