#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "ScaffoldSet.h"
#include "GapRead.h"

using namespace std;

int main(int argc, char *argv[]){
    /*
    char * fileName = (char *)malloc(sizeof(char)*30);
    strcpy(fileName, "gapReadSet.fa");
    ofstream ocout;
    ocout.open(fileName);
    */
    ScaffoldSetHead * scaffoldSetHead = GetScaffoldSetFromScaffoldFile(argv[1]);
    
    long int minGapDistance = atoi(argv[3]);
    
    SplitScaffoldToContig(scaffoldSetHead, minGapDistance, NULL);
    //cout<<"aa"<<endl;
   
    GapRegionReadSet * gapRegionReadSet = GetReadSetInGapRegion(scaffoldSetHead, argv[2], argv[4], argv[5], atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), atoi(argv[9]));
    //cout<<"aa--end"<<endl;
    //WriteGapRegionReadSet(gapRegionReadSet, fileName);
    FillSingleGapByRead(scaffoldSetHead, gapRegionReadSet, atoi(argv[10]), argv[11]);
    
    //SearchSharedKmerRead(scaffoldSetHead, gapRegionReadSet);
      
}

