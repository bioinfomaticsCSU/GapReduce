#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "getcontigcoordinate.h"

using namespace std;

int main(int argc, char *argv[]){
    
    GapContigIndexHead * gapContigIndexHead = GetGapContigIndexHead(argv[1]);
    
    cout<<"aa"<<endl;
    
    SubContigCoordinateSetHead * subContigCoordinateSetHead = GetSubContigCoordinateSetHeadFromCoordinateFile(argv[2]);
    long int contigNumber = 0;
    cout<<"bb"<<endl;
    ContigSet * contigSet = GetContigSet(argv[3], contigNumber, subContigCoordinateSetHead);
    
    ReferenceCoordinateSetHead * referenceCoordinateSetHead = ConstructReferenceCoordinate(argv[4], subContigCoordinateSetHead);
    cout<<"cc"<<endl;
    OutPutGapInScaffold(gapContigIndexHead, referenceCoordinateSetHead, subContigCoordinateSetHead, contigSet, contigNumber, argv[5]);

}

