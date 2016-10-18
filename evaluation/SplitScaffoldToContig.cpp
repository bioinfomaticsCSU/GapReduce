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
    
    ScaffoldSetHead * scaffoldSetHead = GetScaffoldSetFromScaffoldFile(argv[1]);
    long int minGapDistance = atoi(argv[2]);
    SplitScaffoldToContig(scaffoldSetHead, minGapDistance, argv[3]);
    WriteContigSetFromScffoldSet(scaffoldSetHead,argv[4]); 
}

