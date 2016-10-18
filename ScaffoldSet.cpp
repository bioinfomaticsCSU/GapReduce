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

long int GetRealMappingPositionFromEndContig(long int position, long int realContigLength, long int endLength, long int intervalLength, long int softLength){
    
    if(realContigLength <= 2*endLength){
        return position - softLength;
    }
    
    if(position < endLength){
        return position - softLength;
        
    }
    if(position >= endLength){
        return realContigLength - 2*endLength - intervalLength + position;
    }
    
}

MapRegionGapIndex * GetGapReadFromReadPositionOfContig(ScaffoldSetHead * scaffoldSetHead, long int readLength, long int contigIndex, long int readPosition, long int maxInsertSize, long int minInsertSize, bool isPairedRead, bool isReverse){
    
    MapRegionGapIndex * head = NULL;
    MapRegionGapIndex * first = NULL;
    
    ScaffoldSet * scaffoldSet = scaffoldSetHead->scaffoldSet;
    long int scaffoldCount = scaffoldSetHead->scaffoldCount;
    //cout<<"count:"<<scaffoldCount<<endl;
    long int scaffoldLength = 0;
    long int scaffoldIndex;
    long int contigStartPosition = 0;
    long int contigEndPosition = 0;
    
    long int startIndex = 0;
    for(long int i = 0; i < scaffoldCount; i++){
        //cout<<"contigCount:"<<scaffoldSet[i].contigCount<<"--"<<contigIndex<<endl;
        for(long int j = 0; j < scaffoldSet[i].contigCount; j++){
            if(scaffoldSet[i].contigIndex[j] == contigIndex){
                //cout<<contigStartPosition<<"--"<<contigEndPosition<<"--"<<scaffoldSet[i].contigIndex[j]<<endl;
                //cout<<scaffoldSet[i].contigStartCoordinate[j]<<endl;
                contigStartPosition = scaffoldSet[i].contigStartCoordinate[j];
                //cout<<scaffoldSet[i].contigEndCoordinate[j]<<endl;
                contigEndPosition = scaffoldSet[i].contigEndCoordinate[j];
                scaffoldLength = scaffoldSet[i].scaffoldLength;
                scaffoldIndex = i;
                startIndex = j;
                //cout<<contigStartPosition<<"--"<<contigEndPosition<<endl;
                break;
            }
        }
    }
    //cout<<"ss"<<endl;
    long int startMapRegion = 0;
    long int endMapRegion = 0;
    long int mapRegionCount = 0;
    
    if((!isReverse && isPairedRead) ||(isReverse && !isPairedRead)){
        
        startMapRegion = contigStartPosition + readPosition + minInsertSize;
        endMapRegion = contigStartPosition + readPosition + maxInsertSize;
        //cout<<"tt--"<<startMapRegion<<"--"<<endMapRegion<<endl;
        for(long int j = startIndex; j < scaffoldSet[scaffoldIndex].contigCount - 1; j++){
            
            if(endMapRegion < scaffoldSet[scaffoldIndex].contigEndCoordinate[j] || startMapRegion > scaffoldSet[scaffoldIndex].contigStartCoordinate[j+1]){
                
            }else{
                MapRegionGapIndex * temp = (MapRegionGapIndex *)malloc(sizeof(MapRegionGapIndex));
                temp->gapIndex = mapRegionCount;
                temp->distance = scaffoldSet[scaffoldIndex].contigEndCoordinate[j] - scaffoldSet[scaffoldIndex].contigStartCoordinate[startIndex] + 1 - readPosition;
                temp->next = NULL;
                if(head == NULL){
                    head = temp;
                    first = head;
                }else{
                    head->next = temp;
                    head = temp;
                }
            }
            mapRegionCount++;
            
        }
        
    }else{
    
        endMapRegion = contigStartPosition + readPosition + readLength - minInsertSize;
        startMapRegion = contigStartPosition + readPosition + readLength - maxInsertSize;
        
        for(long int j = startIndex; j > 0; j--){
            
            if(endMapRegion < scaffoldSet[scaffoldIndex].contigEndCoordinate[j-1] || startMapRegion > scaffoldSet[scaffoldIndex].contigStartCoordinate[j]){

            }else{
                MapRegionGapIndex * temp = (MapRegionGapIndex *)malloc(sizeof(MapRegionGapIndex));
                temp->gapIndex = mapRegionCount;
                temp->distance = scaffoldSet[scaffoldIndex].contigStartCoordinate[startIndex] - scaffoldSet[scaffoldIndex].contigStartCoordinate[j] + readPosition + readLength;
                temp->next = NULL;
                if(head == NULL){
                    head = temp;
                    first = head;
                }else{
                    head->next = temp;
                    head = temp;
                }
            }
            mapRegionCount++;
            
        }
    }
    
    return first;
    
}

ScaffoldSetHead * GetScaffoldSetFromScaffoldFile(char * scaffoldFileName){
    
    long int i = 0;
    long int j = 0;
    
    long int maxSize = 10000;
    char * scaffold = NULL;
    if(NULL == (scaffold = (char*)malloc(sizeof(char)*maxSize))){
        perror("malloc error!");
        exit(1);
    }
    
    long int scaffoldCount = 0;
    
    FILE * fp;
    
    if((fp = fopen(scaffoldFileName, "r")) == NULL){
        printf("%s, does not exist!", scaffoldFileName);
        exit(0);
    }
    
    while((fgets(scaffold, maxSize, fp)) != NULL){ 
       if(scaffold[0] == '>'){  
           scaffoldCount++; 
       }  
    }  
    
    fclose(fp);
    
    ScaffoldSetHead * scaffoldSetHead = NULL;
    if(NULL == (scaffoldSetHead = (ScaffoldSetHead*)malloc(sizeof(ScaffoldSetHead)))){
        perror("ScaffoldSetHead malloc error!");
        exit(1);
    }
    scaffoldSetHead->scaffoldCount = scaffoldCount;
    if(NULL == (scaffoldSetHead->scaffoldSet = (ScaffoldSet*)malloc(sizeof(ScaffoldSet)*scaffoldCount))){
        perror("ScaffoldSet malloc error!");
        exit(1);
    }
    scaffoldSetHead->contigToGapIndex = NULL;
    scaffoldSetHead->contigCount = 0;
    scaffoldSetHead->gapCount = 0;
    
    for(i = 0; i < scaffoldCount; i++){ 
        scaffoldSetHead->scaffoldSet[i].scaffold = NULL;
        scaffoldSetHead->scaffoldSet[i].scaffoldName = NULL;
        scaffoldSetHead->scaffoldSet[i].scaffoldLength = 0;
        scaffoldSetHead->scaffoldSet[i].contigStartCoordinate = NULL;
        scaffoldSetHead->scaffoldSet[i].contigEndCoordinate = NULL;
        scaffoldSetHead->scaffoldSet[i].contigIndex = NULL;
        scaffoldSetHead->scaffoldSet[i].gapDistance = NULL;
        scaffoldSetHead->scaffoldSet[i].contigCount = 0;
        scaffoldSetHead->scaffoldSet[i].gapCount = 0;
    }
    
    if((fp = fopen(scaffoldFileName, "r")) == NULL){
        printf("%s, does not exist!", scaffoldFileName);
        exit(0);
    }
    
    long int scaffoldIndex = -1;
    
    while((fgets(scaffold, maxSize, fp)) != NULL){ 
       
       if(scaffold[0] == '>'){  
           
           if(strlen(scaffold) == maxSize-1){              
               while((fgets(scaffold, maxSize, fp)) != NULL){
                   if(strlen(scaffold) != maxSize-1){
                       break;
                   }
               }        
           }
           scaffoldIndex++;
           scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldName = (char *)malloc(sizeof(char)*20);
           sprintf(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldName, ">scaffold_%ld", scaffoldIndex);
           continue;
           
       }
       
       long int extendLength = strlen(scaffold);
       if(scaffold[extendLength-1] == '\n'){
           extendLength--;
       }
       long int scaffoldLength = 0;
       char * tempScaffold = NULL;
       if(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold != NULL){
           scaffoldLength = scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldLength;
           tempScaffold = (char *)malloc(sizeof(char)*(scaffoldLength+1));
           strncpy(tempScaffold, scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold, scaffoldLength);
           free(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold);
               
           scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold = (char *)malloc(sizeof(char)*(scaffoldLength + extendLength + 1));
           strncpy(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold, tempScaffold, scaffoldLength);
                   
           strncpy(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold + scaffoldLength, scaffold, extendLength);
           scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold[scaffoldLength + extendLength] = '\0';
           
           scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldLength = scaffoldLength + extendLength;
           
           free(tempScaffold);
       }else{
           scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold = (char *)malloc(sizeof(char)*(extendLength+1));
           strncpy(scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold, scaffold, extendLength);
           scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffold[extendLength] = '\0';
           scaffoldSetHead->scaffoldSet[scaffoldIndex].scaffoldLength = extendLength;
       }    
    }  
    
    fclose(fp);
    /*
    for(i = 0; i < scaffoldCount; i++){ 
        cout<<scaffoldSetHead->scaffoldSet[i].scaffoldName<<"--"<< scaffoldSetHead->scaffoldSet[i].scaffoldLength<<endl;
        cout<<scaffoldSetHead->scaffoldSet[i].scaffold<<endl;
    }
    */
    
    return scaffoldSetHead;
    
}

void SplitScaffoldToContig(ScaffoldSetHead * scaffoldSetHead, long int minGapDistance, char * gapInformationFile){
    
    ScaffoldSet * scaffoldSet = scaffoldSetHead->scaffoldSet;
    long int scaffoldCount = scaffoldSetHead->scaffoldCount;
    
    for(long int i = 0; i < scaffoldCount; i++){
        cout<<"tt--"<<i<<endl;
        long int gapDistance = 0;
        bool token = true;
        cout<<scaffoldSet[i].scaffold<<endl;
        for(long int j = 0; j < scaffoldSet[i].scaffoldLength; j++){
            
            gapDistance = 0;
            if(scaffoldSet[i].scaffold[j] != 'N' && scaffoldSet[i].scaffold[j] != 'n' && token !=false){
                scaffoldSet[i].contigCount++;
                token = false;
            }
            while(j < scaffoldSet[i].scaffoldLength && (scaffoldSet[i].scaffold[j] == 'N' || scaffoldSet[i].scaffold[j] == 'n')){
                gapDistance++;
                j++;
            }
            
            if(gapDistance > minGapDistance){
                cout<<scaffoldSet[i].gapCount<<"--"<<j<<"--"<<gapDistance<<endl;
                token = true;
                scaffoldSet[i].gapCount++;
                j--;
            }
            
        }
        
        //scaffoldSet[i].gapCount = scaffoldSet[i].contigCount;
        //scaffoldSet[i].contigCount++;
        
        cout<<"scaffoldIndex:"<<i<<",gapCount:"<<scaffoldSet[i].gapCount<<endl;
        
        if(scaffoldSet[i].contigCount<=1){
            scaffoldSetHead->contigCount++;
            scaffoldSet[i].contigStartCoordinate = (long int *)malloc(sizeof(long int)*scaffoldSet[i].contigCount);
            scaffoldSet[i].contigEndCoordinate = (long int *)malloc(sizeof(long int)*scaffoldSet[i].contigCount);
            scaffoldSet[i].contigStartCoordinate[0] = -1;
            scaffoldSet[i].contigEndCoordinate[0] = - 1;
            for(long int j = 0; j < scaffoldSet[i].scaffoldLength; j++){
                if(scaffoldSet[i].contigStartCoordinate[0] == -1 && (scaffoldSet[i].scaffold[j] != 'N' && scaffoldSet[i].scaffold[j] != 'n')){
                    scaffoldSet[i].contigStartCoordinate[0] = j;
                }
                if(scaffoldSet[i].contigStartCoordinate[0] != -1 && (scaffoldSet[i].scaffold[j] == 'N' || scaffoldSet[i].scaffold[j] == 'n')){
                    scaffoldSet[i].contigEndCoordinate[0] = j-1;
                }
            }
            if(scaffoldSet[i].contigEndCoordinate[0] == -1){
                scaffoldSet[i].contigEndCoordinate[0] = scaffoldSet[i].scaffoldLength-1;
            }
            
            
            continue;
        }
        
        scaffoldSet[i].contigStartCoordinate = (long int *)malloc(sizeof(long int)*scaffoldSet[i].contigCount);
        scaffoldSet[i].contigEndCoordinate = (long int *)malloc(sizeof(long int)*scaffoldSet[i].contigCount);
        scaffoldSet[i].gapDistance = (long int *)malloc(sizeof(long int)*(scaffoldSet[i].gapCount));
        
        long int contigIndex = 0;
        long int gapIndex = 0;
        scaffoldSet[i].contigStartCoordinate[0] = 0;
        scaffoldSet[i].contigEndCoordinate[scaffoldSet[i].contigCount-1] = scaffoldSet[i].scaffoldLength-1;
        
        for(long int j = 0; j < scaffoldSet[i].scaffoldLength; j++){
            gapDistance = 0;
            while(j < scaffoldSet[i].scaffoldLength && (scaffoldSet[i].scaffold[j] == 'N' || scaffoldSet[i].scaffold[j] == 'n')){
                gapDistance++;
                j++;
            }
            
            if(gapDistance > minGapDistance){
                if(j != gapDistance){
                    scaffoldSet[i].contigEndCoordinate[contigIndex] = j-gapDistance-1;
                    scaffoldSet[i].contigStartCoordinate[contigIndex+1] = j;
                    contigIndex++;
                }else{
                    scaffoldSet[i].contigStartCoordinate[0] = j;
                }
                //cout<<"index:"<<index<<"--"<<scaffoldSet[i].contigEndCoordinate[index]<<"--"<<scaffoldSet[i].contigStartCoordinate[index]<<endl;
                scaffoldSet[i].gapDistance[gapIndex] = gapDistance;
                gapIndex++;
                j--;
            }
        }
        
        if(scaffoldSet[i].scaffold[scaffoldSet[i].scaffoldLength-1] == 'N' || scaffoldSet[i].scaffold[scaffoldSet[i].scaffoldLength-1] == 'n'){
            scaffoldSet[i].contigEndCoordinate[scaffoldSet[i].contigCount-1] = scaffoldSet[i].scaffoldLength - gapDistance - 1;
        }

        scaffoldSetHead->contigCount = scaffoldSetHead->contigCount + scaffoldSet[i].contigCount;
        scaffoldSetHead->gapCount = scaffoldSetHead->gapCount + scaffoldSet[i].gapCount;
  
    }

    scaffoldSetHead->contigToGapIndex = (ContigToGapIndex *)malloc(sizeof(ContigToGapIndex)*scaffoldSetHead->contigCount);
    scaffoldSetHead->gapToContigIndex = (GapToContigIndex *)malloc(sizeof(GapToContigIndex)*scaffoldSetHead->gapCount);
    for(long int i = 0; i < scaffoldSetHead->contigCount; i++){
        scaffoldSetHead->contigToGapIndex[i].nextGapContigIndex = -1;
        scaffoldSetHead->contigToGapIndex[i].leftGapIndex = -1;
        scaffoldSetHead->contigToGapIndex[i].rightGapIndex = -1;
    }
    for(long int i = 0; i < scaffoldSetHead->gapCount; i++){
        scaffoldSetHead->gapToContigIndex[i].leftContigIndex = -1;
        scaffoldSetHead->gapToContigIndex[i].rightContigIndex = -1;
        scaffoldSetHead->gapToContigIndex[i].gapDistance = 0;
    }
    
    long int nextGapIndex = 0;
    long int tempGapCount = 0;
    long int tempContigCount = 0;
    
    for(long int i = 0; i < scaffoldCount; i++){    
        long int p = 0;
        if(scaffoldSet[i].contigStartCoordinate[0] > 0){
            p = 1;
        }
        for(long int j = 0; j < scaffoldSet[i].gapCount; j++){
            
            if(p == 1 && j == 0){
                scaffoldSetHead->gapToContigIndex[tempGapCount].leftContigIndex = -1;
            }else{
                scaffoldSetHead->gapToContigIndex[tempGapCount].leftContigIndex = tempContigCount+j-p;
            }
            scaffoldSetHead->gapToContigIndex[tempGapCount].rightContigIndex = tempContigCount +j-p+ 1;
            scaffoldSetHead->gapToContigIndex[tempGapCount].gapDistance = scaffoldSet[i].gapDistance[j];
            tempGapCount++;
        }
        tempContigCount = tempContigCount + scaffoldSet[i].contigCount;
        if(scaffoldSet[i].scaffold[scaffoldSet[i].scaffoldLength -1] == 'N' || scaffoldSet[i].scaffold[scaffoldSet[i].scaffoldLength -1] == 'n'){
            scaffoldSetHead->gapToContigIndex[tempGapCount - 1].rightContigIndex = -1;
        }
    }
    
    
    tempGapCount = 0;
    tempContigCount = 0;
    for(long int i = 0; i < scaffoldCount; i++){
        scaffoldSet[i].contigIndex = (long int *)malloc(sizeof(long int)*scaffoldSet[i].contigCount);    
        long int p = 0;
        long int p1 = 0;
        if(scaffoldSet[i].contigStartCoordinate[0] > 0){
            p = 1;
        }
        if(scaffoldSet[i].scaffold[scaffoldSet[i].scaffoldLength -1] == 'N' || scaffoldSet[i].scaffold[scaffoldSet[i].scaffoldLength -1] == 'n'){
            p1 = 1;
        }
        
        long int j = 0;
        for(j = 0; j < scaffoldSet[i].contigCount; j++){               
            scaffoldSetHead->contigToGapIndex[tempContigCount].nextGapContigIndex = tempContigCount + 1;
            if(p==0&&j==0){
                scaffoldSetHead->contigToGapIndex[tempContigCount].leftGapIndex = -1;
            }else{
                scaffoldSetHead->contigToGapIndex[tempContigCount].leftGapIndex = tempGapCount + j-1+p;
            }    
            scaffoldSetHead->contigToGapIndex[tempContigCount].rightGapIndex = tempGapCount+ j+p;
            scaffoldSet[i].contigIndex[j] = tempContigCount;
            tempContigCount++;
        }
        scaffoldSetHead->contigToGapIndex[tempContigCount-1].nextGapContigIndex = -1;
        if(p1!=1){
            scaffoldSetHead->contigToGapIndex[tempContigCount-1].rightGapIndex = -1;
        }
        tempGapCount = tempGapCount + scaffoldSet[i].contigCount + p + p1 - 1;
    }
    
    if(gapInformationFile != NULL){    
        FILE * fp;
        if((fp = fopen(gapInformationFile, "w+")) == NULL){
            printf("%s, does not exist!", gapInformationFile);
            exit(0);
        }
        
        fprintf(fp, "%ld\n", tempGapCount);
        for(long int i = 0; i < tempGapCount; i++){
            fprintf(fp, "%ld\t%ld\t%ld\n", scaffoldSetHead->gapToContigIndex[i].gapDistance, scaffoldSetHead->gapToContigIndex[i].leftContigIndex, scaffoldSetHead->gapToContigIndex[i].rightContigIndex);
        }  
    }
    
    
    
    
    for(long int i = 0; i < scaffoldSetHead->contigCount; i++){
        cout<<scaffoldSetHead->contigToGapIndex[i].nextGapContigIndex<<"--"<<scaffoldSetHead->contigToGapIndex[i].leftGapIndex<<"--"<<scaffoldSetHead->contigToGapIndex[i].rightGapIndex<<endl;
    }
    
    /*
    for(long int i = 0; i < tempGapCount; i++){
        cout<<scaffoldSetHead->gapToContigIndex[i].leftContigIndex<<"--"<<scaffoldSetHead->gapToContigIndex[i].rightContigIndex<<endl;
    }
    */
    /*
    for(long int i = 0; i < tempGapCount; i++){
        cout<<scaffoldSetHead->gapToContigIndex[i].leftContigIndex<<"--"<<scaffoldSetHead->gapToContigIndex[i].rightContigIndex<<"--"<<scaffoldSetHead->gapToContigIndex[i].gapDistance<<endl;
    }
    */
}

void WriteContigSetFromScffoldSet(ScaffoldSetHead * scaffoldSetHead, char * contigSetFileName){
    
    FILE * fpEndContigIndex;
    char endContigIndexFileName[35];
    strcpy(endContigIndexFileName, "scaffold_end_contig_index.fa");
    fpEndContigIndex = fopen(endContigIndexFileName, "w+");
    
    
    FILE * fp;
    if((fp = fopen(contigSetFileName, "w+")) == NULL){
        printf("%s, does not exist!", contigSetFileName);
        exit(0);
    }
    
    ScaffoldSet * scaffoldSet = scaffoldSetHead->scaffoldSet;
    long int scaffoldCount = scaffoldSetHead->scaffoldCount;
    
    long int contigIndex = 0;
    
    for(long int i = 0; i < scaffoldCount; i++){
        long int coordinateIndex = 0;
        for(long int j = 0; j < scaffoldSet[i].contigCount; j++){
            fprintf(fp, ">contig_%ld\n", contigIndex);
            //cout<<"contig:"<<contigIndex<<endl;
            contigIndex++;
            if(scaffoldSet[i].contigCount == 1){
                fprintf(fp, "%s\n", scaffoldSet[i].scaffold); 
                break;
            }

            long int contigLength = scaffoldSet[i].contigEndCoordinate[j] - scaffoldSet[i].contigStartCoordinate[j] + 1;
            if(contigLength < 1){
                printf("N in scaffold not right");
                exit(1);
            }
            //cout<<scaffoldSet[i].contigEndCoordinate[j]<<"--"<<scaffoldSet[i].contigStartCoordinate[j]<<"--"<<contigLength<<endl;
            char * tempContig = (char *)malloc(sizeof(char)*(contigLength+1));
            strncpy(tempContig, scaffoldSet[i].scaffold + scaffoldSet[i].contigStartCoordinate[j], contigLength);
            tempContig[contigLength] = '\0';
            fprintf(fp, "%s\n", tempContig);    
            free(tempContig);        
        }
        long int startN = 0;
        long int endN = 0;
        if(scaffoldSet[i].contigStartCoordinate[0] > 0){
            startN = 1;
        }
        if(scaffoldSet[i].scaffold[scaffoldSet[i].scaffoldLength -1] == 'N' || scaffoldSet[i].scaffold[scaffoldSet[i].scaffoldLength -1] == 'n'){
            endN = 1;
        }
        fprintf(fpEndContigIndex, "%ld\t%ld\t%ld\n", contigIndex - 1, startN, endN);
    }

}

void WriteEndContigFromScffoldSet(ScaffoldSetHead * scaffoldSetHead, long int subContigLength, long int intervalLength, char * endContigFileName){
    
    
    FILE * fp;
    if((fp = fopen(endContigFileName, "w+")) == NULL){
        printf("%s, does not exist!", endContigFileName);
        exit(0);
    }
    
    char * tempN = (char *)malloc(sizeof(char)*(intervalLength+1));
    for(long int i = 0; i < intervalLength; i++){
        tempN[i] = 'N';
    }
    tempN[intervalLength] = '\0';
    char * tempEndContig = (char *)malloc(sizeof(char)*(subContigLength/2+1));
    char * tempEndContig1 = (char *)malloc(sizeof(char)*(subContigLength/2+1));
    
    ScaffoldSet * scaffoldSet = scaffoldSetHead->scaffoldSet;
    long int scaffoldCount = scaffoldSetHead->scaffoldCount;
    
    long int contigIndex = 0;
    
    for(long int i = 0; i < scaffoldCount; i++){
        long int coordinateIndex = 0;
        for(long int j = 0; j < scaffoldSet[i].contigCount; j++){
            fprintf(fp, ">end_contig_%ld\n", contigIndex);
            contigIndex++;
            if(scaffoldSet[i].contigCount == 1){
                if(scaffoldSet[i].scaffoldLength > subContigLength){
                    strncpy(tempEndContig, scaffoldSet[i].scaffold, subContigLength/2);
                    strncpy(tempEndContig1, scaffoldSet[i].scaffold + scaffoldSet[i].scaffoldLength - subContigLength/2, subContigLength/2);
                    tempEndContig[subContigLength/2] = '\0';
                    tempEndContig1[subContigLength/2] = '\0';
                    fprintf(fp, "%s%s%s\n", tempEndContig,tempN,tempEndContig1);
                }else{
                    fprintf(fp, "%s\n", scaffoldSet[i].scaffold);
                }
                break;
            }

            long int contigLength = scaffoldSet[i].contigEndCoordinate[j] - scaffoldSet[i].contigStartCoordinate[j] + 1;
            if(contigLength < 1){
                printf("N in scaffold not right");
                exit(1);
            }
            char * tempContig = (char *)malloc(sizeof(char)*(contigLength+1));
            strncpy(tempContig, scaffoldSet[i].scaffold + scaffoldSet[i].contigStartCoordinate[j], contigLength);
            tempContig[contigLength] = '\0';
            
            if(contigLength > subContigLength){
                strncpy(tempEndContig, tempContig, subContigLength/2);
                strncpy(tempEndContig1, tempContig + contigLength - subContigLength/2, subContigLength/2);
                tempEndContig[subContigLength/2] = '\0';
                tempEndContig1[subContigLength/2] = '\0';
                fprintf(fp, "%s%s%s\n", tempEndContig,tempN,tempEndContig1);
            }else{
                fprintf(fp, "%s\n", tempContig);
            }
        
            free(tempContig);        
        }
    }
    free(tempEndContig);
    free(tempEndContig1);
    
}

ContigSetHead * GetContigSet(char * contigSetFile){
    
    ContigSetHead * contigSetHead = (ContigSetHead *)malloc(sizeof(ContigSetHead));
    contigSetHead->contigSet = NULL;
    contigSetHead->contigCount = 0;
    
    long int maxSize = 10000;
    char * contig = NULL;
    if(NULL == (contig = (char*)malloc(sizeof(char)*maxSize))){
        perror("malloc error!");
        exit(1);
    }
    
    FILE * fp; 
    if((fp = fopen(contigSetFile, "r")) == NULL){
        printf("%s, does not exist!", contigSetFile);
        exit(0);
    }
    while((fgets(contig, maxSize, fp)) != NULL){ 
       if(contig[0] == '>'){  
           contigSetHead->contigCount++; 
       }  
    }  
    fclose(fp);
    
    contigSetHead->contigSet = (ContigSet *)malloc(sizeof(ContigSet)*contigSetHead->contigCount);
    for(long int i = 0; i < contigSetHead->contigCount; i++){
        contigSetHead->contigSet[i].contig = NULL;
        contigSetHead->contigSet[i].contigLength = 0;
    }
    
    if((fp = fopen(contigSetFile, "r")) == NULL){
        printf("%s, does not exist!", contigSetFile);
        exit(0);
    }
    
    long int contigIndex = -1;
    while((fgets(contig, maxSize, fp)) != NULL){ 
       
       if(contig[0] == '>'){  
           
           if(strlen(contig) == maxSize-1){              
               while((fgets(contig, maxSize, fp)) != NULL){
                   if(strlen(contig) != maxSize-1){
                       break;
                   }
               }        
           }
           contigIndex++;
           continue;
           
       }
       
       long int extendLength = strlen(contig);
       if(contig[extendLength-1] == '\n'){
           extendLength--;
       }
       long int contigLength = 0;
       char * tempContig = NULL;
       if(contigSetHead->contigSet[contigIndex].contig != NULL){
           contigLength = contigSetHead->contigSet[contigIndex].contigLength;
           tempContig = (char *)malloc(sizeof(char)*(contigLength+1));
           strncpy(tempContig, contigSetHead->contigSet[contigIndex].contig, contigLength);
           free(contigSetHead->contigSet[contigIndex].contig);
               
           contigSetHead->contigSet[contigIndex].contig = (char *)malloc(sizeof(char)*(contigLength + extendLength + 1));
           strncpy(contigSetHead->contigSet[contigIndex].contig, tempContig, contigLength);
                   
           strncpy(contigSetHead->contigSet[contigIndex].contig + contigLength, contig, extendLength);
           contigSetHead->contigSet[contigIndex].contig[contigLength + extendLength] = '\0';
           
           contigSetHead->contigSet[contigIndex].contigLength = contigLength + extendLength;
           
           free(tempContig);
       }else{
           contigSetHead->contigSet[contigIndex].contig = (char *)malloc(sizeof(char)*(extendLength+1));
           strncpy(contigSetHead->contigSet[contigIndex].contig, contig, extendLength);
           contigSetHead->contigSet[contigIndex].contig[extendLength] = '\0';
           contigSetHead->contigSet[contigIndex].contigLength = extendLength;
       }    
    }  
    
    fclose(fp);
    
    return contigSetHead;
}



