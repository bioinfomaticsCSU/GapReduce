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

GapContigIndexHead * GetGapContigIndexHead(char * gapContigIndexFile){
    
    FILE * fp;
    
    if((fp = fopen(gapContigIndexFile, "r")) == NULL){
        printf("%s, does not exist!", gapContigIndexFile);
        exit(0);
    }
    
    long int maxSize = 1000;
    char * line = NULL;
    if(NULL == (line = (char*)malloc(sizeof(char)*maxSize))){
        perror("malloc error!");
        exit(1);
    }
    
    long int gapCount = 0;
    while((fgets(line, maxSize, fp)) != NULL){ 
        if(line[strlen(line) - 1] != '\n'){
            continue;
        }
        gapCount++; 
    }  
    fclose(fp);
    
    GapContigIndexHead * gapContigIndexHead = (GapContigIndexHead *)malloc(sizeof(GapContigIndexHead));
    gapContigIndexHead->gapCount = gapCount;
    gapContigIndexHead->gapContigIndex = (GapContigIndex *)malloc(sizeof(GapContigIndex)*gapCount);
    for(long int i = 0; i < gapCount; i++){
        gapContigIndexHead->gapContigIndex[i].gapDistance = 0;
        gapContigIndexHead->gapContigIndex[i].leftContigIndex = -1;
        gapContigIndexHead->gapContigIndex[i].rightContigIndex = -1;
    }
    
    if((fp = fopen(gapContigIndexFile, "r")) == NULL){
        printf("%s, does not exist!", gapContigIndexFile);
        exit(0);
    }
    
    long int gapIndex = 0;
    while((fgets(line, maxSize, fp)) != NULL){ 
        if(line[strlen(line) - 1] != '\n'){
            printf("line too long!");
            exit(0);
        }
        char * p = strtok(line, "\t");
        sscanf(p, "%ld", &gapContigIndexHead->gapContigIndex[gapIndex].gapDistance);
        p = strtok(NULL, "\t");
        sscanf(p, "%ld", &gapContigIndexHead->gapContigIndex[gapIndex].leftContigIndex);
        p = strtok(NULL, "\t");
        sscanf(p, "%ld", &gapContigIndexHead->gapContigIndex[gapIndex].rightContigIndex);
        gapIndex++;
    }  
    fclose(fp);
    
    return gapContigIndexHead;
    
}


ContigSet * GetContigSet(char * contigFileName, long int & contigNumber, SubContigCoordinateSetHead * subContigCoordinateSetHead){

    long int maxSize = 10000;
    char * contig = NULL;
    if(NULL == (contig = (char*)malloc(sizeof(char)*maxSize))){
        perror("malloc error!");
        exit(1);
    }
    
    contigNumber = 0;
    
    FILE * fp;
    
    if((fp = fopen(contigFileName, "r")) == NULL){
        printf("%s, does not exist!", contigFileName);
        exit(0);
    }
    
    while((fgets(contig, maxSize, fp)) != NULL){ 
       if(contig[0] == '>'){  
           contigNumber++; 
       }  
    }  
    
    fclose(fp);
    
    ContigSet * contigSet = NULL;
    if(NULL == (contigSet = (ContigSet*)malloc(sizeof(ContigSet)*contigNumber))){
        perror("malloc error!");
        exit(1);
    }
    long int i = 0;
    for(i=0;i<contigNumber;i++){
        contigSet[i].contig = NULL;
        contigSet[i].contigName = NULL;
        contigSet[i].contigLength = 0;
        contigSet[i].subContigCoordinateIndexSet = NULL;
        contigSet[i].subContigNumber = 0;
    }
    
    if((fp = fopen(contigFileName, "r")) == NULL){
        printf("%s, does not exist!", contigFileName);
        exit(0);
    }
    i = 0;
    while((fgets(contig, maxSize, fp)) != NULL){ 
       cout<<"b--"<<i<<endl;
       if(contig[0] == '>'){  
           if(NULL == (contigSet[i].contigName = (char*)malloc(sizeof(char)*(strlen(contig)-1)))){
               perror("malloc error!");
               exit(1);
           }
           strncpy(contigSet[i].contigName, contig+1, strlen(contig)-2); 
           contigSet[i].contigName[strlen(contig)-2] = '\0';
           i++;
       }else{
           long int length = 0;
           if(contigSet[i-1].contig != NULL){
               length = strlen(contigSet[i-1].contig);
           }
           long int length1 = strlen(contig);
           if(contig[length1 - 1] == '\n'){
               length1--;
           }
           char * tempContig = NULL;
           if(NULL == (tempContig = (char*)malloc(sizeof(char)*(length + length1 + 1)))){
               perror("malloc error!");
           }
           if(length != 0){
               strncpy(tempContig, contigSet[i-1].contig, length);
           }
           strncpy(tempContig + length, contig, length1);
           free(contigSet[i-1].contig);
           contigSet[i-1].contig = tempContig;
           contigSet[i-1].contig[length + length1] = '\0';
           contigSet[i-1].contigLength = length + length1;
       }  
    }
    
    fclose(fp);
    /*
    for(i=0;i<contigNumber;i++){
        printf(">%s\n", contigSet[i].contigName);
        printf("%s\n", contigSet[i].contig);
    }
    */
    free(contig);
    cout<<"ee"<<endl;
    for(i = 0; i < subContigCoordinateSetHead->subContigNumber; i++){
        long int index = SearchContigIndexFromContigName(contigSet, contigNumber, subContigCoordinateSetHead->subContigCoordinateSet[i].contigName);
        if(index < 0){
            continue;
        }
        subContigCoordinateSetHead->subContigCoordinateSet[i].contigIndex = index;
        contigSet[index].subContigNumber++;
    }
    cout<<"ff"<<endl;
    for(i = 0; i < contigNumber; i++){
        if(contigSet[i].subContigNumber <= 0){
            continue;
        }
        
        contigSet[i].subContigCoordinateIndexSet = (long int *)malloc(sizeof(long int)*contigSet[i].subContigNumber);
        for(long int j = 0; j<contigSet[i].subContigNumber; j++){
            contigSet[i].subContigCoordinateIndexSet[j] = -1;
        }
    }
    cout<<"gg"<<endl;
    for(i = 0; i < subContigCoordinateSetHead->subContigNumber; i++){
        long int index = SearchContigIndexFromContigName(contigSet, contigNumber, subContigCoordinateSetHead->subContigCoordinateSet[i].contigName);
        if(index < 0){
            continue;
        }
        for(long int j = 0; j < contigSet[index].subContigNumber; j++){
            if(contigSet[index].subContigCoordinateIndexSet[j] == -1){
                contigSet[index].subContigCoordinateIndexSet[j] = i;
                break;
            }
        }    
    }
    cout<<"hh"<<endl;
    return contigSet;
}

void OutputSubContigSet(ContigSet * contigSet, long int contigNumber, SubContigCoordinateSetHead * subContigCoordinateSetHead){
    
    char subContigSetFileName[] = "subContig.fa";
    FILE * fp;
    if((fp = fopen(subContigSetFileName, "w")) == NULL){
        printf("%s, does not exist!", subContigSetFileName);
        exit(0);
    }
    
    char * read = NULL;
    long int readLength = 60;
    if(NULL == (read = (char*)malloc(sizeof(char)*(readLength + 1)))){
        perror("malloc error!");
        exit(1);
    }
    
    for(long int i = 0; i < subContigCoordinateSetHead->subContigNumber; i++){
        long int index = subContigCoordinateSetHead->subContigCoordinateSet[i].contigIndex;
        
        fprintf(fp, ">subContig.%ld\n", i);
        
        
        long int tempContigLength = labs(subContigCoordinateSetHead->subContigCoordinateSet[i].contigEnd - subContigCoordinateSetHead->subContigCoordinateSet[i].contigStart) + 2;
        char * tempContig = NULL;
        if(NULL == (tempContig = (char*)malloc(sizeof(char)*tempContigLength))){
            perror("malloc error!");
            exit(1);
        }
            
        if(subContigCoordinateSetHead->subContigCoordinateSet[i].isReverse == false){
            strncpy(tempContig, contigSet[index].contig + subContigCoordinateSetHead->subContigCoordinateSet[i].contigStart - 1, tempContigLength - 1);
        }else{
            strncpy(tempContig, contigSet[index].contig + subContigCoordinateSetHead->subContigCoordinateSet[i].contigEnd - 1, tempContigLength - 1);
        }
        tempContig[tempContigLength -1] = '\0';
        
        long int tempLength = readLength;
        while(tempLength <= tempContigLength - 1){
            strncpy(read, tempContig + tempLength - readLength, readLength);
            read[readLength] = '\0';
            fprintf(fp, "%s\n", read);
            //cout<<read<<endl;
            tempLength = tempLength + readLength;
        }
        //fprintf(fp, "aa--%s\n", read);
        if(tempLength - readLength < tempContigLength - 1){
            //fprintf(fp, "%ld--%ld\n", tempLength - readLength, tempContigLength-1);
            strncpy(read, tempContig + tempLength - readLength, readLength - (tempLength - tempContigLength + 1));
            read[readLength - (tempLength - tempContigLength + 1)] = '\0';
            //cout<<read<<endl;
            fprintf(fp, "%s\n", read);
        }
        
        free(tempContig);      
    }
    
}


long int SearchContigIndexFromContigName(ContigSet * contigSet, long int contigNumber, char * contigName){
    //setvbuf(stdout,NULL,_IONBF,0);  

    for(long int i = 0; i < contigNumber; i++){
        if(strlen(contigSet[i].contigName) < strlen(contigName)){
            continue;
        }
        if(strlen(contigSet[i].contigName) == strlen(contigName)){
            if(strcmp(contigSet[i].contigName, contigName) == 0){
                return i;
            }
        }
        if(strlen(contigSet[i].contigName) > strlen(contigName)){
            long int length = strlen(contigName);
            long int j = 0;
            for(j = 0; j < length; j++){
                if(contigSet[i].contigName[j] != contigName[j]){
                    break;
                }
            }
            if(j != length){
                continue;
            }else{
                return i;
            }
        }
    }
    
    return -1;
    
}

long int SearchReferenceIndexFromContigName(ReferenceCoordinateSetHead * referenceCoordinateSetHead, char * referenceName){
    //setvbuf(stdout,NULL,_IONBF,0);  

    for(long int i = 0; i < referenceCoordinateSetHead->referenceNumber; i++){
        if(strlen(referenceCoordinateSetHead->referenceCoordinateSet[i].referenceName) < strlen(referenceName)){
            continue;
        }
        if(strlen(referenceCoordinateSetHead->referenceCoordinateSet[i].referenceName) == strlen(referenceName)){
            if(strcmp(referenceCoordinateSetHead->referenceCoordinateSet[i].referenceName, referenceName) == 0){
                return i;
            }
        }
        if(strlen(referenceCoordinateSetHead->referenceCoordinateSet[i].referenceName) > strlen(referenceName)){
            long int length = strlen(referenceName);
            long int j = 0;
            for(j = 0; j < length; j++){
                if(referenceCoordinateSetHead->referenceCoordinateSet[i].referenceName[j] != referenceName[j]){
                    break;
                }
            }
            if(j != length){
                continue;
            }else{
                return i;
            }
        }
    }
    
    return -1;
    
}

SubContigCoordinateSetHead * GetSubContigCoordinateSetHeadFromCoordinateFile(char * coordinateFileName){
    
    long int i = 0;
    
    FILE * fp;
    
    if((fp = fopen(coordinateFileName, "r")) == NULL){
        printf("%s, does not exist!", coordinateFileName);
        exit(0);
    }
    
    long int subContigNumber = 0;
    char * line = (char *)malloc(sizeof(char)*10000);
    
    while(!feof(fp)){
        fgets(line, 10000, fp);
        if(feof(fp))break;
        if(i<4){
            i++;
            continue;
        }
        subContigNumber++;
    }
    
    SubContigCoordinateSetHead * subContigCoordinateSetHead = (SubContigCoordinateSetHead *)malloc(sizeof(SubContigCoordinateSetHead));
    subContigCoordinateSetHead->subContigNumber = subContigNumber;
    subContigCoordinateSetHead->subContigCoordinateSet = (SubContigCoordinate*)malloc(sizeof(SubContigCoordinate)*subContigNumber);
    
    for(i=0;i<subContigNumber;i++){
        subContigCoordinateSetHead->subContigCoordinateSet[i].referenceName = NULL;
        subContigCoordinateSetHead->subContigCoordinateSet[i].referenceIndex = -1;
        subContigCoordinateSetHead->subContigCoordinateSet[i].contigName = NULL;
        subContigCoordinateSetHead->subContigCoordinateSet[i].contigIndex = -1;
        subContigCoordinateSetHead->subContigCoordinateSet[i].referenceStart = 0;
        subContigCoordinateSetHead->subContigCoordinateSet[i].referenceEnd = 0;
        subContigCoordinateSetHead->subContigCoordinateSet[i].contigStart = 0;
        subContigCoordinateSetHead->subContigCoordinateSet[i].contigEnd = 0;
        subContigCoordinateSetHead->subContigCoordinateSet[i].referenceOverlapLength = 0;
        subContigCoordinateSetHead->subContigCoordinateSet[i].contigOverlapLength = 0;
        subContigCoordinateSetHead->subContigCoordinateSet[i].identityPercent = 0;
        subContigCoordinateSetHead->subContigCoordinateSet[i].isReverse = false;
    }
    
    fclose(fp);
    
    if((fp = fopen(coordinateFileName, "r")) == NULL){
        printf("%s, does not exist!", coordinateFileName);
        exit(0);
    }
    
    char * p = NULL;
    char * referenceName = NULL;
    
    i = 0;
    long int index = 0;
    
    while(!feof(fp)){
        fgets(line, 10000, fp);
        if(feof(fp))break;
        if(i<4){
            i++;
            continue;
        }
        
        for(long int j = 0; j<13; j++){
            if(j==0){
                p = strtok(line, "\t");
            }else{
                p = strtok(NULL, "\t");
            }
            if(j==0){
                sscanf(p, "%ld", &subContigCoordinateSetHead->subContigCoordinateSet[index].referenceStart);
            }
            if(j==1){
                sscanf(p, "%ld", &subContigCoordinateSetHead->subContigCoordinateSet[index].referenceEnd);
            }
            if(j==2){
                sscanf(p, "%ld", &subContigCoordinateSetHead->subContigCoordinateSet[index].contigStart);
            }
            if(j==3){
                sscanf(p, "%ld", &subContigCoordinateSetHead->subContigCoordinateSet[index].contigEnd);
            }
            if(j==4){
                sscanf(p, "%ld", &subContigCoordinateSetHead->subContigCoordinateSet[index].referenceOverlapLength);
            }
            if(j==5){
                sscanf(p, "%ld", &subContigCoordinateSetHead->subContigCoordinateSet[index].contigOverlapLength);
            }
            if(j==6){
                sscanf(p, "%lf", &subContigCoordinateSetHead->subContigCoordinateSet[index].identityPercent);
            }
            if(j==10){
                int t = 0;
                sscanf(p, "%d", &t);
                if(t == -1){
                    subContigCoordinateSetHead->subContigCoordinateSet[index].isReverse = true;
                }else{
                    subContigCoordinateSetHead->subContigCoordinateSet[index].isReverse = false;
                }
            }
            if(j==11){
                subContigCoordinateSetHead->subContigCoordinateSet[index].referenceName = new char[strlen(p)+1];
                strcpy( subContigCoordinateSetHead->subContigCoordinateSet[index].referenceName, p);
                subContigCoordinateSetHead->subContigCoordinateSet[index].referenceName[strlen(p)] = '\0';
            }
            if(j==12){
                subContigCoordinateSetHead->subContigCoordinateSet[index].contigName = new char[strlen(p)+1];
                strcpy(subContigCoordinateSetHead->subContigCoordinateSet[index].contigName, p);
                subContigCoordinateSetHead->subContigCoordinateSet[index].contigName[strlen(p)] = '\0';
            }
            
        }
        index++;
    }
    
    return subContigCoordinateSetHead;
    
    
}

ReferenceCoordinateSetHead * ConstructReferenceCoordinate(char * referenceFileName, SubContigCoordinateSetHead * subContigCoordinateSetHead){
    
    long int i = 0;
    
    long int maxSize = 10000000;
    char * contig = NULL;
    if(NULL == (contig = (char*)malloc(sizeof(char)*(maxSize+1)))){
        perror("malloc error!");
        exit(1);
    }
    
    long int referenceNumber = 0;
    
    FILE * fp;
    
    if((fp = fopen(referenceFileName, "r")) == NULL){
        printf("%s, does not exist!", referenceFileName);
        exit(0);
    }
    
    while((fgets(contig, maxSize, fp)) != NULL){ 
       if(contig[0] == '>'){  
           referenceNumber++; 
       }  
    }  
    
    ReferenceCoordinateSetHead * referenceCoordinateSetHead = (ReferenceCoordinateSetHead*)malloc(sizeof(ReferenceCoordinateSetHead));
    referenceCoordinateSetHead->referenceNumber = referenceNumber;
    referenceCoordinateSetHead->referenceCoordinateSet = (ReferenceCoordinate *)malloc(sizeof(ReferenceCoordinate)*referenceNumber);
    
    for(i = 0; i < referenceNumber; i++){
        referenceCoordinateSetHead->referenceCoordinateSet[i].referenceName = NULL;
        referenceCoordinateSetHead->referenceCoordinateSet[i].reference = NULL;
        referenceCoordinateSetHead->referenceCoordinateSet[i].referenceLength = 0;
        referenceCoordinateSetHead->referenceCoordinateSet[i].subContigNumber = 0;
        referenceCoordinateSetHead->referenceCoordinateSet[i].subContigCoordinateIndexSet = NULL;
    }
    
    fclose(fp);
    
    if((fp = fopen(referenceFileName, "r")) == NULL){
        printf("%s, does not exist!", referenceFileName);
        exit(0);
    }
    
    i = 0;
    long int index = 0;
    
    while((fgets(contig, maxSize, fp)) != NULL){ 
       cout<<"bb--"<<i<<endl;
       if(contig[0] == '>'){  
           if(NULL == (referenceCoordinateSetHead->referenceCoordinateSet[i].referenceName = (char*)malloc(sizeof(char)*(strlen(contig)-1)))){
               perror("malloc error!");
               exit(1);
           }
           strncpy(referenceCoordinateSetHead->referenceCoordinateSet[i].referenceName, contig+1, strlen(contig)-2); 
           referenceCoordinateSetHead->referenceCoordinateSet[i].referenceName[strlen(contig)-2] = '\0';
           i++;
           index = i - 1;
       }else{
           long int length = 0;
           if(referenceCoordinateSetHead->referenceCoordinateSet[index].reference != NULL){
               length = strlen(referenceCoordinateSetHead->referenceCoordinateSet[index].reference);
           }
           long int length1 = strlen(contig);
           if(contig[length1 - 1] == '\n'){
               length1--;
           }
           char * tempContig = NULL;
           bool extend = false;
           while(length + length1 > maxSize){
               extend = true;
               maxSize = maxSize + maxSize;
           }
           if(referenceCoordinateSetHead->referenceCoordinateSet[index].reference == NULL){
               tempContig = (char *)malloc(sizeof(char)*(maxSize+1));
               strncpy(tempContig, contig, length1);
               tempContig[length1] = '\0';
               referenceCoordinateSetHead->referenceCoordinateSet[index].reference = tempContig;
               referenceCoordinateSetHead->referenceCoordinateSet[index].referenceLength = length1;
           }else{
               if(extend == false){
                   strncpy(referenceCoordinateSetHead->referenceCoordinateSet[index].reference+length, contig, length1);
               }else{
                   tempContig = (char *)malloc(sizeof(char)*maxSize);
                   strncpy(tempContig, referenceCoordinateSetHead->referenceCoordinateSet[index].reference, length);
                   strncpy(tempContig + length, contig, length1);
                   free(referenceCoordinateSetHead->referenceCoordinateSet[index].reference);
                   referenceCoordinateSetHead->referenceCoordinateSet[index].reference = tempContig;
               }
               referenceCoordinateSetHead->referenceCoordinateSet[index].reference[length+length1] = '\0';
               referenceCoordinateSetHead->referenceCoordinateSet[index].referenceLength = length + length1;
           }
       }  
    }
    cout<<"eee"<<endl;
    long int subContigNumber = subContigCoordinateSetHead->subContigNumber;
    for(i = 0; i < subContigNumber; i++){
        
        long int index = SearchReferenceIndexFromContigName(referenceCoordinateSetHead, subContigCoordinateSetHead->subContigCoordinateSet[i].referenceName);
        if(index < 0){
            continue;
        }
        subContigCoordinateSetHead->subContigCoordinateSet[i].referenceIndex = index;
        referenceCoordinateSetHead->referenceCoordinateSet[index].subContigNumber++;
        
    }
    cout<<"fff"<<endl;
    for(i = 0; i < referenceNumber; i++){
        
        referenceCoordinateSetHead->referenceCoordinateSet[i].subContigCoordinateIndexSet = (long int *)malloc(sizeof(long int)*referenceCoordinateSetHead->referenceCoordinateSet[i].subContigNumber);
        
        for(long int j = 0; j < referenceCoordinateSetHead->referenceCoordinateSet[i].subContigNumber; j++){
            referenceCoordinateSetHead->referenceCoordinateSet[i].subContigCoordinateIndexSet[j] = -1;
        }
        
    }
    cout<<"ggg"<<endl;
    for(i = 0; i < subContigNumber; i++){
        long int index = subContigCoordinateSetHead->subContigCoordinateSet[i].referenceIndex;
        if(index < 0){
            continue;
        }
        for(long int j = 0; j < referenceCoordinateSetHead->referenceCoordinateSet[index].subContigNumber; j++){
            if(referenceCoordinateSetHead->referenceCoordinateSet[index].subContigCoordinateIndexSet[j] == -1){
                referenceCoordinateSetHead->referenceCoordinateSet[index].subContigCoordinateIndexSet[j] = i;
                break;
            }
        }
        
    }
    
    
    cout<<"hhh"<<endl;
    //OutputReferenceContigCoordinate(referenceCoordinateHead);
    
    return referenceCoordinateSetHead;

}

bool ReverseComplement(char * temp1, char * temp2){
    long int len = strlen(temp1);
    long int i = 0;
    long int j = 0;
    for(i=0;i<len;i++){
        if(temp1[i]=='A'){
            temp2[len-1-i]='T';
        }else if(temp1[i]=='T'){
            temp2[len-1-i]='A';
        }else if(temp1[i]=='G'){
            temp2[len-1-i]='C';
        }else if(temp1[i]=='C'){
            temp2[len-1-i]='G';
        }else{
             return false;
        }
    }
    temp2[len]='\0';
    return true;
}

void GetMinGapDistanceIndex(SubContigCoordinate * subContigCoordinateSet, long int * leftIndexSet, long int leftCount, long int * rightIndexSet, long int rightCount, long int & leftIndex, long int & rightIndex, long int gapDistance){
    long int d = -1;
    for(long int i = 0; i < leftCount; i++){
        for(long int j = 0; j < rightCount; j++){
            long int tempDistance = labs(subContigCoordinateSet[rightIndexSet[j]].referenceStart - subContigCoordinateSet[leftIndexSet[i]].referenceEnd - 1);
            //cout<<"tempDistance:"<<i<<"--"<<j<<"--"<<tempDistance<<endl;
            if(d == -1){
                d = tempDistance;
                leftIndex = leftIndexSet[i];
                rightIndex = rightIndexSet[j];
                continue;
            }
            if(tempDistance < d && d!=-1){
                d = tempDistance;
                leftIndex = leftIndexSet[i];
                rightIndex = rightIndexSet[j];
            }
        }
    }
}

void OutPutGapInScaffold(GapContigIndexHead * gapContigIndexHead, ReferenceCoordinateSetHead * referenceCoordinateSetHead, SubContigCoordinateSetHead * subContigCoordinateSetHead, ContigSet * contigSet, long int contigNumber, char * gapFileName){
    
    long int referenceNumber = referenceCoordinateSetHead->referenceNumber;
    
    FILE * fp;
    
    if((fp = fopen(gapFileName, "w+")) == NULL){
        printf("%s, does not exist!", gapFileName);
        exit(0);
    }
    
    long int gapIndex = 0;
    
    for(long int i = 0; i < gapContigIndexHead->gapCount; i++){
        
        cout<<"gapIndex:"<<i<<endl;
        long int leftContigIndex = gapContigIndexHead->gapContigIndex[i].leftContigIndex;
        long int rightContigIndex = gapContigIndexHead->gapContigIndex[i].rightContigIndex;
        
        long int leftMapContigIndex = -1;
        long int rightMapContigIndex = -1;
        if(leftContigIndex != -1){
            if(contigSet[leftContigIndex].subContigCoordinateIndexSet!=NULL){
                leftMapContigIndex = contigSet[leftContigIndex].subContigCoordinateIndexSet[0];
            }
        }
        if(rightContigIndex != -1){
            if(contigSet[rightContigIndex].subContigCoordinateIndexSet!=NULL){
                rightMapContigIndex = contigSet[rightContigIndex].subContigCoordinateIndexSet[0];
            }
        }
        
        fprintf(fp, ">gap_index_%ld\n", i);
        if((leftMapContigIndex == -1 || rightMapContigIndex == -1) && leftContigIndex != -1 && rightContigIndex != -1){
            char * gapRegion = (char *)malloc(sizeof(char)*(gapContigIndexHead->gapContigIndex[i].gapDistance + 1));
            for(long int j = 0; j < gapContigIndexHead->gapContigIndex[i].gapDistance; j++){
                gapRegion[j] = 'N';
            }
            gapRegion[gapContigIndexHead->gapContigIndex[i].gapDistance] = '\0';
            fprintf(fp, "%s\n", gapRegion);
            free(gapRegion);
        }else{
            if(leftContigIndex != -1 && rightContigIndex != -1){
                GetMinGapDistanceIndex(subContigCoordinateSetHead->subContigCoordinateSet, contigSet[leftContigIndex].subContigCoordinateIndexSet, contigSet[leftContigIndex].subContigNumber, contigSet[rightContigIndex].subContigCoordinateIndexSet, contigSet[rightContigIndex].subContigNumber, leftMapContigIndex, rightMapContigIndex, gapContigIndexHead->gapContigIndex[i].gapDistance);
            }
            long int referenceEnd = 0;
            long int leftReferenceIndex = -1;
            if(leftContigIndex != -1){
                referenceEnd = subContigCoordinateSetHead->subContigCoordinateSet[leftMapContigIndex].referenceEnd;
                leftReferenceIndex = subContigCoordinateSetHead->subContigCoordinateSet[leftMapContigIndex].referenceIndex;
            }
            long int referenceStart = 0;
            long int rightReferenceIndex = -1;
            if(rightContigIndex != -1){
                referenceStart = subContigCoordinateSetHead->subContigCoordinateSet[rightMapContigIndex].referenceStart;
                rightReferenceIndex = subContigCoordinateSetHead->subContigCoordinateSet[rightMapContigIndex].referenceIndex;
            }else{
                referenceStart = strlen(referenceCoordinateSetHead->referenceCoordinateSet[leftReferenceIndex].reference);
            }
            //fprintf(fp, "\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n", referenceEnd, referenceStart, leftReferenceIndex, rightReferenceIndex, contigSet[leftContigIndex].subContigNumber, contigSet[rightContigIndex].subContigNumber);
            //cout<<referenceEnd<<"--"<<referenceStart<<"--"<<leftReferenceIndex<<"--"<<rightReferenceIndex<<endl;
            if(referenceEnd < referenceStart && leftReferenceIndex == rightReferenceIndex){
                char * reference = referenceCoordinateSetHead->referenceCoordinateSet[leftReferenceIndex].reference;
                char * gapRegion = (char *)malloc(sizeof(char)*(referenceStart - referenceEnd));
                strncpy(gapRegion, reference + referenceEnd, referenceStart - referenceEnd - 1);
                gapRegion[referenceStart - referenceEnd - 1] = '\0';
                fprintf(fp, "%s\n", gapRegion);
                free(gapRegion);
            }else{
                char * gapRegion = (char *)malloc(sizeof(char)*(gapContigIndexHead->gapContigIndex[i].gapDistance + 1));
                for(long int j = 0; j < gapContigIndexHead->gapContigIndex[i].gapDistance; j++){
                    gapRegion[j] = 'N';
                }
                gapRegion[gapContigIndexHead->gapContigIndex[i].gapDistance] = '\0';
                fprintf(fp, "%s\n", gapRegion);
                free(gapRegion); 
            }
        }
    }
    
    
    
}

void OutputContigSetCoordinate(ContigSet * contigSet, SubContigCoordinateSetHead * subContigCoordinateSetHead, long int contigNumber){
    
    long int i = 0;
    long int j = 0;
    
    long int multipleSubContig = 0;
    long int noMappedContig = 0;
    long int singleMappedContig = 0;
    
    for(i = 0; i < contigNumber; i++){
        
        if(contigSet[i].subContigNumber<=0){
            noMappedContig++;
        }
        
        if(contigSet[i].subContigNumber<=1){
            //continue;
        }else{
            multipleSubContig++;
        }
        
        if(contigSet[i].subContigNumber==1){
            singleMappedContig++;
        }
        /*
        printf("contigLength:%ld\n", contigSet[i].contigLength);
        
        for(j=0; j<contigSet[i].subContigNumber;j++){
            long int index = contigSet[i].subContigCoordinateIndexSet[j];
            printf("contigIndex:%ld\treferenceIndex:%ld\t%ld\t%ld\t%ld\t%ld\n", i, 
                subContigCoordinateSetHead->subContigCoordinateSet[index].referenceIndex, 
                subContigCoordinateSetHead->subContigCoordinateSet[index].contigStart, 
                subContigCoordinateSetHead->subContigCoordinateSet[index].contigEnd, 
                subContigCoordinateSetHead->subContigCoordinateSet[index].referenceStart, 
                subContigCoordinateSetHead->subContigCoordinateSet[index].referenceEnd);
            
            
        }
        
        printf("-------------------------\n");
        */

    }
    
    printf("contigNumber:%ld\tmultipleSubContig:%ld\tsingleMappedContig:%ld\tnoMappedContig:%ld\n", contigNumber, multipleSubContig, singleMappedContig, noMappedContig);
    /*
    for(i = 0; i < contigNumber; i++){
        
        if(contigSet[i].subContigNumber>0){
            continue;
        }
        
        printf("contigName:%s\tcontigLength:%ld\n", contigSet[i].contigName, contigSet[i].contigLength);
        printf("%s\n", contigSet[i].contig);

    }
    
    for(i = 0; i < contigNumber; i++){
        
        //printf("contigName:%s\tcontigLength:%ld\t%f\n", contigSet[i].contigName, contigSet[i].contigLength, subContigCoordinateSetHead->subContigCoordinateSet[contigSet[i].subContigCoordinateIndexSet[0]].identityPercent);
        printf("%ld--%ld\n", contigSet[i].contigLength, subContigCoordinateSetHead->subContigCoordinateSet[contigSet[i].subContigCoordinateIndexSet[0]].contigOverlapLength);
        
        
        if(contigSet[i].contigLength != subContigCoordinateSetHead->subContigCoordinateSet[contigSet[i].subContigCoordinateIndexSet[0]].contigOverlapLength){
            printf("unequal!\n");
        }
    }
    */

}


void OutputReferenceCoordinateSet(ReferenceCoordinateSetHead * referenceCoordinateSetHead, SubContigCoordinateSetHead * subContigCoordinateSetHead){
    
    long int referenceNumber = referenceCoordinateSetHead->referenceNumber;
    long int * subContigCoordinateIndexSet = NULL;
    
    for(long int i = 0; i < referenceNumber; i++){
        
        printf("index:%ld\treferenceName:%s\n", i, referenceCoordinateSetHead->referenceCoordinateSet[i].referenceName);
        long int subContigNumber = referenceCoordinateSetHead->referenceCoordinateSet[i].subContigNumber;
        subContigCoordinateIndexSet = referenceCoordinateSetHead->referenceCoordinateSet[i].subContigCoordinateIndexSet;
        
        for(long int j = 0; j < subContigNumber; j++){
            
            long int subContigIndex = subContigCoordinateIndexSet[j];
            
            printf("referenceIndex:%ld\tcontigIndex:%ld\t%ld\t%ld\t%ld\t%ld\t%d\n", 
                subContigCoordinateSetHead->subContigCoordinateSet[subContigIndex].referenceIndex,
                subContigCoordinateSetHead->subContigCoordinateSet[subContigIndex].contigIndex,
                subContigCoordinateSetHead->subContigCoordinateSet[subContigIndex].referenceStart, 
                subContigCoordinateSetHead->subContigCoordinateSet[subContigIndex].referenceEnd,
                subContigCoordinateSetHead->subContigCoordinateSet[subContigIndex].contigStart, 
                subContigCoordinateSetHead->subContigCoordinateSet[subContigIndex].contigEnd, 
                subContigCoordinateSetHead->subContigCoordinateSet[subContigIndex].isReverse);
            
            
        }
        
    }
    

}
























