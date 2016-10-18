#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

using namespace std;

typedef struct ContigSet{
    char * contig;
    long int contigLength;
    char * contigName;
    ContigSet(){
        contig = NULL;
        contigLength = 0;
        contigName = NULL;
    }

}ContigSet;

ContigSet * GetContigSet(char * contigFileName, long int & contigNumber){

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
    }
    
    if((fp = fopen(contigFileName, "r")) == NULL){
        printf("%s, does not exist!", contigFileName);
        exit(0);
    }
    i = 0;
    while((fgets(contig, maxSize, fp)) != NULL){ 
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
    free(contig);
    return contigSet;
}

long int SearchReferenceIndexFromName(ContigSet * reference, long int referenceNumber, char * referenceName){ 

    for(long int i = 0; i < referenceNumber; i++){
        if(strlen(reference[i].contigName) < strlen(referenceName)){
            continue;
        }
        if(strlen(reference[i].contigName) == strlen(referenceName)){
            if(strcmp(reference[i].contigName, referenceName) == 0){
                return i;
            }
        }
        if(strlen(reference[i].contigName) > strlen(referenceName)){
            long int length = strlen(referenceName);
            long int j = 0;
            for(j = 0; j < length; j++){
                if(reference[i].contigName[j] != referenceName[j]){
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


void OutPutPerfectScaffold(ContigSet * reference, long int referenceNumber, char * coordinateFile, char * scaffoldFileName, char * gapFileName){
    
    FILE * fp;
    if((fp = fopen(coordinateFile, "r")) == NULL){
        printf("%s, does not exist!", coordinateFile);
        exit(0);
    }
    
    FILE * fp1;
    if((fp1 = fopen(scaffoldFileName, "w+")) == NULL){
        printf("%s, does not exist!", scaffoldFileName);
        exit(0);
    }
    
    FILE * fp2;
    if((fp2 = fopen(gapFileName, "w+")) == NULL){
        printf("%s, does not exist!", gapFileName);
        exit(0);
    }
    
    long int gapIndex = 0;
    long int tempEnd = -1;
    long int previousIndex = -1;
    long int scaffoldCount = 0;
    long int maxSize = 10000;
    char * line = (char *)malloc(sizeof(char)*maxSize);
    fgets(line, maxSize, fp);
    fgets(line, maxSize, fp);
    fgets(line, maxSize, fp);
    fgets(line, maxSize, fp);
    while(fgets(line, maxSize, fp)){
        char * p = NULL;
        p = strtok(line, "\t");
        long int referenceStart = -1;
        long int referenceEnd = -1;
        sscanf(p, "%ld", &referenceStart);
        p = strtok(NULL, "\t");
        sscanf(p, "%ld", &referenceEnd);
        referenceStart--;
        referenceEnd--;
        long int i = 0;
        while(i<10){
            p = strtok(NULL, "\t");
            i++;
        }
        char * referenceName = new char[strlen(p)+1];
        strcpy(referenceName, p);
        referenceName[strlen(p)] = '\0';
        long int index = SearchReferenceIndexFromName(reference, referenceNumber, referenceName);
        if(previousIndex == -1 || previousIndex != index){
            if(previousIndex != -1){
                fprintf(fp1, "\n");
            }
            fprintf(fp1, ">scaffold:%ld\n", scaffoldCount);
            scaffoldCount++;
            tempEnd = -1;
        }
            
        long int overlapLength = 0;
        if(tempEnd >= referenceEnd){
            cout<<"overlap"<<endl;
            previousIndex = index;
            continue;
        }
            
        if(tempEnd > referenceStart){
            overlapLength = tempEnd - referenceStart + 1;
        }
        
        long int lengthN = referenceStart - tempEnd - 1;
        //cout<<referenceStart<<"--"<<tempEnd<<"--"<<lengthN<<"--"<<gapIndex<<endl;
        if(lengthN > 0 && tempEnd != -1){
            for(long int t = 0; t< lengthN; t++){
                fprintf(fp1, "N");
            }
            fprintf(fp2, ">scaffold_%ld\tgap_%ld\t%ld\n", scaffoldCount-1, gapIndex, lengthN);
            gapIndex++;
            char * tempGap = (char *)malloc(sizeof(char)*(lengthN+1));
            strncpy(tempGap, reference[index].contig+tempEnd+1, lengthN);
            tempGap[lengthN] = '\0';
            fprintf(fp2, "%s\n", tempGap);
            free(tempGap);
        }
            
        char * tempSubContig = (char *)malloc(sizeof(char)*(referenceEnd-referenceStart+2 - overlapLength));
        strncpy(tempSubContig, reference[index].contig+referenceStart, referenceEnd-referenceStart+1-overlapLength);
        tempSubContig[referenceEnd-referenceStart+1-overlapLength] = '\0';
        fprintf(fp1, "%s", tempSubContig);
        tempEnd = referenceEnd;
        
        for(long int j = 0; j < referenceEnd-referenceStart+2-overlapLength; j++){
            long int gapDistance = 0;
            while(j < referenceEnd-referenceStart+2-overlapLength && (tempSubContig[j] == 'N' || tempSubContig[j] == 'n')){
                gapDistance++;
                j++;
            }
            
            if(gapDistance > 0){
                fprintf(fp2, ">scaffold_%ld\tgap_%ld\t%ld\n", scaffoldCount-1, gapIndex, gapDistance);
                for(long int t = 0; t< gapDistance; t++){
                    fprintf(fp2, "N");
                }
                fprintf(fp2, "\n");
                gapIndex++;
                j--;
            }
        }
        
        previousIndex = index;
        
    }
}


int main(int argc, char *argv[]){
    
    long int referenceNumber = 0;
    ContigSet * reference = GetContigSet(argv[1], referenceNumber);
    OutPutPerfectScaffold(reference, referenceNumber,argv[2],argv[3],argv[4]);
    
}

