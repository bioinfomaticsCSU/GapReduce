#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fillGap.h"

using namespace std;


char * GetContigFromContigSet(char * contigSetFile, long int contigIndex){
    
    char * resultContig = NULL;
    
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
    
    long int contigCount = -1;
    
    while((fgets(contig, maxSize, fp)) != NULL){ 
       
       if(contig[0] == '>'){  
           contigCount++;
           if(contigCount > contigIndex){
               break;
           }
           continue;
           
       }
       if(contigCount == contigIndex){
           
           long int extendLength = strlen(contig);
           if(contig[extendLength-1] == '\n'){
               extendLength--;
           }
           long int contigLength = 0;
           char * tempContig = NULL;
           if(resultContig != NULL){
               contigLength = strlen(resultContig);
               tempContig = (char *)malloc(sizeof(char)*(contigLength+1));
               strncpy(tempContig, resultContig, contigLength);
               free(resultContig);
                   
               resultContig = (char *)malloc(sizeof(char)*(contigLength + extendLength + 1));
               strncpy(resultContig, tempContig, contigLength);
                       
               strncpy(resultContig + contigLength, contig, extendLength);
               resultContig[contigLength + extendLength] = '\0';
               free(tempContig);
           }else{
               resultContig = (char *)malloc(sizeof(char)*(extendLength+1));
               strncpy(resultContig, contig, extendLength);
               resultContig[extendLength] = '\0';
           }       
       }
    }  
    
    fclose(fp);
    
    return resultContig;
}

double ScoreTtest(KmerReadIndex * readIndex, ReadSet * readSet, long int mean, long int std, long int gap, bool isLeft){
    
    if(readIndex == NULL){
        return -1;
    }
    
    long int avg = 0;
    long int num = 0;
    
    double tempStdScore = 0;
    
    KmerReadIndex * temp = readIndex;
    while(temp != NULL){
        avg = avg + readSet[temp->index].insertSize + gap;
        tempStdScore = tempStdScore + double(labs(readSet[temp->index].insertSize + gap - mean))/double(std);
        temp = temp->next;
        num++;
    }
    
    return tempStdScore/num;
    /*
    if(num == 1){
        return -1;
    }
    
    avg = avg/num;
    
    temp = readIndex;
    while(temp != NULL){
        std = std + pow(readSet[temp->index].insertSize + gap -avg, 2);
        temp = temp->next;
    }
    std = sqrt(std/num);
    
    double t = double(labs(mean - avg))/double(std/num);
    return t;
    */
}

double ScoreDistribution(long int * distance, long int num, long int mean, long int std){
    
    double score = 0;
    for(long int i = 0; i < num; i++){
        score = score + double(labs(distance[i] - mean))/double(std);
    }
    score = score/num;
    return score;
}



int ScoreKmerCount(long int kmerCount[][2], double score[][2], long int rowCount){
    
    for(int i = 0; i < rowCount; i++){
        if(kmerCount[i][0] == 0 || kmerCount[i][1] == 0){
            //kmerCount[i][0] = 0;
            //kmerCount[i][1] = 0;
        }else{
            //kmerCount[i][0] = kmerCount[i][0] + kmerCount[i][1];
        }
        kmerCount[i][0] = kmerCount[i][0] + kmerCount[i][1];
    }
    
    long int maxCount = -1;
    long int secondMaxCount = -1;
    long int maxRowIndex = -1;
    long int secondMaxRowIndex = -1;
    
    for(int i = 0; i < rowCount; i++){
        if(kmerCount[i][0] > maxCount){
            maxCount = kmerCount[i][0];
            maxRowIndex = i;
        }
    }
    for(int i = 0; i < rowCount; i++){
        if(kmerCount[i][0] <= maxCount && kmerCount[i][0] > secondMaxCount && i != maxRowIndex){
            secondMaxCount = kmerCount[i][0];
            secondMaxRowIndex = i;
        }
    }
    
    if(maxCount != 0 && secondMaxCount ==0){
        return maxRowIndex;
    }else if(secondMaxCount != 0 && double(maxCount)/double(secondMaxCount) >= 2 && score[maxRowIndex][0] < 2 && score[maxRowIndex][1] < 2){
        return maxRowIndex;
    }else if(secondMaxCount != 0 && double(maxCount)/double(secondMaxCount) <= 1000 && score[secondMaxRowIndex][0] < 2 && score[secondMaxRowIndex][1] < 2){
        //return secondMaxRowIndex;
    }
    return -1;
    
}

bool SingleKmerReadConsencusContig(char * read, char * kmer, char * contig, bool isOut){
    long int contigLength = strlen(contig);
    long int kmerLength = strlen(kmer);
    long int readLength = strlen(read);
    long int position = KMPIndexOfContig(read, kmer);
    long int j = 0;
    if(isOut == true){
        while(position - 1 >=0){
            if(read[position -1] != contig[contigLength - kmerLength - j]){
                return false;
            }
            j++;
            position--;
        }
    }else{
        while(position + kmerLength <= readLength - 1){
            if(read[position + kmerLength] != contig[kmerLength - 1 + j]){
                return false;    
            }
            j++;
            position++;
        }
    }
    return true;
}

bool KmerConsencus(ReadSet * readSet, long int readLength, char * kmer, KmerReadIndex * temp, char * contig, bool isOut){
    long int * readCount = (long int *)malloc(sizeof(long int)*readLength);
    long int * baseMatchCount = (long int *)malloc(sizeof(long int)*readLength);
    for(long int i = 0; i < readLength; i++){
        readCount[i] = 0;
        baseMatchCount[i] = 0;
    }
    long int contigLength = strlen(contig);
    long int kmerLength = strlen(kmer);
    long int allReadCount = 0;
    while(temp != NULL){ 
        long int position = KMPIndexOfContig(readSet[temp->index].read, kmer);
        long int j = 0;
        if(isOut == true){
            while(position - 1 >=0){
                if(readSet[temp->index].read[position -1] == contig[contigLength - kmerLength - j]){
                    baseMatchCount[j]++;
                }
                readCount[j]++;
                j++;
                position--;
            }
        }else{
            while(position + kmerLength <= readLength - 1){
                if(readSet[temp->index].read[position + kmerLength] == contig[kmerLength - 1 + j]){
                    baseMatchCount[j]++;
                }
                readCount[j]++;
                j++;
                position++;
            }
        }
        allReadCount++;
        temp = temp->next;
    }
    bool result = true;
    for(long int i = 0; i < readLength && readCount[i] >0; i++){
        if((double)(readCount[i]) / (double)(allReadCount) > 0.8 && (double)(baseMatchCount[i]) / (double)(readCount[i]) < 0.2){
            result = false;
            break;
        }
    }
    free(readCount);
    free(baseMatchCount);
    return result;
    
}

char * TranverseGraphFromNode(DBGraphHead * deBruijnGraphHead, ReadSetHead * readSetHead, KmerSetHead * largeKmerSetHead, KmerSetHead * kmerSetHead, char * startContig, long int startNodeIndex, long int gapDistance, long int & extendLength, bool isOut){
    
    if(startNodeIndex < 0){
        cout<<"startNode NULL"<<endl;
        return NULL;
    }
    cout<<"aa--"<<startNodeIndex<<endl;
    DBGraph * deBruijnGraph = deBruijnGraphHead->deBruijnGraph;
    long int graphNodeCount = deBruijnGraphHead->nodeNumber;
    long int kmerLength = kmerSetHead->kmerLength;
    long int readLength = readSetHead->readLength;
    long int largeKmerLength = largeKmerSetHead->kmerLength;
    
    char * kmer = (char *)malloc(sizeof(char)*(kmerLength + 1));
    char * largeKmer = (char *)malloc(sizeof(char)*(largeKmerLength + 1));
    cout<<"bb"<<endl;
    int * visited = (int *)malloc(sizeof(int)*graphNodeCount);
    for(long int i = 0; i < graphNodeCount; i++){
        visited[i] = 0;
    }
    cout<<"cc"<<endl;
    long int index = 1;
    long int nodeIndex = startNodeIndex;
    cout<<"outNode:"<<nodeIndex<<"--";
    
    long int contigLength = 0;
    char * contig = NULL;
    if(startContig == NULL){
        contigLength = strlen(deBruijnGraph[nodeIndex].contig);
        cout<<"dd--"<<contigLength<<endl;
        contig = (char *)malloc(sizeof(char)*(contigLength + 1));
        cout<<"ff--"<<contigLength<<endl;
        strncpy(contig, deBruijnGraph[nodeIndex].contig, contigLength);
    }else{
        contigLength = strlen(startContig);
        cout<<"dd--"<<contigLength<<endl;
        contig = (char *)malloc(sizeof(char)*(contigLength + 1));
        cout<<"ff--"<<contigLength<<endl;
        strncpy(contig, startContig, contigLength);
    }
    contig[contigLength] = '\0';
    
    cout<<"gg--"<<contig<<endl;
    
    //cout<<"outNode:"<<nodeIndex;
    while(index == 1){
        cout<<"index"<<endl;
        if(visited[nodeIndex] > 2){
            break;
        }
        visited[nodeIndex]++;
         
        long int nodeCount = 0;
        if(isOut == true){
            nodeCount = NodeCount(deBruijnGraph[nodeIndex].outNode);
        }else{
            nodeCount = NodeCount(deBruijnGraph[nodeIndex].inNode);
        }
       
        if(nodeCount == 0){
            break;
        }
        if(nodeCount == 1){
            long int previousNodeIndex = nodeIndex;
            if(isOut == true){
                nodeIndex = deBruijnGraph[nodeIndex].outNode->index;
            }else{
                nodeIndex = deBruijnGraph[nodeIndex].inNode->index;
            }
            
            if(visited[nodeIndex] > 2){
                break;
            }

            long int nodeLength = strlen(deBruijnGraph[nodeIndex].contig);
            long int tempContigLength = contigLength + nodeLength - kmerLength + 1;
            extendLength = extendLength + nodeLength - kmerLength + 1;
            char * tempContig = new char[tempContigLength + 1];
            if(isOut == true){
                strncpy(tempContig, contig, contigLength);
                strncpy(tempContig + contigLength, deBruijnGraph[nodeIndex].contig + kmerLength - 1, nodeLength - kmerLength + 1);
            }else{
                strncpy(tempContig, deBruijnGraph[nodeIndex].contig, nodeLength - kmerLength + 1);
                strncpy(tempContig + nodeLength - kmerLength + 1, contig, contigLength);
            }
            
            tempContig[tempContigLength] = '\0';
            free(contig);
            contig = tempContig;
            tempContig = NULL;
            contigLength = tempContigLength;
            continue;
        }
        cout<<"nodeIndex:"<<nodeIndex<<"--"<<nodeCount<<endl;
        double scoreTtestOfLargeKmer[nodeCount][2];
        double scoreTtest[nodeCount][2];
        long int largeKmerCount[nodeCount][2];
        long int kmerCount[nodeCount][2];
        long int * tempNodeIndex = new long int[nodeCount];
        GraphNode * tempGraphNode = NULL;
        if(isOut == true){
            tempGraphNode = deBruijnGraph[nodeIndex].outNode;
        }else{
            tempGraphNode = deBruijnGraph[nodeIndex].inNode;
        }
        for(long int i = 0; i < nodeCount; i++){
            scoreTtestOfLargeKmer[i][0] = 0;
            scoreTtest[i][1] = 0;
            largeKmerCount[i][0] = 0;
            largeKmerCount[i][1] = 0;
            kmerCount[i][0] = 0;
            kmerCount[i][1] = 0;
            tempNodeIndex[i] = tempGraphNode->index;
            
            if(contigLength >= largeKmerLength -1){
                if(isOut == true){
                    strncpy(largeKmer, contig + contigLength - largeKmerLength + 1, largeKmerLength-1);
                    largeKmer[largeKmerLength - 1] = deBruijnGraph[tempGraphNode->index].contig[kmerLength - 1];
                }else{
                    largeKmer[0] = deBruijnGraph[tempGraphNode->index].contig[strlen(deBruijnGraph[tempGraphNode->index].contig) - kmerLength];
                    strncpy(largeKmer+1, contig, largeKmerLength-1);
                }
                
                largeKmer[largeKmerLength] = '\0';
                cout<<"dfdfd--"<<largeKmer<<endl;
                long int largeHashKey = SearchKmerOfKmerSet(largeKmer, largeKmerSetHead->leftKmerSet, largeKmerLength, largeKmerSetHead->leftKmerCount);
                cout<<"dfdfd--"<<largeKmer<<endl;
                long int largeHashKey1 = SearchKmerOfKmerSet(largeKmer, largeKmerSetHead->rightKmerSet, largeKmerLength, largeKmerSetHead->rightKmerCount);
                if(gapDistance - extendLength > readSetHead->insertSize + 2*readSetHead->std){
                    if(isOut == true){
                        largeHashKey1 = -1;
                    }else{
                        largeHashKey = -1;
                    }
                }
                if(largeHashKey != -1){
                    largeKmerCount[i][0] = largeKmerSetHead->leftKmerSet[largeHashKey].kmerCount;
                    KmerReadIndex * temp = largeKmerSetHead->leftKmerSet[largeHashKey].readIndex;
                    cout<<"largeDistance:";
                    long int distance = 0;
                    long int * tempInsertSize = (long int *)malloc(sizeof(long int)*largeKmerCount[i][0]);
                    long int p = 0;
                    while(temp != NULL){
                        
                        if(false == SingleKmerReadConsencusContig(readSetHead->leftReadSet[temp->index].read, largeKmer, contig, isOut)){
                            temp = temp->next;
                            continue;
                        }
                        
                        long int position = KMPIndexOfContig(readSetHead->leftReadSet[temp->index].read, largeKmer);
                        if(isOut == true){
                            tempInsertSize[p] = readSetHead->leftReadSet[temp->index].insertSize + extendLength + readSetHead->readLength - position - largeKmerLength + 1;
                        }else{
                            tempInsertSize[p] = readSetHead->leftReadSet[temp->index].insertSize + gapDistance - extendLength + readSetHead->readLength - position - largeKmerLength + 1;
                        }
                        
                        cout<<tempInsertSize[p]<<"--";
                        distance = distance +  tempInsertSize[p]; 
                        temp = temp->next;
                        p++;
                    }
                    /*
                    cout<<p<<"--ss--"<<largeKmerCount[i][0]<<endl; 
                    if(double(p)/double(largeKmerCount[i][0]) < 0.2){
                        largeKmerCount[i][0] = 0;
                    }else{
                        largeKmerCount[i][0] = p;
                    }
                    */
                    largeKmerCount[i][0] = p;
                    //cout<<distance/largeKmerCount[i][0]<<endl;
                    if(largeKmerCount[i][0] >= 1){
                        scoreTtestOfLargeKmer[i][0] = ScoreDistribution(tempInsertSize, largeKmerCount[i][0], readSetHead->insertSize, readSetHead->std);
                    }
                    free(tempInsertSize);
                    cout<<"Ttest:"<<scoreTtestOfLargeKmer[i][0]<<endl;
                }
                if(largeHashKey1 != -1){
                    largeKmerCount[i][1] = largeKmerSetHead->rightKmerSet[largeHashKey1].kmerCount;
                    KmerReadIndex * temp = largeKmerSetHead->rightKmerSet[largeHashKey1].readIndex;
                    cout<<"largeDistance:";
                    long int distance = 0;
                    long int * tempInsertSize = (long int *)malloc(sizeof(long int)*largeKmerCount[i][1]);
                    long int p = 0;
                    while(temp != NULL){
                        
                        if(false == SingleKmerReadConsencusContig(readSetHead->rightReadSet[temp->index].read, largeKmer, contig, isOut)){
                            temp = temp->next;
                            continue;
                        }
                        
                        long int position = KMPIndexOfContig(readSetHead->rightReadSet[temp->index].read, largeKmer);
                        if(isOut == true){
                            tempInsertSize[p] = readSetHead->rightReadSet[temp->index].insertSize + gapDistance - extendLength + position + largeKmerLength;
                        }else{
                            tempInsertSize[p] = readSetHead->rightReadSet[temp->index].insertSize + extendLength + position + 1;
                        }
                        
                        cout<<tempInsertSize[p]<<"--";
                        distance = distance + tempInsertSize[p]; 
                        temp = temp->next;
                        p++;
                    }
                    /*
                    if(double(p)/double(largeKmerCount[i][1]) < 0.2){
                        largeKmerCount[i][1] = 0;
                    }else{
                        largeKmerCount[i][1] = p;
                    }
                    */
                    cout<<p<<"--s11s--"<<largeKmerCount[i][1]<<endl; 
                    //cout<<distance/largeKmerCount[i][1]<<endl;
                    largeKmerCount[i][1] = p;
                    if(largeKmerCount[i][1] >= 1){
                        scoreTtestOfLargeKmer[i][1] = ScoreDistribution(tempInsertSize, largeKmerCount[i][1], readSetHead->insertSize, readSetHead->std);
                    }
                    free(tempInsertSize);
                    cout<<"Ttest:"<<scoreTtestOfLargeKmer[i][1]<<endl;
                }
                cout<<largeKmer<<endl;
            }
            
            if(isOut == true){
                strncpy(kmer, deBruijnGraph[nodeIndex].contig + strlen(deBruijnGraph[nodeIndex].contig) - kmerLength + 1, kmerLength-1);
                kmer[kmerLength - 1] = deBruijnGraph[tempGraphNode->index].contig[kmerLength - 1];
            }else{
                kmer[0] = deBruijnGraph[tempGraphNode->index].contig[strlen(deBruijnGraph[tempGraphNode->index].contig) - kmerLength];
                strncpy(kmer+1, contig, kmerLength-1);
            }
            kmer[kmerLength] = '\0';
            cout<<"kmer:"<<kmer<<endl;
            long int hashKey = SearchKmerOfKmerSet(kmer, kmerSetHead->leftKmerSet, kmerLength, kmerSetHead->leftKmerCount);
            long int hashKey1 = SearchKmerOfKmerSet(kmer, kmerSetHead->rightKmerSet, kmerLength, kmerSetHead->rightKmerCount);
            if(gapDistance - extendLength > readSetHead->insertSize + 2*readSetHead->std){
                if(isOut == true){
                    hashKey1 = -1;
                }else{
                    hashKey = -1;
                }
            }
            if(hashKey != -1){
                kmerCount[i][0] = kmerSetHead->leftKmerSet[hashKey].kmerCount;
                KmerReadIndex * temp = kmerSetHead->leftKmerSet[hashKey].readIndex;
                cout<<"Distance:";
                long int distance = 0;
                long int * tempInsertSize = (long int *)malloc(sizeof(long int)*kmerCount[i][0]);
                long int p = 0;
                while(temp != NULL){
                    
                    if(false == SingleKmerReadConsencusContig(readSetHead->leftReadSet[temp->index].read, kmer, contig, isOut)){
                        temp = temp->next;
                        continue;
                    }
                    
                    long int position = KMPIndexOfContig(readSetHead->leftReadSet[temp->index].read, kmer);
                    if(isOut == true){
                        tempInsertSize[p] = readSetHead->leftReadSet[temp->index].insertSize + extendLength + readSetHead->readLength - position - kmerLength + 1;
                    }else{
                        tempInsertSize[p] = readSetHead->leftReadSet[temp->index].insertSize + gapDistance - extendLength + readSetHead->readLength - position - kmerLength + 1;
                    }
                    cout<<tempInsertSize[p]<<"--";
                    distance = distance +  tempInsertSize[p]; 
                    temp = temp->next;
                    p++;
                }
                
                if(double(p)/double(kmerCount[i][0]) < 0.2){
                    kmerCount[i][0] = 0;
                }else{
                    kmerCount[i][0] = p;
                }
                
                cout<<p<<"--s22s--"<<kmerCount[i][0]<<endl; 
                //cout<<distance/kmerCount[i][0]<<endl;
                kmerCount[i][0] = p;
                if(kmerCount[i][0] >= 1){
                    scoreTtest[i][0] = ScoreDistribution(tempInsertSize, kmerCount[i][0], readSetHead->insertSize, readSetHead->std);
                }else{
                    //scoreTtest[i][0] = 100;
                }
                free(tempInsertSize);
                cout<<"Ttest:"<<scoreTtest[i][0]<<endl;
            }
            if(hashKey1 != -1){
                kmerCount[i][1] = kmerSetHead->rightKmerSet[hashKey1].kmerCount;
                KmerReadIndex * temp = kmerSetHead->rightKmerSet[hashKey1].readIndex;
                cout<<"Distance:";
                long int distance = 0;
                long int * tempInsertSize = (long int *)malloc(sizeof(long int)*kmerCount[i][1]);
                long int p = 0;
                while(temp != NULL){
                    
                    if(false == SingleKmerReadConsencusContig(readSetHead->rightReadSet[temp->index].read, kmer, contig, isOut)){
                        temp = temp->next;
                        continue;
                    }
                    
                    long int position = KMPIndexOfContig(readSetHead->rightReadSet[temp->index].read, kmer);
                    if(isOut == true){
                        tempInsertSize[p] = readSetHead->rightReadSet[temp->index].insertSize + gapDistance - extendLength + position + kmerLength;
                    }else{
                        tempInsertSize[p] = readSetHead->rightReadSet[temp->index].insertSize + extendLength + position + 1;
                    }
                    cout<<tempInsertSize[p]<<"--";
                    distance = distance +  tempInsertSize[p]; 
                    temp = temp->next;
                    p++;
                }
                /*
                if(double(p)/double(kmerCount[i][1]) < 0.2){
                    kmerCount[i][1] = 0;
                }else{
                    kmerCount[i][1] = p;
                }
                */
                cout<<p<<"--s33s--"<<kmerCount[i][1]<<endl; 
                //cout<<distance/kmerCount[i][1]<<endl;
                kmerCount[i][1] = p;
                if(kmerCount[i][1] >= 1){
                    scoreTtest[i][1] = ScoreDistribution(tempInsertSize, kmerCount[i][1], readSetHead->insertSize, readSetHead->std);
                }else{
                    //scoreTtest[i][1] = 100;
                }
                free(tempInsertSize);
                cout<<"Ttest:"<<scoreTtest[i][1]<<endl;
            }
            cout<<kmer<<endl; 
            tempGraphNode = tempGraphNode->next;
        }
        
        for(long int i = 0; i < nodeCount; i++){
            cout<<largeKmerCount[i][0]<<"--"<<largeKmerCount[i][1]<<endl;
            if(scoreTtestOfLargeKmer[i][0] > 2 || scoreTtestOfLargeKmer[i][1] > 2){
                //largeKmerCount[i][0] = 0;
                //largeKmerCount[i][1] = 0;
            }
            
        }
        int maxIndexOfLargeKmer = ScoreKmerCount(largeKmerCount, scoreTtestOfLargeKmer, nodeCount);
        cout<<"largeMax--"<<maxIndexOfLargeKmer<<endl;
        for(long int i = 0; i < nodeCount; i++){
            cout<<kmerCount[i][0]<<"--"<<kmerCount[i][1]<<endl;
            if(scoreTtest[i][0] > 2 || scoreTtest[i][1] > 2){
                //kmerCount[i][0] = 0;
                //kmerCount[i][1] = 0;
            }
        }
        int maxIndex = ScoreKmerCount(kmerCount, scoreTtest, nodeCount);
        cout<<"max--"<<maxIndex<<endl;
        
        long int tempMaxIndex = -1;
        if(maxIndexOfLargeKmer != -1){
            tempMaxIndex = maxIndexOfLargeKmer;
        }else if(maxIndex != -1){
            tempMaxIndex = maxIndex;
        }
        
        
        if(tempMaxIndex != -1){
            nodeIndex = tempNodeIndex[tempMaxIndex];
            if(visited[nodeIndex] > 2){
                break;
            }
            
            long int nodeLength = strlen(deBruijnGraph[nodeIndex].contig);
            long int tempContigLength = contigLength + nodeLength - kmerLength + 1;
            extendLength = extendLength + nodeLength - kmerLength + 1;
            char * tempContig = new char[tempContigLength + 1];
            if(isOut == true){
                strncpy(tempContig, contig, contigLength);
                strncpy(tempContig + contigLength, deBruijnGraph[nodeIndex].contig + kmerLength - 1, nodeLength - kmerLength + 1);
            }else{
                strncpy(tempContig, deBruijnGraph[nodeIndex].contig, nodeLength - kmerLength + 1);
                strncpy(tempContig + nodeLength - kmerLength + 1, contig, contigLength);
            }
            tempContig[tempContigLength] = '\0';
            free(contig);
            contig = tempContig;
            tempContig = NULL;
            contigLength = tempContigLength;
            
            //cout<<"--"<<nodeIndex;
            //cout<<"m"<<nodeIndex<<"--";
            continue;
        }else{
            //cout<<endl;
            break;
        }
 
    }
    cout<<"path_end"<<endl;
    free(kmer);
    free(largeKmer);
    free(visited);
    return contig;
}

char * GetGapRegionFromSingle(char * gap, long int extendLength, long int gapDistance, bool isLeft){
    
    if(gap == NULL || extendLength <= 0){
        return NULL;
    }
    long int gapLength = strlen(gap);
    char * gapRegion = (char *)malloc(sizeof(char)*(gapDistance + 1));
    
    if(extendLength <= gapDistance){
        if(isLeft == true){
            strncpy(gapRegion, gap + gapLength - extendLength, extendLength);
            for(long int i = 0; i < gapDistance - extendLength; i++){
                *(gapRegion + extendLength + i) = 'N';
            } 
        }else{
            for(long int i = 0; i < gapDistance - extendLength; i++){
                *(gapRegion + i) = 'N';
            } 
            strncpy(gapRegion + gapDistance - extendLength, gap, extendLength);
        }
        gapRegion[gapDistance] = '\0';
        
    }else{
        if(isLeft == true){
            strncpy(gapRegion, gap + gapLength - extendLength, gapDistance);
        }else{
            strncpy(gapRegion, gap + extendLength - gapDistance, gapDistance);
        }
        gapRegion[gapDistance] = '\0';
    }
    return gapRegion;

}

char * GetGapRegion(char * leftGap, long int leftExtendLength, char * rightGap, long int rightExtendLength, long int gapDistance, long int cutLength){
    
    long int leftGapLength = strlen(leftGap);
    long int rightGapLength = strlen(rightGap);
    
    if(leftExtendLength <= 0){
        leftExtendLength = 0;
    }
    if(rightExtendLength <= 0){
        rightExtendLength = 0;
    }
    
    char * gapRegion = (char *)malloc(sizeof(char)*(gapDistance + 1));
    
    if(leftExtendLength + rightExtendLength <= gapDistance){
        strncpy(gapRegion, leftGap + leftGapLength - leftExtendLength, leftExtendLength);
        for(long int i = 0; i < gapDistance - leftExtendLength - rightExtendLength; i++){
            *(gapRegion + leftExtendLength + i) = 'N';
        }
        strncpy(gapRegion + gapDistance - rightExtendLength, rightGap, rightExtendLength);
        gapRegion[gapDistance] = '\0';
        cout<<leftGap<<endl;
        cout<<rightGap<<endl;
        cout<<"fill--1"<<endl;
        return gapRegion;
    }
    if(leftExtendLength >= gapDistance && rightExtendLength >= gapDistance){
        /*
        if(leftExtendLength > rightExtendLength){
            strncpy(gapRegion, leftGap + leftGapLength - leftExtendLength, gapDistance);
            gapRegion[gapDistance] = '\0';
            return gapRegion;
        }else{
            strncpy(gapRegion, rightGap + rightExtendLength - gapDistance, gapDistance);
            gapRegion[gapDistance] = '\0';
            return gapRegion;
        }
        */
        int t = 0;
        if(gapDistance%2!=0){
            t = 1;
        }
        strncpy(gapRegion, leftGap + leftGapLength - leftExtendLength, (gapDistance/2)+t);
        strncpy(gapRegion, rightGap + rightExtendLength - gapDistance/2, gapDistance/2);
        
        
    }
    if(leftExtendLength >= gapDistance){
        strncpy(gapRegion, leftGap + leftGapLength - leftExtendLength, gapDistance);
        for(long int i = 0; i < gapDistance/2 && i < rightExtendLength; i++){
            gapRegion[gapDistance - 1 - i] = rightGap[rightExtendLength - 1 - i];
        }
        gapRegion[gapDistance] = '\0';
        return gapRegion;
    }
    
    if(rightExtendLength >= gapDistance){
        strncpy(gapRegion, rightGap + rightExtendLength - gapDistance, gapDistance);
        for(long int i = 0; i < gapDistance/2 && i < leftExtendLength; i++){
            gapRegion[i] = leftGap[leftGapLength - leftExtendLength + i];
        }
        gapRegion[gapDistance] = '\0';
        return gapRegion;
    }
    
    if(leftExtendLength + rightExtendLength > gapDistance){
        long int overlapLength = leftExtendLength + rightExtendLength - gapDistance;
        strncpy(gapRegion, leftGap + leftGapLength - leftExtendLength, leftExtendLength);
        strncpy(gapRegion + leftExtendLength, rightGap + overlapLength, rightExtendLength - overlapLength);
        gapRegion[gapDistance] = '\0';
        return gapRegion;
    }
    
    for(long int i = 0; i < gapDistance; i++){
        gapRegion[i] = 'N';
    }
    gapRegion[gapDistance] = '\0';
    return gapRegion;
}

char * ExtractPathFromGraph(char * leftContig, char * rightContig, long int gapDistance, DBGraphHead * deBruijnGraphHead, ReadSetHead * readSetHead, KmerSetHead * kmerSetHead, KmerSetHead * largeKmerSetHead){
    
    long int kmerLength = kmerSetHead->kmerLength;
    long int cutLength = 5;
    if(leftContig == NULL){
        leftContig = (char *)malloc(sizeof(char)*(kmerLength*2+1));
        for(long int i = 0; i < 2*kmerLength+1;i++){
            leftContig[i] = 'G';
        }
        leftContig[2*kmerLength] = '\0';
    }
    if(rightContig == NULL){
        rightContig = (char *)malloc(sizeof(char)*(kmerLength*2+1));
        for(long int i = 0; i < 2*kmerLength+1;i++){
            rightContig[i] = 'G';
        }
        rightContig[2*kmerLength] = '\0';
    }
    
    long int leftContigLength = strlen(leftContig);
    long int rightContigLength = strlen(rightContig);
    
    long int tempLeftNodeIndex = -1;
    long int tempRightNodeIndex = -1;
    
    char * tempLeftKmer = new char[kmerLength + 1];
    char * tempRightKmer = new char[kmerLength + 1];
    long int tempLeftKmerMapPosition = -1;
    long int tempRightKmerMapPosition = -1;
    char * leftStartContig = NULL;
    char * rightStartContig = NULL;
    strncpy(tempLeftKmer, leftContig + leftContigLength - kmerLength - cutLength, kmerLength);
    strncpy(tempRightKmer, rightContig + cutLength, kmerLength);
    tempLeftKmer[kmerLength] = '\0';
    tempRightKmer[kmerLength] = '\0';
    
    cout<<tempLeftKmer<<"--"<<tempRightKmer<<endl;
    tempLeftNodeIndex = FindKmerInDeBruijnGraph(tempLeftKmer, kmerLength, deBruijnGraphHead->deBruijnGraph, deBruijnGraphHead->nodeNumber, tempLeftKmerMapPosition);
    tempRightNodeIndex = FindKmerInDeBruijnGraph(tempRightKmer, kmerLength, deBruijnGraphHead->deBruijnGraph, deBruijnGraphHead->nodeNumber, tempRightKmerMapPosition);
    cout<<tempLeftNodeIndex<<"--"<<tempRightNodeIndex<<endl;
    long int leftExtendLength = 0;
    if(tempLeftNodeIndex >= 0){
        leftExtendLength = strlen(deBruijnGraphHead->deBruijnGraph[tempLeftNodeIndex].contig) - tempLeftKmerMapPosition - kmerLength - cutLength;
        cout<<"leftNodeLength:"<<strlen(deBruijnGraphHead->deBruijnGraph[tempLeftNodeIndex].contig)<<endl;
        if(strlen(deBruijnGraphHead->deBruijnGraph[tempLeftNodeIndex].contig) < readSetHead->readLength){
            leftStartContig = (char *)malloc(sizeof(char)*(readSetHead->readLength + 1));
            strncpy(leftStartContig, leftContig + leftContigLength - cutLength - kmerLength - tempLeftKmerMapPosition - readSetHead->readLength + strlen(deBruijnGraphHead->deBruijnGraph[tempLeftNodeIndex].contig), readSetHead->readLength - strlen(deBruijnGraphHead->deBruijnGraph[tempLeftNodeIndex].contig));
            strncpy(leftStartContig + readSetHead->readLength - strlen(deBruijnGraphHead->deBruijnGraph[tempLeftNodeIndex].contig), deBruijnGraphHead->deBruijnGraph[tempLeftNodeIndex].contig, strlen(deBruijnGraphHead->deBruijnGraph[tempLeftNodeIndex].contig));
            leftStartContig[readSetHead->readLength] = '\0';
        }
    }
    long int rightExtendLength = 0;
    if(tempRightNodeIndex >= 0){
        rightExtendLength = tempRightKmerMapPosition - cutLength;
        cout<<"rightNodeLength:"<<strlen(deBruijnGraphHead->deBruijnGraph[tempRightNodeIndex].contig)<<endl;
        if(strlen(deBruijnGraphHead->deBruijnGraph[tempRightNodeIndex].contig) < readSetHead->readLength){
            rightStartContig = (char *)malloc(sizeof(char)*(readSetHead->readLength + 1));
            strncpy(rightStartContig, deBruijnGraphHead->deBruijnGraph[tempRightNodeIndex].contig, strlen(deBruijnGraphHead->deBruijnGraph[tempRightNodeIndex].contig));
            strncpy(rightStartContig + strlen(deBruijnGraphHead->deBruijnGraph[tempRightNodeIndex].contig), rightContig + cutLength + strlen(deBruijnGraphHead->deBruijnGraph[tempRightNodeIndex].contig) - tempRightKmerMapPosition, readSetHead->readLength - strlen(deBruijnGraphHead->deBruijnGraph[tempRightNodeIndex].contig));
            rightStartContig[readSetHead->readLength] = '\0';
        }
    }
    if(leftStartContig != NULL){
        //cout<<leftContig<<endl;
        cout<<"leftStartContig"<<endl;
        cout<<deBruijnGraphHead->deBruijnGraph[tempLeftNodeIndex].contig<<endl;
        cout<<leftStartContig<<endl;
    }
    if(rightStartContig != NULL){
        //cout<<rightContig<<endl;
        cout<<"rightStartContig"<<endl;
        cout<<deBruijnGraphHead->deBruijnGraph[tempRightNodeIndex].contig<<endl;
        cout<<rightStartContig<<endl;
    }
    char * leftGap = NULL;
    char * rightGap = NULL;
    cout<<"extendLength:"<<leftExtendLength<<"--"<<rightExtendLength<<endl;
    leftGap = TranverseGraphFromNode(deBruijnGraphHead, readSetHead, largeKmerSetHead, kmerSetHead, leftStartContig, tempLeftNodeIndex, gapDistance, leftExtendLength, 1);
    cout<<"tt"<<endl;
    if(leftGap != NULL){
        cout<<"gapRegionLeft:"<<leftGap<<endl;
        
    }
    cout<<"ff"<<endl;
    rightGap = TranverseGraphFromNode(deBruijnGraphHead, readSetHead, largeKmerSetHead, kmerSetHead, rightStartContig, tempRightNodeIndex, gapDistance, rightExtendLength, 0);
    
    
    if(rightGap != NULL){
        cout<<"gapRegionRight:"<<rightGap<<endl;
        
    }
    cout<<"extendLength:"<<leftExtendLength<<"--"<<rightExtendLength<<endl;
    char * gapRegion = NULL;
    if(leftGap != NULL && rightGap != NULL){
        gapRegion = GetGapRegion(leftGap, leftExtendLength, rightGap, rightExtendLength, gapDistance, cutLength);
    }else if(leftGap != NULL && rightGap == NULL){
        gapRegion = GetGapRegionFromSingle(leftGap, leftExtendLength, gapDistance, 1);
    }else if(leftGap == NULL && rightGap != NULL){
        gapRegion = GetGapRegionFromSingle(rightGap, rightExtendLength, gapDistance, 0);
    }
    if(gapRegion == NULL){
        gapRegion = (char *)malloc(sizeof(char)*(gapDistance + 1));
        for(long int i = 0; i < gapDistance; i++){
            gapRegion[i] = 'N';
        }
        gapRegion[gapDistance] = '\0';
    }

    cout<<"finalGapRegion:"<<gapRegion<<endl;
    cout<<"extract_path_end!"<<endl;
    if(leftGap != NULL){
        cout<<"aa"<<endl;
        free(leftGap);
    }
    if(rightGap != NULL){
        cout<<"bb"<<endl;
        free(rightGap);
    }
    
    cout<<"extract_path_end!"<<endl;
    return gapRegion;
}




