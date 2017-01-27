#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "GapRead.h"

using namespace std;
using namespace BamTools;


GapRegionReadSet * GetReadSetInGapRegion(ScaffoldSetHead * scaffoldSetHead, char * endContigFile, char * endContigBamLeft, char * endContigBamRight, long int maxInsertSize, long int minInsertSize, long int intervalLength, bool isPairedRead){
    

    ContigSetHead * endContigSetHead = GetContigSet(endContigFile);
    ContigSet * endContigSet = endContigSetHead->contigSet;

    ScaffoldSet * scaffoldSet = scaffoldSetHead->scaffoldSet;
    GapToContigIndex * gapToContigIndex = scaffoldSetHead->gapToContigIndex;
    ContigToGapIndex * contigToGapIndex = scaffoldSetHead->contigToGapIndex;
    long int scaffoldCount = scaffoldSetHead->scaffoldCount;
    long int contigCount = scaffoldSetHead->contigCount;
    long int gapCount = scaffoldSetHead->gapCount;
    
    GapRegionReadSet * gapRegionReadSet = (GapRegionReadSet *)malloc(sizeof(GapRegionReadSet));
    gapRegionReadSet->gapCount = gapCount;
    gapRegionReadSet->singleMapReadSet = (GapRegionSingleMapReadSet *)malloc(sizeof(GapRegionSingleMapReadSet)*gapCount);
    gapRegionReadSet->pairedMapReadSet = (GapRegionPairedMapReadSet *)malloc(sizeof(GapRegionPairedMapReadSet)*gapCount);
    //cout<<"ss"<<endl;
    for(long int i = 0; i < gapCount; i++){
        gapRegionReadSet->singleMapReadSet[i].leftContigSingleReadSet = NULL;
        gapRegionReadSet->singleMapReadSet[i].rightContigSingleReadSet = NULL;
        gapRegionReadSet->singleMapReadSet[i].leftContigSingleReadCount = 0;
        gapRegionReadSet->singleMapReadSet[i].rightContigSingleReadCount = 0;
        gapRegionReadSet->pairedMapReadSet[i].pairedReadSet = NULL;
        gapRegionReadSet->pairedMapReadSet[i].pairedReadCount = 0;
        
    }
    
    //cout<<"tt"<<endl;
    
    BamReader bamReaderLeft;
    BamReader bamReaderRight;
    string bamFileNameLeft = endContigBamLeft;
    string bamFileNameRight = endContigBamRight;
    bamReaderLeft.Open(bamFileNameLeft);
    bamReaderRight.Open(bamFileNameRight);
    
    BamAlignment alignmentLeft;
    BamAlignment alignmentRight;
    long int t = 0;
    while(bamReaderLeft.GetNextAlignmentCore(alignmentLeft) && bamReaderRight.GetNextAlignmentCore(alignmentRight)){
        //cout<<"--"<<t++<<endl;
        while((alignmentLeft.AlignmentFlag & 0x900) != 0){
            bamReaderLeft.GetNextAlignmentCore(alignmentLeft);
            continue;
        } 
        while((alignmentRight.AlignmentFlag & 0x900) != 0){
            bamReaderRight.GetNextAlignmentCore(alignmentRight);
            continue;
        }
        
        bool leftSoftClip = false;
        bool rightSoftClip = false;
        long int leftSoftLength = 0;
        long int rightSoftLength = 0;
        for(int t = 0; t < 2; t++){
            std::vector< int > clipSizes;
            std::vector< int > readPositions;
            std::vector< int > genomePositions;
            BamAlignment tempAlignment;
            if(t == 0){
                tempAlignment = alignmentLeft;
            }else{
                tempAlignment = alignmentRight;
            }
            bool tempSoftClip = tempAlignment.GetSoftClips(clipSizes, readPositions, genomePositions);
            if(tempSoftClip != false){
                long int leftPosition = GetRealMappingPositionFromEndContig(tempAlignment.Position, endContigSet[tempAlignment.RefID].contigLength, maxInsertSize, intervalLength, 0);
                if(leftPosition ==1 && readPositions[0] == clipSizes[0]){
                    long int gapIndex = contigToGapIndex[tempAlignment.RefID].leftGapIndex;
                    if(gapIndex >=0){
                        gapRegionReadSet->singleMapReadSet[gapIndex].rightContigSingleReadCount++;
                        leftSoftLength = clipSizes[0];
                    }else{
                        tempSoftClip = false;
                    }
                    
                    
                }else if(leftPosition > endContigSet[tempAlignment.RefID].contigLength - tempAlignment.Length 
                      && readPositions[0] + clipSizes[0] == tempAlignment.Length){
                    
                    long int gapIndex = contigToGapIndex[tempAlignment.RefID].rightGapIndex;
                    if(gapIndex >= 0){
                        gapRegionReadSet->singleMapReadSet[gapIndex].leftContigSingleReadCount++;
                        rightSoftLength = clipSizes[0];
                    }else{
                        //tempSoftClip = false;
                    }
                    
                    
                }else{
                    //tempSoftClip = false;
                }
                if(t == 0){
                    leftSoftClip = tempSoftClip;
                }else{
                    rightSoftClip = tempSoftClip;
                }
            }
            vector<int>().swap(clipSizes);
            vector<int>().swap(readPositions);
            vector<int>().swap(genomePositions);
        }
        
        if((alignmentLeft.IsMapped() && alignmentRight.IsMapped()) ){                                    
            //cout<<"aa--"<<endl;
            //continue;
            if(abs(alignmentLeft.RefID - alignmentRight.RefID) != 1 || abs(alignmentLeft.RefID - alignmentRight.RefID) == 1){
                //cout<<"ff"<<endl;
                long int leftPosition = GetRealMappingPositionFromEndContig(alignmentLeft.Position, endContigSet[alignmentLeft.RefID].contigLength, maxInsertSize, intervalLength, leftSoftLength);
                MapRegionGapIndex * temp = GetGapReadFromReadPositionOfContig(scaffoldSetHead, alignmentLeft.Length, alignmentLeft.RefID, leftPosition, maxInsertSize, minInsertSize, isPairedRead, alignmentLeft.IsReverseStrand());
                //cout<<"ff--end"<<endl;
                if((!alignmentLeft.IsReverseStrand() && isPairedRead) ||(alignmentLeft.IsReverseStrand() && !isPairedRead)){
                    
                    while(temp!=NULL){
                        //cout<<"ttt--"<<temp->gapIndex<<endl;
                        gapRegionReadSet->singleMapReadSet[contigToGapIndex[alignmentLeft.RefID].rightGapIndex + temp->gapIndex].leftContigSingleReadCount++;
                        temp = temp->next;
                    }
                    /*
                    if(contigToGapIndex[alignmentLeft.RefID].rightGapIndex != -1 
                        && alignmentLeft.Position > endContigSet[alignmentLeft.RefID].contigLength - maxInsertSize
                        && alignmentLeft.Position < endContigSet[alignmentLeft.RefID].contigLength - minInsertSize + gapToContigIndex[contigToGapIndex[alignmentLeft.RefID].rightGapIndex].gapDistance
                        ){
                        cout<<"xxx"<<endl;
                        gapRegionReadSet->singleMapReadSet[contigToGapIndex[alignmentLeft.RefID].rightGapIndex].leftContigSingleReadCount++;
                    }
                    */
                }else{
                    
                    while(temp!=NULL){
                        //cout<<"ttt--"<<temp->gapIndex<<endl;
                        gapRegionReadSet->singleMapReadSet[contigToGapIndex[alignmentLeft.RefID].leftGapIndex - temp->gapIndex].rightContigSingleReadCount++;
                        temp = temp->next;
                    }
                    /*
                    if(contigToGapIndex[alignmentLeft.RefID].leftGapIndex != -1 
                        && alignmentLeft.Position < maxInsertSize
                        && alignmentLeft.Position > minInsertSize - gapToContigIndex[contigToGapIndex[alignmentLeft.RefID].leftGapIndex].gapDistance
                        ){
                        gapRegionReadSet->singleMapReadSet[contigToGapIndex[alignmentLeft.RefID].leftGapIndex].rightContigSingleReadCount++;
                    }
                    */       
                }
                
                //cout<<"gg"<<endl;
                long int rightPosition = GetRealMappingPositionFromEndContig(alignmentRight.Position, endContigSet[alignmentRight.RefID].contigLength, maxInsertSize, intervalLength, rightSoftLength);
                temp = GetGapReadFromReadPositionOfContig(scaffoldSetHead, alignmentRight.Length, alignmentRight.RefID, rightPosition, maxInsertSize, minInsertSize, isPairedRead, alignmentRight.IsReverseStrand());
                //cout<<"gg--end"<<endl;
                if((!alignmentRight.IsReverseStrand() && isPairedRead) ||(alignmentRight.IsReverseStrand() && !isPairedRead)){
                    while(temp!=NULL){
                        gapRegionReadSet->singleMapReadSet[contigToGapIndex[alignmentRight.RefID].rightGapIndex  + temp->gapIndex].leftContigSingleReadCount++;
                        temp = temp->next;
                    }
                    /*
                    if(contigToGapIndex[alignmentRight.RefID].rightGapIndex != -1 
                        && alignmentRight.Position > endContigSet[alignmentRight.RefID].contigLength - maxInsertSize
                        && alignmentRight.Position < endContigSet[alignmentRight.RefID].contigLength - minInsertSize + gapToContigIndex[contigToGapIndex[alignmentRight.RefID].rightGapIndex].gapDistance
                        ){
                        gapRegionReadSet->singleMapReadSet[contigToGapIndex[alignmentRight.RefID].rightGapIndex].leftContigSingleReadCount++;
                    } 
                    */
                }else{
                    while(temp!=NULL){
                        gapRegionReadSet->singleMapReadSet[contigToGapIndex[alignmentRight.RefID].leftGapIndex - temp->gapIndex].rightContigSingleReadCount++;
                        temp = temp->next;
                    }
                    
                    /*
                    if(contigToGapIndex[alignmentRight.RefID].leftGapIndex != -1 
                        && alignmentRight.Position < maxInsertSize
                        && alignmentRight.Position > minInsertSize - gapToContigIndex[contigToGapIndex[alignmentRight.RefID].leftGapIndex].gapDistance
                        ){
                        gapRegionReadSet->singleMapReadSet[contigToGapIndex[alignmentRight.RefID].leftGapIndex].rightContigSingleReadCount++;
                    } 
                    */     
                }
                
                
                continue;
            }
        
            
        }
        if((alignmentLeft.IsMapped() && !alignmentRight.IsMapped()) ){
            //cout<<"bb--"<<endl;
            long int leftPosition = GetRealMappingPositionFromEndContig(alignmentLeft.Position, endContigSet[alignmentLeft.RefID].contigLength, maxInsertSize, intervalLength, leftSoftLength);
            MapRegionGapIndex * temp = GetGapReadFromReadPositionOfContig(scaffoldSetHead, alignmentLeft.Length, alignmentLeft.RefID, leftPosition, maxInsertSize, minInsertSize, isPairedRead, alignmentLeft.IsReverseStrand());
            if((!alignmentLeft.IsReverseStrand() && isPairedRead) ||(alignmentLeft.IsReverseStrand() && !isPairedRead)){
                
                while(temp!=NULL){
                    gapRegionReadSet->singleMapReadSet[contigToGapIndex[alignmentLeft.RefID].rightGapIndex + temp->gapIndex].leftContigSingleReadCount++;
                    temp = temp->next;
                }
                
                /*
                if(contigToGapIndex[alignmentLeft.RefID].rightGapIndex != -1 
                    && alignmentLeft.Position > endContigSet[alignmentLeft.RefID].contigLength - maxInsertSize
                    && alignmentLeft.Position < endContigSet[alignmentLeft.RefID].contigLength - minInsertSize + gapToContigIndex[contigToGapIndex[alignmentLeft.RefID].rightGapIndex].gapDistance
                    ){
                    gapRegionReadSet->singleMapReadSet[contigToGapIndex[alignmentLeft.RefID].rightGapIndex].leftContigSingleReadCount++;
                }
                */
            }else{
                while(temp!=NULL){
                    gapRegionReadSet->singleMapReadSet[contigToGapIndex[alignmentLeft.RefID].leftGapIndex - temp->gapIndex].rightContigSingleReadCount++;
                    temp = temp->next;
                }
                /*
                if(contigToGapIndex[alignmentLeft.RefID].leftGapIndex != -1 
                    && alignmentLeft.Position < maxInsertSize
                    && alignmentLeft.Position > minInsertSize - gapToContigIndex[contigToGapIndex[alignmentLeft.RefID].leftGapIndex].gapDistance
                    ){
                    gapRegionReadSet->singleMapReadSet[contigToGapIndex[alignmentLeft.RefID].leftGapIndex].rightContigSingleReadCount++;
                } 
                */      
            }
        }
        
        if((!alignmentLeft.IsMapped() && alignmentRight.IsMapped())){
            //cout<<"cc--"<<endl;
            long int rightPosition = GetRealMappingPositionFromEndContig(alignmentRight.Position, endContigSet[alignmentRight.RefID].contigLength, maxInsertSize, intervalLength, rightSoftLength);
            MapRegionGapIndex * temp = GetGapReadFromReadPositionOfContig(scaffoldSetHead, alignmentRight.Length, alignmentRight.RefID, rightPosition, maxInsertSize, minInsertSize, isPairedRead, alignmentRight.IsReverseStrand());
            
            if((!alignmentRight.IsReverseStrand() && isPairedRead) ||(alignmentRight.IsReverseStrand() && !isPairedRead)){
                
                while(temp!=NULL){
                    gapRegionReadSet->singleMapReadSet[contigToGapIndex[alignmentRight.RefID].rightGapIndex  + temp->gapIndex].leftContigSingleReadCount++;
                    temp = temp->next;
                }
                
                /*
                if(contigToGapIndex[alignmentRight.RefID].rightGapIndex != -1 
                    && alignmentRight.Position > endContigSet[alignmentRight.RefID].contigLength - maxInsertSize
                    && alignmentRight.Position < endContigSet[alignmentRight.RefID].contigLength - minInsertSize + gapToContigIndex[contigToGapIndex[alignmentRight.RefID].rightGapIndex].gapDistance
                    ){
                    gapRegionReadSet->singleMapReadSet[contigToGapIndex[alignmentRight.RefID].rightGapIndex].leftContigSingleReadCount++;
                } 
                */
            }else{
                while(temp!=NULL){
                    gapRegionReadSet->singleMapReadSet[contigToGapIndex[alignmentRight.RefID].leftGapIndex - temp->gapIndex].rightContigSingleReadCount++;
                    temp = temp->next;
                }
                /*
                if(contigToGapIndex[alignmentRight.RefID].leftGapIndex != -1 
                    && alignmentRight.Position < maxInsertSize
                    && alignmentRight.Position > minInsertSize - gapToContigIndex[contigToGapIndex[alignmentRight.RefID].leftGapIndex].gapDistance
                    ){
                    gapRegionReadSet->singleMapReadSet[contigToGapIndex[alignmentRight.RefID].leftGapIndex].rightContigSingleReadCount++;
                } 
                */     
            }
        }
    }
    cout<<"ff"<<endl;
    for(long int i = 0; i < gapCount; i++){
        gapRegionReadSet->singleMapReadSet[i].leftContigSingleReadSet = (SingleRead *)malloc(sizeof(SingleRead)*gapRegionReadSet->singleMapReadSet[i].leftContigSingleReadCount);
        gapRegionReadSet->singleMapReadSet[i].rightContigSingleReadSet = (SingleRead *)malloc(sizeof(SingleRead)*gapRegionReadSet->singleMapReadSet[i].rightContigSingleReadCount);
        gapRegionReadSet->pairedMapReadSet[i].pairedReadSet = (PairedRead *)malloc(sizeof(PairedRead)*gapRegionReadSet->pairedMapReadSet[i].pairedReadCount);
        for(long int j = 0; j < gapRegionReadSet->singleMapReadSet[i].leftContigSingleReadCount; j++){
            gapRegionReadSet->singleMapReadSet[i].leftContigSingleReadSet[j].readIndex = -1;
            gapRegionReadSet->singleMapReadSet[i].leftContigSingleReadSet[j].read = NULL;
            gapRegionReadSet->singleMapReadSet[i].leftContigSingleReadSet[j].orientation = false;
            gapRegionReadSet->singleMapReadSet[i].leftContigSingleReadSet[j].contigPosition = -1;
            gapRegionReadSet->singleMapReadSet[i].leftContigSingleReadSet[j].contigIndex = -1;
        }
        for(long int j = 0; j < gapRegionReadSet->singleMapReadSet[i].rightContigSingleReadCount; j++){
            gapRegionReadSet->singleMapReadSet[i].rightContigSingleReadSet[j].readIndex = -1;
            gapRegionReadSet->singleMapReadSet[i].rightContigSingleReadSet[j].read = NULL;
            gapRegionReadSet->singleMapReadSet[i].rightContigSingleReadSet[j].orientation = false;
            gapRegionReadSet->singleMapReadSet[i].rightContigSingleReadSet[j].contigPosition = -1;
            gapRegionReadSet->singleMapReadSet[i].rightContigSingleReadSet[j].contigIndex = -1;
        }
        for(long int j = 0; j < gapRegionReadSet->pairedMapReadSet[i].pairedReadCount; j++){
            gapRegionReadSet->pairedMapReadSet[i].pairedReadSet[j].readIndex = -1;
            gapRegionReadSet->pairedMapReadSet[i].pairedReadSet[j].leftRead = NULL;
            gapRegionReadSet->pairedMapReadSet[i].pairedReadSet[j].leftOrientation = false;
            gapRegionReadSet->pairedMapReadSet[i].pairedReadSet[j].leftContigPosition = -1;
            gapRegionReadSet->pairedMapReadSet[i].pairedReadSet[j].leftContigIndex = -1;
            gapRegionReadSet->pairedMapReadSet[i].pairedReadSet[j].rightRead = NULL;
            gapRegionReadSet->pairedMapReadSet[i].pairedReadSet[j].rightOrientation = false;
            gapRegionReadSet->pairedMapReadSet[i].pairedReadSet[j].rightContigPosition = -1;
            gapRegionReadSet->pairedMapReadSet[i].pairedReadSet[j].rightContigIndex = -1;
        }

        gapRegionReadSet->singleMapReadSet[i].leftContigSingleReadCount = 0;
        gapRegionReadSet->singleMapReadSet[i].rightContigSingleReadCount = 0;
        gapRegionReadSet->pairedMapReadSet[i].pairedReadCount = 0;
    }
    cout<<"ee"<<endl;
    bamReaderLeft.Close();
    bamReaderRight.Close();
    
    bamReaderLeft.Open(bamFileNameLeft);
    bamReaderRight.Open(bamFileNameRight);
    
    long int readIndex = 0;
    long int gapIndex = 0;
    
    long int readCount = 0;
    long int pairedMapCount = 0;
    long int singleMapCoun = 0;
    long int noMapCount = 0;
    long int softCount = 0;
    long int softEndCount = 0;
    long int tt = 0;
    cout<<"tt"<<endl;
    while(bamReaderLeft.GetNextAlignment(alignmentLeft) && bamReaderRight.GetNextAlignment(alignmentRight)){
        //cout<<"--"<<readIndex<<endl;
        readCount++;
        while((alignmentLeft.AlignmentFlag & 0x900) != 0){
            bamReaderLeft.GetNextAlignment(alignmentLeft);
            continue;
        } 
        while((alignmentRight.AlignmentFlag & 0x900) != 0){
            bamReaderRight.GetNextAlignment(alignmentRight);
            continue;
        }
        bool leftSoftClip = false;
        bool rightSoftClip = false;
        long int leftSoftLength = 0;
        long int rightSoftLength = 0;
        //cout<<"aa"<<endl;
        for(int t = 0; t < 2; t++){
            std::vector< int > clipSizes;
            std::vector< int > readPositions;
            std::vector< int > genomePositions;
            BamAlignment tempAlignment;
            if(t == 0){
                tempAlignment = alignmentLeft;
            }else{
                tempAlignment = alignmentRight;
            }
            bool tempSoftClip = tempAlignment.GetSoftClips(clipSizes, readPositions, genomePositions);
            if(tempSoftClip != false){
                softCount++;
            }
            if(tempSoftClip != false){
                long int leftPosition = GetRealMappingPositionFromEndContig(tempAlignment.Position, endContigSet[tempAlignment.RefID].contigLength, maxInsertSize, intervalLength, 0);
                if((leftPosition ==1 && readPositions[0] == clipSizes[0]) ){
                    
                    gapIndex = contigToGapIndex[tempAlignment.RefID].leftGapIndex;
                    if(gapIndex >= 0){
                        long int index = gapRegionReadSet->singleMapReadSet[gapIndex].rightContigSingleReadCount++;
                        gapRegionReadSet->singleMapReadSet[gapIndex].rightContigSingleReadSet[index].readIndex = readIndex;
                        
                        gapRegionReadSet->singleMapReadSet[gapIndex].rightContigSingleReadSet[index].read = (char *)malloc(sizeof(char)*(tempAlignment.Length +1));
                        tempAlignment.QueryBases.copy(gapRegionReadSet->singleMapReadSet[gapIndex].rightContigSingleReadSet[index].read, tempAlignment.Length, 0);
                        gapRegionReadSet->singleMapReadSet[gapIndex].rightContigSingleReadSet[index].read[tempAlignment.Length] = '\0';
                                    
                        if(!isPairedRead){
                            char * tempRead = (char *)malloc(tempAlignment.Length +1);
                            ReverseComplement(gapRegionReadSet->singleMapReadSet[gapIndex].rightContigSingleReadSet[index].read, tempRead);
                            free(gapRegionReadSet->singleMapReadSet[gapIndex].rightContigSingleReadSet[index].read);
                            gapRegionReadSet->singleMapReadSet[gapIndex].rightContigSingleReadSet[index].read = tempRead;
                            tempRead = NULL;
                        }
                        
                        gapRegionReadSet->singleMapReadSet[gapIndex].rightContigSingleReadSet[index].orientation = tempAlignment.IsReverseStrand();
                        gapRegionReadSet->singleMapReadSet[gapIndex].rightContigSingleReadSet[index].contigIndex = tempAlignment.RefID;
                        gapRegionReadSet->singleMapReadSet[gapIndex].rightContigSingleReadSet[index].contigPosition = (maxInsertSize + minInsertSize)/2;
                        
                        leftSoftLength = clipSizes[0];
                        
                    }else{
                        tempSoftClip = false;
                    }
                    
                }else if(leftPosition > endContigSet[tempAlignment.RefID].contigLength - tempAlignment.Length 
                      && readPositions[0] + clipSizes[0] == tempAlignment.Length){
                    /*
                    char * tempRead = new char[tempAlignment.Length +1];
                    char * tempRead1 = new char[tempAlignment.Length +1];
                    strcpy(tempRead1, "AGCCCGCATATTTGAAGAAGCCGATGGCCGCGAGATTGGCGGCCACCCCCGCAAGCAGGAGCCACCGGCGCGGCGTCTCGACCAGCATGATCCCGAGCAGG");
                    alignmentLeft.QueryBases.copy(tempRead, alignmentLeft.Length, 0);
                    tempRead[tempAlignment.Length] = '\0';
                    tempRead1[tempAlignment.Length] = '\0';
                    if(strcmp(tempRead, tempRead1) == 0){
                        cout<<"tempRead:"<<endl;
                        cout<<tempAlignment.Length<<endl;
                        cout<<endContigSet[tempAlignment.RefID].contigLength<<endl;
                        cout<<leftPosition<<endl;
                        tempAlignment.QueryBases.copy(tempRead, alignmentLeft.Length, 0);
                        cout<<tempRead<<endl;
                        //exit(0);
                    }
                    */
                    
                    gapIndex = contigToGapIndex[tempAlignment.RefID].rightGapIndex;
                    
                    if(gapIndex >=0){
                        long int index = gapRegionReadSet->singleMapReadSet[gapIndex].leftContigSingleReadCount++;
                        gapRegionReadSet->singleMapReadSet[gapIndex].leftContigSingleReadSet[index].readIndex = readIndex;
                        
                        gapRegionReadSet->singleMapReadSet[gapIndex].leftContigSingleReadSet[index].read = (char *)malloc(sizeof(char)*(tempAlignment.Length +1));
                        tempAlignment.QueryBases.copy(gapRegionReadSet->singleMapReadSet[gapIndex].leftContigSingleReadSet[index].read, tempAlignment.Length, 0);
                        gapRegionReadSet->singleMapReadSet[gapIndex].leftContigSingleReadSet[index].read[tempAlignment.Length] = '\0';
                                    
                        if(isPairedRead){
                            char * tempRead = (char *)malloc(tempAlignment.Length +1);
                            ReverseComplement(gapRegionReadSet->singleMapReadSet[gapIndex].leftContigSingleReadSet[index].read, tempRead);
                            free(gapRegionReadSet->singleMapReadSet[gapIndex].leftContigSingleReadSet[index].read);
                            gapRegionReadSet->singleMapReadSet[gapIndex].leftContigSingleReadSet[index].read = tempRead;
                            tempRead = NULL;
                        }
                        
                        gapRegionReadSet->singleMapReadSet[gapIndex].leftContigSingleReadSet[index].orientation = tempAlignment.IsReverseStrand();
                        gapRegionReadSet->singleMapReadSet[gapIndex].leftContigSingleReadSet[index].contigIndex = tempAlignment.RefID;
                        gapRegionReadSet->singleMapReadSet[gapIndex].leftContigSingleReadSet[index].contigPosition = (maxInsertSize + minInsertSize)/2;
                        
                        rightSoftLength = clipSizes[0];
                        
                    }else{
                        //tempSoftClip = false;
                    }
                    
                    
                }else{
                    //tempSoftClip = false;
                }
                if(t == 0){
                    leftSoftClip = tempSoftClip;
                }else{
                    rightSoftClip = tempSoftClip;
                }
                if(tempSoftClip == true){
                    softEndCount++;
                }
            }
            vector<int>().swap(clipSizes);
            vector<int>().swap(readPositions);
            vector<int>().swap(genomePositions);
        }
        //cout<<"bb"<<endl;
        
        if((alignmentLeft.IsMapped() && alignmentRight.IsMapped())){                                    
            //continue;
            pairedMapCount++;
            if(abs(alignmentLeft.RefID - alignmentRight.RefID) != 1 || abs(alignmentLeft.RefID - alignmentRight.RefID) == 1){
                
                BamAlignment tempAlignment;
                BamAlignment tempNoAlignment;
                long int leftPosition = GetRealMappingPositionFromEndContig(alignmentLeft.Position, endContigSet[alignmentLeft.RefID].contigLength, maxInsertSize, intervalLength, leftSoftLength);
                long int rightPosition = GetRealMappingPositionFromEndContig(alignmentRight.Position, endContigSet[alignmentRight.RefID].contigLength, maxInsertSize, intervalLength, rightSoftLength);
                MapRegionGapIndex * tempLeft = GetGapReadFromReadPositionOfContig(scaffoldSetHead, alignmentLeft.Length, alignmentLeft.RefID, leftPosition, maxInsertSize, minInsertSize, isPairedRead, alignmentLeft.IsReverseStrand());
                MapRegionGapIndex * tempRight = GetGapReadFromReadPositionOfContig(scaffoldSetHead, alignmentRight.Length, alignmentRight.RefID, rightPosition, maxInsertSize, minInsertSize, isPairedRead, alignmentRight.IsReverseStrand());
                
                /*
                char * tempRead = new char[alignmentLeft.Length +1];
                char * tempRead1 = new char[alignmentLeft.Length +1];
                strcpy(tempRead1, "AGCCCGCATATTTGAAGAAGCCGATGGCCGCGAGATTGGCGGCCACCCCCGCAAGCAGGAGCCACCGGCGCGGCGTCTCGACCAGCATGATCCCGAGCAGG");
                alignmentLeft.QueryBases.copy(tempRead, alignmentLeft.Length, 0);
                tempRead[alignmentLeft.Length] = '\0';
                tempRead1[alignmentLeft.Length] = '\0';
                if(strcmp(tempRead, tempRead1) == 0){
                    cout<<"tempRead00:"<<endl;
                    cout<<alignmentLeft.Length<<endl;
                    cout<<endContigSet[alignmentLeft.RefID].contigLength<<endl;
                    cout<<leftPosition<<endl;
                    
                    if(leftSoftClip == true){
                        cout<<"leftTrue"<<endl;
                    }
                    if(rightSoftClip == true){
                        cout<<"rightTrue"<<endl;
                    }
                    
                    
                    if(tempLeft != NULL){
                        
                        cout<<tempLeft->gapIndex<<endl;
                    }
                    exit(0);
                }
                */
                
                tempAlignment = alignmentLeft;
                tempNoAlignment = alignmentRight;
                MapRegionGapIndex * temp = tempLeft;
                long int softLength = leftSoftLength;
                for(int p = 0; p < 2; p++){
                    if(p == 1){
                        tempAlignment = alignmentRight;
                        tempNoAlignment = alignmentLeft;
                        temp = tempRight;
                        softLength = rightSoftLength;
                    }
                    
                    if((!tempAlignment.IsReverseStrand() && isPairedRead) ||(tempAlignment.IsReverseStrand() && !isPairedRead)){
                        while(temp != NULL){
                            gapIndex = contigToGapIndex[tempAlignment.RefID].rightGapIndex;
                            long int index = gapRegionReadSet->singleMapReadSet[gapIndex + temp->gapIndex].leftContigSingleReadCount++;
                            gapRegionReadSet->singleMapReadSet[gapIndex + temp->gapIndex].leftContigSingleReadSet[index].readIndex = readIndex;
                            
                            gapRegionReadSet->singleMapReadSet[gapIndex + temp->gapIndex].leftContigSingleReadSet[index].read = (char *)malloc(sizeof(char)*(tempNoAlignment.Length +1));
                            tempNoAlignment.QueryBases.copy(gapRegionReadSet->singleMapReadSet[gapIndex + temp->gapIndex].leftContigSingleReadSet[index].read, tempNoAlignment.Length, 0);
                            gapRegionReadSet->singleMapReadSet[gapIndex + temp->gapIndex].leftContigSingleReadSet[index].read[tempNoAlignment.Length] = '\0';
                            
                            if(tempNoAlignment.IsReverseStrand()){
                                char * tempRead = (char *)malloc(tempNoAlignment.Length +1);
                                ReverseComplement(gapRegionReadSet->singleMapReadSet[gapIndex + temp->gapIndex].leftContigSingleReadSet[index].read, tempRead);
                                free(gapRegionReadSet->singleMapReadSet[gapIndex + temp->gapIndex].leftContigSingleReadSet[index].read);
                                gapRegionReadSet->singleMapReadSet[gapIndex + temp->gapIndex].leftContigSingleReadSet[index].read = tempRead;
                                tempRead = NULL;
                            }
                            
                            
                            /*
                            cout<<"readIndex:"<<readIndex<<endl;
                            char * tempRead = new char[tempNoAlignment.Length +1];
                            tempAlignment.QueryBases.copy(tempRead, tempNoAlignment.Length, 0);
                            tempRead[tempNoAlignment.Length] = '\0';
                            cout<<tempRead<<"--"<<gapRegionReadSet->singleMapReadSet[gapIndex + temp->gapIndex].leftContigSingleReadSet[index].read<<endl;
                            exit(0);
                            */
                            gapRegionReadSet->singleMapReadSet[gapIndex + temp->gapIndex].leftContigSingleReadSet[index].orientation = tempAlignment.IsReverseStrand();
                            gapRegionReadSet->singleMapReadSet[gapIndex + temp->gapIndex].leftContigSingleReadSet[index].contigIndex = tempAlignment.RefID;
                            gapRegionReadSet->singleMapReadSet[gapIndex + temp->gapIndex].leftContigSingleReadSet[index].contigPosition = temp->distance;
                            
                            temp = temp->next;
                        
                        }
                    }else{
                        while(temp != NULL){
                            gapIndex = contigToGapIndex[tempAlignment.RefID].leftGapIndex;
                            long int index = gapRegionReadSet->singleMapReadSet[gapIndex - temp->gapIndex].rightContigSingleReadCount++;
                            gapRegionReadSet->singleMapReadSet[gapIndex - temp->gapIndex].rightContigSingleReadSet[index].readIndex = readIndex;
                            
                            gapRegionReadSet->singleMapReadSet[gapIndex - temp->gapIndex].rightContigSingleReadSet[index].read = (char *)malloc(sizeof(char)*(tempNoAlignment.Length +1));
                            tempNoAlignment.QueryBases.copy(gapRegionReadSet->singleMapReadSet[gapIndex - temp->gapIndex].rightContigSingleReadSet[index].read, tempNoAlignment.Length, 0);
                            gapRegionReadSet->singleMapReadSet[gapIndex - temp->gapIndex].rightContigSingleReadSet[index].read[tempNoAlignment.Length] = '\0';
                            
                            if(tempNoAlignment.IsReverseStrand()){
                                char * tempRead = (char *)malloc(tempNoAlignment.Length +1);
                                ReverseComplement(gapRegionReadSet->singleMapReadSet[gapIndex - temp->gapIndex].rightContigSingleReadSet[index].read, tempRead);
                                free(gapRegionReadSet->singleMapReadSet[gapIndex - temp->gapIndex].rightContigSingleReadSet[index].read);
                                gapRegionReadSet->singleMapReadSet[gapIndex - temp->gapIndex].rightContigSingleReadSet[index].read = tempRead;
                                tempRead = NULL;
                            }
                            
                            gapRegionReadSet->singleMapReadSet[gapIndex - temp->gapIndex].rightContigSingleReadSet[index].orientation = tempAlignment.IsReverseStrand();
                            gapRegionReadSet->singleMapReadSet[gapIndex - temp->gapIndex].rightContigSingleReadSet[index].contigIndex = tempAlignment.RefID;
                            gapRegionReadSet->singleMapReadSet[gapIndex - temp->gapIndex].rightContigSingleReadSet[index].contigPosition = temp->distance - softLength;
                        
                            temp = temp->next;
                        }       
                    }
                    
                }
                readIndex++;
                continue;
            }    
            
        }
        
        if((alignmentLeft.IsMapped() && !alignmentRight.IsMapped() ) || (!alignmentLeft.IsMapped() && alignmentRight.IsMapped())){
            
            singleMapCoun++;
            MapRegionGapIndex * temp = NULL;    
            BamAlignment tempAlignment;
            BamAlignment tempNoAlignment;
            long int softLength = 0;
            if((alignmentLeft.IsMapped() && !alignmentRight.IsMapped()) || (leftSoftClip == true && rightSoftClip == false)){
                tempAlignment = alignmentLeft;
                tempNoAlignment = alignmentRight;
                long int leftPosition = GetRealMappingPositionFromEndContig(alignmentLeft.Position, endContigSet[alignmentLeft.RefID].contigLength, maxInsertSize, intervalLength, leftSoftLength);  
                temp = GetGapReadFromReadPositionOfContig(scaffoldSetHead, alignmentLeft.Length, alignmentLeft.RefID, leftPosition, maxInsertSize, minInsertSize, isPairedRead, alignmentLeft.IsReverseStrand());
                softLength = leftSoftLength;
                
            }else{
                tempAlignment = alignmentRight;
                tempNoAlignment = alignmentLeft;
                long int rightPosition = GetRealMappingPositionFromEndContig(alignmentRight.Position, endContigSet[alignmentRight.RefID].contigLength, maxInsertSize, intervalLength, rightSoftLength);
                temp = GetGapReadFromReadPositionOfContig(scaffoldSetHead, alignmentRight.Length, alignmentRight.RefID, rightPosition, maxInsertSize, minInsertSize, isPairedRead, alignmentRight.IsReverseStrand());
                softLength = rightSoftLength;
            }
            if((!tempAlignment.IsReverseStrand() && isPairedRead) ||(tempAlignment.IsReverseStrand() && !isPairedRead)){
                while(temp != NULL){
                    gapIndex = contigToGapIndex[tempAlignment.RefID].rightGapIndex;
                    long int index = gapRegionReadSet->singleMapReadSet[gapIndex + temp->gapIndex].leftContigSingleReadCount++;
                    //cout<<"33index--"<<index<<endl;
                    gapRegionReadSet->singleMapReadSet[gapIndex + temp->gapIndex].leftContigSingleReadSet[index].readIndex = readIndex;
                    
                    //gapRegionReadSet->singleMapReadSet[gapIndex].leftContigSingleReadSet[index].read = tempNoAlignment.QueryBases;
                    gapRegionReadSet->singleMapReadSet[gapIndex + temp->gapIndex].leftContigSingleReadSet[index].read = (char *)malloc(sizeof(char)*(tempNoAlignment.Length +1));
                    tempNoAlignment.QueryBases.copy(gapRegionReadSet->singleMapReadSet[gapIndex + temp->gapIndex].leftContigSingleReadSet[index].read, tempNoAlignment.Length, 0);
                    gapRegionReadSet->singleMapReadSet[gapIndex + temp->gapIndex].leftContigSingleReadSet[index].read[tempNoAlignment.Length] = '\0';
                    
                    gapRegionReadSet->singleMapReadSet[gapIndex + temp->gapIndex].leftContigSingleReadSet[index].orientation = tempAlignment.IsReverseStrand();
                    gapRegionReadSet->singleMapReadSet[gapIndex + temp->gapIndex].leftContigSingleReadSet[index].contigIndex = tempAlignment.RefID;
                    gapRegionReadSet->singleMapReadSet[gapIndex + temp->gapIndex].leftContigSingleReadSet[index].contigPosition = temp->distance;
                
                    temp = temp->next;
                }
            }else{
                while(temp != NULL){
                    gapIndex = contigToGapIndex[tempAlignment.RefID].leftGapIndex;
                    long int index = gapRegionReadSet->singleMapReadSet[gapIndex - temp->gapIndex].rightContigSingleReadCount++;
                    //cout<<"44index--"<<index<<endl;
                    gapRegionReadSet->singleMapReadSet[gapIndex - temp->gapIndex].rightContigSingleReadSet[index].readIndex = readIndex;
                    
                    //gapRegionReadSet->singleMapReadSet[gapIndex].rightContigSingleReadSet[index].read = tempNoAlignment.QueryBases;
                    gapRegionReadSet->singleMapReadSet[gapIndex - temp->gapIndex].rightContigSingleReadSet[index].read = (char *)malloc(sizeof(char)*(tempNoAlignment.Length +1));
                    tempNoAlignment.QueryBases.copy(gapRegionReadSet->singleMapReadSet[gapIndex - temp->gapIndex].rightContigSingleReadSet[index].read, tempNoAlignment.Length, 0);
                    gapRegionReadSet->singleMapReadSet[gapIndex - temp->gapIndex].rightContigSingleReadSet[index].read[tempNoAlignment.Length] = '\0';
                    /*
                    if(gapIndex == 1){
                        cout<<">gap"<<endl;
                        cout<<tempNoAlignment.QueryBases<<endl;
                    }
                    */
                    //cout<<tempNoAlignment.Length<<endl;
                    //cout<<gapRegionReadSet->singleMapReadSet[gapIndex].rightContigSingleReadSet[index].read<<endl;
                    /*
                    char * tempRead = new char[alignmentLeft.Length +1];
                    char * tempRead1 = new char[alignmentLeft.Length +1];
                    strcpy(tempRead1, "AGCCCGCATATTTGAAGAAGCCGATGGCCGCGAGATTGGCGGCCACCCCCGCAAGCAGGAGCCACCGGCGCGGCGTCTCGACCAGCATGATCCCGAGCAGG");
                    alignmentLeft.QueryBases.copy(tempRead, alignmentLeft.Length, 0);
                    tempRead[alignmentLeft.Length] = '\0';
                    tempRead1[alignmentLeft.Length] = '\0';
                    if(strcmp(tempRead, tempRead1) == 0){
                        cout<<"tempRead11:"<<endl;
                        cout<<alignmentLeft.Length<<endl;
                        cout<<endContigSet[alignmentLeft.RefID].contigLength<<endl;
                        //cout<<leftPosition<<endl;
                        
                        exit(0);
                    }
                    */
                    
                    
                    gapRegionReadSet->singleMapReadSet[gapIndex - temp->gapIndex].rightContigSingleReadSet[index].orientation = tempAlignment.IsReverseStrand();
                    gapRegionReadSet->singleMapReadSet[gapIndex - temp->gapIndex].rightContigSingleReadSet[index].contigIndex = tempAlignment.RefID;
                    gapRegionReadSet->singleMapReadSet[gapIndex - temp->gapIndex].rightContigSingleReadSet[index].contigPosition = temp->distance - softLength;
                    
                    temp = temp->next;
                }       
            }
        }
        
        if(!alignmentLeft.IsMapped() && !alignmentRight.IsMapped()){
            noMapCount++;
        }
        
        
        readIndex++;
    }
    
    
    cout<<"readCount:"<<readCount<<"--readIndex:"<<readIndex<<endl;
    cout<<"pairedMapCount:"<<pairedMapCount<<endl;
    cout<<"singleMapCoun:"<<singleMapCoun<<endl;
    cout<<"noMapCount:"<<noMapCount<<endl;
    cout<<"softCount:"<<softCount<<endl;
    cout<<"softEndCount:"<<softEndCount<<endl;
    
    //exit(0);
    /*
    for(long int j = 0; j < gapRegionReadSet->pairedMapReadSet[1].pairedReadCount; j++){
        cout<<gapRegionReadSet->pairedMapReadSet[1].pairedReadSet[j].leftContigIndex<<"--"<<gapRegionReadSet->pairedMapReadSet[1].pairedReadSet[j].rightContigIndex<<endl;
    }
    */
    return gapRegionReadSet;
    
}


int KMPIndexOfContigOfMisMatch(char * contig, char * pattern){
    long int i = 0;
    long int j = 0;
    long int index = 0;
    long int token = 0;
    long int p = 0;
    long int len1 = strlen(contig);
    long int len2 = strlen(pattern);
    if(len1<len2){
        cout<<"IndexOfContig Wrong!"<<endl;
        cout<<pattern<<"--"<<contig<<endl;
        //return -1;
        exit(0);
    }
    for(i=0;i<=len1-len2;i++){
        p = i;
        token = 0;
        index = 0;
        for(j=0;j<len2;j++){
            if(contig[p]!=pattern[j]){
                if(token == 0){
                    token = 1;
                }else{
                    index = -1;
                    break;
                }
            }
            p++;
        }
        if(index==0){
            return i;
        }
    }
    return -1;
}

void SearchSharedKmerRead(ScaffoldSetHead * scaffoldSetHead, GapRegionReadSet * gapRegionReadSet){
    
    long int gapIndex = 1;
    long int gapCount = gapRegionReadSet->gapCount;
    GapRegionSingleMapReadSet * singleMapReadSet = gapRegionReadSet->singleMapReadSet;

    long int kmerLength = 15;
    char * kmer = (char *)malloc(sizeof(char)*kmerLength+1);
    
    long int readIndex = 0;
    
    
    for(long int i = 0; i < singleMapReadSet[gapIndex].leftContigSingleReadCount; i++){
        char * read = singleMapReadSet[gapIndex].leftContigSingleReadSet[i].read;
        long int t = -1;
        for(long int p = 0; p < strlen(read)-kmerLength; p++){
            
            strncpy(kmer, read+p, kmerLength);
            kmer[kmerLength] = '\0';
            for(long int j = 0; j < singleMapReadSet[gapIndex].rightContigSingleReadCount; j++){
                t = KMPIndexOfContigOfMisMatch(singleMapReadSet[gapIndex].rightContigSingleReadSet[j].read, kmer);
                if(t!=-1){
                    break;
                }
            }
            if(t!=-1){
                break;
            }
        }
        if(t!=-1){
            cout<<">read_left_"<<readIndex<<endl;
            cout<<read<<endl;
            readIndex++;
        }
        
    }
    
    for(long int i = 0; i < singleMapReadSet[gapIndex].rightContigSingleReadCount; i++){
        char * read = singleMapReadSet[gapIndex].rightContigSingleReadSet[i].read;
        long int t = -1;
        for(long int p = 0; p < strlen(read)-kmerLength; p++){
            
            strncpy(kmer, read+p, kmerLength);
            kmer[kmerLength] = '\0';
            for(long int j = 0; j < singleMapReadSet[gapIndex].leftContigSingleReadCount; j++){
                t = KMPIndexOfContigOfMisMatch(singleMapReadSet[gapIndex].leftContigSingleReadSet[j].read, kmer);
                if(t!=-1){
                    break;
                }
            }
            if(t!=-1){
                break;
            }
        }
        if(t!=-1){
            cout<<">read_right_"<<readIndex<<endl;
            cout<<read<<endl;
            readIndex++;
        }
        
    }
    
}


void FillSingleGapByRead(ScaffoldSetHead * scaffoldSetHead, GapRegionReadSet * gapRegionReadSet, long int libraryIndex, char * outputName){
    
    long int gapIndex = 0;
    long int gapCount = gapRegionReadSet->gapCount;
    GapRegionSingleMapReadSet * singleMapReadSet = gapRegionReadSet->singleMapReadSet;
    GapRegionPairedMapReadSet * pairedMapReadSet = gapRegionReadSet->pairedMapReadSet;
    GapToContigIndex * gapToContigIndex = scaffoldSetHead->gapToContigIndex;
    
    char * leftReadFileName = (char *)malloc(sizeof(char)*200);
    char * rightReadFileName = (char *)malloc(sizeof(char)*200);
    FILE * fp;
    FILE * fp1;
    
    for(long int i = 0; i < gapCount; i++){
        //cout<<"gapIndex:"<<i<<endl;
        //cout<<"leftPosition------------"<<endl;
    
        sprintf(leftReadFileName, "./%s/gap_read_out_put_%ld/leftReadSet_gapIndex_%ld.fa", outputName, libraryIndex, i);
        sprintf(rightReadFileName, "./%s/gap_read_out_put_%ld/rightReadSet_gapIndex_%ld.fa", outputName, libraryIndex, i);
        fp = fopen(leftReadFileName, "w+");
        fp1 = fopen(rightReadFileName, "w+");
        
        for(long int j = 0; j < singleMapReadSet[i].leftContigSingleReadCount; j++){
            //cout<<singleMapReadSet[i].leftContigSingleReadSet[j].contigPosition<<endl;
            //cout<<singleMapReadSet[i].leftContigSingleReadSet[j].read<<endl;
            //fprintf(fp, ">%ld\n", j);
            fprintf(fp, "%s\t%ld\n", singleMapReadSet[i].leftContigSingleReadSet[j].read, singleMapReadSet[i].leftContigSingleReadSet[j].contigPosition);
        }
        /*
        if(singleMapReadSet[i].leftContigSingleReadCount%2 != 0){
            long int readLength = strlen(singleMapReadSet[i].leftContigSingleReadSet[0].read);
            fprintf(fp, ">false\n");
            for(long int p = 0; p < readLength; p++){
                fprintf(fp, "A");
            }
            fprintf(fp, "\t0\n");
        }
        */
        //cout<<"rightPosition------------"<<endl;
        for(long int j = 0; j < singleMapReadSet[i].rightContigSingleReadCount; j++){
            //cout<<singleMapReadSet[i].rightContigSingleReadSet[j].contigPosition<<endl;
            //cout<<singleMapReadSet[i].rightContigSingleReadSet[j].read<<endl;
            //fprintf(fp1, ">%ld\n", j);
            fprintf(fp1, "%s\t%ld\n", singleMapReadSet[i].rightContigSingleReadSet[j].read, singleMapReadSet[i].rightContigSingleReadSet[j].contigPosition);
        }
        /*
        if(singleMapReadSet[i].rightContigSingleReadCount%2 != 0){
            long int readLength = strlen(singleMapReadSet[i].rightContigSingleReadSet[0].read);
            //fprintf(fp1, ">false\n");
            for(long int p = 0; p < readLength; p++){
                fprintf(fp1, "A");
            }
            fprintf(fp1, "\t0\n");
        }
        */
        fclose(fp);
        fclose(fp1);
        
    }
    
    //cout<<"fillRegionCount:"<<fillRegionCount<<endl;
    
    /*
    for(long int i = 1; i < 2; i++){
        for(long int j = 0; j < pairedMapReadSet[i].pairedReadCount; j++){
            fprintf(fp, ">gapIndex:%ld--left_%ld\n", i, j);
            fprintf(fp, "%s\n", pairedMapReadSet[i].pairedReadSet[j].leftRead);
            fprintf(fp, ">gapIndex:%ld--right_%ld\n", i, j);
            fprintf(fp, "%s\n", pairedMapReadSet[i].pairedReadSet[j].leftRead);
        } 
    }
    */
    
}

void WriteGapRegionReadSet(GapRegionReadSet * gapRegionReadSet, char * fileName){
    
    FILE * fp;
    if((fp = fopen(fileName, "w+")) == NULL){
        printf("%s, does not exist!", fileName);
        exit(0);
    }
    
    long int gapCount = gapRegionReadSet->gapCount;
    GapRegionSingleMapReadSet * singleMapReadSet = gapRegionReadSet->singleMapReadSet;
    GapRegionPairedMapReadSet * pairedMapReadSet = gapRegionReadSet->pairedMapReadSet;
    
    fprintf(fp, "gapCount:%ld\n", gapCount);
    fprintf(fp, "\n-----------------------------------\n");
    fprintf(fp, "pairedMapReadSet\n");
    fprintf(fp, "\n-----------------------------------\n");
    
    //cout<<"--------------"<<endl;
    
    
    for(long int i = 0; i < gapCount; i++){
        long int readCount = pairedMapReadSet[i].pairedReadCount;
        fprintf(fp, "gapIndex:%ld--readCount:%ld\n", i, readCount);
        for(long int j = 0; j < readCount; j++){
            fprintf(fp, "readIndex:%ld--%ld\n", pairedMapReadSet[i].pairedReadSet[j].readIndex, j);
            fprintf(fp, "%ld,%ld,%d\n", pairedMapReadSet[i].pairedReadSet[j].leftContigIndex, pairedMapReadSet[i].pairedReadSet[j].leftContigPosition, pairedMapReadSet[i].pairedReadSet[j].leftOrientation);
            fprintf(fp, "%ld,%ld,%d\n", pairedMapReadSet[i].pairedReadSet[j].rightContigIndex, pairedMapReadSet[i].pairedReadSet[j].rightContigPosition, pairedMapReadSet[i].pairedReadSet[j].rightOrientation);
        }
    }
    fprintf(fp, "\n-----------------------------------\n");
    fprintf(fp, "singleMapReadSet\n");
    fprintf(fp, "\n-----------------------------------\n");
    for(long int i = 0; i < gapCount; i++){
        long int readCount = singleMapReadSet[i].leftContigSingleReadCount;
        for(long int j = 0; j < readCount; j++){
            fprintf(fp, "leftReadIndex:%ld\n", singleMapReadSet[i].leftContigSingleReadSet[j].readIndex);
            fprintf(fp, "%ld,%ld,%d\n", singleMapReadSet[i].leftContigSingleReadSet[j].contigIndex, singleMapReadSet[i].leftContigSingleReadSet[j].contigPosition, singleMapReadSet[i].leftContigSingleReadSet[j].orientation);
            
        }
        readCount = singleMapReadSet[i].rightContigSingleReadCount;
        for(long int j = 0; j < readCount; j++){
            fprintf(fp, "rightReadIndex:%ld\n", singleMapReadSet[i].rightContigSingleReadSet[j].readIndex);
            fprintf(fp, "%ld,%ld,%d\n", singleMapReadSet[i].rightContigSingleReadSet[j].contigIndex, singleMapReadSet[i].rightContigSingleReadSet[j].contigPosition, singleMapReadSet[i].rightContigSingleReadSet[j].orientation);
            
        }
    }
    
    
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
        }else if(temp1[i]=='N'){
            temp2[len-1-i]='N';
        }else{
             return false;
        }
    }
    temp2[len]='\0';
    return true;
}

