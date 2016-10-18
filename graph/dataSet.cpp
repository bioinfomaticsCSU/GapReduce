#include "fstream"
#include <cstring>
#include <ctype.h>
#include <iostream>
#include "stdio.h"
#include "stdlib.h"
#include <pthread.h>

#include "dataSet.h"

using namespace std;

long int file_size(char* filename){  
    FILE *fp=fopen(filename,"r");  
    if(!fp) return -1;  
    fseek(fp,0L,SEEK_END);  
    long int size=ftell(fp);  
    fclose(fp);  
      
    return size;  
} 

double RatioOfN(char * contig){
    if(contig == NULL){
        return 1;
    }
    long int count = 0;
    long int len = strlen(contig);
    for(long int i = 0; i < len; i++){
        if(contig[i] == 'N'){
           count++;
        }
    }
    return double(count)/double(len);
}


void CopyReadofATGC(char * destination, char * source, long int length){
    for(long int i = 0; i < length; i ++){
        *(destination+i) = toupper(*(source+i));
        if(*(destination+i) == 'N'){
            *(destination+i) = 'A';
        } 
    }
    
}

int KMPIndexOfContig(char * contig, char * pattern){
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
                index = -1;
                break;
            }
            p++;
        }
        if(index==0){
            return i;
        }
    }
    return -1;
}

void AppendRight( char * temp, char * temp1, char * temp2,long int kmerLength){
     long int i = 0;
     long int j = 0;
     long int len1 = strlen(temp1);
     long int len2 = strlen(temp2);
     for(i =0;i<len1+len2-kmerLength+1;i++){
         if(i<len1){
             temp[i] = temp1[i];
         }else{
             temp[i] = temp2[kmerLength-1+j];
             j++;
         }
     }
     temp[i] = '\0';
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

bool ReverseComplement(char * temp){
    long int len = strlen(temp);
    char a;
    char b;
    long int m = len/2;
    for(long int i = 0; i < len/2; i++){
        a = temp[i];
        b = temp[len - 1 - i];
        if(a == 'A'){
            temp[len - 1 - i] = 'T';
        }else if(a == 'T'){
            temp[len - 1 - i] = 'A';
        }else if(a == 'G'){
            temp[len - 1 - i] = 'C';
        }else if(a == 'C'){
            temp[len - 1 - i] = 'G';
        }else if(a == 'N'){
            temp[len - 1 - i] = 'N';
        }
        
        if(b == 'A'){
            temp[i] = 'T';
        }else if(b == 'T'){
            temp[i] = 'A';
        }else if(b == 'G'){
            temp[i] = 'C';
        }else if(b == 'C'){
            temp[i] = 'G';
        }else if(b == 'N'){
            temp[i] = 'N';
        }
        
    }
    if(len%2 != 0){
        if(temp[m] == 'A'){
            temp[m] = 'T';
        }else if(temp[m] == 'T'){
            temp[m] = 'A';
        }else if(temp[m] == 'G'){
            temp[m] = 'C';
        }else if(temp[m] == 'C'){
            temp[m] = 'G';
        }else if(temp[m] == 'N'){
            temp[m] = 'N';
        }
    }
    temp[len]='\0';
    return true;
}

bool ReverseComplementReadSet(ReadSet * readSet, long int readCount){
    for(long int i = 0; i < readCount; i++){
        ReverseComplement(readSet[i].read);
    }
}

ReadSet * GetReadSet(char * readAddress, long int & readCount){
    FILE * fp;
    if((fp = fopen(readAddress, "r")) == NULL){
        printf("%s, does not exist!", readAddress);
        exit(0);
    }
    
    long int readLength = 0;
    
    long int maxLength = 1000;
    char * read = (char *)malloc(sizeof(char)*maxLength);
    
    while((fgets(read, maxLength, fp)) != NULL){ 
       readLength = strlen(read);
       if(readLength == maxLength -1){  
           while((fgets(read, maxLength, fp)) != NULL){
               readLength = readLength + strlen(read);
               if(strlen(read) != maxLength -1){
                   break;
               }
           }
       }
       readCount++;  
    }
    fclose(fp);
    
    ReadSet * readSet = (ReadSet *)malloc(sizeof(ReadSet)*readCount);
    if((fp = fopen(readAddress, "r")) == NULL){
        printf("%s, does not exist!", readAddress);
        exit(0);
    }
    char * p = NULL;
    long int readIndex = 0;
    long int tempLength = 0;
    char * tempRead = (char *)malloc(sizeof(char)*(readLength+1));
    while((fgets(read, maxLength, fp)) != NULL){ 
       tempLength = strlen(read);
       CopyReadofATGC(tempRead, read, tempLength);
       if(strlen(read) == maxLength -1){  
           while((fgets(read, maxLength, fp)) != NULL){
               CopyReadofATGC(tempRead+tempLength, read, strlen(read));
               if(strlen(read) != maxLength -1){
                   break;
               }
               tempLength = tempLength + strlen(read);
           }
       }
       p = strtok(tempRead, "\t");
       readSet[readIndex].read = (char *)malloc(sizeof(char)*(strlen(p) + 1));
       strncpy(readSet[readIndex].read, p, strlen(p));
       readSet[readIndex].read[strlen(p)] = '\0';
       p = strtok(NULL, "\t");
       sscanf(p, "%ld", &readSet[readIndex].insertSize);
       readIndex++;
    }
    fclose(fp); 
    free(read);
    free(tempRead);
    return readSet;
}


ReadSetHead * GetReadSetHead(char * leftReadAddress, char * rightReadAddress, bool isPaired){
    
    ReadSetHead * readSetHead = (ReadSetHead*)malloc(sizeof(ReadSetHead));
    readSetHead->readLength = 0;
    readSetHead->leftReadSet = NULL;
    readSetHead->rightReadSet = NULL;
    readSetHead->leftReadCount = 0;
    readSetHead->rightReadCount = 0;
    readSetHead->isPaired = isPaired;
    
    readSetHead->leftReadSet = GetReadSet(leftReadAddress, readSetHead->leftReadCount);
    readSetHead->rightReadSet = GetReadSet(rightReadAddress, readSetHead->rightReadCount);
    
    
    if(readSetHead->isPaired == true){
        ReverseComplementReadSet(readSetHead->leftReadSet, readSetHead->leftReadCount);
    }else{
        ReverseComplementReadSet(readSetHead->rightReadSet, readSetHead->rightReadCount);
    }
    
    
    if(readSetHead->leftReadCount != 0){
        readSetHead->readLength = strlen(readSetHead->leftReadSet[0].read);
    }else if(readSetHead->rightReadCount != 0){
        readSetHead->readLength = strlen(readSetHead->rightReadSet[0].read);
    }
    /*/
    for(long int i = 0; i < readSetHead->rightReadCount; i++){
        cout<<readSetHead->rightReadSet[i].read<<'\t'<<readSetHead->rightReadSet[i].insertSize<<endl;
    }
    */
    return readSetHead;
    
}

unsigned long int Hash(char * str, unsigned int len, unsigned long int max)  
{  
   unsigned long int hash = 0;  
   unsigned long int i = 0;  
  
   for(i = 0; i < len; str++, i++) {  
      hash = (*str) + (hash << 6) + (hash << 16) - hash;  
   }  
  
   return hash % max;  
}

long int SearchKmerOfKmerSet(char * kmer, KmerSet * kmerSet, long int kmerLength, long int kmerCount){
    if(kmerSet == NULL){
        return -1;
    }
    unsigned long int hashKey = Hash(kmer, kmerLength, kmerCount);
    while(kmerSet[hashKey].kmer != NULL){
        if(strcmp(kmer, kmerSet[hashKey].kmer) == 0){
            return hashKey;
        }
        hashKey = (hashKey+1)%kmerCount;
    }
    return -1;
}

void InsertKmerToKmerSet(char * kmer, long int readIndex, KmerSet * kmerSet, long int kmerLength, long int kmerCount){
    unsigned long int hashKey = Hash(kmer, kmerLength, kmerCount);
    while(kmerSet[hashKey].kmer != NULL){
        if(strcmp(kmer, kmerSet[hashKey].kmer) == 0){
            break;
        }
        hashKey = (hashKey+1)%kmerCount;
    }
    kmerSet[hashKey].kmerCount++;
    KmerReadIndex * kmerReadIndex = (KmerReadIndex *)malloc(sizeof(KmerReadIndex));
    kmerReadIndex->index = readIndex;
    kmerReadIndex->next = NULL;
    if(kmerSet[hashKey].kmer != NULL){
        kmerReadIndex->next = kmerSet[hashKey].readIndex->next;
        kmerSet[hashKey].readIndex->next = kmerReadIndex;
    }else{
        kmerSet[hashKey].kmer = (char *)malloc(sizeof(char)*(kmerLength+1));
        strncpy(kmerSet[hashKey].kmer, kmer, kmerLength);
        kmerSet[hashKey].kmer[kmerLength] = '\0';
        kmerSet[hashKey].readIndex = kmerReadIndex;
    }
}

KmerSet * GetKmerSet(ReadSet * readSet, long int readCount, long int readLength, long int kmerLength, long int & kmerCount){
    
    if(readCount <= 0){
        return NULL;
    }
    
    kmerCount = 1.2*(readLength - kmerLength + 1)*readCount;
    KmerSet * kmerSet = (KmerSet *)malloc(sizeof(KmerSet)*kmerCount);
    for(long int i = 0; i < kmerCount; i++){
        kmerSet[i].kmer = NULL;
        kmerSet[i].kmerCount = 0;
        kmerSet[i].readIndex = NULL;
    }
    char * tempKmer = (char *)malloc(sizeof(char)*(kmerLength+1));
    for(long int i = 0; i < readCount; i++){
        for(long int j = 0; j < readLength - kmerLength + 1; j++){
            strncpy(tempKmer, readSet[i].read + j, kmerLength);
            tempKmer[kmerLength] = '\0';
            InsertKmerToKmerSet(tempKmer, i, kmerSet, kmerLength, kmerCount);
        }
    }
    free(tempKmer);
    return kmerSet;
}

KmerSetHead * GetKmerSetHead(ReadSetHead * readSetHead, long int kmerLength, long int minKmerFrequency){
    
    KmerSetHead * kmerSetHead = (KmerSetHead *)malloc(sizeof(KmerSetHead));
    kmerSetHead->kmerLength = kmerLength;
    kmerSetHead->minKmerFrequency = minKmerFrequency;
    kmerSetHead->leftKmerSet = NULL;
    kmerSetHead->rightKmerSet = NULL;
    kmerSetHead->leftKmerCount = 0;
    kmerSetHead->rightKmerCount = 0;
    
    kmerSetHead->leftKmerSet = GetKmerSet(readSetHead->leftReadSet, readSetHead->leftReadCount, readSetHead->readLength, kmerLength, kmerSetHead->leftKmerCount);
    
    kmerSetHead->rightKmerSet = GetKmerSet(readSetHead->rightReadSet, readSetHead->rightReadCount, readSetHead->readLength, kmerLength, kmerSetHead->rightKmerCount);
    
    /*
    for(long int i = 0; i < kmerSetHead->rightKmerCount; i++){
        if(kmerSetHead->rightKmerSet[i].kmer != NULL){
            
            cout<<kmerSetHead->rightKmerSet[i].kmer<<'\t'<<kmerSetHead->rightKmerSet[i].kmerCount<<endl;
            
            
            
            KmerReadIndex * tempKmerReadIndex = kmerSetHead->leftKmerSet[i].readIndex;
            cout<<"index:"<<'\t';
            while(tempKmerReadIndex != NULL){
                cout<<tempKmerReadIndex->index<<'\t';
                tempKmerReadIndex = tempKmerReadIndex->next;
            }
            cout<<endl;
            
        }
        
    }
    */
    
    /*
    char * tempKmer = new char[22];
    strcpy(tempKmer, "AAAAGTTTCATTTGATTATAT");
    tempKmer[21] = '\0';
    long int key = SearchKmerOfKmerSet(tempKmer, kmerSetHead->leftKmerSet, kmerLength, kmerSetHead->leftKmerCount);
    cout<<kmerSetHead->leftKmerSet[key].kmer<<"--"<<kmerSetHead->leftKmerSet[key].kmerCount<<endl;
    */
    
    return kmerSetHead;
}



