#include "fstream"
#include <cstring>
#include <iostream>
#include "stdio.h"
#include "stdlib.h"
#include <pthread.h>

#include "graph.h"


using namespace std;


int NodeCount(GraphNode * temp){
    int i = 0;
    while(temp!=NULL){
        i++;
        temp = temp->next;
    }
    return i;
}


int DBGSearchNode(GraphNode * tempGraphNode, long int num, long int newNum){
    while(tempGraphNode!=NULL){
        if(tempGraphNode->index==num){
            tempGraphNode->index = newNum;
            return 1;
        }
        tempGraphNode = tempGraphNode->next;
    }
    return 0;
}

int DBGSearchNode(GraphNode * tempGraphNode, long int num){
    while(tempGraphNode!=NULL){
        if(tempGraphNode->index==num){
            return 1;
        }
        tempGraphNode = tempGraphNode->next;
    }
    return 0;
}



int DBGRemoveInNode(DBGraph * graph, long int num){
    GraphNode * tempGraphNode = graph->inNode;
    GraphNode * temp = NULL;
    while(tempGraphNode!=NULL){
        if(tempGraphNode->index==num){
            if(temp==NULL){
                graph->inNode = tempGraphNode->next;
                tempGraphNode->next = NULL;
                delete tempGraphNode;
                return 1;
            }
            temp->next = tempGraphNode->next;
            tempGraphNode->next = NULL;
            delete tempGraphNode; 
            return 1;
        }
        temp = tempGraphNode;
        tempGraphNode = tempGraphNode->next;
    }
    return 0;
}

int DBGRemoveOutNode(DBGraph * graph, long int num){
    GraphNode * tempGraphNode = graph->outNode;
    GraphNode * temp = NULL;
    while(tempGraphNode!=NULL){
        if(tempGraphNode->index==num){
            if(temp==NULL){
                graph->outNode = tempGraphNode->next;
                tempGraphNode->next = NULL;
                delete tempGraphNode;
                return 1;
            }
            temp->next = tempGraphNode->next;
            tempGraphNode->next = NULL;
            delete tempGraphNode; 
            return 1;
        }
        temp = tempGraphNode;
        tempGraphNode = tempGraphNode->next;
    }
    return 0;
}

void AdjustDBG(DBGraph *deBruijnGraph, int kmerLength,long int max){
    long int i = 0;
    long int j = 0;
    long int m = 0;
    long int n = 0;
    GraphNode * tempGraphNode = NULL;
    
    while(i<max){
        if(deBruijnGraph[i].contig==NULL){
            i++;
            j++;
            continue;
        }
        if(j>0){
            deBruijnGraph[i-j].contig = deBruijnGraph[i].contig;
            deBruijnGraph[i-j].inNode = deBruijnGraph[i].inNode;
            deBruijnGraph[i-j].outNode = deBruijnGraph[i].outNode;
            deBruijnGraph[i].contig = NULL;
            deBruijnGraph[i].inNode = NULL;
            deBruijnGraph[i].outNode = NULL;
            tempGraphNode = deBruijnGraph[i-j].inNode;
            while(tempGraphNode!=NULL){
                if(DBGSearchNode(deBruijnGraph[tempGraphNode->index].outNode,i,i-j)!=1){
                    //cout<<i<<"--"<<j<<"--"<<deBruijnGraph[i].contig<<"--"<<deBruijnGraph[tempGraphNode->index].contig<<endl;
                    cout<<"1SimplePathMerge Error!"<<endl;
                    exit(0);
                }
                tempGraphNode = tempGraphNode->next;
            }
            tempGraphNode = deBruijnGraph[i-j].outNode;
            while(tempGraphNode!=NULL){
                if(DBGSearchNode(deBruijnGraph[tempGraphNode->index].inNode,i,i-j)!=1){
                    //cout<<i<<"--"<<j<<"--"<<deBruijnGraph[i].contig<<"--"<<deBruijnGraph[tempGraphNode->index].contig<<endl;
                    cout<<"2SimplePathMerge Error!"<<endl;
                    exit(0);
                }
                tempGraphNode = tempGraphNode->next;
            }
        }
        i++;
    }
    //cout<<"AdjustDBG End!"<<endl;
}


void SimplePathMerge(DBGraph *deBruijnGraph, int kmerLength, long int max){
    long int i = 0;
    long int j = 0;
    long int m = 0;
    long int n = 0;
    long int t = 0;
    long int previousIndex = -1;
    GraphNode *tempGraphNode = NULL;
    while(i<max){
        
        if(deBruijnGraph[i].contig==NULL){
            i++;
            continue;
        }
        while(NodeCount(deBruijnGraph[i].inNode)==1){
            previousIndex = deBruijnGraph[i].inNode->index;
            if(NodeCount(deBruijnGraph[previousIndex].outNode)!=1||DBGSearchNode(deBruijnGraph[previousIndex].inNode,i)==1){
                break;
            }
            t = 0;
            n = strlen(deBruijnGraph[previousIndex].contig);
            m = strlen(deBruijnGraph[i].contig);
            //cout<<"m:"<<m<<"--n:"<<n<<"--current:"<<i<<"--pre--"<<previousIndex<<endl;
            //cout<<deBruijnGraph[previousIndex].contig<<endl;
            //cout<<deBruijnGraph[i].contig<<endl;
            char *tempContig = new char[m+n-kmerLength + 2];
            AppendRight(tempContig,deBruijnGraph[previousIndex].contig, deBruijnGraph[i].contig, kmerLength);
            
            delete deBruijnGraph[i].contig;
            deBruijnGraph[i].contig = tempContig;
            delete deBruijnGraph[i].inNode;
            deBruijnGraph[i].inNode = deBruijnGraph[previousIndex].inNode;
            deBruijnGraph[previousIndex].inNode = NULL;
            delete deBruijnGraph[previousIndex].outNode;
            deBruijnGraph[previousIndex].outNode = NULL;
            delete deBruijnGraph[previousIndex].contig;
            deBruijnGraph[previousIndex].contig = NULL;
            //cout<<deBruijnGraph[i].contig<<endl;
            
            tempGraphNode = deBruijnGraph[i].inNode;
            while(tempGraphNode!=NULL){
                if(DBGSearchNode(deBruijnGraph[tempGraphNode->index].outNode,previousIndex,i)!=1){
                    cout<<"0SimplePathMerge Error!"<<endl;
                    cout<<tempGraphNode->index<<"--"<<previousIndex<<"--"<<i<<endl;
                    exit(0);
                }
                tempGraphNode = tempGraphNode->next;
            }
        }
        i++;
    }
         
    //cout<<"SimplePathMerge End!"<<endl;
    
}

void RemoveTip(DBGraph *deBruijnGraph, int kmerLength, long int graphNodeNumber){
    long int i = 0;
    long int j = 0;
    long int inNodeCount = 0;
    long int outNodeCount = 0;
    GraphNode * tempGraphNode = new GraphNode;
    GraphNode * first = tempGraphNode;
    while(i<graphNodeNumber){
        if(deBruijnGraph[i].contig==NULL || strlen(deBruijnGraph[i].contig)>=2*kmerLength){
            i++;
            continue;
        }
        inNodeCount = NodeCount(deBruijnGraph[i].inNode);
        outNodeCount = NodeCount(deBruijnGraph[i].outNode);
        if((inNodeCount==1&&outNodeCount==0)||(inNodeCount==0&&outNodeCount==1)||(inNodeCount==0&&outNodeCount==0)){
            GraphNode * newGraphNode = new GraphNode;
            newGraphNode->index = i;
            tempGraphNode->next = newGraphNode;
            tempGraphNode = newGraphNode;
        }
        i++;
    }
    
    first = first->next;
    while(first!=NULL){
        inNodeCount = NodeCount(deBruijnGraph[first->index].inNode);
        outNodeCount = NodeCount(deBruijnGraph[first->index].outNode);
        if(inNodeCount==1&&outNodeCount==0){
            j = deBruijnGraph[first->index].inNode->index;
            DBGRemoveOutNode(deBruijnGraph+j,first->index);
            delete deBruijnGraph[first->index].inNode;
            deBruijnGraph[first->index].inNode = NULL;
            delete deBruijnGraph[first->index].contig;
            deBruijnGraph[first->index].contig = NULL; 
        }
        if(inNodeCount==0&&outNodeCount==1){
            j = deBruijnGraph[first->index].outNode->index;
            DBGRemoveInNode(deBruijnGraph+j,first->index);
            delete deBruijnGraph[first->index].outNode;
            deBruijnGraph[first->index].outNode = NULL;
            delete deBruijnGraph[first->index].contig;
            deBruijnGraph[first->index].contig = NULL; 
        }
        if(inNodeCount==0&&outNodeCount==0){
            delete deBruijnGraph[first->index].contig;
            deBruijnGraph[first->index].contig = NULL; 
        }
        first = first->next;
    }
    //cout<<"RemoveTip End!"<<endl;
}



void MergeBubble(DBGraph *deBruijnGraph, int kmerLength, long int graphNodeNumber){
    
    long int i = 0;
    long int j = 0;
    
    long int outNodeNumber = 0;
    long int inNodeNumber = 0;
    //cout<<"aaaaaaaaaaaaaaaaaa"<<endl;
    for(i = 0;i<graphNodeNumber;i++){
        //cout<<"fffffffffffffffff--"<<i<<endl;
        if(deBruijnGraph[i].contig==NULL || NodeCount(deBruijnGraph[i].outNode)!=2){
            continue;
        }
        outNodeNumber = NodeCount(deBruijnGraph[i].outNode);
        long int endNode[outNodeNumber];
        long int adjNode[outNodeNumber];
        for(j=0;j<outNodeNumber;j++){
            endNode[j] = -1;
        }
        GraphNode * tempGraphNode = deBruijnGraph[i].outNode;
        j = 0;
        while(tempGraphNode != NULL){
            adjNode[j] = tempGraphNode->index;
            if(NodeCount(deBruijnGraph[tempGraphNode->index].inNode)==1 && NodeCount(deBruijnGraph[tempGraphNode->index].outNode) == 1){
                endNode[j] = deBruijnGraph[tempGraphNode->index].outNode->index;
            }
            j++;
            tempGraphNode = tempGraphNode->next;
        }
        //cout<<"fffffffffffffffff--"<<outNodeNumber<<endl;
        
        if(endNode[0] == -1 || endNode[1] == -1 || (endNode[0] != endNode[1])){
            continue;
        }
        
        if(NodeCount(deBruijnGraph[endNode[0]].inNode) != 2){
            continue;
        }
        
        long int upContigLength = strlen(deBruijnGraph[adjNode[0]].contig);
        long int downContigLength = strlen(deBruijnGraph[adjNode[1]].contig);
                    
        if(upContigLength!=downContigLength){
            continue;
        }
        
        
        for(j=0;j<outNodeNumber-1;j++){
            //cout<<"first--"<<endl;
            //cout<<"first--"<<endNode[j]<<endl;
            if(endNode[j] == -1){
                continue;
            }
            //cout<<"sencond"<<endl;
            for(long int p = j+1; p<outNodeNumber; p++){
                if(deBruijnGraph[adjNode[j]].contig == NULL || endNode[p] == -1){
                    continue;
                }
                if(endNode[j] == endNode[p]){
                    
                    long int misNumber = 0;
                    for(long int t = 0;t<upContigLength;t++){
                        if(deBruijnGraph[adjNode[j]].contig[t]!=deBruijnGraph[adjNode[p]].contig[t]){
                            misNumber++;
                            if(misNumber>1){
                                break;
                            }
                        }
                    }
                    if(misNumber>1){
                        continue;
                    }
                    
                    //cout<<"aa--"<<i<<"--"<<adjNode[j]<<"--"<<endNode[p]<<endl;
                    delete [] deBruijnGraph[adjNode[p]].contig;
                    delete deBruijnGraph[adjNode[p]].outNode;
                    delete deBruijnGraph[adjNode[p]].inNode;
                    deBruijnGraph[adjNode[p]].outNode = NULL;
                    deBruijnGraph[adjNode[p]].inNode = NULL;
                    deBruijnGraph[adjNode[p]].contig = NULL;
                    //cout<<"bb"<<endl;
                    tempGraphNode = deBruijnGraph[i].outNode;
                    GraphNode * previousGraphNode = NULL;
                    while(tempGraphNode != NULL){
                        if(tempGraphNode->index == adjNode[p]){
                            if(previousGraphNode!=NULL){
                               previousGraphNode->next = tempGraphNode->next;
                               //cout<<"sucess--"<<endl;
                            }else{
                               deBruijnGraph[i].outNode = tempGraphNode->next;
                            }
                            tempGraphNode->next = NULL;
                            delete tempGraphNode;
                            break;
                        }
                        previousGraphNode = tempGraphNode;
                        tempGraphNode = tempGraphNode->next;
                    }
                    //cout<<"cc--"<<deBruijnGraph[i].outNode->index<<"--"<<NodeCount(deBruijnGraph[i].outNode)<<endl;
                    tempGraphNode = deBruijnGraph[endNode[p]].inNode;
                    previousGraphNode = NULL;
                    while(tempGraphNode != NULL){
                        if(tempGraphNode->index == adjNode[p]){
                            if(previousGraphNode!=NULL){
                               previousGraphNode->next = tempGraphNode->next;
                               //cout<<"sucess11--"<<endl;
                            }else{
                               deBruijnGraph[endNode[p]].inNode = tempGraphNode->next;
                            }
                            tempGraphNode->next = NULL;
                            delete tempGraphNode;
                            break;
                        }
                        previousGraphNode = tempGraphNode;
                        tempGraphNode = tempGraphNode->next;
                    }
                    //cout<<"dd--"<<endl;
                    //exit(0);
                }
            }
        }
        
        //cout<<"endddddddddddddddddddd--"<<outNodeNumber<<endl;
        
    }
}

void RemoveCycle(DBGraph *deBruijnGraph, int kmerLength, long int graphNodeNumber, KmerSetHead * kmerSetHead){
    long int i = 0;
    long int j = 0;
    long int m = 0;
    long int n = 0;
    long int p = 0;
    long int q = 0;
    long int t = 0;
    int repeat = 2;
    long int inNodeCount = 0;
    long int outNodeCount = 0;
    while(i<graphNodeNumber){
        if(deBruijnGraph[i].contig==NULL){
            i++;
            continue;
        }
        inNodeCount = NodeCount(deBruijnGraph[i].inNode);
        outNodeCount = NodeCount(deBruijnGraph[i].outNode);
        if(inNodeCount==1&&outNodeCount==1){
                                                      
            j = deBruijnGraph[i].inNode->index;
            if(deBruijnGraph[i].inNode->index == deBruijnGraph[i].outNode->index){
                m = strlen(deBruijnGraph[i].contig);
                n = strlen(deBruijnGraph[j].contig);
                repeat = 2;
                /*
                double avg = 0;
                for(t=0; t<kmerSetHashTableHead->setNumber;t++){
                    double avgKmerFrequency = GetAvgKmerFrequency(deBruijnGraph[i].contig, t, 0, m-kmerSetHashTableHead->kmerLength+1, kmerSetHashTableHead);
                    avgKmerFrequency = avgKmerFrequency/kmerSetHashTableHead->avgKmerFrequency[t];
                    avg = avg + avgKmerFrequency;
                }
                avg = avg/kmerSetHashTableHead->setNumber;
                if(avg<3){
                    repeat = 1;
                }
                */
                char *tempContig = NULL;
                if(repeat == 1){
                    tempContig = new char[2*n+m-2*kmerLength + 3];
                    q = 0;
                    t = 0;
                    for(p = 0;p<2*n+m-2*kmerLength + 2;p++){
                        if(p<n){
                            tempContig[p] = deBruijnGraph[j].contig[p];
                        }else if(p<m+n-kmerLength+1){
                            tempContig[p] = deBruijnGraph[i].contig[kmerLength-1+q];
                            q++;
                        }else{
                            tempContig[p] = deBruijnGraph[j].contig[kmerLength-1+t];
                            t++;
                        }
                    }
                    tempContig[p] = '\0';
                }else{
                    
                    tempContig = new char[3*n+2*m-4*kmerLength + 5];
                    AppendRight(tempContig,deBruijnGraph[j].contig,deBruijnGraph[i].contig,kmerLength);
                    AppendRight(tempContig,tempContig,deBruijnGraph[j].contig,kmerLength);
                    AppendRight(tempContig,tempContig,deBruijnGraph[i].contig,kmerLength);
                    AppendRight(tempContig,tempContig,deBruijnGraph[j].contig,kmerLength);
                }
                
                delete deBruijnGraph[i].contig;
                delete deBruijnGraph[i].outNode;
                delete deBruijnGraph[i].inNode;
                deBruijnGraph[i].contig = NULL;
                deBruijnGraph[i].outNode = NULL;
                deBruijnGraph[i].inNode = NULL;
                delete deBruijnGraph[j].contig;
                deBruijnGraph[j].contig = tempContig;
                DBGRemoveInNode(deBruijnGraph+j,i);
                DBGRemoveOutNode(deBruijnGraph+j,i);
                
            }
        }
        i++;
    }
    //cout<<"RemoveCycle End!"<<endl;
}

void RemoveCycle1(DBGraph *deBruijnGraph, int kmerLength, long int graphNodeNumber){
    long int i = 0;
    long int j = 0;
    long int m = 0;
    long int n = 0;
    long int p = 0;
    long int q = 0;
    long int t = 0;
    long int inNodeCount = 0;
    long int outNodeCount = 0;

    while(i<graphNodeNumber){
        if(deBruijnGraph[i].contig==NULL){
            i++;
            continue;
        }
        inNodeCount = NodeCount(deBruijnGraph[i].inNode);
        outNodeCount = NodeCount(deBruijnGraph[i].outNode);
        if(inNodeCount==2&&outNodeCount==1){
            j = deBruijnGraph[i].outNode->index;
            if(!(NodeCount(deBruijnGraph[j].inNode)==1&&NodeCount(deBruijnGraph[j].outNode)==2)){
                i++;
                continue;
            }
            if(DBGSearchNode(deBruijnGraph[i].inNode,j)==0||DBGSearchNode(deBruijnGraph[j].outNode,i)==0){
                i++;
                continue;
            }
           
            m = strlen(deBruijnGraph[i].contig);
            n = strlen(deBruijnGraph[j].contig);
            char *tempContig = new char[n+m-kmerLength + 2];
            AppendRight(tempContig,deBruijnGraph[i].contig,deBruijnGraph[j].contig,kmerLength);
            char * tempContig1 = new char[2*strlen(tempContig)-kmerLength + 2];
            AppendRight(tempContig1,tempContig,tempContig,kmerLength);
            delete []deBruijnGraph[i].contig;
            delete tempContig;
            deBruijnGraph[i].contig = tempContig1;
            
            delete deBruijnGraph[i].outNode;
            deBruijnGraph[i].outNode = NULL;
            DBGRemoveInNode(deBruijnGraph + i,j);
            delete deBruijnGraph[j].inNode;
            deBruijnGraph[j].inNode = NULL;
            DBGRemoveOutNode(deBruijnGraph + j,i);
            deBruijnGraph[i].outNode = deBruijnGraph[j].outNode;
            DBGSearchNode(deBruijnGraph[deBruijnGraph[i].outNode->index].inNode,j,i);
            deBruijnGraph[j].outNode = NULL;
            delete []deBruijnGraph[j].contig;
            deBruijnGraph[j].contig = NULL;
                        
        }
        i++;
    }
    //cout<<"RemoveCycle1 End!"<<endl;
}

DBGraph * OptimizeDBGraphSpace(DBGraph * deBruijnGraph, unsigned long int & graphNodeCount){
    long int i = 0;
    long int j = 0;
    graphNodeCount = 0;
    while(deBruijnGraph[graphNodeCount].contig!=NULL){
        graphNodeCount++;
    }

    DBGraph * newDeBruijnGraph = new DBGraph[graphNodeCount+1];
    for(i=0;i<graphNodeCount;i++){
        newDeBruijnGraph[i].contig = deBruijnGraph[i].contig;
        newDeBruijnGraph[i].inNode = deBruijnGraph[i].inNode;
        newDeBruijnGraph[i].outNode = deBruijnGraph[i].outNode;
        deBruijnGraph[i].contig = NULL;
        deBruijnGraph[i].inNode = NULL;
        deBruijnGraph[i].outNode = NULL;
    }
    delete [] deBruijnGraph;
    return newDeBruijnGraph;
}



void GetDBGraphFromAddress(DBGraphHead * deBruijnGraphHead, char * str){
    long int i = 0;
    long int j = 0;
    long int m = 0;
    long int n = 0;
    long int p = 0;
    long int q = 0;
    long int graphNodeCount = 0;
    char *temp = new char[10000000];   	
	ifstream icin;
	icin.open(str);
	while(icin.getline(temp,10000000)){
		i++;
	}
	icin.close();
	ifstream icin1;
	icin1.open(str);
	graphNodeCount = (long int)(i/2);
	deBruijnGraphHead->nodeNumber = graphNodeCount;
    DBGraph * deBruijnGraph = new DBGraph[graphNodeCount+1];
    deBruijnGraphHead->deBruijnGraph = deBruijnGraph;
    i = 0;
    while(icin1.getline(temp,10000000)){
        int tempLength = strlen(temp);
		if(i%2!=0){
            deBruijnGraph[j].contig = new char[tempLength+1];
            strcpy(deBruijnGraph[j].contig,temp);
            deBruijnGraph[j].contig[tempLength] = '\0';
            j++;
            
        }else{
            n = 0;
            for(m=0;m<tempLength;m++){
                if(temp[m]==':'&&n==0){
                    n = 1;
                    continue;
                }
                if(temp[m]==':'&&n==1){
                   if(temp[m+1]!='0'){
                       m = m + 5;
                       while(temp[m]!='o'){
                           q = 0;
                           char * tempNum = new char[10];
                           while(temp[m]!=','){
                               tempNum[q] = temp[m];
                               m++;
                               q++;
                           }
                           tempNum[q] = '\0';
                           long int num = atoi(tempNum);
                           InsertLeftMateNode(deBruijnGraph,j,num);
                           m++;
                           delete tempNum;
                       }
                       
                   }
                   n++;
                   continue;
                }
                if(temp[m]==':'&&n==2){
                   if(temp[m+1]!='0'){
                       m = m + 5;
                       while(temp[m]!='l'){
                           q = 0;
                           char * tempNum = new char[10];
                           while(temp[m]!=','){
                               tempNum[q] = temp[m];
                               m++;
                               q++;
                           }
                           tempNum[q] = '\0';
                           long int num = atoi(tempNum);
                           InsertRightMateNode(deBruijnGraph,j,num);
                           m++;
                           delete tempNum;
                       }
                       
                   }
                   n++;
                }
            }
        }
        i++;
	}
	delete [] temp;	

}




long int FindKmerInDeBruijnGraph(char * kmer, long int kmerLength, DBGraph *deBruijnGraph, long int graphNodeCount, long int & mapPosition){
    
    for(long int i = 0; i < graphNodeCount; i++){    
        if(deBruijnGraph[i].contig != NULL){
            long int position = KMPIndexOfContig(deBruijnGraph[i].contig, kmer);
            if(position != -1){
                mapPosition = position;
                return i;
            }
        }   
    }
    return -1;
     
}



void InsertLeftMateNode(DBGraph * deBruijnGraph, long int right, long int left){
    //cout<<"ee"<<endl;
    GraphNode * temp = deBruijnGraph[right].inNode;
    while(temp!=NULL){
        if(temp->index == left){
            return;
        }
        temp = temp->next;
    }
    //cout<<"ee1"<<endl;

    temp = new GraphNode;
    
    temp->index = left;
    //cout<<right<<endl;
    temp->next = deBruijnGraph[right].inNode;
    //cout<<"ee2:"<<left<<endl;
    deBruijnGraph[right].inNode = temp;
}

void InsertRightMateNode(DBGraph * deBruijnGraph, long int left, long int right){
    GraphNode * temp = deBruijnGraph[left].outNode;
    while(temp!=NULL){
        if(temp->index == right){
            return;
        }
        temp = temp->next;
    }
    temp = new GraphNode;
    temp->index = right;
    temp->next = deBruijnGraph[left].outNode;
    deBruijnGraph[left].outNode = temp;
}

long int DBGraphInsertNode(long int previousPos, char * kmer, int kmerLength, DBGraph * deBruijnGraph, unsigned long int graphHashCount){
    unsigned long int i = 0;
    i = Hash(kmer, kmerLength, graphHashCount);

    while(deBruijnGraph[i].contig!=NULL){
        if(strcmp(deBruijnGraph[i].contig, kmer)==0){         
            if(!(previousPos == -1 || previousPos == i)){
                InsertLeftMateNode(deBruijnGraph, i, previousPos);
                InsertRightMateNode(deBruijnGraph, previousPos, i); 
            }
            return i;
        }
        i = (i + 1)%graphHashCount;
    }
    
    deBruijnGraph[i].contig = new char[kmerLength + 1];
    strncpy(deBruijnGraph[i].contig, kmer, kmerLength);
    deBruijnGraph[i].contig[kmerLength] = '\0';
    if(previousPos != -1){
        InsertLeftMateNode(deBruijnGraph, i, previousPos);
        InsertRightMateNode(deBruijnGraph, previousPos, i);
    }
    
    return i;
}

long int CreateDBGraphFromReadSet(DBGraphHead * deBruijnGraphHead, ReadSetHead * readSetHead, KmerSetHead * kmerSetHead, bool isLeftRead){
    
    KmerSet * kmerSet = NULL;
    long int kmerCount = 0;
    ReadSet * readSet = NULL;
    long int readCount = 0;
    bool isPaired = readSetHead->isPaired;
    long int readLength = readSetHead->readLength;
    if(isLeftRead == true){
        readSet = readSetHead->leftReadSet;
        readCount = readSetHead->leftReadCount;
        kmerSet = kmerSetHead->leftKmerSet;
        kmerCount = kmerSetHead->leftKmerCount;
    }else{
        readSet = readSetHead->rightReadSet;
        readCount = readSetHead->rightReadCount;
        kmerSet = kmerSetHead->rightKmerSet;
        kmerCount = kmerSetHead->rightKmerCount;
    }
    
    char * kmer = (char *)malloc(sizeof(char)*(kmerSetHead->kmerLength+1));
    long int previousNodeIndex = -1;
    for(long int i = 0; i < readCount; i++){
        previousNodeIndex = -1;
        for(long int j = 0; j < readLength - kmerSetHead->kmerLength + 1; j++){
            strncpy(kmer, readSet[i].read + j, kmerSetHead->kmerLength);
            kmer[kmerSetHead->kmerLength] = '\0';
            long int hashKey = -1;
            hashKey = SearchKmerOfKmerSet(kmer, kmerSet, kmerSetHead->kmerLength, kmerCount);
            if(hashKey !=-1){
                if(kmerSet[hashKey].kmerCount > kmerSetHead->minKmerFrequency){
                    previousNodeIndex = DBGraphInsertNode(previousNodeIndex, kmer, kmerSetHead->kmerLength, deBruijnGraphHead->deBruijnGraph, deBruijnGraphHead->nodeNumber);
                    continue;
                }
            }
            previousNodeIndex = -1;
        }
    }
    free(kmer);
}


DBGraphHead * CreateDBGraphHead(ReadSetHead * readSetHead, KmerSetHead * kmerSetHead){

    DBGraphHead * deBruijnGraphHead = (DBGraphHead *)malloc(sizeof(DBGraphHead));
    deBruijnGraphHead->deBruijnGraph = NULL;
    deBruijnGraphHead->nodeNumber = 0;
    
    char * kmer = (char *)malloc(sizeof(char)*(kmerSetHead->kmerLength+1));
    char * tempKmerRS = (char *)malloc(sizeof(char)*(kmerSetHead->kmerLength+1));
    char * tempReadRS = (char *)malloc(sizeof(char)*(readSetHead->readLength+1));
    char * read = NULL;
    for(long int i = 0; i < kmerSetHead->leftKmerCount; i++){
        if(kmerSetHead->leftKmerSet[i].kmer != NULL){
            if(kmerSetHead->leftKmerSet[i].kmerCount > kmerSetHead->minKmerFrequency){
                deBruijnGraphHead->nodeNumber++;
            }
        }
    }
    
    for(long int i = 0; i < kmerSetHead->rightKmerCount; i++){
        if(kmerSetHead->rightKmerSet[i].kmer != NULL){
            if(kmerSetHead->rightKmerSet[i].kmerCount > kmerSetHead->minKmerFrequency){
                deBruijnGraphHead->nodeNumber++;
            }
        }
    }
    
    deBruijnGraphHead->nodeNumber = deBruijnGraphHead->nodeNumber*1.2;
    cout<<"graphNodeNumber:"<<deBruijnGraphHead->nodeNumber<<endl;
    if(deBruijnGraphHead->nodeNumber == 0){
        return deBruijnGraphHead;
    }
    
    deBruijnGraphHead->deBruijnGraph = (DBGraph *)malloc(sizeof(DBGraph)*deBruijnGraphHead->nodeNumber);
    for(long int i = 0; i < deBruijnGraphHead->nodeNumber; i++){
        deBruijnGraphHead->deBruijnGraph[i].contig = NULL;
        deBruijnGraphHead->deBruijnGraph[i].outNode = NULL;
        deBruijnGraphHead->deBruijnGraph[i].inNode = NULL;
    }  
    CreateDBGraphFromReadSet(deBruijnGraphHead, readSetHead, kmerSetHead, 1);
    CreateDBGraphFromReadSet(deBruijnGraphHead, readSetHead, kmerSetHead, 0);
    
    char * graphAddress = (char *)malloc(sizeof(char)*20);
    strcpy(graphAddress, "graph.fa");
    //WriteDBGraph(deBruijnGraphHead->deBruijnGraph, deBruijnGraphHead->nodeNumber, graphAddress);
    
    SimplePathMerge(deBruijnGraphHead->deBruijnGraph,kmerSetHead->kmerLength,deBruijnGraphHead->nodeNumber);
    
    //strcpy(graphAddress, "graph0.fa");
    //WriteDBGraph(deBruijnGraphHead->deBruijnGraph, deBruijnGraphHead->nodeNumber, graphAddress);
     
    AdjustDBG(deBruijnGraphHead->deBruijnGraph,kmerSetHead->kmerLength,deBruijnGraphHead->nodeNumber);
    //strcpy(graphAddress, "graph1.fa");
    //WriteDBGraph(deBruijnGraphHead->deBruijnGraph, deBruijnGraphHead->nodeNumber, graphAddress);
    
    for(long int i = 0; i < 3;i++){
        RemoveTip(deBruijnGraphHead->deBruijnGraph,kmerSetHead->kmerLength,deBruijnGraphHead->nodeNumber);
        SimplePathMerge(deBruijnGraphHead->deBruijnGraph,kmerSetHead->kmerLength,deBruijnGraphHead->nodeNumber);
        AdjustDBG(deBruijnGraphHead->deBruijnGraph,kmerSetHead->kmerLength,deBruijnGraphHead->nodeNumber);
    }
    
    //strcpy(graphAddress, "graph2.fa");
    //WriteDBGraph(deBruijnGraphHead->deBruijnGraph, deBruijnGraphHead->nodeNumber, graphAddress);
    
    RemoveCycle(deBruijnGraphHead->deBruijnGraph,kmerSetHead->kmerLength,deBruijnGraphHead->nodeNumber, kmerSetHead);
    SimplePathMerge(deBruijnGraphHead->deBruijnGraph,kmerSetHead->kmerLength,deBruijnGraphHead->nodeNumber);
    AdjustDBG(deBruijnGraphHead->deBruijnGraph,kmerSetHead->kmerLength,deBruijnGraphHead->nodeNumber);
    RemoveCycle1(deBruijnGraphHead->deBruijnGraph,kmerSetHead->kmerLength,deBruijnGraphHead->nodeNumber);
    SimplePathMerge(deBruijnGraphHead->deBruijnGraph,kmerSetHead->kmerLength,deBruijnGraphHead->nodeNumber);
    AdjustDBG(deBruijnGraphHead->deBruijnGraph,kmerSetHead->kmerLength,deBruijnGraphHead->nodeNumber);
    
    //strcpy(graphAddress, "graph3.fa");
    //WriteDBGraph(deBruijnGraphHead->deBruijnGraph, deBruijnGraphHead->nodeNumber, graphAddress);
    
    //MergeBubble(deBruijnGraphHead->deBruijnGraph,kmerSetHead->kmerLength,deBruijnGraphHead->nodeNumber);
    //SimplePathMerge(deBruijnGraphHead->deBruijnGraph,kmerSetHead->kmerLength,deBruijnGraphHead->nodeNumber);
    //AdjustDBG(deBruijnGraphHead->deBruijnGraph,kmerSetHead->kmerLength,deBruijnGraphHead->nodeNumber);
    
    deBruijnGraphHead->deBruijnGraph = OptimizeDBGraphSpace(deBruijnGraphHead->deBruijnGraph, deBruijnGraphHead->nodeNumber);
    
    WriteDBGraph(deBruijnGraphHead->deBruijnGraph, deBruijnGraphHead->nodeNumber, graphAddress); 
    
    /*
    char * tempGapRegion = new char[709];
    strcpy(tempGapRegion, "TTTACAAAAGAAGATAAATTAAAAGTTTCATTTGATTATATTGATTGGAAGAATACAGAATTTGATCAATTAGGCCAAGAAAACTACTATAACTATAGAAAATTCGGAATTATACCAGAAATGGAATATGAAATGGAAGAGGTTAAACAAATCGAGCAATATATTAAAGAGCAAGAAGAAGCTGAACAATAGAGGCGATAACATGATTTTCGAAGAAAAACTAAATGAAATGTACAACGAGATTGCGAATAAAATTAGTAGCATGATACCAGTAGAATGGGAAAAGGTATATACAATGGCTTATATAGATGATGGAGGAGGTGAAGTATTCTTTAATTATACTAAACAGCGATGAATTGAATTATTACACCGATATACCTAAGGAGTATAACATTTCTGTGCAAGTATTTGATGATTTATGGATGGATTTATATGATTTGTTTGAGGAATTAAGAAATTTATTTAAAGAAGAAGGACTAGAACCATGGACATCATGCGAATTTGATTTTACAAGAGAAGGTAAATTAAAAGTTTCATTTGATTATATTGATTGGATAAATTCAGAATTTGGTCAAGTAGGTCGACAAAATTACTATAAGTATAGAAAATTTGGAATTTTACCAGAAACGGAATATGAAATTAATAAAGTTAAAGAAATCGAGCAATATATTAAAGAGCAAGATGAAGCTGAACTATAGGGGCGATAAT");
    tempGapRegion[708] = '\0';
    long int kmerLength = kmerSetHead->kmerLength;
    char * tempLeftKmer = new char[kmerLength + 1];
    for(long int i = 0; i < 708 - kmerLength + 1; i++){ 
        strncpy(tempLeftKmer, tempGapRegion + i, kmerLength);
        tempLeftKmer[kmerLength] = '\0';
        long int tempIndex = FindKmerInDeBruijnGraph(tempLeftKmer, kmerLength, deBruijnGraphHead->deBruijnGraph, deBruijnGraphHead->nodeNumber);    
        cout<<tempLeftKmer<<"--"<<tempIndex<<"--"<<i + kmerLength - 1<<endl;       
    }
    */
    /*
    char * tempGapRegion = new char[151];
    strcpy(tempGapRegion, "TTCATCATCGTCACTTCCTTTAGTATTCTTCTGGTAAAAGCATCACATAATAAAAAGCGTCCACGTCATCTTCACGAATGACGTAGACTTTCTTAGGTAATGCATTTTGATTTTTTTCATAGTTTGTATAGTGATATTCCAATTTGTATG");
    tempGapRegion[150] = '\0';
    long int kmerLength = kmerSetHead->kmerLength;
    char * tempLeftKmer = new char[kmerLength + 1];
    for(long int i = 0; i < 150 - kmerLength + 1; i++){ 
        strncpy(tempLeftKmer, tempGapRegion + i, kmerLength);
        tempLeftKmer[kmerLength] = '\0';
        long int tempIndex = FindKmerInDeBruijnGraph(tempLeftKmer, kmerLength, deBruijnGraphHead->deBruijnGraph, deBruijnGraphHead->nodeNumber);    
        cout<<tempLeftKmer<<"--"<<tempIndex<<"--"<<i + kmerLength - 1<<endl;       
    }
    */
    
    return deBruijnGraphHead;
    
}



void WriteDBGraph(DBGraph * deBruijnGraph, long int count, char * str){
    long int i = 0;
    long int j = 0;
    long int m = 0;
    long int n = 0;
    GraphNode * tempGraphNode;
    ofstream ocout;
    ocout.open(str);
    while(i<count){
        if(deBruijnGraph[i].contig == NULL){
            i++;
            continue;
        }
        tempGraphNode = deBruijnGraph[i].inNode;
        ocout<<">"<<i<<"-contig:inNodeCount:"<<NodeCount(tempGraphNode)<<"---";
        while(tempGraphNode!=NULL){
            ocout<<tempGraphNode->index<<",";
            tempGraphNode = tempGraphNode->next;
        }
        tempGraphNode = deBruijnGraph[i].outNode;
        ocout<<"outNodeCount:"<<NodeCount(tempGraphNode)<<"---";
        while(tempGraphNode!=NULL){
            ocout<<tempGraphNode->index<<",";
            tempGraphNode = tempGraphNode->next;
        }
        ocout<<"len:"<<strlen(deBruijnGraph[i].contig)<<endl;
        ocout<<deBruijnGraph[i].contig<<endl;
        i++;   
    }
}
