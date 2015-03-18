/*
 * Transcriptome.c
 * 
 * Copyright (c) 2011-2013 BGI-Shenzhen <soap at genomics dot org dot cn>. 
 *
 * This file is part of SOAPdenovo-Trans.
 *
 * SOAPdenovo-Trans is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SOAPdenovo-Trans is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SOAPdenovo-Trans.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include "stdinc.h" 
#include "unistd.h"
#include "stdlib.h"
#include "newhash.h"
#include "extfunc.h" 
#include "extvab.h" 
#include "dfibHeap.h"
#include "fibHeap.h"
#include "darray.h"

static DARRAY * scaf;
static DARRAY * gap;
unsigned int scaffIndex=1;
#define LINEAR 1
#define FORK 2
#define BUBBLE 3
#define COMPLEX 4
static int max_loci=0;
static unsigned int *discover ;
static unsigned int *finish ;
static unsigned int *repeat;
static char * orig;
static unsigned int * id;
static unsigned int * score;
static unsigned int  heavyContig;
static unsigned int num_unnecessary_connection;
static unsigned int num_simply_linear;
int max_step=5;
static int debug=0;
static int * flag;
static unsigned int loci_id=0;
static unsigned int loci_count=0;
static char *curr_type;
typedef struct LociHeap
{
	unsigned int contigID;
	struct LociHeap *next;
}LOCIHEAP;

static LOCIHEAP * lHeap;

void makeNewHeap(unsigned int contigID)
{
	lHeap = (LOCIHEAP*)malloc(sizeof(LOCIHEAP));
	lHeap->contigID=contigID;
	lHeap->next=NULL;
}
void InserHeap(unsigned int contigID)
{
	if(lHeap)
	{
		LOCIHEAP *oneHeap = (LOCIHEAP*)malloc(sizeof(LOCIHEAP));
		oneHeap->contigID=contigID;
		oneHeap->next=lHeap;
		lHeap=oneHeap;
	}
	else
	{
		makeNewHeap(contigID);
	}
}
unsigned int outputHeap()
{
	unsigned int oneContigID=0;
	if(lHeap)
	{
		LOCIHEAP *oneHeap = lHeap->next;
		oneContigID=lHeap->contigID;
		free(lHeap);
		lHeap=oneHeap;
	}
	return oneContigID;
}
void setUniqueContig(int cutoff)
{
	unsigned int index;
	unsigned int uniqueCount=0,other=0;
	printf("mask the short contigs, which are shorter than %d\n",cutoff);
	for(index=1;index<=num_ctg;index++)
	{
		if(contig_array[index].length >= cutoff)
		{
			contig_array[index].unique=1;	
			uniqueCount++;
		}
		else 
		{
			contig_array[index].unique=0;
			other++;
		}
	}
//	printf("The num of the unique contig:%d\nThe num of the masked contig:%d\n",uniqueCount,other);
	
}
static void ArcToConnect2()
{
	unsigned int index;
	preARC* parc=NULL;
	unsigned int to_contig,bal_contig,bal_to_contig;
	unsigned int weight;
	CONNECT * connection;
	for(index=1;index<=num_ctg;index++)
	{
		parc=contig_array[index].arcs;
		while(parc)
		{			
			weight=parc->multiplicity;
			to_contig=parc->to_ed;
			if(debug)
				printf("%d->%d\n",index,to_contig);
			if(to_contig <=0)
			{
				parc=parc->next;
				continue;
			}
			bal_contig=getTwinCtg(index);
			bal_to_contig=getTwinCtg(to_contig);
			connection=getCntBetween(index,to_contig);
			if(!connection)
			{
				connection=add1Connect(index,to_contig, 0,weight ,0);	
				if(connection)
					connection->SECount+=weight;
				connection=add1Connect(bal_to_contig,bal_contig,0,weight,0);
				if(connection)
					connection->SECount+=weight;
			}
			parc=parc->next;
		}
	}	
}
void updatedGapLength()
{
	unsigned int index;
	unsigned int uniqueCount=0,other=0;
	int count=0;
	CONNECT * connection;
	for(index=1;index<=num_ctg;index++)
	{
		connection=contig_array[index].downwardConnect;
		while(connection)
		{
			if(connection->SECount >=3 )
			{
				if(connection->gapLen < ins_size_var && connection->gapLen > -ins_size_var)
				{
					connection->gapLen=0;
					count++;
				}
				else
					printf("connection:%d->%d\tSECount:%d\tPECount:%d\tgap:%d\n",index,connection->contigID,connection->SECount,connection->PECount,connection->gapLen);
			}
			connection=connection->next;
		}		
	}
	printf("update %d gaplen of connection\n",count);
}

static  int  validConnect (unsigned int ctg,unsigned int *next_ctg )
{
	CONNECT *cn_temp;
	int count = 0;

	if (!contig_array[ctg].downwardConnect)
	{
		return count;
	}

	cn_temp = contig_array[ctg].downwardConnect;

	while (cn_temp)
	{
		if (!cn_temp->deleted && !cn_temp->mask && contig_array[cn_temp->contigID].unique)
		{
			*next_ctg=cn_temp->contigID;
			count++;
		}

		cn_temp = cn_temp->next;
	}
	return count;
		
}
static CONNECT *getNextContig (unsigned int ctg)
{
	CONNECT *cn_temp, *retCNT = NULL;
	int count = 0, valid_in;
	unsigned int nextCtg, bal_ctg;
	unsigned int next_ctg;

	bal_ctg = getTwinCtg (ctg);
	valid_in = validConnect (bal_ctg, &next_ctg);

	if (valid_in > 1)
		return NULL;

	if (!contig_array[ctg].downwardConnect)
		return NULL;

	cn_temp = contig_array[ctg].downwardConnect;
	while (cn_temp)
	{
		if (cn_temp->mask || cn_temp->deleted)
		{
			cn_temp = cn_temp->next;
			continue;
		}
		count++;
		if (count == 1)
			retCNT = cn_temp;
		else if (count == 2)
			return NULL;
		cn_temp = cn_temp->next;
	}
	return retCNT;
}
static int setConnectDelete (unsigned int from_c, unsigned int to_c, char flag)
{
	CONNECT *cn_temp, *cn_bal;
	unsigned int bal_fc = getTwinCtg (from_c);
	unsigned int bal_tc = getTwinCtg (to_c);

	cn_temp = getCntBetween (from_c, to_c);
	cn_bal = getCntBetween (bal_tc, bal_fc);

	if (!cn_temp || !cn_bal)
	{
		return 0;
	}

	cn_temp->deleted = flag;
	cn_bal->deleted = flag;	
	return 1;
}
void singleRead2connection(char * filelink)
{
	char name[1024];
	sprintf (name,"%s.ctg2Read",filelink);
	FILE * fp = NULL;
	fp=ckopen(name,"r");
	if(fp == NULL)
	{
		printf("can't open the file %s\n",name);
		exit(0);
	}
	char line[1024];
	
	unsigned long long pre_readno,readno,pre_contigID,contigID;
	int pre_pos,pos;
	char orig;
	unsigned newIndex=0;
	int flag=0;
	int  gapLen=0;
	fgets(line,sizeof(line),fp);
	pre_readno=0;
	pre_contigID=0;
	long Count=0;
	CONNECT * connection=NULL;
	while(fgets(line,sizeof(line),fp) != NULL){
		sscanf(line,"%llu\t%llu\t%d\t%c\n",&readno, &contigID,&pos,&orig);
		newIndex = index_array[contigID];
		contigID = newIndex;

		if(!contig_array[contigID].unique)
			continue;
		if(isSameAsTwin(contigID))
			continue;		

		if(pre_readno == readno && pre_contigID != contigID)
		{
			gapLen = pos - pre_pos -(int)( contig_array[pre_contigID].length);
			if(gapLen < 0 )
			{
				continue;
			}
			connection=add1Connect(pre_contigID, contigID,  gapLen, 1,0);
//			if(connection)
				connection->SECount +=1;
			connection=add1Connect(getTwinCtg(contigID), getTwinCtg(pre_contigID),  gapLen, 1,0);
//			if(connection)
				connection->SECount +=1;
			Count++;
		}

		pre_readno=readno;
		pre_contigID=contigID;
		pre_pos=pos;
	}
	printf("%ld links built by the informations of single read\n",Count);	
}
void setContigFlag(unsigned int flag)
{
	unsigned int i;
	
	for(i=1;i<=num_ctg;i++)
	{
		contig_array[i].flag=flag;
	}
}
void setOneContigFlag(unsigned int contigID,unsigned int flag)
{
        contig_array[contigID].flag=flag;
	contig_array[getTwinCtg(contigID)].flag=flag;
}

void propagateComponent(unsigned int contigID)
{
        CONNECT * connection;
        if(contig_array[contigID].flag || !contig_array[contigID].unique)
                return;
        setOneContigFlag(contigID,1);

        for(connection =contig_array[contigID].downwardConnect;
		connection != NULL; 
		connection = connection->next
		)
	{
		if(!connection->deleted)
			propagateComponent(connection->contigID);
	}
	for(connection =contig_array[getTwinCtg(contigID)].downwardConnect;
		connection != NULL; 	
		connection = connection->next
		)
	{
		if(!connection->deleted)
			propagateComponent(connection->contigID);
	}
}

unsigned int getLociCount()
{
        unsigned int lociCount=0;
        unsigned int i;

        setContigFlag(0);

        for(i=1;i<=num_ctg;i++)
        {		
                if(!contig_array[i].flag && contig_array[i].unique)
                {
                        lociCount++;
                        propagateComponent(i);
                }
        }
	
	return lociCount;
}

void fillUpComponent(unsigned int contigID,unsigned int * contigCount,int direction)
{
	CONNECT * connection;
	setOneContigFlag(contigID,1);
	*(unsigned int *)darrayPut(scaf,(* contigCount)++)=direction>0?contigID:getTwinCtg(contigID);

	for(connection = contig_array[contigID].downwardConnect;
		connection != NULL;
		connection = connection->next
		)
	{
		if( !connection->deleted 
			&& !contig_array[connection->contigID].flag 
			&& contig_array[connection->contigID].unique
			)
		{
			fillUpComponent(connection->contigID,contigCount,direction);
		}
	}
	
	for(connection = contig_array[getTwinCtg(contigID)].downwardConnect;
                connection != NULL;
                connection = connection->next
                )
        {
                if( !connection->deleted 
			&& !contig_array[connection->contigID].flag 
			&& contig_array[connection->contigID].unique
			)
                {
                        fillUpComponent(connection->contigID,contigCount,direction*(-1));
                }
        }
}

LOCUS * getLoci(unsigned int count)
{
	unsigned int i;
	unsigned int contigCount;
	unsigned int lociIndex=0;
	setContigFlag(0);

	LOCUS * loci = (LOCUS * ) malloc (count * sizeof(LOCUS));
	for(i=1;i<=num_ctg;i++)
	{
		contigCount=0;
		if(!contig_array[i].flag && contig_array[i].unique)
		{
			fillUpComponent(i,&contigCount,1);	
			loci[lociIndex].contigCount= contigCount;
			loci[lociIndex].contigID = 
				(unsigned int *)malloc((contigCount) * sizeof(unsigned int));
			for(contigCount=0;contigCount<loci[lociIndex].contigCount;contigCount++)
			{
				loci[lociIndex].contigID[contigCount] = 
					*(unsigned int *) darrayGet(scaf,contigCount);
			}
			lociIndex++;
		}
	}
	return loci;
}
static CONNECT *checkConnect(unsigned int from_c,unsigned int to_c)
{
	CONNECT *cn_temp=getCntBetween(from_c,to_c);
	if(!cn_temp)
		return NULL;
	if(!cn_temp->mask&&!cn_temp->deleted)
		return cn_temp;
	return NULL;
}
static void checkCircle()
{
	unsigned int i,ctg;
	CONNECT *cn_temp1,*cn_temp2;
	int counter=0;

	for(i=1;i<=num_ctg;i++){
		cn_temp1 = contig_array[i].downwardConnect;
		while(cn_temp1){
			if(cn_temp1->weak||cn_temp1->deleted){
				cn_temp1 = cn_temp1->next;
				continue;
			}
			ctg = cn_temp1->contigID;
			cn_temp2=getCntBetween(ctg,i);
			if(cn_temp2){
				counter++;
				if(cn_temp1->weight > cn_temp2->weight)
					cn_temp2->deleted=1;
				else
					cn_temp1->deleted=1;
			}
			cn_temp1 = cn_temp1->next;
		}

	}
	printf("%d simple loops found and removed in graph\n",counter);
}
static void deleteWeakCnt(int cut_off)
{
	unsigned int i;
	CONNECT *cn_temp1;
	int weaks=0,counter=0;

	for(i=1;i<=num_ctg;i++){
		cn_temp1 = contig_array[i].downwardConnect;
		while(cn_temp1){
			if(!cn_temp1->mask&&!cn_temp1->deleted){
				counter++;
			}
			if(cn_temp1->weak&&cn_temp1->deleted&&cn_temp1->weight>=cut_off){
				cn_temp1->deleted = 0;
				cn_temp1->weak = 0;
			}
			else if(!cn_temp1->deleted&&cn_temp1->weight>0&&cn_temp1->weight<cut_off)
			{
				cn_temp1->deleted = 1;
				cn_temp1->weak = 1;
				if(!cn_temp1->mask)
					weaks++;
			}
			cn_temp1 = cn_temp1->next;
		}

	}
	printf("%d weak connections(<%d) removed (there were %d active cnnects))\n",weaks,cut_off,counter);
	checkCircle();
}
static void deleteInconsistent(unsigned int lociCount,LOCUS * loci)
{
	unsigned int i,j,contigID,bal_contigID,dest_contigID;
	CONNECT * connection,*bal_connection,*temp;
	LOCUS * curr_loci=NULL;

	for(i=0;i<lociCount;i++)
	{
		curr_loci=&(loci[i]);
		if(curr_loci->contigCount ==1)
			continue;
		for(j=0;j<curr_loci->contigCount;j++)
		{
			contigID=curr_loci->contigID[j];
			bal_contigID=getTwinCtg(contigID);
			contig_array[contigID].flag=0;
			contig_array[bal_contigID].flag=1;
		}

		for(j=0;j<curr_loci->contigCount;j++)
		{
			contigID=curr_loci->contigID[j];
			bal_contigID=getTwinCtg(contigID);
			connection=contig_array[contigID].downwardConnect;
			while(connection)
			{
				dest_contigID=connection->contigID;
				if(!connection->deleted && contig_array[dest_contigID].unique)
				{
					if(contig_array[dest_contigID].flag == 1)					
					{
						connection->deleted=1;
						temp=getCntBetween(getTwinCtg(dest_contigID),bal_contigID);
//						if(temp!=NULL)
						{					
							temp->deleted=1;
						}
					}
				}
				connection=connection->next;
			}
			connection=contig_array[bal_contigID].downwardConnect;
			while(connection)
			{
				dest_contigID=connection->contigID;
				if(!connection->deleted && contig_array[dest_contigID].unique)
				{
					if(contig_array[dest_contigID].flag == 0)					
					{
						connection->deleted=1;
						temp=getCntBetween(getTwinCtg(dest_contigID),contigID);
//						if(temp!=NULL)
						{					
							temp->deleted=1;
						}
					}
				}
				connection=connection->next;
			}
		}
	}
}
void traceAlongConnection(unsigned int destContig,unsigned int currContig,CONNECT * curr_connection,int max_steps,
	int min,int max,int index,int len,int gap_length,int *num_route)
{
	int pos=index,i;
	unsigned int * array;
	if(pos > max_steps || len > max || *num_route >= max_step-1)
		return ;
	if (currContig == destContig && pos == 0)
	{
		printf ("traceAlongConnection: start and destination are the same\n");
		return;
	}
	if(pos >0)
		so_far[pos-1]=currContig;
	if(currContig == destContig && len>=min)
	{				
		array = found_routes[*num_route];
		for(i=0;i<pos;i++)
			array[i]=so_far[i];
		if(pos<max_steps)
			array[pos]=0;
		(*num_route)++;
	}
	if(pos==max_steps || len >= max)
		return ;
	if(pos++ >0)
		len += contig_array[currContig].length+gap_length;
	CONNECT * connection=contig_array[currContig].downwardConnect;
	while(connection)
	{
		if(!(connection->deleted)&& contig_array[connection->contigID].unique && connection!=curr_connection)
			traceAlongConnection(destContig,connection->contigID,curr_connection,max_steps,
				min,max,pos,len,connection->gapLen,num_route);
		connection=connection->next;		
	}
	
}
static int linearC2C(unsigned int starter, CONNECT * cnt2c1, unsigned int c2, int min_dis, int max_dis)
{
	int out_num,in_num;
	CONNECT *curr_cn;
	int len=0;
	unsigned int ctg,c1,bal_c1,next_c1,bal_c2,prev_c2,bal_start;

	c1=cnt2c1->contigID;
	
	if (c1 == c2)
	{
		printf ("linearC2C: c1(%d) and c2(%d) are the same contig\n", c1, c2);
		return -1;
	}

	bal_c1=getTwinCtg(c1);
	in_num=validConnect(bal_c1, &next_c1);

	if(in_num >1)
		return 0;

	bal_start=getTwinCtg(starter);
	curr_cn=cnt2c1;

	while((curr_cn=getNextContig(c1))!=NULL)
	{
		c1=curr_cn->contigID;
		len += curr_cn->gapLen + contig_array[c1].length;
		if(c1==c2)
			return 1;
		if(len > max_dis || c1==starter || c1==bal_start)
			return 0;
	}

	bal_c2=getTwinCtg(c2);
	curr_cn=NULL;
	bal_c1=getTwinCtg(c1);
	ctg=bal_c2;

	while((curr_cn=getNextContig(ctg))!=NULL)
	{
		ctg=curr_cn->contigID;
		len += curr_cn->gapLen + contig_array[c1].length;
		if(ctg==bal_c1)
			return 1;
		if(len > max_dis || c1==starter || c1==bal_start)
			return 0;
	}

	c2=getTwinCtg(ctg);
	min_dis -= len;
	max_dis -=len;
	if(max_dis <0)
		return 0;
	
	curr_cn=getCntBetween(c1,c2);
	if(curr_cn)
	{
		setConnectDelete(c1,c2, 0);
		return 1;
	}

	len = (min_dis + max_dis)/2 >0 ?(min_dis + max_dis)/2:0;
	curr_cn = allocateCN(c2,len);
	curr_cn->weight=0;
	curr_cn->next=contig_array[c1].downwardConnect;
	contig_array[c1].downwardConnect=curr_cn;
	
	bal_c1=getTwinCtg(c1);
	bal_c2=getTwinCtg(c2);
	curr_cn=allocateCN(bal_c2,len);
	curr_cn->weight=0;
	curr_cn->next=contig_array[c1].downwardConnect;
	contig_array[c1].downwardConnect=curr_cn;
	return 1;	
}
static void bal_simply_linear(LOCUS * locus)
{
	int index,locusCount=locus->contigCount;
	CONNECT *connection,*cn1,*cn2,*temp_cn;
	int out_num;
	int min,max;
	int linear;
	for(index=0;index<locusCount;index++)
	{
		connection=contig_array[getTwinCtg(locus->contigID[index])].downwardConnect;
		out_num=0;
		while(connection)
		{
			if(!(connection->deleted) && contig_array[connection->contigID].unique)
			{
				out_num++;
				if(out_num == 1)
				{
					cn1=connection;
				}
				else if(out_num == 2)
				{
					cn2=connection;
				}
				else
					break;				
			}
			connection=connection->next;
		}
		if(out_num != 2)
			continue;
		if(cn1->gapLen > cn2->gapLen)
		{
			temp_cn = cn1;
			cn1=cn2;
			cn2=temp_cn;
		}

		min=cn2->gapLen - cn1->gapLen - contig_array[cn1->contigID].length - ins_size_var/2;
		max=cn2->gapLen - cn1->gapLen - contig_array[cn1->contigID].length + ins_size_var/2;

		if(max < 0)
			continue;
		setConnectDelete(getTwinCtg(locus->contigID[index]),cn2->contigID,1);
		linear=linearC2C(getTwinCtg(locus->contigID[index]),cn1,cn2->contigID,min,max);
		if(!linear)
			setConnectDelete(getTwinCtg(locus->contigID[index]),cn2->contigID,0);
		else
			num_simply_linear++;
	}
}
static void simply_linear(LOCUS * locus)
{
	int index,locusCount=locus->contigCount;
	CONNECT *connection,*cn1,*cn2,*temp_cn;
	int out_num;
	int min,max;
	int linear;
	for(index=0;index<locusCount;index++)
	{
		connection=contig_array[locus->contigID[index]].downwardConnect;
		out_num=0;
		while(connection)
		{
			if(!(connection->deleted) && contig_array[connection->contigID].unique)
			{
				out_num++;
				if(out_num == 1)
				{
					cn1=connection;
				}
				else if(out_num == 2)
				{
					cn2=connection;
				}
				else
					break;				
			}
			connection=connection->next;
		}
		if(out_num != 2)
			continue;
		if(cn1->gapLen > cn2->gapLen)
		{
			temp_cn = cn1;
			cn1=cn2;
			cn2=temp_cn;
		}

		min=cn2->gapLen - cn1->gapLen - contig_array[cn1->contigID].length - ins_size_var/2;
		max=cn2->gapLen - cn1->gapLen - contig_array[cn1->contigID].length + ins_size_var/2;

		if(max < 0)
			continue;
		setConnectDelete(locus->contigID[index],cn2->contigID,1);
		linear=linearC2C(locus->contigID[index],cn1,cn2->contigID,min,max);
		if(!linear)
			setConnectDelete(locus->contigID[index],cn2->contigID,0);
		else
			num_simply_linear++;
	}
}
void deleteUnnecessary(LOCUS * locus)
{
	int i,	 locusCount=locus->contigCount;
	CONNECT * bal_connection,*connection;
	int flag;
	int num_route=0;
	
	for(i=0;i<locusCount;i++)
	{
		connection=contig_array[locus->contigID[i]].downwardConnect;
		while(connection)
		{
			flag=0;
			num_route=0;
			if(!(connection->deleted) && contig_array[connection->contigID].unique)
			{
				if(!connection->SECount && connection->PECount)
				{
					traceAlongConnection(connection->contigID,locus->contigID[i],connection,max_step,
						connection->gapLen-2*ins_size_var,connection->gapLen+2*ins_size_var,
						0,0,0,&num_route);	
					if(num_route)
					{
						num_unnecessary_connection ++ ;
						setConnectDelete(locus->contigID[i],connection->contigID,1);
					}
				}
			}
			connection=connection->next;
		}
	}	
}
void linearization(LOCUS * loci,unsigned int lociCount)
{
	unsigned int i,j;
	LOCUS * locus;
	setContigFlag(0);
	so_far = (unsigned int *) ckalloc (2*max_n_routes * sizeof (unsigned int));
	found_routes = (unsigned int **) ckalloc (2*max_step * sizeof (unsigned int *));

	for (j = 0; j < max_step; j++)
	{
		found_routes[j] = (unsigned int *) ckalloc (max_n_routes * sizeof (unsigned int));
	}

	for(i=0;i<lociCount;i++)
	{
		if(debug)
			printf("%d\t%d\n",i,lociCount);
		locus=&loci[i];
		if(locus->contigCount <=2)
			continue;
//		simply_linear(locus);
//		bal_simply_linear( locus);
		deleteUnnecessary(locus);
	}

	printf("Delete the unnecessary connections : %d+%d\n",num_unnecessary_connection,num_simply_linear);
	
	free(so_far);
	for (j = 0; j < max_step; j++)
	{
		free(found_routes[j]);
	}
	free(found_routes);
}
void tourLoci(unsigned int contigID,unsigned int * order,unsigned int *id)
{
	if(orig[contigID]=='-')
		return;
	id[*order] = contigID;
	discover[getTwinCtg(contigID)] = *order;
	discover[contigID]=	(*order)++;
	
	CONNECT * connection=NULL ;
	
	connection = contig_array[contigID].downwardConnect;
	while(connection)
	{
		if(!(connection->deleted) &&  contig_array[connection->contigID].unique)
		{
			if(discover[connection->contigID] == 0)	
			{
				tourLoci(connection->contigID, order, id);
			}
		}
		connection = connection->next;
	}
	id[*order] = contigID;
	finish[getTwinCtg(contigID)]= *order;
	finish[contigID]=(*order)++;

}
void found_repeat(unsigned int contigID,int *total)
{
	if(orig[contigID]=='-')
		return;
	repeat[(*total)++]=contigID;
	discover[contigID]=0;
	discover[getTwinCtg(contigID)]=0;
	finish[contigID]=0;
	finish[getTwinCtg(contigID)]=0;
	unsigned int  twin = getTwinCtg(contigID);
	CONNECT * connection = contig_array[twin].downwardConnect;

	while(connection)
	{
		if(!(connection->deleted) &&  contig_array[connection->contigID].unique )
		{
			twin = getTwinCtg(connection->contigID);
			if(discover[twin] >0)
				found_repeat(twin, total);
		}
		connection = connection->next;
	}
}
unsigned int repeatCounter=0;
void avoidoneLociLoop(LOCUS * loci)
{	
	unsigned int contigCount = loci->contigCount;	
	loci->repeatMark=0;
	unsigned int i,order=1,contigID;
	CONNECT * connection=NULL ;
	for(i=0;i<contigCount;i++)
	{
		id[i]=0;
		id[i+contigCount]=0;
	}
	for(i=0;i<contigCount;i++)
	{
		contigID = loci->contigID[i];
		orig[contigID]='+';
	}
	for(i=0;i<contigCount;i++)
	{
		contigID = loci->contigID[i];
		if(discover[contigID] == 0)		
		{					
			tourLoci(contigID,&order,id);
		}		
	}
	if(debug)
	{
		for(i=1;i<=2*contigCount;i++)
		{
			printf("id[%d]:%d\n",i,id[i]);
		}
	}

	int total,j,min_j=-1;
	int minWeight;
	int m,n,curr_in_time,curr_out_time,time;
	for(i=2*contigCount;i>0;i--)
	{
		total=0;
		time=0;
		minWeight=9999;
		if(discover[id[i]] != 0 )
		{
			found_repeat(id[i], &total);
			min_j=-1;
			if(total > 1)
			{
				repeatCounter++;
				for(m = 0;m<total;m++)
				{
					curr_in_time=0;
					curr_out_time=0;
					for(n =0;n<total;n++)
					{
						if(m == n)
							continue;
						
						connection=getCntBetween(repeat[n],repeat[m]);
						if(connection && (connection->deleted !=1))
							curr_in_time++;
						
						connection=getCntBetween(repeat[m],repeat[n]);
						if(connection && (connection->deleted !=1))
							curr_out_time++;
					}
					if(curr_in_time >1)
					{
						time++;
						for(n=0;n<total;n++)
						{
							if(m==n)
								continue;
							connection=getCntBetween(repeat[n],repeat[m]);
							if(connection)
								connection->deleted=1;
							connection=getCntBetween(getTwinCtg(repeat[m]),getTwinCtg(repeat[n]));
							if(connection)
								connection->deleted=1;
						}
					}
					if(curr_out_time>1)
					{
						time++;
						for(n=0;n<total;n++)
						{
							if(m==n)
								continue;
							connection=getCntBetween(repeat[m],repeat[n]);
							if(connection)
								connection->deleted=1;
							connection=getCntBetween(getTwinCtg(repeat[n]),getTwinCtg(repeat[m]));
							if(connection)
								connection->deleted=1;
						}
					}
				}	
				for(j=total-1;j>=0;j--)
				{			
					if(j == 0)
					{
						connection = getCntBetween(repeat[0],repeat[total-1]);		
						if(connection && connection->deleted != 1)
						{
							if(minWeight > connection->weight)
							{
								minWeight = connection->weight;
								min_j = j;
							}
						}
						else
						{
							min_j=-2;
							break;
						}
					}
					else
					{
						connection = getCntBetween(repeat[j],repeat[j-1]);
						if(connection  && connection->deleted != 1)
						{
							if(minWeight > connection->weight)
							{
								minWeight = connection->weight;
								min_j = j;
							}
						}
						else
						{
							min_j=-2;
							break;
						}
					}
				}
				if(min_j == -1)
				{
					printf("error : repeat question\n");
				}
				else if(min_j == 0)
				{
					connection = getCntBetween(repeat[0],repeat[total-1]);
					connection->deleted=1;
					connection=getCntBetween(getTwinCtg(repeat[total-1]),getTwinCtg(repeat[0]));
					connection->deleted=1;
				}
				else if(min_j >0)
				{
					connection = getCntBetween(repeat[min_j],repeat[min_j-1]);
					connection->deleted=1;
					connection=getCntBetween(getTwinCtg(repeat[min_j-1]),getTwinCtg(repeat[min_j]));
					connection->deleted=1;
				}									
			}			
		}
	}
}
void avoidLoop(unsigned int lociCount,LOCUS * loci)
{
	unsigned int i;
	discover = (unsigned int *)malloc((num_ctg+1)* sizeof(unsigned int));
	finish = (unsigned int *) malloc((num_ctg+1)* sizeof (unsigned int ));
	repeat = (unsigned  int *)malloc((num_ctg+1) * sizeof(unsigned int ));
	id = (unsigned  int *)malloc(2*(num_ctg+1) * sizeof(unsigned int ));
	orig=(char*)malloc((num_ctg+1) * sizeof(char));
	for(i=1;i<=num_ctg;i++)
	{
		discover[i]=0;
		finish[i]=0;
		repeat[i]=0;
		orig[i]='-';
	}
	for(i=0;i<lociCount;i++)
	{
		if(debug)
			printf("%d\t%d\n",i,lociCount);
		if(loci[i].contigCount >1 )
			avoidoneLociLoop(&(loci[i]));
	}
	free(discover);
	free(finish);
	free(repeat);
	free(id);
	
	printf("%d complex loops found and removed in graph\n",repeatCounter);
}
static int num_connection[256];


static unsigned int getLocusKind(LOCUS * locus)
{
	int distribution[3]={0,0,0};
	unsigned int i ;
	int connectionCount;
	CONNECT * connection;
	if(locus->contigCount <= 2)
		return LINEAR;
	
	for(i=0;i<locus->contigCount;i++)
	{
		connectionCount =0 ;		
		connection = contig_array[locus->contigID[i]].downwardConnect;
		while(connection)
		{
			if(!connection->deleted && contig_array[connection->contigID].unique)
			{
				connectionCount++;
			}
			connection = connection->next;
		}
		if(connectionCount == 0)
		{
			distribution[0]++;
		}
		else if(connectionCount == 2)
		{
			distribution[1]++;
		}
		else if(connectionCount >= 3)
		{
			distribution[2]++;
		}

		if(connectionCount<256)
			num_connection[connectionCount]++;
		connectionCount =0 ;
                connection = contig_array[getTwinCtg(locus->contigID[i])].downwardConnect;
                while(connection)
                {
                     if(!connection->deleted && contig_array[connection->contigID].unique)
			{
				connectionCount++;
			}
                     connection = connection->next;
                }
                if(connectionCount == 0)
                {
                        distribution[0]++;
                }
                else if(connectionCount == 2)
                {
                        distribution[1]++;
                }
                else if(connectionCount >= 3)
                {
                        distribution[2]++;
                }	

		if(connectionCount<256)
			num_connection[connectionCount]++;
	}		
	if(distribution[0] == 2 && distribution[1] == 0 && distribution[2] == 0)
	{
		return LINEAR;
	}
	else if(distribution[0] == 3 && distribution[1] == 1 && distribution[2] == 0)
	{
		return FORK;
	}
	else if(distribution[0] == 2 && distribution[1] == 2 && distribution[2] == 0)
	{
		return BUBBLE;
	}
	else
		return COMPLEX;
}

void outputOneTranscriptome(FILE * fq,FILE * fo,unsigned int count,unsigned int len)
{
	if(count==1)
		return;
	unsigned int prev_ctg,index , max_steps=5,currentContigID;	
	int num_route,gap_c,j;
//	boolean weak;	
	short gapLen=0;	

	fprintf(fq,">scaffold%d %d %d Locus_%d_%d %s\n",scaffIndex,count,len,loci_id,loci_count,curr_type);
	fprintf(fo,">scaffold%d %d %d Locus_%d_%d %s\n",scaffIndex,count,len,loci_id,loci_count,curr_type);
	scaffIndex++;
	loci_count++;
	len = prev_ctg = 0;
	so_far = (unsigned int *)ckalloc(max_n_routes*sizeof(unsigned int));
	found_routes = (unsigned int **)ckalloc(max_steps*sizeof(unsigned int *));
	for(j=0;j<max_steps;j++)
		found_routes[j] = (unsigned int *)ckalloc(2*max_n_routes*sizeof(unsigned int));
	for(index=0;index<count;index++)
	{
		currentContigID = *(unsigned int *)darrayGet(scaf,index);	
		if(!isLargerThanTwin(currentContigID))
		{			
			fprintf(fq,"%-10d %-10d +   %d \n",
				index_array[currentContigID],
				len,
				contig_array[currentContigID].length+overlaplen);
//			weak = printCnts(fq,currentContigID);
		}
		else
		{
			fprintf(fq,"%-10d %-10d -   %d \n",
				index_array[getTwinCtg(currentContigID)],
				len,	
				contig_array[currentContigID].length+overlaplen);
//			 weak = printCnts(fq,currentContigID);
		}
		if(prev_ctg)
		{
			num_route = num_trace = 0;
			traceAlongArc(currentContigID,prev_ctg,max_steps,
				gapLen-ins_size_var,gapLen+ins_size_var,0,0,&num_route);
			if(num_route == 1)
			{
				output1gap(fo,max_steps);
				gap_c++;
			}		
		}
		fprintf(fo,"%-10d %-10d\n",
			currentContigID,
			len
			);
		len += contig_array[currentContigID].length + *(unsigned int *)darrayGet(gap,index);
		prev_ctg = currentContigID;
		gapLen = *(unsigned int *)darrayGet(gap,index)>0?*(unsigned int *)darrayGet(gap,index):0;
	}
	for(j=0;j<max_steps;j++)
		free((void *)found_routes[j]);
	free((void *)found_routes);
	free((void *)so_far);
}

void outputLinearTranscriptome(LOCUS * locus ,FILE * fq,FILE * fo)
{
	unsigned int index;
	unsigned int leftContig=0;
	unsigned int twin;
	CONNECT * connection;
	unsigned int out;
	unsigned int len=0;
	if(locus->contigCount == 1)
	{
		return;
	}
	for(index = 0; index < locus->contigCount ; index++)
	{
		twin = getTwinCtg(locus->contigID[index]);
		connection = contig_array[twin].downwardConnect;
		out=0;
		while(connection)
		{
			if(!connection->deleted && contig_array[connection->contigID].unique)
			{
				out++;
			}
			connection = connection->next;
		}
		if(!out)
		{
			leftContig = locus->contigID[index];
			break;	
		}
	}
	index = 0;
	len = contig_array[leftContig].length;
	*(unsigned int*)darrayPut(scaf,index++) = leftContig;	
	while(leftContig)
	{		
		connection = contig_array[leftContig].downwardConnect;					
		leftContig = 0;		
		while(connection)
		{
			if(!connection->deleted && contig_array[connection->contigID].unique)
			{
				leftContig = connection->contigID;
				*(int *)darrayPut(gap,index-1) = connection->gapLen;
				*(unsigned int*)darrayPut(scaf,index++) = leftContig;
				len += contig_array[leftContig].length + connection->gapLen;
				break;
			}	
			connection = connection ->next;
		}		
	}
	outputOneTranscriptome(fq,fo,index,len);
}
void outputForkTranscriptome(LOCUS * locus ,FILE * fq,FILE * fo)
{
	unsigned int index;
	unsigned int leftContig1=0,leftContig2=0;
	unsigned int twin;
	CONNECT * connection;
	unsigned int out;
	unsigned int len=0;
	for(index = 0; index < locus->contigCount ; index++)
	{
		twin = getTwinCtg(locus->contigID[index]);
		connection = contig_array[twin].downwardConnect;
		out=0;
		while(connection)
		{
			if(!connection->deleted && contig_array[connection->contigID].unique)
			{
				out++;
			}
			connection = connection->next;
		}
		if(!out)
		{
			if(!leftContig1)
				leftContig1= locus->contigID[index];
			else
				leftContig2= locus->contigID[index];				
		}
		if(leftContig2)
			break;
	}
	if(leftContig2)
	{
		index = 0;
		len = contig_array[leftContig1].length;
		*(unsigned int*)darrayPut(scaf,index++) = leftContig1;
		while(leftContig1)
		{
			connection = contig_array[leftContig1].downwardConnect;						
			leftContig1= 0;			
			while(connection)
			{
				if(!connection->deleted && contig_array[connection->contigID].unique)
				{
					leftContig1= connection->contigID;
					*(int *)darrayPut(gap,index-1) = connection->gapLen;
					*(unsigned int*)darrayPut(scaf,index++) = leftContig1;
					len += contig_array[leftContig1].length + connection->gapLen;
					break;
				}	
				connection = connection->next;
			}			
		}	
		outputOneTranscriptome(fq,fo,index,len);
		index = 0;
		len = contig_array[leftContig2].length;
		*(unsigned int*)darrayPut(scaf,index++) = leftContig2;
		while(leftContig2)
		{
			connection = contig_array[leftContig2].downwardConnect;							
			leftContig2= 0;			
			while(connection)
			{
				if(!connection->deleted && contig_array[connection->contigID].unique)
				{
					leftContig2= connection->contigID;
					*(int *)darrayPut(gap,index-1) = connection->gapLen;
					*(unsigned int*)darrayPut(scaf,index++) = leftContig2;
					len += contig_array[leftContig2].length + connection->gapLen;
					break;
				}	
				connection = connection->next;
			}			
		}	
		outputOneTranscriptome(fq,fo,index,len);
		
	}
	else
	{
		unsigned int leftContig = leftContig1;
		unsigned int forkContig1=0;
		unsigned int forkContig2=0;
		while(leftContig)
		{			
			connection = contig_array[leftContig].downwardConnect;
			leftContig = 0;
			while(connection)
			{
				if(!connection->deleted && contig_array[connection->contigID].unique)
				{
					leftContig = connection->contigID;
					if(!forkContig1)
						forkContig1=leftContig;
					else
						forkContig2=leftContig;
				}	
				connection = connection->next;
			}
			if(forkContig2)
				break;
			else
				forkContig1 = 0;
		}
		index = 0;		
		leftContig = leftContig1;
		len = contig_array[leftContig].length;
		*(unsigned int*)darrayPut(scaf,index++) = leftContig;	
		while(leftContig)
		{
			connection = contig_array[leftContig].downwardConnect;							
			leftContig = 0;		
			while(connection)
			{
				if(!connection->deleted && contig_array[connection->contigID].unique)
				{
					leftContig = connection->contigID;
					*(int *)darrayPut(gap,index-1) = connection->gapLen;
					*(unsigned int*)darrayPut(scaf,index++) = leftContig;	
					len += contig_array[leftContig].length + connection->gapLen;
					break;
				}	
				connection = connection->next;
			}
			
		}	
		outputOneTranscriptome(fq,fo,index,len);
		index = 0;
		leftContig = leftContig1;
		len = contig_array[leftContig].length;
		*(unsigned int*)darrayPut(scaf,index++) = leftContig;
		while(leftContig)
		{
			connection = contig_array[leftContig].downwardConnect;						
			leftContig = 0;			
			while(connection)
			{
				if(!connection->deleted && contig_array[connection->contigID].unique)
				{				
					leftContig = connection->contigID;
					if(leftContig != forkContig1)
					{
						*(int *)darrayPut(gap,index-1) = connection->gapLen;
						*(unsigned int*)darrayPut(scaf,index++) = leftContig;
						len += contig_array[leftContig].length + connection->gapLen;
						break;
					}
				}	
				connection = connection->next;
			}			
		}	
		outputOneTranscriptome(fq,fo,index,len);
	}
}
void outputBubbleTranscriptome(LOCUS * locus ,FILE * fq,FILE * fo)
{
	unsigned int index;
	unsigned int leftContig=0;
	unsigned int twin;
	CONNECT * connection;
	unsigned int out;
	unsigned int len=0;
	for(index = 0; index < locus->contigCount ; index++)
	{
		twin = getTwinCtg(locus->contigID[index]);
		connection = contig_array[twin].downwardConnect;
		out=0;
		while(connection)
		{
			if(!connection->deleted && contig_array[connection->contigID].unique)
			{
				out++;
			}
			connection = connection->next;
		}
		if(!out)
		{
			leftContig = locus->contigID[index];
			break;	
		}
	}
	if(debug)
		printf("leftcontig:\t%d\n",leftContig);
	
	unsigned int leftTemp = leftContig;
	unsigned int bubbleContig1=0;
	unsigned int bubbleContig2=0;
	unsigned int tempContig;
	while(leftContig)
	{		
		connection = contig_array[leftContig].downwardConnect;
		tempContig = leftContig;
		leftContig = 0;
		while(connection)
		{
			if(!connection->deleted && contig_array[connection->contigID].unique)
			{
				if(debug)
					printf("%d->%d\n",tempContig,connection->contigID);
				leftContig = connection->contigID;
				if(!bubbleContig1)
					bubbleContig1=leftContig;
				else
					bubbleContig2=leftContig;
			}	
			connection = connection->next;
		}
		if(bubbleContig2)
			break;
		else
			bubbleContig1=0;
		
	}
	if(debug)
	{
		printf("bubbleContig1:\t%d\n",bubbleContig1);
		printf("bubbleContig2:\t%d\n",bubbleContig2);
	}

	index = 0;
	leftContig = leftTemp;
	len = contig_array[leftContig].length;
	*(unsigned int*)darrayPut(scaf,index++) = leftContig;	
	while(leftContig)
	{
		connection = contig_array[leftContig].downwardConnect;				
		leftContig = 0;
		while(connection)
		{
			if(!connection->deleted && contig_array[connection->contigID].unique)
			{
				leftContig = connection->contigID;
				*(int *)darrayPut(gap,index-1) = connection->gapLen;
				*(unsigned int*)darrayPut(scaf,index++) = leftContig;
				len += contig_array[leftContig].length + connection->gapLen;
				break;
			}	
			connection = connection->next;
		}		
	}	
	outputOneTranscriptome(fq,fo,index,len);
	
	index = 0;
	leftContig = leftTemp;
	len = contig_array[leftContig].length;
	*(unsigned int*)darrayPut(scaf,index++) = leftContig;
	int first=1;
	while(leftContig)
	{
		connection = contig_array[leftContig].downwardConnect;		
		leftContig = 0;
		while(connection)
		{
			if(!connection->deleted && contig_array[connection->contigID].unique)
			{				
				leftContig = connection->contigID;
				if(leftContig != bubbleContig1 || !first)
				{
					*(int *)darrayPut(gap,index-1) = connection->gapLen;
					*(unsigned int*)darrayPut(scaf,index++) = leftContig;
					len += contig_array[leftContig].length + connection->gapLen;
					break;
				}
				else
					first=0;
			}	
			connection = connection->next;
		}
		
	}	
	outputOneTranscriptome(fq,fo,index,len);
}
void computeScore(unsigned int contigID)
{
	CONNECT * connection=NULL;
	double maxWeight=0,maxScore=0;
	unsigned int  pre_contigID=0;

	connection=contig_array[getTwinCtg(contigID)].downwardConnect;
	while(connection)
	{
		if(!connection->deleted && contig_array[connection->contigID].unique){
			if(heavyContig && getTwinCtg(connection->contigID) == heavyContig){
				maxWeight = connection->weight;
				pre_contigID = getTwinCtg(connection->contigID);
				maxScore=10000*connection->weight+score[pre_contigID];
				break;
				//mao 2011 10 18
			}else if(heavyContig && (connection->weight + score[getTwinCtg(connection->contigID)]> maxScore))
			{
				
				maxWeight = connection->weight;
				pre_contigID=getTwinCtg(connection->contigID);
				maxScore=connection->weight+score[pre_contigID];
			}
			else if(connection->weight  > maxWeight)
			{
				maxWeight = connection->weight;
				pre_contigID=getTwinCtg(connection->contigID);
			}
		}
		connection=connection->next;
	}
	if(pre_contigID == 0)
		return;
	if(heavyContig && (contigID == heavyContig || pre_contigID == heavyContig))
		score[contigID]=10000*maxWeight+score[pre_contigID] >score[contigID] ? 10000*maxWeight+score[pre_contigID] : score[contigID]  ;
	else
		score[contigID]=maxWeight+score[pre_contigID] >score[contigID] ?  maxWeight+score[pre_contigID] : score[contigID];
}
	

void propagate(unsigned int contigID,LOCUS * locus)
{
	CONNECT* connection=NULL;
	int i=0;

	connection=contig_array[contigID].downwardConnect;
	while(connection)
	{
		if(!connection->deleted && contig_array[connection->contigID].unique)
		{
			InserHeap(connection->contigID);
		}
		connection=connection->next;
	}
}
void getBestWay(unsigned int contigID,unsigned int *count)
{	
	if(contigID ==0)
	{
		*count=0;
		return;
	}
	InserHeap(contigID);
	(*count)++;
	unsigned int twin = getTwinCtg(contigID);
	CONNECT * connection = contig_array[twin].downwardConnect;
	double maxWeight=0,maxScore=0;
	unsigned int  pre_contigID=0;
	
	while(connection)
	{
		if(!connection->deleted && contig_array[connection->contigID].unique)
		{

			if(contig_array[getTwinCtg(connection->contigID)].flag == 1)
			{
				connection=connection->next;
				continue;
			}
			else if(heavyContig && getTwinCtg(connection->contigID) == heavyContig)
			{
				maxWeight = connection->weight;
				pre_contigID=getTwinCtg(connection->contigID);
				break;
			}
			//mao 2011 10 18
			else if(heavyContig && connection->weight + score[getTwinCtg(connection->contigID)] >=  score[contigID])
			{
				maxWeight = connection->weight;
				pre_contigID=getTwinCtg(connection->contigID);	
				break;
			}
			else if(connection->weight > maxWeight)
			{
				maxWeight = connection->weight;
				pre_contigID=getTwinCtg(connection->contigID);
			}
		}
		connection=connection->next;
	}
	if(pre_contigID)
	{
		setOneContigFlag(pre_contigID,1);
		getBestWay(pre_contigID,count);
	}
}
void findHeavyUsedContig(LOCUS * locus)
{
	unsigned int curr_contigID;
	unsigned int maxCovContig=0;
	double maxCov=0;
	unsigned int index=0;
	
	for(;index<locus->contigCount;index++)
	{
		curr_contigID=locus->contigID[index];
		if(!flag[curr_contigID])
		{
			if(contig_array[curr_contigID].cvg>maxCov)
			{
				maxCov=contig_array[curr_contigID].cvg;
				maxCovContig=curr_contigID;
			}
		}
	}
	heavyContig=maxCovContig;
}
typedef struct allTranscript
{
	unsigned int * contigID ;
	int contigCount;
	int length;
	struct allTranscript * next;
}ALLT;
static  int getIndex(unsigned int contigID,LOCUS * locus)
{
	int contigCount = locus->contigCount;
	int index;
	for(index=0;index<contigCount;index++)
	{
		if(locus->contigID[index] == contigID)
			return index;
	}
	return -1;
}
static ALLT * ALLT_ckalloc(ALLT * old, unsigned int newContig,int len,int max)
{
	ALLT * temp=NULL,*newPath=NULL;	
	if(old == NULL)
	{	
		newPath = (ALLT*) ckalloc(sizeof(ALLT));
		newPath->next = NULL;
		newPath->contigCount=1;
		newPath->contigID=(unsigned int *)ckalloc(sizeof(unsigned int )*2);
		newPath->contigID[0]=newContig;
		newPath->length = contig_array[newContig].length;
	}
	else
	{
		while(old)
		{
			if(old->contigCount + 1 < max)
			{
				newPath = (ALLT*) ckalloc(sizeof(ALLT));
				newPath->next = temp;
				newPath->contigCount=1+old->contigCount;
				newPath->contigID=(unsigned int *)ckalloc(sizeof(unsigned int )*(newPath->contigCount+1));
				int index;
				for(index=0;index<newPath->contigCount-1;index++)
				{
					newPath->contigID[index]=old->contigID[index];					
				}
				newPath->contigID[index]=newContig;
				newPath->length += contig_array[newContig].length + len;

				temp = newPath;
			}
			old = old->next;
		}

	}
	return newPath;
}
static void ALLT_free(ALLT * one)
{
	if(one==NULL)
		return;
	ALLT * temp;
	while(one)
	{
		temp=one->next;
		one->next=NULL;
		free(one->contigID);
		free(one);
		one=temp;
	}
	one=NULL;
}
static ALLT * getAllPath(LOCUS * locus,ALLT * curr_path,ALLT **allPath,int max)
{
	if(curr_path == NULL)
		return NULL;
	unsigned int lastContig = curr_path->contigID[curr_path->contigCount-1];
	CONNECT * cnt = contig_array[lastContig].downwardConnect;
	int change=0;
	int contig_index;
	ALLT * newPath;
	while(cnt)
	{
		if((cnt->deleted) || contig_array[cnt->contigID].unique != 1)
		{
			cnt = cnt->next;
			continue;
		}
		change++;
		newPath = ALLT_ckalloc(curr_path,cnt->contigID,cnt->gapLen,max);		
		
		newPath = getAllPath(locus,newPath,allPath,max);
		if(newPath)
		{
			contig_index = getIndex(cnt->contigID,  locus);
			if(contig_index <0)
			{
				printf("getIndex:\tthere is no %d\n",cnt->contigID);
				exit(0);
			}
			newPath->next=allPath[contig_index];
			allPath[contig_index]=newPath;
		}		
		cnt = cnt->next;			
	}
	if(change>0)
	{
		ALLT_free(curr_path);	
		curr_path=NULL;
	}
	return curr_path;
	
}
void outputPath(FILE * fq,FILE * fo,ALLT* path)
{
	if(path->contigCount==1)
		return;
	unsigned int prev_ctg,index , max_steps=5,currentContigID;	
	int num_route,gap_c,j;
//	boolean weak;	
	short gapLen=0;		

	fprintf(fq,">scaffold%d Locus_%d_%d %d %d %s\n",scaffIndex,loci_id,loci_count,path->contigCount,path->length,curr_type);
	fprintf(fo,">scaffold%d Locus_%d_%d %d %d %s\n",scaffIndex,loci_id,loci_count,path->contigCount,path->length,curr_type);

	scaffIndex++;
	loci_count++;
	int len=0;
	int gaplen=0;
	CONNECT * cnt;
	for(index=0;index<path->contigCount;index++)
	{
		currentContigID = path->contigID[index];	
		if(!isLargerThanTwin(currentContigID))
		{			
			fprintf(fq,"%-10d %-10d +   %d \n",
				index_array[currentContigID],
				len,
				contig_array[currentContigID].length+overlaplen);
//			weak = printCnts(fq,currentContigID);
			printf("%-10d %-10d +   %d \n",
                                index_array[currentContigID],
                                len,
                                contig_array[currentContigID].length+overlaplen);
		}
		else
		{
			printf("%-10d %-10d -   %d \n",
				index_array[getTwinCtg(currentContigID)],
				len,	
				contig_array[currentContigID].length+overlaplen);
//			 weak = printCnts(fq,currentContigID);
		}
		if(index < path->contigCount -1)
		{
			cnt = getCntBetween(currentContigID,path->contigID[index+1]);
			gaplen = cnt->gapLen>0 ? cnt->gapLen : 0;
			len += contig_array[currentContigID].length + gaplen;
		}
	}
	
}
static void  allPath(LOCUS * locus ,FILE * fq,FILE * fo,unsigned int *leftContigID,int count)
{
	int contigCount = locus->contigCount;
	ALLT ** allpath = (ALLT**) ckalloc(sizeof(ALLT*)*contigCount);
	int index;
	for(index=0;index<contigCount;index++)
	{
		allpath[index]=NULL;
	}
	int curr_index ;
	for(index=0;index<count;index++)
	{
		curr_index= getIndex(leftContigID[index],  locus);
		allpath[curr_index]=ALLT_ckalloc(NULL, leftContigID[index],0,contigCount);
		allpath[curr_index] = getAllPath(locus,allpath[curr_index],allpath,contigCount);		
	}

	for(index=0;index<contigCount;index++)
	{
		if(allpath[index] == NULL)
			continue;
		ALLT * temp = allpath[index];
		while(temp)
		{
			outputPath(fq,fo,temp);
			temp = temp->next;
		}	
		ALLT_free(allpath[index]);

	}
	free(allpath);

}
void outputOneComplexTranscriptome(LOCUS * locus ,FILE * fq,FILE * fo)
{	
	unsigned int index;
	unsigned int contigCount = locus->contigCount;

	unsigned int *leftContig = (unsigned int *)malloc(contigCount * sizeof(unsigned int));
	int leftContigCount=0;	
	unsigned int contigID,twin,maxScoreContigID,nextContigID;
	CONNECT * connection;
	unsigned int out,count=0;
	int len=0;	
	unsigned int  totalScore = 0;

	for(index = 0; index < locus->contigCount ; index++)
	{
		twin = getTwinCtg(locus->contigID[index]);
		connection = contig_array[twin].downwardConnect;
		out=0;
		while(connection)
		{
			if(!connection->deleted && contig_array[connection->contigID].unique)
			{
				out++;
			}
			connection = connection->next;
		}
		if(!out)
		{
			leftContig[leftContigCount++]= locus->contigID[index];
		}
	}	

	unsigned int current_contig;
	maxScoreContigID=0;
	//mao 2011-10-14
	int * score_time = (int *) ckalloc (sizeof(int)*(num_ctg+1));
	for(index = 0; index < locus->contigCount ; index++)
	{
		contigID = locus->contigID[index];
		score_time [contigID]=0;
		score[contigID]=0;
	}
	for(index=0;index < leftContigCount;index++)
	{
		InserHeap(leftContig[index]);
		while(contigID= outputHeap())
		{
			computeScore(contigID);
			if(score[contigID] > totalScore)
			{
				totalScore = score[contigID];
				maxScoreContigID=contigID;
			}

			if(score_time[contigID] <1000)
			{
				propagate(contigID,locus);
				score_time[contigID]++;
			}
		}
	}
	free(leftContig);
	free(score_time);

        for(index = 0; index < locus->contigCount ; index++)
        {
        	setOneContigFlag(locus->contigID[index],0);
        }
//	setOneContigFlag(maxScoreContigID,1);	
	getBestWay(maxScoreContigID,&count);	
	int total_weight=0;
	if(count>=2)	
	{
		index=0;
		contigID=outputHeap();			
		len=contig_array[contigID].length;
		*(unsigned int*)darrayPut(scaf,index++) = contigID;
//		setOneContigFlag(contigID, 1);
		flag[contigID]=1;
		flag[getTwinCtg(contigID)]=1;
		while(nextContigID=outputHeap())
		{
			connection=getCntBetween(contigID,nextContigID);
			if(connection)
			{
				total_weight+=connection->weight;
				*(int *)darrayPut(gap,index-1) = connection->gapLen;
				*(unsigned int*)darrayPut(scaf,index++) = nextContigID;
				len += contig_array[nextContigID].length + connection->gapLen;
//				setOneContigFlag(nextContigID,1);
				flag[nextContigID]=1;
				flag[getTwinCtg(nextContigID)]=1;
			}
			else
				printf("error:it doesn't exist %d->%d\n",contigID,nextContigID);
			contigID=nextContigID;
		}
//		printf(">scaffold%d\t%d\n",scaffIndex,total_weight);
		outputOneTranscriptome( fq,  fo, index, len + overlaplen);
	}
	findHeavyUsedContig(locus);	

}
unsigned int getUsedContigCount(LOCUS * locus)
{
	unsigned int count=0,index=0;
	for(;index<locus->contigCount;index++)
	{
		if(flag[locus->contigID[index]] == 1)
			count++;
	}
	return count;	
}
static int num_light=0;
static void deleteLightCnt(LOCUS * locus)
{
        int index;
        CONNECT * cnt;
        int sum;
        int count;

	for(index=0;index<locus->contigCount;index++)
	{
		unsigned int contigID = locus->contigID[index];
		sum=0;
		count=0;
		cnt = contig_array[contigID].downwardConnect;
		while(cnt)
		{
			if(!(cnt->deleted || cnt->weak))
			{
				sum += cnt->weight;
				count++;
			}
			cnt = cnt->next;
		}

		if(sum == 0 || count <=1)
			continue;
		cnt = contig_array[contigID].downwardConnect;
		while(cnt)
		{
			if(!(cnt->deleted || cnt->weak))
			{
				int weight = (int)(cnt ->weight);
				printf("test_connection:\t%4.4f\t%d\t%d\n",(double)weight/sum,weight,sum);
				if(weight < (double)sum * 0.05)
				{
					cnt->deleted =1;
					cnt->weak=1;
					
					CONNECT * temp = getCntBetween(getTwinCtg(cnt->contigID),getTwinCtg(contigID));
					temp->deleted=1;
					temp->weak=1;
					
					num_light++;
				}
			}
			cnt = cnt->next;
		}
	}
}
static unsigned int getSubLociCount(LOCUS * locus)
{
        unsigned int lociCount=0;
        unsigned int i;

        setContigFlag(0);

        for(i=0;i<locus->contigCount;i++)
        {		
        	unsigned int contigID= locus->contigID[i];
                if(!contig_array[contigID].flag && contig_array[contigID].unique)
                {
                        lociCount++;
                        propagateComponent(contigID);
                }
        }
	
	return lociCount;
}
static LOCUS * getSubLoci(LOCUS * locus,unsigned int count)
{
	unsigned int i;
	unsigned int contigCount;
	unsigned int lociIndex=0;
	setContigFlag(0);

	LOCUS * loci = (LOCUS * ) malloc (count * sizeof(LOCUS));
	for(i=0;i<locus->contigCount;i++)
       {		
        	unsigned int contigID= locus->contigID[i];
		contigCount=0;
		if(!contig_array[contigID].flag && contig_array[contigID].unique)
		{
			fillUpComponent(contigID,&contigCount,1);	
			loci[lociIndex].contigCount= contigCount;
			loci[lociIndex].contigID = 
				(unsigned int *)malloc((contigCount) * sizeof(unsigned int));
			for(contigCount=0;contigCount<loci[lociIndex].contigCount;contigCount++)
			{
//				printf("",contigCount,contigCount<=loci[lociIndex].contigCount);
				loci[lociIndex].contigID[contigCount] = 
					*(unsigned int *) darrayGet(scaf,contigCount);
			}
			lociIndex++;
		}
	}

	return loci;
}
void outputComplexTranscriptome(LOCUS * locus ,FILE * fq,FILE * fo)
{
	int counter=0;
	heavyContig=0;
	unsigned int index=0;
	if(locus->contigCount >max_loci)
		max_loci=locus->contigCount;
	for(;index<locus->contigCount;index++)
	{
		flag[locus->contigID[index]]=0;
		flag[getTwinCtg(locus->contigID[index])]=0;
	}
	while(counter++ < max_num 
		&&  getUsedContigCount(locus)< locus->contigCount)
		outputOneComplexTranscriptome(locus, fq,  fo);

	if(debug)
		printf("COMPLEX : %d\n",counter);
}
void outputTranscriptome(LOCUS * locus,unsigned int kinds ,FILE * fq,FILE * fo)
{		
	if(kinds == LINEAR)
	{
		curr_type="LINEAR";
		outputLinearTranscriptome(locus,fq,fo);
	}
	else if(kinds == FORK)
	{
		curr_type="FORK";
		outputForkTranscriptome(locus,fq,fo);
	}
	else if(kinds == BUBBLE)
	{
		curr_type="BUBBLE";
		outputBubbleTranscriptome(locus,fq,fo);
	}
	else if(kinds == COMPLEX){ 
		curr_type="COMPLEX";
		outputComplexTranscriptome(locus,fq,fo);
	}
}

void transcript(unsigned int lociCount,LOCUS * loci,FILE * fq,FILE * fo)
{
	unsigned int i;
	LOCUS * locus;
	unsigned int kinds=0;
	setContigFlag(0);
	unsigned int kind[4]={0,0,0,0};
	for(i = 0 ; i <lociCount; i++)
	{	
		if(debug)
			printf("loci:\t%d\n",i);
		locus = &(loci[i]);
		if(locus->contigCount == 1)
			continue;
		kinds=getLocusKind(locus);
		kind[kinds-1]++;
		if(debug)
			printf("kind:%d\tcontigCount:%d\n",kinds,locus->contigCount);
		loci_id=i;
		loci_count=0;
		outputTranscriptome(locus,kinds,fq,fo);		
	}
	printf("There are %d locis, which contain more than 2 contigs\n",kind[0]+kind[1]+kind[2]+kind[3]);
	printf("The loci can be classified to four kinds:\n");
	for(i=0;i<4;i++)
	{
		if(i==0)
			printf("LINEAR:%d\n",kind[i]);
		else if(i==1)
			printf("FORK:%d\n",kind[i]);
		else if(i==2)
			printf("BUBBLE:%d\n",kind[i]);
		else if(i==3)
			printf("COMPLEX:%d\n",kind[i]);
	}
//	printf("delete light connection %d\n",num_light);
}
static void removeUnnecessaryConnection(unsigned int contigID,int cut_off)
{
	if(cut_off ==0 ||cut_off>10 )
		return;
	int cov[10];
	int i,j;
	for(i=0;i<10;i++)
	{
		cov[i]=0;
	}
	CONNECT * connection = contig_array[contigID].downwardConnect;
	while(connection)
	{
		if(!connection->deleted && contig_array[connection->contigID].unique)
		{
			for(i=0;i<10;i++)
			{
				if(connection->weight > cov[i])
				{
					for(j=9;j>i;j--)
					{
						cov[j]=cov[i];					
					}
					cov[i]=connection->weight;
				}
			}
		}
		connection=connection->next;
	}
	connection = contig_array[contigID].downwardConnect;
	while(connection)
	{
		if(!connection->deleted && contig_array[connection->contigID].unique)
		{
			if(connection->weight < cov[cut_off-1])
			{
				connection->deleted=1;
				connection->weak=1;

				CONNECT * cnt = getCntBetween(getTwinCtg(connection->contigID),getTwinCtg(contigID));
				cnt->deleted=1;
				cnt->weak=1;			
			}
		}
		connection=connection->next;
	}
}
void deleteUnlikelyCnt(int cut_off)
{
	if(cut_off == 0 || cut_off >10)
		return ;
	unsigned int index;
	for(index=1;index<=num_ctg;index++)
	{		
		if(contig_array[index].unique)
			continue;
		CONNECT * cnt= contig_array[index].downwardConnect;
		int output=0;
		while(cnt)
		{
			if(!cnt->deleted && contig_array[cnt->contigID].unique)
				output++;
			cnt=cnt->next;				
		}
		if(output>cut_off)
			removeUnnecessaryConnection(index,cut_off);
	}
}
void transcriptome(char * outfile)
{
	unsigned int lociCount=0;
	FILE * flc;	
	LOCUS * loci;
	unsigned int i,j;
	char  name[256];
	CONNECT *connection,*bal_connection;

	setUniqueContig(ctg_mask);
	singleRead2connection(outfile);	
//	ArcToConnect2();
//	updatedGapLength();	
	deleteWeakCnt(3);


	lociCount = getLociCount();
//	printf("Divided into %d loci by the connections\n",lociCount);	
	scaf = (DARRAY * )createDarray(1000,sizeof(unsigned int));
	gap = (DARRAY *) createDarray(1000,sizeof(unsigned int));

	loci = getLoci(lociCount);	
	printf("Begin to linearize the graph\n");
	linearization(loci,lociCount);	

	printf("Remove the error connections\n");
	deleteInconsistent(lociCount,loci);
	avoidLoop(lociCount,loci);

	printf("The second time to linearize the graph \n");
        linearization(loci,lociCount);
//	output_loci_graph(outfile,lociCount,loci);	
// output transciptome
	deleteUnlikelyCnt(max_cnt);	

	
	for(i=0;i<lociCount;i++)
	{
		free(loci[i].contigID);
	}
	free(loci);		
	lociCount = getLociCount();
	loci = getLoci(lociCount);
	printf("Divided into %d loci by the connections\n\n",lociCount);	
	
	sprintf(name,"%s.scaf",outfile);
	FILE * fq = ckopen (name,"w");
	sprintf(name,"%s.scaf_gap",outfile);	
	FILE * fo = ckopen (name,"w");

	score = (unsigned int *)malloc((num_ctg+1)*sizeof(unsigned int ));
	flag=(int *)malloc((num_ctg+1)*sizeof( int ));

	if(debug) 	//	some information are output and useful to debug
	{
		sprintf(name,"%s.lociInformation",outfile);	
		flc = ckopen (name,"w");
		for(i = 0 ; i < lociCount ; i++)
		{
			fprintf(flc,">loci[%d]\t%d\n",i,loci[i].contigCount);
			for(j=0;j<loci[i].contigCount;j++)
			{
				connection=contig_array[loci[i].contigID[j]].downwardConnect;			
				fprintf(flc,"%d\tlength:%d\n",loci[i].contigID[j],contig_array[loci[i].contigID[j]].length);
				while(connection)
				{
					if(!(connection->deleted) && contig_array[connection->contigID].unique == 1)
					{
						fprintf(flc,"%d->%d\tgap:%d\tweight:%d\tSE:%d\tPE:%d\n",loci[i].contigID[j],connection->contigID,connection->gapLen,connection->weight,connection->SECount,connection->PECount);	
					}
					connection=connection->next;
				}
			}
		}
		fclose(flc);
	}

	if(debug)		//	some information saved and are useful to debug
	{
		sprintf(name,"%s.CntGap",outfile);	
		flc = ckopen (name,"w");
		for(i=1;i<=num_ctg;i++)
		{
			if(!contig_array[i].unique)
				continue;
			CONNECT * cnt  = contig_array[i].downwardConnect;
			while(cnt)
			{
				if(!(cnt->deleted) && contig_array[cnt->contigID].unique == 1)
					fprintf(flc,"%d %d -> %d\n",cnt->gapLen,i,cnt->contigID);
				cnt=cnt->next;
			}
		}
		fclose(flc);
	}

	curr_type=(char *)ckalloc(sizeof(char)*100);
	transcript(lociCount,loci,fq,fo);
	
	if(debug==1)
	{
		for(i=0;i<256;i++)
		{
			if(num_connection[i] >0)
				printf("connection%d:\t%d\n",i,num_connection[i]);
		}
	}
	if(debug)
		printf("the number of the max loci:\t%d\n",max_loci);
	free(score);

//free
	fclose(fq);
	fclose(fo);
	freeDarray(scaf);
	freeDarray(gap);
	free(flag);
	for(i=0;i<lociCount;i++)
	{
		free(loci[i].contigID);
	}
	free(loci);
}
