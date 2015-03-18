/*
 * ReadTrace.c
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


typedef struct readOnContig
{
	unsigned long long  readID;
	int read_pos,contig_pos;
	int alignLength;
	char orign;
	struct readOnContig * next;
}READONCONTIG;
void getReadOnScaf(char *file)
{
	char name[1024];
	char line[1024];
	FILE * fp1=NULL,*fp2=NULL,*output=NULL;

	sprintf(name,"%s.readInformation",file);
	fp1=ckopen(name, "r");
	sprintf(name,"%s.contigPosInscaff",file);
	fp2=ckopen(name, "r");
	sprintf(name,"%s.readOnScaf",file);
	output=ckopen(name,"w");

	unsigned  long long   readID;
	unsigned int contigID;
	int read_pos,contig_pos;
	int alignLength;
	char orign;

	READONCONTIG ** contig2read = (READONCONTIG **) ckalloc ( sizeof(READONCONTIG *) * (num_ctg + 1));
	READONCONTIG *temp=NULL;
	unsigned int index;
	for(index=0;index<=num_ctg;index++)
	{
		contig2read[index]=NULL;
		contig_array[index].flag=0;
	}
	while(fgets(line,sizeof(line),fp1))
	{
		sscanf(line,"%llu %d %d %d %d %c",&readID,&read_pos,&contigID,&contig_pos,&alignLength,&orign);
		temp=(READONCONTIG *) ckalloc ( sizeof(READONCONTIG ));
		temp->next=NULL;
		temp->readID=readID;
		temp->read_pos=read_pos;
		temp->contig_pos=contig_pos;
		temp->alignLength=alignLength;
		temp->orign=orign;

		temp->next = contig2read[contigID];
		contig2read[contigID] = temp;		
	}
	int isFirst=0;
	int scafPos,scafAlignLength;
	char scafOrign;
	while(fgets(line,sizeof(line),fp2))
	{
		if(line[0]=='>')
		{
			fprintf(output,"%s",line);
			isFirst=1;
		}
		else
		{
			sscanf(line ,"%d %d %c %d",&contigID,&contig_pos,&orign,&alignLength);

			temp = contig2read[contigID] ;
			contig_array[contigID].flag=1;
			contig_array[getTwinCtg(contigID)].flag=1;
			while(temp)
			{
				if(isFirst)
				{
					scafPos=contig_pos + temp->contig_pos;
					scafAlignLength = temp->alignLength;
				}
				else
				{
					scafPos=contig_pos + temp->contig_pos - overlaplen;
					if( temp->contig_pos < overlaplen)
						scafAlignLength = temp->alignLength - overlaplen + temp->contig_pos;
					else
						scafAlignLength = temp->alignLength;
				}
				
				if(orign ==  temp->orign)
				{
					scafOrign = '+';
				}
				else
				{
					scafOrign = '-';
				}
				fprintf(output , "%llu\t%d\t%d\t%c\t%d\n" , temp->readID,temp->read_pos,scafPos,scafOrign,scafAlignLength);
				temp = temp->next;
			}
			isFirst=0;
		}
	}
	for(index=0;index<=num_ctg;index++)
	{
		if ((contig_array[index].length + overlaplen) < 100 || contig_array[index].flag)
		{
			continue;
		}
		fprintf(output,">C%d\n", index);
		contig_array[index].flag=1;
		contig_array[getTwinCtg(index)].flag=1;
		temp = contig2read[index] ;
		while(temp)
		{
			fprintf(output , "%llu\t%d\t%d\t%c\t%d\n" , temp->readID,temp->read_pos,temp->contig_pos,temp->orign,temp->alignLength);
			temp=temp->next;
		}
	}
	fclose(fp1);
	fclose(fp2);
	fclose(output);

	for(index=0;index<=num_ctg;index++)
	{		
		while(contig2read[index])
		{
			temp=contig2read[index]->next;
			contig2read[index]->next = NULL;
			free(contig2read[index]);
			contig2read[index]=temp;
		}
	}
	free(contig2read);
}
