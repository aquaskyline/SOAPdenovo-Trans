/*
 * loadReadPath.c
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
#include "newhash.h"
#include "extfunc.h"
#include "extvab.h"
#include "zlib.h"
#include <string.h>


static int split(unsigned int *arr, char *str,int contig_count)//, const char *del) 
{ 
	char *s =NULL; 
	char *del = "\t";
	int index=0;
	unsigned int contigID;
	unsigned int pre_contigID;
	s=strtok(str,del); 
	while(s != NULL) 
	{ 	
		sscanf (s ,"%d", &contigID);
		if(contigID== pre_contigID && index == contig_count)
			break;
		arr[index++] = contigID; 
		pre_contigID=contigID;
		s = strtok(NULL,del); 
	} 
/*
	int i;
	for(i=0;i<index;i++)
	{
		printf("%d\t",arr[i]);	
	}
	printf("\nindex:\t%d\n",index);
*/
	return index;
} 
void loadReadPath (char * graphfile)
{
	char name[256],line[10240];
	unsigned int index;
	FILE *fp;
	contig_path_array = (CONTIG_PATHID *)ckalloc((num_ctg+1) * sizeof(CONTIG_PATHID));
	for(index=0;index<=num_ctg;index++)
	{
		contig_path_array[index].path_count=0;
	}

	sprintf(name,"%s.contigOnpath",graphfile);	
	fp=ckopen( name, "r");
	unsigned int contigID,pathID;
	while (fgets (line, sizeof (line), fp) != NULL)
	{
		sscanf (line ,"%d %d", &contigID, &pathID);
		contig_path_array[contigID].pathID = (unsigned int *)ckalloc (pathID* sizeof(unsigned int ));
	}
	fclose(fp);

	sprintf(name,"%s.read_path",graphfile);	
	fp=ckopen( name, "r");
	unsigned int path_count;
	fgets (line, sizeof (line), fp);
	sscanf (line ,"%d",&path_count);
	path_contig_array = (PATH_CONTIGID *)ckalloc (path_count * sizeof(PATH_CONTIGID));
	
	int contig_count,cov;
	unsigned int path_index=0;
	int i;
	while (fgets (line, sizeof (line), fp) != NULL)
	{
		if(line[0] == '>')
		{
			sscanf (line ,"%s %d %d ", &name, &contig_count,&cov);
//			printf("check:\t%s\t%d\t%d\n",name,contig_count,path_count);
			path_contig_array[path_index].contigID = (unsigned int *)ckalloc (contig_count* sizeof(unsigned int ));
			path_contig_array[path_index].contig_count = contig_count;
			path_contig_array[path_index].coverage=cov;
		}
		else
		{
			contig_count = split(&(path_contig_array[path_index].contigID[0]) , &(line[0]),path_contig_array[path_index].contig_count);
			if(contig_count >path_contig_array[path_index].contig_count)
			{				
				printf("ERROR:\tinput %d\tfault:\t%d\n",contig_count,path_contig_array[path_index].contig_count);
				exit(0);
			}
			else
			{
				for(i=0;i<path_contig_array[path_index].contig_count;i++)
				{
					unsigned int contigID = path_contig_array[path_index].contigID[i];
					contig_path_array[contigID].pathID[contig_path_array[contigID].path_count++]=path_index;
				}
			}
			path_index++;			
		}
	}
	fclose(fp);	
}
