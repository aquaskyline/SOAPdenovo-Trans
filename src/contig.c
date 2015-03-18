/*
 * Contig.c
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

static void initenv (int argc, char **argv);
static void display_contig_usage ();
char shortrdsfile[256], graphfile[256];
static boolean repeatSolve;
static int M = 1;

FILE * edge_info=NULL;
int file_num=0;
void output_edge_info(unsigned int edgeid)
{
	unsigned int in_cov=0,out_cov=0;
	ARC * arc;
	unsigned int twin = getTwinEdge(edgeid);
	
	arc=edge_array[edgeid].arcs;
	int in_time=0,out_time=0;
	while(arc)
	{
		out_time++;
		out_cov += edge_array[arc->to_ed].cvg;
		arc = arc->next;
	}

	arc=edge_array[twin].arcs;
	while(arc)
	{
		in_time++;
		in_cov += edge_array[arc->to_ed].cvg;
		arc = arc->next;
	}

	fprintf(edge_info,"%d\t%d\t%d\t%d\t%d\n",in_cov,in_time,edge_array[edgeid].cvg,out_cov,out_time);
	
}
void find_edge(unsigned int edgeid,int poolid,double avg_cov)
{
	if(pool[edgeid] !=0)// || (edge_array[edgeid].length<=48 && edge_array[edgeid].cvg <=avg_cov))
		return ;
	if(edge_array[edgeid].length<=48 && edge_array[edgeid].cvg <=avg_cov)
	{
		if(edge_info != NULL)
			output_edge_info(edgeid);
		return ;
	}
	pool[edgeid]=poolid;
	unsigned int twin_ed=getTwinEdge(edgeid);
	pool[twin_ed]=poolid;
	
	ARC * arc;
	arc=edge_array[edgeid].arcs;
	while(arc)
	{
		find_edge(arc->to_ed,poolid,avg_cov);
		arc=arc->next;
	}

	if(twin_ed == edgeid)
		return ;
	
	arc=edge_array[twin_ed].arcs;
	while(arc)
	{
		find_edge(arc->to_ed,poolid,avg_cov);
		arc=arc->next;
	}
}

void resetZeroEdge()
{
	unsigned int index;
	for(index=1;index<=num_ed;index++)
	{
		if(edge_array[index].cvg > 0)
			continue;
		unsigned short min_cov=0;
		ARC * arc;

		arc = edge_array[index].arcs;
		while(arc)
		{
			if(min_cov > edge_array[arc->to_ed].cvg && edge_array[arc->to_ed].cvg >0)
				min_cov=edge_array[arc->to_ed].cvg;
			if(min_cov == 0)
				min_cov=edge_array[arc->to_ed].cvg;
			arc=arc->next;				
		}

		arc = edge_array[getTwinEdge(index)].arcs;
		while(arc)
		{
			if(min_cov > edge_array[arc->to_ed].cvg && edge_array[arc->to_ed].cvg>0)
				min_cov=edge_array[arc->to_ed].cvg;
			if(min_cov == 0)
				min_cov=edge_array[arc->to_ed].cvg;
			arc=arc->next;				
		}
		edge_array[index].cvg=min_cov;		
		edge_array[getTwinEdge(index)].cvg=min_cov;
		if(index == getTwinEdge(index))
			index++;		
	}
}
void divideComponent(char * outfile)
{
	pool=(int*)malloc(sizeof(int )*(num_ed+1));	
	unsigned int index,twin_ed;
	int pool_id=1;

	for(index=1;index<=num_ed;index++)
	{
		pool[index]=0;
	}
	for(index=1;index<=num_ed;index++)
	{
		if(pool[index] != 0 )
			continue;
		
		twin_ed=getTwinEdge(index);

		find_edge(index,pool_id,0.0);
		if(index != twin_ed)
		{
			find_edge(twin_ed,pool_id,0.0);
			index++;
		}		
		pool_id++;
	}	
	int pool_count=pool_id-1;
	int *contigCount=(int*)malloc((1+pool_count)*sizeof(int));	
	int *length=( int* ) malloc ((1+pool_count)*sizeof( int ));
	double * avg_cov=( double* ) malloc ((1+pool_count)*sizeof( double ));
	int  * pool_backup =  (int*)malloc(sizeof(int )*(num_ed+1));	
	for(pool_id=1;pool_id<=pool_count;pool_id++)
	{
		length[pool_id]=0;
		contigCount[pool_id]=0;
		avg_cov[pool_id]=0.0;
	}
	unsigned int twin;
	for(index=1;index<=num_ed;index++)
	{
		twin=getTwinEdge(index);		
		contigCount[pool[index]]++;
		length[pool[index]]+=edge_array[index].length;
		avg_cov[pool[index]]+=edge_array[index].cvg*edge_array[index].length;
		if(twin!=index)
			index++;
	}
	
	for(pool_id=1;pool_id<=pool_count;pool_id++)
	{
		avg_cov[pool_id] /= length[pool_id];		
	}
	for(index=1;index<=num_ed;index++)
	{
		pool_backup[index]=pool[index];
		pool[index]=0;
	}
	char  fileName[1024];
	sprintf (fileName, "%s.edgeInfo_%d", graphfile,file_num);	
	edge_info = fopen(fileName,"w");
	fprintf(edge_info,"in\tin_time\tself\tout\tout_time\n");
	pool_id=1;
	for(index=1;index<=num_ed;index++)
	{
		if(pool[index] != 0 )//|| (edge_array[index].length<=48 && edge_array[index].cvg <= avg_cov[pool_backup[index]]))
			continue;
		if(edge_array[index].length<=48 && edge_array[index].cvg <= avg_cov[pool_backup[index]])
		{
			if(edge_info != NULL)
				output_edge_info(index);
			continue;
		}
		
		twin_ed=getTwinEdge(index);

		find_edge(index,pool_id,avg_cov[pool_backup[index]]);
		if(index != twin_ed)
		{
			find_edge(twin_ed,pool_id,avg_cov[pool_backup[index]]);
			index++;
		}		
		pool_id++;
	}	
	fclose(edge_info);
	free(avg_cov);
	free(length);
	free(contigCount);
	free(pool_backup);

	sprintf(fileName,"%s.poolid",graphfile);
	
	output_pool(fileName);

	free(pool);	
	file_num++;
}
int call_heavygraph (int argc, char **argv) 
{
	time_t start_t, stop_t, time_bef, time_aft;
	time (&start_t);
	boolean ret;
	initenv (argc, argv);
	loadVertex (graphfile);
	loadEdge (graphfile);

	char name[256];

//	swapedge();
//	sortedge();
//	freshArc();
	
	//0531
	if (M > 0)
	{
		time (&time_bef);
		bubblePinch (0.90, graphfile, M);
		time (&time_aft);
		printf ("time spent on bubblePinch: %ds\n\n", (int) (time_aft - time_bef));
	}
	deleteWeakEdge(de);
	
	cutTipsInGraph (0, 0);

//	resetCov();
	deleteUnlikeArc();
	delowHighArc(delowArc);
	int change=1;
	int changed_time=1;
	while(change){	
		printf("%d time to delete light arcs \n",changed_time++);
		deleteSimpleLoop();
		change=deleteLightArc();
		if(change)
		{
			linearConcatenate();
			compactEdgeArray ();
		}
		printf("\n");
	}		
	//output_graph(graphfile);
//	deleteLightContig();
	deleteShortContig(cut_length);
	linearConcatenate();
	compactEdgeArray ();
	printf("\n");
	
	output_contig (edge_array, num_ed, graphfile, overlaplen + 1);

//	sprintf(name,"%s.AllComponent_back",graphfile);
//	divideComponent(name);
	
	output_updated_edges (graphfile);
	output_heavyArcs (graphfile);
	if (vt_array)
	{
		free ((void *) vt_array);
		vt_array = NULL;
	}
	if (edge_array)
	{
		free_edge_array (edge_array, num_ed_limit);
		edge_array = NULL;
	}
	destroyArcMem ();
	time (&stop_t);
	printf ("time elapsed: %dm\n\n", (int) (stop_t - start_t) / 60);
	return 0;
}


/*****************************************************************************
 * Parse command line switches
 *****************************************************************************/ 
void initenv (int argc, char **argv) 
{
	int copt;
	int inpseq, outseq;
	extern char *optarg;
	char temp[100];

	inpseq = outseq = repeatSolve = 0;
	optind = 1;
	while ((copt = getopt (argc, argv, "g:e:q:Q:H:S:M:")) != EOF)	
	{
		switch (copt)	
		{
		
			case 'M':
				sscanf (optarg, "%s", temp);	//
				M = atoi (temp);
				break;
	/*		case 'D':
				sscanf (optarg, "%s", temp);
				de= atoi (temp) >= 0 ? atoi (temp) : 0;
				de *=10;
				break;*/
			case 'g':
				inGraph = 1;
				sscanf (optarg, "%s", graphfile);	//
				break;
/*			case 'R':
				repeatSolve = 1;
				break;				
*/
			case 'S':
				sscanf (optarg, "%s", temp);
				cut_length= atoi (temp) >= 0 ? atoi (temp) : 0;
				break;
			case 'e':
				sscanf (optarg, "%s", temp);
				de= atoi (temp) >= 0 ? atoi (temp) : 0;
				de *=10;
				break;
			case 'q':
				sscanf (optarg, "%s", temp);
				da= atoi (temp) >= 0 ? atoi (temp) : 0;
				break;
			case 'Q':
				sscanf (optarg, "%s", temp);
				dA= atoi (temp) >= 0 ? atoi (temp) : 0;
				break;
			case 'H':
				sscanf (optarg, "%s", temp);
				delowArc= atoi (temp) >= 100 ? atoi (temp) : 200;
				break;
			default:
			if (inGraph == 0)	//
			{
				display_contig_usage ();
				exit (-1);
			}
		}
	}
	if (inGraph == 0)	//
	{
		display_contig_usage ();
		exit (-1);
	}
}
static void display_contig_usage () 
{
	printf ("\ncontig -g inputGraph [-e EdgeCovCutoff -M mergeLevel\n");
	printf ("  -g\t<string>\tinputGraph: prefix of input graph file name\n");
	printf ("  -e\t<int>\t\tEdgeCovCutoff: edges with coverage no larger than EdgeCovCutoff will be deleted, [2]\n");	
	printf ("  -M\t<int>\t\tmergeLevel(min 0, max 3): the strength of merging similar sequences during contiging, [1]\n");	
} 
