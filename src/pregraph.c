/*
 * Pregraph.c
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
static char shortrdsfile[256];
static char graphfile[256];
static int cutTips = 1;
static void display_pregraph_usage ();
int call_pregraph (int argc, char **argv) 
{
	time_t start_t, stop_t, time_bef, time_aft;
	time (&start_t);
	initenv (argc, argv);
	if (overlaplen % 2 == 0)	
	{
		overlaplen++;
		printf ("K should be an odd number\n");
	}
	if (overlaplen < 13)	
	{
		overlaplen = 13;
		printf ("K should not be less than 13\n");
	}
#ifdef MER127
	else if (overlaplen > 127)	
		overlaplen = 127;
#endif
#ifdef MER63
	else if(overlaplen > 63)
		overlaplen=63;
#endif
#ifdef MER31
	else if(overlaplen>31)
		overlaplen = 31;
#endif
	
	time (&time_bef);

	prlRead2HashTable (shortrdsfile, graphfile);
	time (&time_aft);
	printf ("time spent on pre-graph construction: %ds\n\n", (int) (time_aft - time_bef));
	printf ("deLowKmer %d, deLowEdge %d\n", deLowKmer, deLowEdge);
	
	time (&time_bef);
	removeMinorOut();
	time (&time_aft);
	printf ("time spent on cut kmer: %ds\n\n", (int) (time_aft - time_bef));
	
	//analyzeTips(hash_table, graphfile);

	if (!deLowKmer && cutTips)	
	{
		time (&time_bef);
		removeSingleTips ();
		removeMinorTips ();
		time (&time_aft);
		printf ("time spent on cutTipe: %ds\n\n", (int) (time_aft - time_bef));
	}
	else	
	{
		time (&time_bef);
		removeMinorTips ();
		time (&time_aft);
		printf ("time spent on cutTipe: %ds\n\n", (int) (time_aft - time_bef));
	}
 

	initKmerSetSize = 0;

		
	//combine each linear part to an edge
	time (&time_bef);
	kmer2edges (graphfile);
	time (&time_aft);
	printf ("time spent on making edges: %ds\n\n", (int) (time_aft - time_bef));
	
	//map read to edge one by one
	time (&time_bef);
	prlRead2edge (shortrdsfile, graphfile);
	time (&time_aft);
	printf ("time spent on mapping reads: %ds\n\n", (int) (time_aft - time_bef));
	output_vertex (graphfile);
	free_Sets (KmerSets, thrd_num);
	free_Sets (KmerSetsPatch, thrd_num);
	time (&stop_t);
	printf ("overall time for lightgraph: %dm\n\n", (int) (stop_t - start_t) / 60);
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

	optind = 1;
	inpseq = outseq = 0;
	while ((copt = getopt (argc, argv, "a:s:o:K:p:d:Di:n")) != EOF)	
	{
		
		//printf("get option\n");
		switch (copt)			
		{
			case 's':
				inpseq = 1;
				sscanf (optarg, "%s", shortrdsfile);
				break;
			case 'o':
				outseq = 1;
				sscanf (optarg, "%s", graphfile);	//
				break;
			case 'K':
				sscanf (optarg, "%s", temp);	//
				overlaplen = atoi (temp);
				break;
			case 'p':
				sscanf (optarg, "%s", temp);	//
				thrd_num = atoi (temp);
				break;
/*			case 'R':
				repsTie = 1;
				break;*/
			case 'd':
				sscanf (optarg, "%s", temp);	//
				deLowKmer = atoi (temp) >= 0 ? atoi (temp) : 0;
				break;
/*			case 'D':
				deLowEdge = 1;
				break;
*/
			case 'a':
				initKmerSetSize = atoi (optarg);
				break;
			case 'i'://mao 9-21
				sscanf (optarg, "%s", temp);
				dd=atoi (temp) >= 0 ? atoi (temp) : 0;
				break;
			case 'n'://mao 2011-10-13
				N_kmer=1;
				break;
			default:
			if (inpseq == 0 || outseq == 0)	//
			{
				display_pregraph_usage ();
				exit (-1);
			}
		}
	}
	if (inpseq == 0 || outseq == 0)	//
	{
		
		//printf("need more\n");
		display_pregraph_usage ();
		exit (-1);
	}
}

static void display_pregraph_usage () 
{	
	printf ("\npregraph -s configFile -o outputGraph [-K kmer -p n_cpu -d kmerFreqCutoff]\n");
	printf ("  -s\t<string>\tconfigFile: the config file of reads\n");
	printf ("  -o\t<string>\toutputGraph: prefix of output graph file name\n");		
#ifdef MER127
	printf ("  -K\t<int>\t\tkmer(min 13, max 127): kmer size, [23]\n");
#endif
#ifdef MER63
	printf ("  -K\t<int>\t\tkmer(min 13, max 63): kmer size, [23]\n");
#endif
#ifdef MER31
	printf ("  -K\t<int>\t\tkmer(min 13, max 31): kmer size, [23]\n");
#endif
	printf ("  -p\t<int>\t\tn_cpu: number of cpu for use, [8]\n");	
	printf ("  -d\t<int>\t\tkmerFreqCutoff: kmers with frequency no larger than KmerFreqCutoff will be deleted, [0]\n");	
	
} 
