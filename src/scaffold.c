/*
 * Scaffold.c
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
#include "transcriptome.c"

static void initenv (int argc, char **argv);
static void display_scaff_usage ();

static boolean LINK, SCAFF;
static char graphfile[256];

int call_scaffold (int argc, char **argv)
{
	time_t start_t, stop_t, time_bef, time_aft;

	time (&start_t);
	initenv (argc, argv);
	loadPEgrads (graphfile);
	time (&time_bef);
	loadUpdatedEdges (graphfile);
	time (&time_aft);
	printf ("time spent on loading edges %ds\n", (int) (time_aft - time_bef));

	if (!SCAFF)
	{
		time (&time_bef);
		PE2Links (graphfile);
		time (&time_aft);
		printf ("time spent on loading pair end info %ds\n\n", (int) (time_aft - time_bef));
		
		time (&time_bef);
		Links2Scaf (graphfile);
		time (&time_aft);
		printf ("time spent on creating scaffolds %ds\n\n", (int) (time_aft - time_bef));
		
		time(&time_bef);
		transcriptome(graphfile);
		time(&time_aft);
		printf("time spent on creating transcriptome %ds\n",(int)(time_aft-time_bef));
//		scaffolding (100, graphfile);
	}

	prlReadsCloseGap (graphfile);
	//  locateReadOnScaf(graphfile);
	ScafStat (100, graphfile);
	if(read_trace)
	{
		getReadOnScaf(graphfile);
		if(RPKM)   //Must add  '-R'  parameter    RPKMStat(graphfile);
			RPKMStat(graphfile);
	}
	
	free_pe_mem ();

	if (index_array)
	{
		free ((void *) index_array);
	}

	freeContig_array ();
	destroyPreArcMem ();
	destroyConnectMem ();
	deleteCntLookupTable ();
	time (&stop_t);
	printf ("time elapsed: %dm\n", (int) (stop_t - start_t) / 60);
	return 0;
}

/*****************************************************************************
 * Parse command line switches
 *****************************************************************************/

void initenv (int argc, char **argv)
{
	int copt;
	int inpseq;
	extern char *optarg;
	char temp[256];

	inpseq = 0;
	LINK = 0;
	SCAFF = 0;
	optind = 1;

	while ((copt = getopt (argc, argv, "g:L:p:G:N:FuSt:c:rR")) != EOF)
	{
		switch (copt)
		{
		case 'g':
			inGraph = 1;
			sscanf (optarg, "%s", graphfile);	//
			break;
		case 'G':
			sscanf (optarg, "%s", temp);	//
			GLDiff = atoi (temp);
			break;
		case 'L':
			sscanf (optarg, "%s", temp);
			ctg_mask = atoi (temp);
			break;
		case 'N':                                   
			sscanf ( optarg, "%s", temp );         
			known_genome_size = atoi ( temp );    
			break;	
		case 'F':
			fillGap = 1;
			break;
		case 'S':
			SCAFF = 1;
			break;
		case 'u':
			maskRep = 0;
			break;
		case 'p':
			sscanf (optarg, "%s", temp);	//
			thrd_num = atoi (temp);
			break;
		case 't':
			sscanf (optarg, "%s", temp);    //
			max_num = atoi (temp)>0? atoi (temp):5;
			break;
		case 'c':
			sscanf (optarg, "%s", temp);    //
			max_cnt  = atoi (temp)>=0? atoi (temp):0;
			break;		
		case 'r':			
			read_trace= 1;
			break;
		case 'R':                                   
			RPKM=1;
			read_trace= 1;
			break;	
		default:
			if (inGraph == 0)	//
			{
				display_scaff_usage ();
				exit (-1);
			}
		}
	}

	if (inGraph == 0)	//
	{
		display_scaff_usage ();
		exit (-1);
	}
}

static void display_scaff_usage ()
{
	printf ("\nscaff -g inputGraph [-R -S -F] [-p n_cpu -L minContigLen -t locusMaxOutput -G gapLenDiff]\n");
	printf ("  -g\t<string>\tinputGraph: prefix of input graph file name\n");
//	printf ("  -r\t(optional)\toutput the information between read and scaffold, [NO]\n");
	printf ("  -R\t\t(optional) output assembly RPKM statistics, [NO]\n");	
	printf ("  -S\t(optional)\tscaffold structure exists, [NO]\n");		
	printf ("  -F\t(optional)\tfill gaps in scaffolds, [NO]\n");	
	printf ("  -p\t<int>\t\tn_cpu: number of cpu for use, [8]\n");
	printf ("  -L\t<int>\t\tminContigLen: shortest contig for scaffolding, [100]\n");
	printf ("  -t\t<int>\t\tlocusMaxOutput: output the number of transcripts no more than locusMaxOutput in one locus, [5]\n");	
	printf ("  -G\t<int>\t\tgapLenDiff: allowed length difference between estimated and filled gap, [50]\n");	
}
