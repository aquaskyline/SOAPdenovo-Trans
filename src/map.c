/*
 * Map.c
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

static void display_map_usage ();

static int getMinOverlap (char *gfile)
{
	char name[256], ch;
	FILE *fp;
	int num_kmer, overlaplen = 23;
	char line[1024];

	sprintf (name, "%s.preGraphBasic", gfile);
	fp = fopen (name, "r");

	if (!fp)
	{
		return overlaplen;
	}

	while (fgets (line, sizeof (line), fp) != NULL)
	{
		if (line[0] == 'V')
		{
			sscanf (line + 6, "%d %c %d", &num_kmer, &ch, &overlaplen);
		}
		else if (line[0] == 'M')
		{
			sscanf (line, "MaxReadLen %d MinReadLen %d MaxNameLen %d", &maxReadLen, &minReadLen, &maxNameLen);
		}
	}

	fclose (fp);
	return overlaplen;
}

int call_align (int argc, char **argv)
{
	time_t start_t, stop_t, time_bef, time_aft;
	N_kmer=0;
	time (&start_t);
	if(strlen(graphfile) == 0)
	{
		initenv (argc, argv);
	}
	overlaplen = getMinOverlap (graphfile);
	initenv (argc, argv);
	
	printf ("K = %d\n", overlaplen);
	time (&time_bef);
	ctg_short = overlaplen + 2;
	printf ("contig len cutoff: %d\n", ctg_short);
	prlContig2nodes (graphfile, ctg_short);
	time (&time_aft);
	printf ("time spent on De bruijn graph construction: %ds\n\n", (int) (time_aft - time_bef));
	
	//map long read to edge one by one
	time (&time_bef);
//	prlLongRead2Ctg (shortrdsfile, graphfile);
	time (&time_aft);
//	printf ("time spent on mapping long reads: %ds\n\n", (int) (time_aft - time_bef));
	
	//map read to edge one by one
	time (&time_bef);
	prlRead2Ctg (shortrdsfile, graphfile);
	time (&time_aft);
	printf ("time spent on mapping reads: %ds\n\n", (int) (time_aft - time_bef));
	free_Sets (KmerSets, thrd_num);
	time (&stop_t);
	printf ("overall time for alignment: %dm\n\n", (int) (stop_t - start_t) / 60);
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

	while ((copt = getopt (argc, argv, "s:g:K:p:rfR")) != EOF)
	{
		//printf("get option\n");
		switch (copt)
		{
		case 's':
			inpseq = 1;
			sscanf (optarg, "%s", shortrdsfile);
			break;
		case 'g':
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
		case 'r':
			read_trace=1;
			break;
		case 'R':                                   
			RPKM=1;
			read_trace= 1;
			break;	
		case 'f':
			fill=1;
			break;			
		default:
			if (inpseq == 0 || outseq == 0)	//
			{
				display_map_usage ();
				exit (1);
			}
		}
	}

	if (inpseq == 0 || outseq == 0)	//
	{
		//printf("need more\n");
		display_map_usage ();
		exit (1);
	}
}

static void display_map_usage ()
{
	printf ("\nmap -s configFile -g inputGraph [-f -R] [-p n_cpu]\n");
	printf ("  -s\t<string>\tconfigFile: the config file of reads\n");
	printf ("  -g\t<string>\tinputGraph: prefix of input graph file name\n");	
//	printf ("  -r\t(optional)\toutput the information between read and contig, [NO]\n");	
	printf ("  -f\t(optional)\toutput gap related reads for SRkgf to fill gap, [NO]\n");
	printf ("  -R\t(optional)\toutput assembly RPKM statistics, [NO]\n");	
	printf ("  -p\t<int>\t\tn_cpu: number of cpu for use, [8]\n");	
}
