/*
 * main.c
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
#include "global.h"

extern int call_pregraph (int arc, char **argv);
extern int call_heavygraph (int arc, char **argv);
extern int call_map2contig (int arc, char **argv);
//extern call_transcriptome(int arc, char **argv);
extern int call_scaffold (int arc, char **argv);
extern int call_align (int arc, char **argv);

static void display_usage ();
static void display_all_usage ();
static void pipeline (int argc, char **argv);

int main (int argc, char **argv)
{
	printf ("\nVersion 1.04\n\n");
	argc--;
	argv++;

	/*
	   __uint128_t temp;
	   ubyte8 long1=0x5ad6c7ef8a;
	   ubyte8 long2=0x87a3c27a2b;
	   temp = long1;
	   temp <<= 64;
	   temp |=long2;
	   long2 = (ubyte8)temp;
	   long1 = (ubyte8)(temp>>64);

	   printf("%p,%p,%p\n",long1,long2,temp);
	 */
	if (argc == 0)
	{
		display_usage ();
		return 0;
	}

	if (strcmp ("pregraph", argv[0]) == 0)
	{
		call_pregraph (argc, argv);
	}
	else if (strcmp ("contig", argv[0]) == 0)
	{
		call_heavygraph (argc, argv);
	}
	else if (strcmp ("map", argv[0]) == 0)
	{
		call_align (argc, argv);
	}
	//call_map2contig(argc,argv);
/*
	else if (strcmp ("trans", argv[0]) == 0)
	{
		call_transcriptome (argc, argv);
	}
*/
	else if (strcmp ("scaff", argv[0]) == 0)
	{
		call_scaffold (argc, argv);
	}
	else if (strcmp ("all", argv[0]) == 0)
	{
		pipeline (argc, argv);
	}
	else
	{
		display_usage ();
	}

	return 0;
}

static void display_usage ()
{
	printf ("\nUsage: SOAPdenovo-Trans <command> [option]\n");
	printf ("    pregraph     construction kmer-graph\n");
	printf ("    contig       eliminate errors and output contigs\n");
	printf ("    map          map reads to contigs\n");
//	printf ("    trans          get all path\n");
	printf ("    scaff        scaffolding\n");
	printf ("    all          doing all the above in turn\n");
}

static void pipeline (int argc, char **argv)
{
	char *options[16];
	unsigned char getK, getRfile, getOfile, getD, getDD, getL, getR, getP, getF, getf, getk;
	char readfile[256], outfile[256];
	char temp[128];
	char *name;
	int kmer = 0, cutoff_len = 0, ncpu = 0, lowK = 0, lowC = 0, kmer_small = 0;
	char kmer_s[16], len_s[16], ncpu_s[16], M_s[16], lowK_s[16], lowC_s[16], kmer_small_s[16];
	int i, copt, index, M = 1;
	extern char *optarg;
	time_t start_t, stop_t;

	time (&start_t);
	getK = getRfile = getOfile = getD = getDD = getL = getR = getP = getF = getf = getk = 0;

//	while ((copt = getopt (argc, argv, "a:s:o:K:M:L:p:G:d:D:RuFk:f")) != EOF)
	while ((copt = getopt (argc, argv, "a:s:o:K:L:p:G:d:D:uFk:fi:ne:q:Q:H:rt:c:M:R")) != EOF)
	{
		switch (copt)
		{
		case 's':
			getRfile = 1;
			sscanf (optarg, "%s", readfile);
			break;
		case 'o':
			getOfile = 1;
			sscanf (optarg, "%s", outfile);	//
			break;
		case 'K':
			getK = 1;
			sscanf (optarg, "%s", temp);	//
			kmer = atoi (temp);
			break;
		case 'G':
			sscanf (optarg, "%s", temp);	//
			GLDiff = atoi (temp);
			break;
		case 'M':
			sscanf (optarg, "%s", temp);	//
			M = atoi (temp);
			break;
		case 'p':
			getP = 1;
			sscanf (optarg, "%s", temp);	//
			ncpu = atoi (temp);
			break;
		case 'L':
			getL = 1;
			sscanf (optarg, "%s", temp);	//
			ctg_mask = atoi (temp);
			break;
/*		case 'R':
			getR = 1;
			break;*/
		case 'u':
			maskRep = 0;
			break;
		case 'd':
			getD = 1;
			sscanf (optarg, "%s", temp);
			lowK = atoi (temp);
			break;
		case 'D':
			getDD = 1;
			sscanf (optarg, "%s", temp);
			lowC = atoi (temp);
			break;
		case 'a':
			initKmerSetSize = atoi (optarg);
			break;
		case 'F':
			getF = 1;
			break;
		case 'k':
			getk = 1;
			sscanf(optarg,"%s", temp);   
			kmer_small = atoi(temp);
			break;
		case 'f':
			getf = 1;
			break;
		case 'i'://mao 2011-10-21
			sscanf (optarg, "%s", temp);
			dd=atoi (temp) >= 0 ? atoi (temp) : 0;
			break;
		case 'n'://mao 2011-10-21
			N_kmer=1;
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
		case 'r':
			read_trace=1;
			break;
		case 't':
			sscanf (optarg, "%s", temp);    //
			max_num = atoi (temp)>0? atoi (temp):5;
			break;
		case 'c':
			sscanf (optarg, "%s", temp);    //
			max_cnt  = atoi (temp)>=0? atoi (temp):0;
			break;
		case 'R':                                   
			RPKM=1;
			read_trace= 1;
			break;	
		default:

			if (getRfile == 0 || getOfile == 0)	//
			{
				display_all_usage ();
				exit (-1);
			}
		}
	}

	if (getRfile == 0 || getOfile == 0)	//
	{
		display_all_usage ();
		exit (-1);
	}

	if (thrd_num < 1)
	{
		thrd_num = 1;
	}

	// getK = getRfile = getOfile = getD = getL = getR = 0;
	name = "pregraph";
	index = 0;
	options[index++] = name;
	options[index++] = "-s";
	options[index++] = readfile;

	if (getK)
	{
		options[index++] = "-K";
		sprintf (kmer_s, "%d", kmer);
		options[index++] = kmer_s;
	}

	if (getP)
	{
		options[index++] = "-p";
		sprintf (ncpu_s, "%d", ncpu);
		options[index++] = ncpu_s;
	}

	if (getD)
	{
		options[index++] = "-d";
		sprintf (lowK_s, "%d", lowK);
		options[index++] = lowK_s;
	}

	if (getR)
	{
		options[index++] = "-R";
	}

	options[index++] = "-o";
	options[index++] = outfile;

	for (i = 0; i < index; i++)
	{
		printf ("%s ", options[i]);
	}

	printf ("\n");
	call_pregraph (index, options);

	name = "contig";
	index = 0;
	options[index++] = name;
	options[index++] = "-g";
	options[index++] = outfile;
	options[index++] = "-M";
	sprintf (M_s, "%d", M);
	options[index++] = M_s;
/*
	if (getR)
	{
		options[index++] = "-R";
	}

	if (getDD)
	{
		options[index++] = "-D";
		sprintf (lowC_s, "%d", lowC);
		options[index++] = lowC_s;
	}*/

	for (i = 0; i < index; i++)
	{
		printf ("%s ", options[i]);
	}

	printf ("\n");
	call_heavygraph (index, options);
	
	name = "map";
	index = 0;
	options[index++] = name;
	options[index++] = "-s";
	options[index++] = readfile;
	options[index++] = "-g";
	options[index++] = outfile;

	if (getP)
	{
		options[index++] = "-p";
		sprintf (ncpu_s, "%d", ncpu);
		options[index++] = ncpu_s;
	}

	if (getK)
	{
	    options[index++] = "-K";
	    sprintf (kmer_s, "%d", kmer);
	    options[index++] = kmer_s;
	}
			
	
	if(getk){
		options[index++] = "-k";
		sprintf(kmer_small_s,"%d",kmer_small);
		options[index++] = kmer_small_s;
	}
	if (getf)
	{
		options[index++] = "-f";
	}

	for (i = 0; i < index; i++)
	{
		printf ("%s ", options[i]);	
	}

	printf ("\n");
	call_align (index, options);
	
	name = "scaff";
	index = 0;
	options[index++] = name;
	options[index++] = "-g";
	options[index++] = outfile;

	if (getF)
	{
		options[index++] = "-F";
	}

	if (getP)
	{
		options[index++] = "-p";
		sprintf (ncpu_s, "%d", ncpu);
		options[index++] = ncpu_s;
	}
/*
	if (getL)
	{
		options[index++] = "-L";
		sprintf (len_s, "%d", ctg_mask);
		options[index++] = len_s;
	}
	*/

	for (i = 0; i < index; i++)
	{
		printf ("%s ", options[i]);
	}

	printf ("\n");
	call_scaffold (index, options);
	time (&stop_t);
	printf ("time for the whole pipeline: %dm\n", (int) (stop_t - start_t) / 60);
}

static void display_all_usage ()
{
	printf ("\nSOAPdenovo-Trans all -s configFile -o outputGraph [-R -f -S -F] [-K kmer -p n_cpu -d kmerFreqCutoff -e EdgeCovCutoff -M mergeLevel -L minContigLen -t locusMaxOutput -G gapLenDiff]\n");

	printf ("  -s\t<string>\tconfigFile: the config file of reads\n");
	printf ("  -o\t<string>\toutputGraph: prefix of output graph file name\n");
	
//	printf ("  -r\t(optional)\toutput the information between read and scaffold, [NO]\n");
	printf ("  -R\t(optional)\toutput assembly RPKM statistics\n");
	printf ("  -f\t(optional)\toutput gap related reads for SRkgf to fill gap, [NO]\n");
	printf ("  -S\t(optional)\tscaffold structure exists, [NO]\n");
	printf ("  -F\t(optional)\tfill gaps in scaffolds, [NO]\n");
	
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
	printf ("  -e\t<int>\t\tEdgeCovCutoff: edges with coverage no larger than EdgeCovCutoff will be deleted, [2]\n");
	printf ("  -M\t<int>\t\tmergeLevel(min 0, max 3): the strength of merging similar sequences during contiging, [1]\n");
	printf ("  -L\t<int>\t\tminContigLen: shortest contig for scaffolding, [100]\n");
	printf ("  -t\t<int>\t\tlocusMaxOutput: output the number of transcripts no more than locusMaxOutput in one locus, [5]\n");
	printf ("  -G\t<int>\t\tgapLenDiff: allowed length difference between estimated and filled gap, [50]\n");	
}
