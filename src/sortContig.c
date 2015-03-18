#include "stdinc.h"
#include "newhash.h"
#include "extfunc.h"
#include "extvab.h"

ARC * sort_arc ( ARC * list )
{
	if ( !list )
		{ return list; }

	ARC * head = ( ARC * ) malloc ( sizeof ( ARC ) );
	head->next = list;
	list->prev = head;
	ARC * curr = list;
	ARC * temp = list;
	ARC * temp1 = NULL;

	while ( curr )
	{
		temp = curr;

		if ( temp )
		{
			temp1 = temp->next;

			while ( temp1 )
			{
				if ( temp->to_ed > temp1->to_ed )
					{ temp = temp1; }

				temp1 = temp1->next;
			}
		}

		if ( temp && temp != curr )
		{
			if ( temp->next )
			{
				temp->prev->next = temp->next;
				temp->next->prev = temp->prev;
			}
			else
			{
				temp->prev->next = NULL;
			}

			temp->next = curr;
			temp->prev = curr->prev;
			curr->prev->next = temp;
			curr->prev = temp;
		}
		else
		{
			curr = curr->next;
		}
	}

	list = head->next;
	list->prev = NULL;
	head->next = NULL;
	free ( head );
	return list;
};
static void copyOneEdge ( EDGE * target , EDGE * source )
{
	target->from_vt = source->from_vt;
	target->to_vt = source->to_vt;
	target->length = source->length;
	target->cvg = source->cvg;
	target->multi = source->multi;
	target->flag = source->flag;
	target->bal_edge = source->bal_edge;
	target->seq = source->seq;
	source->seq = NULL;
	target->arcs = source->arcs;
	source->arcs = NULL ;
	target->markers = source->markers;
	source->markers = NULL;
	target->deleted = source->deleted;
};
static void updateArcToEd ( unsigned int ed_index )
{
	ARC * arc = edge_array[ed_index].arcs;

	while ( arc )
	{
		arc->to_ed = index_array[arc->to_ed];
		arc = arc->next;
	}
};
inline void delete1Edge ( unsigned int index )
{
	edge_array[index].deleted = 1;
};
//Copy edge from source to target.
void copy1Edge ( EDGE * source, EDGE * target )
{
	target->from_vt = source->from_vt;
	target->to_vt = source->to_vt;
	target->length = source->length;
	target->cvg = source->cvg;
	target->multi = source->multi;

	if ( target->seq )
	{
		free ( ( void * ) target->seq );
	}

	target->seq = source->seq;
	source->seq = NULL;
	target->arcs = source->arcs;
	source->arcs = NULL;
	target->deleted = source->deleted;
};

//Check whether two bases are equal.
int BaseEqual ( char ch1, char ch2 )
{
	if ( ch1 == ch2 )
		{ return 0; }
	else if ( ch1 > ch2 )
		{ return 1; }
	else
		{ return -1; }
};

int EdgeEqual ( unsigned int prev, unsigned int next )
{
	int i = 0;
	int length = edge_array[prev].length;
	char ch1, ch2;
	int equal = 0;

	for ( i = 0; i < length; ++i )
	{
		ch1 = int2base ( ( int ) getCharInTightString ( edge_array[prev].seq, i ) );
		ch2 = int2base ( ( int ) getCharInTightString ( edge_array[next].seq, i ) );

		if ( ( equal = BaseEqual ( ch1, ch2 ) ) )
		{
			return equal;
		}
	}

	return 0;
};
void swapedge()
{
	unsigned int i;
	ARC * arc, *bal_arc, *temp_arc;
	int count_swap = 0, count_equal = 0;

	for ( i = 1; i <= num_ed; ++i )
	{
		if ( edge_array[i].deleted || EdSameAsTwin ( i ) )
			{ continue; }

		if ( EdSmallerThanTwin ( i ) )
		{
			if ( KmerLarger ( vt_array[edge_array[i].from_vt].kmer, vt_array[edge_array[i + 1].from_vt].kmer ) )
			{
				count_swap++;
				copyEdge ( i, num_ed + 1 + 1 );
				copyEdge ( i + 1, num_ed + 1 );
				copyEdge ( num_ed + 1, i );
				copyEdge ( num_ed + 1 + 1, i + 1 );
				edge_array[i].bal_edge = 2;
				edge_array[i + 1].bal_edge = 0;
				//take care of the arcs
				arc = edge_array[i].arcs;

				while ( arc )
				{
					arc->bal_arc->to_ed = i + 1;
					arc = arc->next;
				}

				arc = edge_array[i + 1].arcs;

				while ( arc )
				{
					arc->bal_arc->to_ed = i;
					arc = arc->next;
				}
			}
			else if ( KmerEqual ( vt_array[edge_array[i].from_vt].kmer, vt_array[edge_array[i + 1].from_vt].kmer ) )
			{
				int temp = EdgeEqual ( i, i + 1 );

				if ( temp == 0 )
				{
					count_equal++;
					edge_array[i].bal_edge = 1;
					delete1Edge ( i + 1 );
					//take care of the arcs
					arc = edge_array[i].arcs;

					while ( arc )
					{
						arc->bal_arc->to_ed = i;
						arc = arc->next;
					}

					bal_arc = edge_array[i + 1].arcs;
					edge_array[i + 1].arcs = NULL;

					while ( bal_arc )
					{
						temp_arc = bal_arc;
						bal_arc = bal_arc->next;

						if ( edge_array[i].arcs )
							{ edge_array[i].arcs->prev = temp_arc; }

						temp_arc->next = edge_array[i].arcs;
						edge_array[i].arcs = temp_arc;
					}
				}
				else if ( temp > 0 )
				{
					count_swap++;
					copyEdge ( i, num_ed + 1 + 1 );
					copyEdge ( i + 1, num_ed + 1 );
					copyEdge ( num_ed + 1, i );
					copyEdge ( num_ed + 1 + 1, i + 1 );
					edge_array[i].bal_edge = 2;
					edge_array[i + 1].bal_edge = 0;
					//take care of the arcs
					arc = edge_array[i].arcs;

					while ( arc )
					{
						arc->bal_arc->to_ed = i + 1;
						arc = arc->next;
					}

					arc = edge_array[i + 1].arcs;

					while ( arc )
					{
						arc->bal_arc->to_ed = i;
						arc = arc->next;
					}
				}
			}

			++i;
		}
		else
		{
			delete1Edge ( i );
			printf( "Warning : Front edge %d is larger than %d.\n", i, i + 1 );
		}
	}

	printf( "%d none-palindrome edge(s) swapped, %d palindrome edge(s) processed.\n", count_swap, count_equal );
};
static int cmp_seq ( const void * a, const void * b )
{
	EDGE * A, *B;
	A = ( EDGE * ) a;
	B = ( EDGE * ) b;

	if ( KmerLarger ( vt_array[A->from_vt].kmer, vt_array[B->from_vt].kmer ) )
	{
		return 1;
	}
	else if ( KmerSmaller ( vt_array[A->from_vt].kmer , vt_array[B->from_vt].kmer ) )
	{
		return -1;
	}
	else
	{
		if ( A->seq[0] > B->seq[0] )
		{
			return 1;
		}
		else if ( A->seq[0] == B->seq[0] )
		{
			int i = 0;

			for ( i = 1; i < A->length && i < B->length; i++ )
			{
				if ( getCharInTightString ( A->seq, i ) > getCharInTightString ( B->seq, i ) )
					{ return 1; }
				else if ( getCharInTightString ( A->seq, i ) < getCharInTightString ( B->seq, i ) )
					{ return -1; }
			}

			if ( i == A->length && i < B->length )
				{ return -1; }
			else if ( i < A->length && i ==  B->length )
				{ return 1; }
			else
			{
				printKmerSeq ( stderr , vt_array[A->from_vt].kmer );
				fprintf ( stderr , "\n" );
				printKmerSeq ( stderr , vt_array[B->from_vt].kmer );
				fprintf ( stderr , "\n" );

				for ( i = 0; i < A->length; i++ )
				{
					printf( "%c", int2base ( ( int ) getCharInTightString ( A->seq, i ) ) );
				}

				printf( "\n" );

				for ( i = 0; i < B->length; i++ )
				{
					printf( "%c", int2base ( ( int ) getCharInTightString ( B->seq, i ) ) );
				}

				printf( "\n" );
				printf( "cmp_seq:\terr\n" );
				exit ( 0 );
				return 0;
			}
		}
		else
		{
			return -1;
		}
	}
};
void sortedge()
{
	unsigned int index ;
	EDGE * sort_edge , * backup_edge ;
	sort_edge = ( EDGE * ) ckalloc ( sizeof ( EDGE ) * ( num_ed + 1 ) );
	backup_edge = ( EDGE * ) ckalloc ( sizeof ( EDGE ) * ( num_ed + 1 ) );
	unsigned int i = 1;

	for ( index = 1 ; index <= num_ed ; index ++ )
	{
		sort_edge[i].from_vt = edge_array[index].from_vt;
		sort_edge[i].seq = edge_array[index].seq;
		sort_edge[i].to_vt = index; // record old id
		sort_edge[i].length = edge_array[index].length;
		i++;
		copyOneEdge ( & ( backup_edge[index] ) , & ( edge_array[index] ) );

		if ( !EdSameAsTwin ( index ) )
		{
			index++;
			copyOneEdge ( & ( backup_edge[index] ) , & ( edge_array[index] ) );
		}
	}

	qsort ( & ( sort_edge[1] ), i - 1, sizeof ( sort_edge[1] ), cmp_seq );
	index_array = ( unsigned int * ) ckalloc ( sizeof ( unsigned int ) * ( num_ed + 1 ) ); // used to record new id
	unsigned int new_index = 1, old_index;

	for ( index = 1; index <= i - 1; index++ )
	{
		old_index = sort_edge[index].to_vt; // old id
		sort_edge[index].seq = NULL;
		index_array[old_index] = new_index++;// old id -> new id

		if ( !EdSameAsTwin ( old_index ) )
		{
			index_array[old_index + 1] = new_index++; // old id -> new id
		}
	}

	for ( index = 1; index <= num_ed; index++ )
	{
		new_index = index_array[index];
		copyOneEdge ( & ( edge_array[new_index] ), & ( backup_edge[index] ) );
		updateArcToEd ( new_index );
	}

	free ( index_array );
	free ( sort_edge );
	free ( backup_edge );
};
void freshArc()
{
	unsigned int i;
	ARC * arc_temp, *parc;

	for ( i = 1; i <= num_ed; ++i )
	{
		if ( edge_array[i].deleted )
			{ continue; }

		edge_array[i].arcs = sort_arc ( edge_array[i].arcs );
	}
};
