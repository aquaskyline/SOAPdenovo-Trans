/*
 * Check.c
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

#include <stdinc.h>
void *ckalloc (unsigned long long amount);
FILE *ckopen (char *name, char *mode);
FILE *ckopen (char *name, char *mode)
{
	FILE *fp;

	if ((fp = fopen (name, mode)) == NULL)
	{
		printf ("Cannot open %s. Now exit to system...\n", name);
		exit (-1);
	}
	return (fp);
}

/* ckalloc - allocate space; check for success */
void *ckalloc (unsigned long long amount)
{
	void *p;

	if ((p = (void *) calloc (1, (unsigned long long) amount)) == NULL && amount != 0)
	{
		printf ("Ran out of memory while applying %lldbytes\n", amount);
		printf ("There may be errors as follows:\n");
		printf ("1) Not enough memory.\n");
		printf ("2) The ARRAY may be overrode.\n");
		printf ("3) The wild pointers.\n");
		fflush (stdout);
		exit (-1);
	}
	return (p);
}

/* reallocate memory */
void *ckrealloc (void *p, size_t new_size, size_t old_size)
{
	void *q;

	q = realloc ((void *) p, new_size);
	if (new_size == 0 || q != (void *) 0)
	{
		return q;
	}
	/* manually reallocate space */
	q = ckalloc (new_size);

	/* move old memory to new space */
	bcopy (p, q, old_size);
	free (p);
	return q;
}
