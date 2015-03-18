/*
 * inc/Sam_header.h
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
#ifndef __SAM_HEADER_H__
#define __SAM_HEADER_H__

#ifdef __cplusplus
extern "C" {
#endif

	void *sam_header_parse2(const char *headerText);
	void *sam_header_merge(int n, const void **dicts);
	void sam_header_free(void *header);
	char *sam_header_write(const void *headerDict);   // returns a newly allocated string

	char **sam_header2list(const void *_dict, char type[2], char key_tag[2], int *_n);

	void *sam_header2tbl(const void *dict, char type[2], char key_tag[2], char value_tag[2]);
	const char *sam_tbl_get(void *h, const char *key);
	int sam_tbl_size(void *h);
	void sam_tbl_destroy(void *h);

#ifdef __cplusplus
}
#endif

#endif
