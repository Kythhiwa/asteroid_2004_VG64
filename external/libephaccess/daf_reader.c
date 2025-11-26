/* 
 * Файл         : daf_reader.c
 * Исполнитель  : ИПА РАН
 * Заказчик     : ОАО НПК СПП
 * Проект       : КР ФЭЛП, ОКР "Эфемериды"
 * Год          : 2015
 */

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "daf_reader.h"

static int _readExactDouble (FILE *f, int *result)
{
  double value;
  
  if (fread(&value, 8, 1, f) != 1)
    return 0;
  
  *result = (int)value;
  if (fabs(*result - value) > value * 1e-14)
    return 0;
  return 1;
}

static int _readExactDoubleAndCheckAgainst (FILE *f, int expected)
{
  double value;
  
  if (fread(&value, 8, 1, f) != 1)
    return 0;
  
  if (fabs(expected - value) > value * 1e-14)
    return 0;
  return 1;
}

int dafReadFile (FILE *f, DAF *daf)
{
  char id[8];
  int first_summary, last_summary;
  int prev_summary;
  int summary, isummary;
  
  memset(daf, 0, sizeof(DAF));
  
  // read 8-byte ID
  if (fread(id, 8, 1, f) != 1)
    return DAF_ERROR_BAD_FORMAT;

  if (memcmp(id, "DAF/SPK", 7) == 0 ||
      memcmp(id, "NAIF/DAF", 8) == 0)
    daf->filetype = DAF_FILE_SPK;
  else if (memcmp(id, "DAF/PCK", 7) == 0)
    daf->filetype = DAF_FILE_PCK;
  else
    return DAF_ERROR_UNSUPPORTED;
  
  // read ND and NI
  if (fread(&daf->nd, 4, 1, f) != 1)
    return DAF_ERROR_BAD_FORMAT;
  if (fread(&daf->ni, 4, 1, f) != 1)
    return DAF_ERROR_BAD_FORMAT;
  
  // there must be at least two integer parameters: numbers of the first
  // and the last block of data
  if (daf->ni < 2)
    return DAF_ERROR_BAD_FORMAT;
  
  // skip 60-bytes name
  if (fseek(f, 60, SEEK_CUR) != 0)
    return DAF_ERROR_BAD_FORMAT;
  
  if (fread(&first_summary, 4, 1, f) != 1)
    return DAF_ERROR_BAD_FORMAT;
  
  if (fread(&last_summary, 4, 1, f) != 1)
    return DAF_ERROR_BAD_FORMAT;
  
  // skip free address
  if (fseek(f, 4, SEEK_CUR) != 0)
    return DAF_ERROR_BAD_FORMAT;
  
  // skip 'reserved' blocks till the beginning of the first summary
#ifdef _WIN32
  if (_fseeki64(f, (__int64)(first_summary - 2) * 1024, SEEK_SET) != 0)
#else
  if (fseek(f, (first_summary - 2) * 1024L, SEEK_SET) != 0)
#endif
    return DAF_ERROR_BAD_FORMAT;
    
  // loop through all blocks of summaries to calculate the number of summaries
  for (summary = first_summary, prev_summary = 0; summary != 0; )
  {
    int next_summary;
    int n_summaries;
    
#ifdef _WIN32
    if (_fseeki64(f, (__int64)(summary - 1) * 1024, SEEK_SET) != 0)
#else
    if (fseek(f, (summary - 1) * 1024L, SEEK_SET) != 0)
#endif
      return DAF_ERROR_BAD_FORMAT;
    
    if (!_readExactDouble(f, &next_summary))
      return DAF_ERROR_BAD_FORMAT;
    if (!_readExactDoubleAndCheckAgainst(f, prev_summary))
      return DAF_ERROR_BAD_FORMAT; // inconsistency
    if (!_readExactDouble(f, &n_summaries))
      return DAF_ERROR_BAD_FORMAT;
  
    daf->nsegments += n_summaries;
    prev_summary = summary;
    summary = next_summary;
  }
  
  if (prev_summary != last_summary)
    return DAF_ERROR_BAD_FORMAT; // inconsistency
  
  daf->all_parameters_d = (double *)malloc(daf->nd * daf->nsegments * sizeof(double));
  if (daf->all_parameters_d == NULL)
    return DAF_ERROR_NO_MEMORY;
  
  daf->all_parameters_i = (int *)malloc((daf->ni - 2) * daf->nsegments * sizeof(int));
  if (daf->all_parameters_i == NULL)
    return DAF_ERROR_NO_MEMORY;
  
  daf->segments = (DafSegment *)malloc(daf->nsegments * sizeof(DafSegment));
  if (daf->segments == NULL)
    return DAF_ERROR_NO_MEMORY;
  
  isummary = 0;
  
  // now loop through all the blocks again and gather summaries
  for (summary = first_summary, prev_summary = 0; summary != 0; )
  {
    int next_summary;
    int n_summaries;
    
    // skip error checking this time
#ifdef _WIN32
    _fseeki64(f, (__int64)(summary - 1) * 1024, SEEK_SET);
#else
    fseek(f, (summary - 1) * 1024L, SEEK_SET);
#endif
    
    _readExactDouble(f, &next_summary);
    fseek(f, 8, SEEK_CUR);
    _readExactDouble(f, &n_summaries);

    while (n_summaries-- > 0)
    {
      int initial_address, final_address;

      daf->segments[isummary].parameters_d = daf->all_parameters_d +
                                             isummary * daf->nd;
      daf->segments[isummary].parameters_i = daf->all_parameters_i +
                                             isummary * (daf->ni - 2);

      if (fread(daf->segments[isummary].parameters_d, sizeof(double),
              daf->nd, f) != daf->nd)
        return DAF_ERROR_BAD_FORMAT;

      if (fread(daf->segments[isummary].parameters_i, sizeof(int),
              daf->ni - 2, f) != daf->ni - 2)
        return DAF_ERROR_BAD_FORMAT;

      if (fread(&initial_address, sizeof(int), 1, f) != 1)
        return DAF_ERROR_BAD_FORMAT;
      if (fread(&final_address, sizeof(int), 1, f) != 1)
        return DAF_ERROR_BAD_FORMAT;
      
      // skip the rest of the summary record that is an integer amount
      // of 8-byte blocks
      if (daf->ni % 2 != 0)
        if (fseek(f, 4, SEEK_CUR) != 0)
          return DAF_ERROR_BAD_FORMAT;

      daf->segments[isummary].offset = (initial_address - 1);
      daf->segments[isummary].length = final_address - initial_address + 1;
      daf->segments[isummary].file = f;
      isummary++;
    }
    prev_summary = summary;
    summary = next_summary;
  }
  
  daf->file = f;
  return 0;
}

int dafSegmentReadRange (DafSegment *segment, int start, int length, double *result)
{
  if (start + length > segment->length)
    return 0;
  if (length < 0 || start < 0)
    return 0;
#ifdef _WIN32
  if (_fseeki64(segment->file, (((__int64)segment->offset) + start) * sizeof(double), SEEK_SET) != 0)
    return 0;
#else
  if (fseek(segment->file, (long)(segment->offset + start) * sizeof(double), SEEK_SET) != 0)
    return 0;
#endif
  if (fread(result, sizeof(double), length, segment->file) != length)
    return 0;
  return 1;
}

void dafFreeData (DAF *daf)
{
  if (daf == NULL)
    return;
  free(daf->all_parameters_d);
  free(daf->all_parameters_i);
  free(daf->segments);
  fclose(daf->file);
}
