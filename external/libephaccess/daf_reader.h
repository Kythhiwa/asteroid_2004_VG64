/* 
 * Файл         : daf_reader.h
 * Исполнитель  : ИПА РАН
 * Заказчик     : ОАО НПК СПП
 * Проект       : КР ФЭЛП, ОКР "Эфемериды"
 * Год          : 2015
 */

#ifndef DAF_READER_H
#define	DAF_READER_H

#include <stdio.h>

#ifdef	__cplusplus
extern "C" {
#endif

typedef struct
{
  double *parameters_d; // double parameters of segment's summary
  int    *parameters_i; // integer parameters of segment's summary
  int     offset;       // offset to segment's data (in 8-byte units)
  int     length;       // length of segment's data (in 8-byte units)
  FILE   *file;         // file that the segment belongs to
} DafSegment;
  
enum
{
  DAF_FILE_SPK,
  DAF_FILE_PCK
};

/* Struct for access to DAF (Double Precision Array Files) */
typedef struct
{
  int  filetype; // DAF_FILE_***
  int  nd;       // number of double parameters in each segment's summary
  int  ni;       // number of integer parameters in each segment's summary
  
  DafSegment *segments;  // array of segments
  int nsegments;         // number of segments

  // array to store all double and integer parameters of all summaries
  // (the summaries themselves store only pointers to these arrays)
  double *all_parameters_d;
  int    *all_parameters_i;
  
  FILE   *file; // to close it when it is over
} DAF;

// error codes returned by dafReadFile function
enum
{
  DAF_ERROR_BAD_FORMAT  = -1,
  DAF_ERROR_UNSUPPORTED = -2,
  DAF_ERROR_NO_MEMORY   = -3
};

int dafReadFile (FILE *f, DAF *daf);

int dafSegmentReadRange (DafSegment *segment, int start, int length, double *result);
  
void dafFreeData (DAF *daf);

#ifdef	__cplusplus
}
#endif

#endif

