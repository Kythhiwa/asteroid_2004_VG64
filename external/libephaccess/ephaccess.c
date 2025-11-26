/* 
 * Файл         : ephaccess.c
 * Исполнитель  : ИПА РАН
 * Заказчик     : ОАО НПК СПП
 * Проект       : КР ФЭЛП, ОКР "Эфемериды"
 * Год          : 2015
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>

#include "ephaccess.h"
#include "daf_reader.h"

#ifdef _MSC_VER
#define strcasecmp _stricmp
#define snprintf _snprintf
#endif

enum
{
  _MODIFIED_DIFFERENCE_ARRAY = 1,
  _POSITION_ONLY = 2,
  _VELOCITY_ONLY = 20
};

enum
{
  // limit for temporary arrays that are quickly allocated on the stack
  _MAX_DEGREE = 20,
  // limits for pre-allocated arrays inside EphAccess
  _MAX_FILES = 100,
};

const double SEC_IN_DAY = (86400.0);

typedef struct
{
  DafSegment *segment; // corresponding segment in a DAF file
  char filename[1024]; // corresponding file name (for diagnostic purposes)
  int filetype;        // copied from DAF
  int representation;  // see _POSITION_ONLY etc
  int object;          // object of interest
  int reference;       // reference object
  double jd0, jd1;     // initial date
  int nintervals;      // number of intervals
  double intlen;       // interval length in days
  int rsize;           // size of interval record
  int degree;          // degree of polynomial

  double dscale;       // distance scaling factor
  double tscale;       // velocity scaling factor
  
  int    cached_interval;
  double cache[(_MAX_DEGREE + 2) * 3];

  // MDA data
  int nrecords;
  double *record_epochs;
  double *record_epochs_each100;
  double initial_epoch;
  double final_epoch;
} _Theory;

#define _NELEM(arr) ((int)(sizeof(arr) / sizeof(*arr)))

struct tagEphAccess
{
  DAF dafs[_MAX_FILES];
  int ndafs;
  
  _Theory *theories;
  int ntheories_allocated;
  int ntheories;
  
  int    distance_units;
  
  double distance_scaling_factor;
  double time_scaling_factor;
  
  // some flags describing form of representation of Earth-Moon system
  // in the loaded ephemerides
  int have_earth_wrt_ssb;
  int have_moon_wrt_earth;
  int have_moon_wrt_emb;
  int have_earth_wrt_emb;
  int have_emb_wrt_ssb;
  
  double leftmost_jd;
  double rightmost_jd;
  
  int last_rc;
  char error[1200];
};

struct tagEphCache
{
  _Theory *theory;
  int      interval;
  double   coefs[(_MAX_DEGREE + 2) * 3]; 
};

int ephObjectByName (const char *name)
{
#define CHECK(id) if (strcasecmp(name, #id) == 0) return EPH_##id
  CHECK(SSB);
  CHECK(MERCURY);
  CHECK(VENUS);
  CHECK(EMB);
  CHECK(MARS_BC);
  CHECK(JUPITER_BC);
  CHECK(SATURN_BC);
  CHECK(URANUS_BC);
  CHECK(NEPTUNE_BC);
  CHECK(PLUTO_BC);
  CHECK(SUN);
  CHECK(MOON);
  CHECK(EARTH);
#undef CHECK 
  return -1;
}

CEXPORT EphAccess * ephCreate (void)
{
  EphAccess *eph = malloc(sizeof(EphAccess));
  
  if (eph != NULL)
  {
    memset(eph, 0, sizeof(EphAccess));
    eph->ntheories_allocated = 0;
    eph->theories = NULL;
    eph->distance_units = EPH_KM;
    eph->time_scaling_factor = SEC_IN_DAY;
    eph->distance_scaling_factor = 1.0;
    eph->leftmost_jd = 1e100;
    eph->rightmost_jd = -1e100;
  }
  
  return eph;
}

CEXPORT const char * ephLastError (EphAccess *eph)
{
  if (eph->last_rc == 0)
    return NULL;
  return eph->error;
}

/* Macros for reporting errors (not variadic thanks to Visual Studio) */
#define ERROR0(eph, rc, format) \
  do \
  {  \
    eph->last_rc = rc; \
    snprintf(eph->error, sizeof(eph->error), format); \
    return rc; \
  } while (0)
#define ERROR1(eph, rc, format, param1) \
  do \
  {  \
    eph->last_rc = rc; \
    snprintf(eph->error, sizeof(eph->error), format, param1); \
    return rc; \
  } while (0)
#define ERROR2(eph, rc, format, param1, param2) \
  do \
  {  \
    eph->last_rc = rc; \
    snprintf(eph->error, sizeof(eph->error), format, param1, param2); \
    return rc; \
  } while (0)
#define ERROR3(eph, rc, format, param1, param2, param3) \
  do \
  {  \
    eph->last_rc = rc; \
    snprintf(eph->error, sizeof(eph->error), format, param1, param2, param3); \
    return rc; \
  } while (0)

#define SUCCESS(eph) return (eph->last_rc = 0)

static int _loadDAF (EphAccess *eph, const char *filename, FILE *f)
{
  DAF *daf;
  int i;

  if (eph->ndafs >= _NELEM(eph->dafs))
    ERROR1(eph, EPH_ERROR_INTERNAL_LIMITS, "limit of DAF files (%d) exceeded",
          _NELEM(eph->dafs));
  
  daf = &eph->dafs[eph->ndafs];
  
  switch (dafReadFile(f, daf))
  {
    case DAF_ERROR_BAD_FORMAT:
    case DAF_ERROR_UNSUPPORTED:
      ERROR1(eph, EPH_ERROR_BAD_FILE_FORMAT, "bad format of file %s", filename);
    case DAF_ERROR_NO_MEMORY:
      ERROR0(eph, EPH_ERROR_OUT_OF_MEMORY, "out of memory");
    case 0:
      break;
    default:
      ERROR0(eph, EPH_ERROR_BAD_FILE_FORMAT, "unknown DAF reader error");;
  }

  if (eph->ntheories_allocated < eph->ntheories + daf->nsegments)
  {
    eph->theories = (_Theory *)realloc(eph->theories,
            sizeof(_Theory) * (eph->ntheories + daf->nsegments));
    if (eph->theories == NULL) 
      ERROR0(eph, EPH_ERROR_OUT_OF_MEMORY, "out of memory");
    
    eph->ntheories_allocated = eph->ntheories + daf->nsegments;
  }
  
  for (i = 0; i < daf->nsegments; i++)
  {
    _Theory *theory = &eph->theories[eph->ntheories + i];
    DafSegment *segment = &daf->segments[i];

    memset(theory, 0, sizeof(_Theory));
    
    if (daf->filetype == DAF_FILE_SPK)
    {
      theory->object = segment->parameters_i[0];
      theory->reference = segment->parameters_i[1];
      theory->representation = segment->parameters_i[3];
      theory->initial_epoch = segment->parameters_d[0];
      theory->final_epoch = segment->parameters_d[1];

      if (theory->object == EPH_MOON && theory->reference == EPH_EARTH)
        eph->have_moon_wrt_earth = 1;
      if (theory->object == EPH_MOON && theory->reference == EPH_EMB)
        eph->have_moon_wrt_emb = 1;
      if (theory->object == EPH_EARTH && theory->reference == EPH_EMB)
        eph->have_earth_wrt_emb = 1;
      if (theory->object == EPH_EARTH && theory->reference == EPH_SSB)
        eph->have_earth_wrt_ssb = 1;
      if (theory->object == EPH_EMB && theory->reference == EPH_SSB)
        eph->have_emb_wrt_ssb = 1;
    }
    else if (daf->filetype == DAF_FILE_PCK)
    {
      theory->object = segment->parameters_i[0];
      theory->representation = segment->parameters_i[2];
    }
    
    if (theory->representation == _POSITION_ONLY)
    {
      double params[4];
      int days;
      
      if (!dafSegmentReadRange(segment, segment->length - 4, 4, params))
        ERROR1(eph, EPH_ERROR_BAD_FILE_FORMAT,
                "can not read theory parameters from %s", filename);

      // params[0] is the starting point of the theory, in seconds since J2000.
      // Convert seconds to integer days (jd0) + fraction of day (jd1)
      days = (int)(params[0] / SEC_IN_DAY);
      theory->jd0 = 2451545.0 + days;
      theory->jd1 = (params[0] - days * SEC_IN_DAY) / SEC_IN_DAY;
      theory->intlen = params[1] / SEC_IN_DAY;
      theory->rsize = (int)params[2];
      theory->nintervals = (int)params[3];

      // check that RSIZE is 3N + 2
      if ((theory->rsize % 3) != 2)
        ERROR2(eph, EPH_ERROR_BAD_FILE_FORMAT,
              "bad RSIZE (%d) in %s", theory->rsize, filename);

      // polynomial degree is N-1
      theory->degree = (theory->rsize - 2) / 3 - 1;

      // DSCALE and TSCALE are not used in this type of ephemeris
      theory->dscale = 1.0;
      theory->tscale = 1.0;
    }
    else if (theory->representation == _VELOCITY_ONLY)
    {
      double params[7];
      
      if (!dafSegmentReadRange(segment, segment->length - 7, 7, params))
        ERROR1(eph, EPH_ERROR_BAD_FILE_FORMAT,
                "can not read theory parameters from %s", filename);

      theory->dscale = params[0];
      theory->tscale = params[1] / SEC_IN_DAY;
      theory->jd0    = (int)params[2];
      theory->jd1    = params[3];
      theory->intlen = params[4];
      theory->rsize  = (int)params[5];
      theory->nintervals = (int)params[6];

      // check that RSIZE is 3N
      if ((theory->rsize % 3) != 0)
        ERROR2(eph, EPH_ERROR_BAD_FILE_FORMAT,
              "bad RSIZE (%d) in %s", theory->rsize, filename);

      // polynomial degree is N-2
      theory->degree = theory->rsize / 3 - 2;
    }
    else if (theory->representation == _MODIFIED_DIFFERENCE_ARRAY)
    {
      double params[1];
      
      if (!dafSegmentReadRange(segment, segment->length - 1, 1, params))
        ERROR1(eph, EPH_ERROR_BAD_FILE_FORMAT,
                "can not read theory parameters from %s", filename);

      theory->nrecords = (int)params[0];
      theory->record_epochs = (double *)malloc(sizeof(double) * theory->nrecords);
      if (theory->record_epochs == NULL) 
	ERROR0(eph, EPH_ERROR_OUT_OF_MEMORY, "out of memory");
      if (!dafSegmentReadRange(segment, segment->length - 1 - theory->nrecords / 100 - theory->nrecords,
                               theory->nrecords, theory->record_epochs))
        ERROR1(eph, EPH_ERROR_BAD_FILE_FORMAT,
	       "can not read theory parameters from %s", filename);
      
      if (theory->nrecords >= 100)
      {
	theory->record_epochs_each100 = (double *)malloc(sizeof(double) * (theory->nrecords / 100));
	if (theory->record_epochs_each100 == NULL) 
	  ERROR0(eph, EPH_ERROR_OUT_OF_MEMORY, "out of memory");
	if (!dafSegmentReadRange(segment, segment->length - 1 - theory->nrecords / 100,
                                 theory->nrecords / 100, theory->record_epochs_each100))
	  ERROR1(eph, EPH_ERROR_BAD_FILE_FORMAT,
		 "can not read theory parameters from %s", filename);
      }
    }
    else
      ERROR2(eph, EPH_ERROR_UNSUPPORTED_FORMAT,
            "unsupported representation (%d) in file %s",
              theory->representation, filename);

    if (theory->degree > _MAX_DEGREE)
      ERROR2(eph, EPH_ERROR_INTERNAL_LIMITS,
              "polynomial degree limit (%d) exceeded in file %s",
              _MAX_DEGREE, filename);
    
    // initialize cache to some non-existing interval number
    theory->cached_interval = -1;
    
    // store the file name for diagnostics
    strncpy(theory->filename, filename, sizeof(theory->filename) - 1);
    
    // store the file type to tell SPK theories from PCK ones
    theory->filetype = daf->filetype;
    
    // store the segment for access
    theory->segment = segment;
    
    if (eph->leftmost_jd > theory->jd0 + theory->jd1)
      eph->leftmost_jd = theory->jd0 + theory->jd1;
    if (eph->rightmost_jd < theory->jd0 + theory->jd1 +
                            theory->intlen * theory->nintervals)
      eph->rightmost_jd = theory->jd0 + theory->jd1 +
                          theory->intlen * theory->nintervals;
  }
  
  eph->ndafs++;
  eph->ntheories += daf->nsegments;

  SUCCESS(eph);
}

CEXPORT int ephLoadFile (EphAccess *eph, const char *filename)
{
  FILE *f;
  char id[8];

  if (eph == NULL || filename == NULL)
    ERROR0(eph, EPH_ERROR_NULL_POINTER, "null pointer given");
  
  f = fopen(filename, "rb");
  
  if (f == NULL)
    ERROR2(eph, EPH_ERROR_OPENING_FILE, "error opening file %s: %s",
          filename, strerror(errno));
  
  // read the format identidier
  if (fread(id, 8, 1, f) != 1)
    ERROR1(eph, EPH_ERROR_READING_FILE, "error reading format ID from file %s",
          filename);

  fseek(f, 0, SEEK_SET);
  
  // check if it is a binary DAF file (PCK or SPK)
  if ((memcmp(id, "DAF/SPK", 7) == 0) ||
      (memcmp(id, "DAF/PCK", 7) == 0) ||
      (memcmp(id, "NAIF/DAF", 8) == 0))
    return _loadDAF(eph, filename, f);
  
  ERROR1(eph, EPH_ERROR_BAD_FILE_FORMAT, "format of file %s not recognized",
        filename);
}

/* Return value is the number of the interval in the theory to which the
 * date (date0+date1) belongs.
 * *t gets a value in the range [-1..1] corresponding to the position of
 * (date0+date1) in the found interval. */
static int _findInterval (_Theory *theory, double date0, double date1, double *t)
{
  // subtract big from big and small from small, to avoid loss of precision
  double diff = (date0 - theory->jd0) + (date1 - theory->jd1);
  int interval = (int)floor(diff / theory->intlen);
  double diff_interval = (date0 - theory->jd0 - interval * theory->intlen) +
                         (date1 - theory->jd1);

  *t = (diff_interval / theory->intlen) * 2.0 - 1.0;
  return interval;
}

/* Calculate values of Chebyshev polynomials.
 * t is in the interval [-1 .. +1] */
static void _calcChebyshevPolynomials (int ncoef, double t, double *polynomials)
{       
  int i;
  
  polynomials[0] = 1.0;	
  polynomials[1] = t; 

  for (i = 2; i <= ncoef; i++) 
    polynomials[i] = 2.0 * polynomials[i - 1] * t - polynomials[i - 2];
}

/* Calculate values of Chebyshev polynomials' antiderivatives.
 * t is in the interval [-1 .. +1] */
static void _calcChebyshevAntiderivatives (int ncoef, double t, double *polynomials,
                                    double *antiderivatives)
{
  double d;
  int i, j, flag;

  antiderivatives[0] = t;
  antiderivatives[1] = (polynomials[2] + polynomials[0]) * 0.25;
  for (i = 2; i < ncoef; i++) 
    antiderivatives[i] = 0.5 * (polynomials[i + 1] / (i + 1) -
                                polynomials[i - 1] / (i - 1));
  
  // Eliminate constant terms
  flag = 0; 
  d = 2.0 * t; 

  for (i = 3, j = 1; i < ncoef; i += 2, j++)
  { 
    d = 0.25 / j + 0.25 / (j + 1); 
    flag = 1 - flag; 
    if (flag)
      d = -d; 
    antiderivatives[i] += d;
  }
}

/* Calculate values of Chebyshev polynomials' derivatives.
 * t is in the interval [-1 .. +1] */
static void _calcChebyshevDerivatives (int ncoef, double t, double *polynomials,
                                double *derivatives)
{
  int i;
  
  derivatives[0] = 0.0;
  derivatives[1] = 1.0;

  for (i = 2; i < ncoef; i++)
    derivatives[i] = t * 2.0 * derivatives[i - 1] +
                         2.0 * polynomials[i - 1] - derivatives[i - 2];
}

// borrowed from calceph
struct SegType1_Record
{
    double tl;                  /*!< final epoch of record */
    double g[15];               /*!< stepsize function vector */
    double refpos[3];           /*!< reference position vector */
    double refvel[3];           /*!< reference velocity vector */
    double dt[15][3];           /*!< Modified divided difference arrays */
    int kqmax1;                 /*!< maximum integration order plus 1 */
    int kq[3];                  /*!< integration order array */
};

static void _interpolate_MDA(const struct SegType1_Record *Precord,
                            double TimeJD0, double Timediff, double *pos, double *vel)
{
    double delta = (TimeJD0 - Precord->tl) + Timediff;

    double tp = delta;

    int mpq2 = Precord->kqmax1 - 2;

    int ks = Precord->kqmax1 - 1;

    double fc[15];

    double wc[14];

    double w[18];

    int i, j, jx;

    fc[0] = 1.E0;
    for (j = 0; j < mpq2; j++)
    {
        fc[j + 1] = tp / Precord->g[j];
        wc[j] = delta / Precord->g[j];
        tp = delta + Precord->g[j];
    }

    /* inverse of 1/(j-0) */
    for (j = 0; j < Precord->kqmax1; j++)
    {
        w[j] = 1.E0 / ((double) (j + 1));
    }

    jx = 0;
    for (; ks >= 2; ks--)
    {
        int ks1 = ks - 1;

        jx++;
        for (j = 0; j < jx; j++)
        {
            w[j + ks] = fc[j + 1] * w[j + ks1] - wc[j] * w[j + ks];
        }
    }

    /* Compute position */
    for (i = 0; i < 3; i++)
    {
        double P_Sum = 0.E0;

        for (j = Precord->kq[i] - 1; j >= 0; j--)
        {
            P_Sum += Precord->dt[j][i] * w[j + ks];
        }
        pos[i] = Precord->refpos[i] + delta * (Precord->refvel[i] + delta * P_Sum);
    }

    /* Update w for velocity */
    if (vel != NULL)
    {
        int ks1 = ks - 1;

        for (j = 0; j < jx; j++)
        {
            w[j + ks] = fc[j + 1] * w[j + ks1] - wc[j] * w[j + ks];
        }
        ks--;

        /* Compute velocity */
        for (i = 0; i < 3; i++)
        {
            double V_Sum = 0.E0;

            for (j = Precord->kq[i] - 1; j >= 0; j--)
            {
                V_Sum += Precord->dt[j][i] * w[j + ks];
            }
            vel[i] = Precord->refvel[i] + delta * V_Sum;
            vel[i] *= 86400;
        }
    }
}


/* Calculate raw ephemeride values by given theory and date (MDA type) */
static int _calculateByTheory_MDA (EphAccess *eph, _Theory *theory,
				   double date0, double date1, double *pos, double *vel,
				   int scale_distance)
{
  double t = ((date0 - 2451545.0) + date1) * 86400;
  int n, n100;

  for (n100 = 0; n100 < (theory->nrecords / 100); n100++)
    if (theory->record_epochs_each100[n100] > t)
      break;

  for (n = n100 * 100; n < (n100 + 1) * 100 && n < theory->nrecords; n++)
    if (theory->record_epochs[n] > t)
      break;

  double data[71];
  if (!dafSegmentReadRange(theory->segment, 71 * n, 71, data))
    ERROR1(eph, EPH_ERROR_BAD_FILE_FORMAT,
           "can not read coefficients from %s", theory->filename);

 
  struct SegType1_Record arecord;
  int k, p;

  // borrowed from calceph
  arecord.tl = data[0];
  for (k = 0; k < 15; k++)
    arecord.g[k] = data[1 + k];
  for (k = 0; k < 3; k++)
    arecord.refpos[k] = data[16 + 2 * k];
  for (k = 0; k < 3; k++)
    arecord.refvel[k] = data[17 + 2 * k];
  for (p = 0; p < 15; p++)
  {
    for (k = 0; k < 3; k++)
      arecord.dt[p][k] = data[22 + k * 15 + p];
  }
  arecord.kqmax1 = (int) data[67];
  for (k = 0; k < 3; k++)
    arecord.kq[k] = (int) data[68 + k];

  _interpolate_MDA(&arecord, (date0 - 2451545.0) * 86400, date1 * 86400, pos, vel);
  SUCCESS(eph);
}


/* Calculate raw ephemeride values by given theory and date */
static int _calculateByTheory (EphAccess *eph, _Theory *theory,
                               double date0, double date1, double *pos, double *vel,
                               int scale_distance)
{
  if (theory->representation == _MODIFIED_DIFFERENCE_ARRAY)
    return _calculateByTheory_MDA(eph, theory, date0, date1, pos, vel, scale_distance);
      
  int degree = theory->degree;
  double t;
  int interval = _findInterval(theory, date0, date1, &t);
  int j;
  double polynomials[_MAX_DEGREE + 1];
  double *coefs = theory->cache;
  int to_read;

  if (interval < 0 || interval >= theory->nintervals)
    ERROR3(eph, EPH_ERROR_DATE_OUT_OF_RANGE,
          "date %lf+%lf is out of range for file %s",
          date0, date1, theory->filename);
  
  to_read = (theory->cached_interval != interval);
  
  if (theory->representation == _POSITION_ONLY)
  {
    double dscale = eph->distance_scaling_factor;
    
    _calcChebyshevPolynomials(degree + 1, t, polynomials);

    if (to_read)
      if (!dafSegmentReadRange(theory->segment, theory->rsize * interval + 2,
     			     theory->rsize - 2, coefs))
        ERROR1(eph, EPH_ERROR_BAD_FILE_FORMAT,
              "can not read coefficients from %s", theory->filename);

    pos[0] = 0.0;
    pos[1] = 0.0;
    pos[2] = 0.0;

    for (j = degree; j >= 0; j--)
    {
      pos[0] += polynomials[j] * coefs[j];
      pos[1] += polynomials[j] * coefs[j + degree + 1];
      pos[2] += polynomials[j] * coefs[j + (degree + 1) * 2];
    }

    if (scale_distance)
    {
      pos[0] = pos[0] * dscale;
      pos[1] = pos[1] * dscale;
      pos[2] = pos[2] * dscale;
    }

    if (vel != NULL)
    {
      double derivatives[_MAX_DEGREE + 1];
      
      _calcChebyshevDerivatives(degree + 1, t, polynomials, derivatives);
      
      vel[0] = 0.0;
      vel[1] = 0.0;
      vel[2] = 0.0;
      
      for (j = degree; j >= 0; j--)
      {
        vel[0] += derivatives[j] * coefs[j];
	      vel[1] += derivatives[j] * coefs[j + degree + 1];
	      vel[2] += derivatives[j] * coefs[j + (degree + 1) * 2];
      }
      vel[0] /= (0.5 * theory->intlen);
      vel[1] /= (0.5 * theory->intlen);
      vel[2] /= (0.5 * theory->intlen);

      if (scale_distance)
      {
        vel[0] = vel[0] * dscale;
        vel[1] = vel[1] * dscale;
        vel[2] = vel[2] * dscale;
      }
    }
  }
  else if (theory->representation == _VELOCITY_ONLY)
  {
    double antiderivatives[_MAX_DEGREE + 1];
    double dscale = theory->dscale * eph->distance_scaling_factor;
    
    if (scale_distance &&
        eph->distance_units == EPH_AU &&
        fabs(theory->dscale - 149597870.7) < 1000.0)
      // theory has its own AU value, keep it as is
      dscale = 1.0;
    
    _calcChebyshevPolynomials(degree + 2, t, polynomials);
    _calcChebyshevAntiderivatives(degree + 1, t, polynomials, antiderivatives);

    if (to_read)
      if (!dafSegmentReadRange(theory->segment, theory->rsize * interval,
              theory->rsize, coefs))
        ERROR1(eph, EPH_ERROR_BAD_FILE_FORMAT,
              "can not read coefficients from %s", theory->filename);

    pos[0] = 0.0;
    pos[1] = 0.0;
    pos[2] = 0.0;
    
    for (j = degree; j >= 0; j--)
    {
      pos[0] += antiderivatives[j] * coefs[j];
      pos[1] += antiderivatives[j] * coefs[j + degree + 2];
      pos[2] += antiderivatives[j] * coefs[j + (degree + 2) * 2];
    }
    pos[0] = 0.5 * theory->intlen * pos[0] + coefs[degree + 1];
    pos[1] = 0.5 * theory->intlen * pos[1] + coefs[degree + 1 + degree + 2];
    pos[2] = 0.5 * theory->intlen * pos[2] + coefs[degree + 1 + (degree + 2) * 2];

    if (scale_distance)
    {
      pos[0] = pos[0] * dscale;
      pos[1] = pos[1] * dscale;
      pos[2] = pos[2] * dscale;
    }
    
    if (vel != NULL)
    {
      vel[0] = 0.0;
      vel[1] = 0.0;
      vel[2] = 0.0;
    
      for (j = theory->degree; j >= 0; j--)
      {
        vel[0] += polynomials[j] * coefs[j];
        vel[1] += polynomials[j] * coefs[j + degree + 2];
        vel[2] += polynomials[j] * coefs[j + (degree + 2) * 2];
      }

      if (scale_distance)
      {
        vel[0] = vel[0] * dscale / theory->tscale;
        vel[1] = vel[1] * dscale / theory->tscale;
        vel[2] = vel[2] * dscale / theory->tscale;
      }
    }
  }
  
  // set cached interval number to re-use coefficients next time
  theory->cached_interval = interval;
  SUCCESS(eph);
}

/* Forward declaration */
static int _calculateRectangular (EphAccess *eph, int object, int reference,
            double date0, double date1, double *xyz, double *vxyz);

/* Linear combination of two ephemerides */
static int _combineTwoEphemerides (EphAccess *eph, int object1, int reference1,
            int object2, int reference2, double coef1, double coef2,
            double date0, double date1, double *xyz, double *vxyz)
{
  double posvel1[6], posvel2[6];
  int i, rc;
  
  rc = _calculateRectangular(eph, object1, reference1, date0, date1, posvel1,
                             (vxyz == NULL) ? NULL : posvel1 + 3);
  if (rc != 0) return rc;

  rc = _calculateRectangular(eph, object2, reference2, date0, date1, posvel2,
                             (vxyz == NULL) ? NULL : posvel2 + 3);
  if (rc != 0) return rc;

  for (i = 0; i < 3; i++)
  {
    xyz[i] = posvel1[i] * coef1 + posvel2[i] * coef2;
    if (vxyz != NULL)
      vxyz[i] = posvel1[i + 3] * coef1 + posvel2[i + 3] * coef2;
  }

  SUCCESS(eph);
}

/* Linear combination of three ephemerides */
static int _combineThreeEphemerides (EphAccess *eph, int object1, int reference1,
            int object2, int reference2, 
            int object3, int reference3, 
            double coef1, double coef2, double coef3,
            double date0, double date1, double *xyz, double *vxyz)
{
  double posvel1[6], posvel2[6], posvel3[6];
  int i, rc;
  
  rc = _calculateRectangular(eph, object1, reference1, date0, date1, posvel1,
                             (vxyz == NULL) ? NULL : posvel1 + 3);
  if (rc != 0) return rc;

  rc = _calculateRectangular(eph, object2, reference2, date0, date1, posvel2,
                             (vxyz == NULL) ? NULL : posvel2 + 3);
  if (rc != 0) return rc;

  rc = _calculateRectangular(eph, object3, reference3, date0, date1, posvel3,
                             (vxyz == NULL) ? NULL : posvel3 + 3);
  if (rc != 0) return rc;
  
  for (i = 0; i < 3; i++)
  {
    xyz[i] = posvel1[i] * coef1 + posvel2[i] * coef2 + posvel3[i] * coef3;
    if (vxyz != NULL)
      vxyz[i] = posvel1[i + 3] * coef1 + posvel2[i + 3] * coef2 +
                posvel3[i + 3] * coef3;
  }

  SUCCESS(eph);
}

/* This function implements the logic of selection of ephemerides depending
 * of object-reference pair.  */
static int _calculateRectangular (EphAccess *eph, int object, int reference,
            double date0, double date1, double *xyz, double *vxyz)
{
  int i;
  int match = 0;

  if (eph == NULL || xyz == NULL)
    ERROR0(eph, EPH_ERROR_NULL_POINTER, "null pointer given");

  if (object == reference)
  {
    // zero result
    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    if (vxyz != NULL)
    {
      vxyz[0] = 0.0;
      vxyz[1] = 0.0;
      vxyz[2] = 0.0;
    }
    SUCCESS(eph);
  }

  double seconds_since_j2000 = ((date0 - 2451545.0) + date1) * 86400.0;
  
  // Search for matching segment
  for (i = 0; i < eph->ntheories; i++)
  {
    if (eph->theories[i].object == object &&
        eph->theories[i].reference == reference)
    {
      match = 1;
      
      if (seconds_since_j2000 >= eph->theories[i].initial_epoch &&
          seconds_since_j2000 <= eph->theories[i].final_epoch)
        break;
    }
  }

  if (i != eph->ntheories)
    return _calculateByTheory(eph, &eph->theories[i], date0, date1, xyz, vxyz, 1);

  if (match)
    ERROR2(eph, EPH_ERROR_DATE_OUT_OF_RANGE,
           "given date is out of bounds for object %d and reference %d", object, reference);
  
  // Theory not found -- try to work around
  
  if (object == EPH_EARTH && reference == EPH_SSB)
  {
    if (eph->have_earth_wrt_emb && eph->have_emb_wrt_ssb)
      return _combineTwoEphemerides(eph, EPH_EARTH, EPH_EMB, EPH_EMB, EPH_SSB,
                                    1.0, 1.0, date0, date1, xyz, vxyz);
  }
  else if (object == EPH_MOON && reference == EPH_SSB)
  {
    if (eph->have_moon_wrt_earth && eph->have_earth_wrt_ssb)
      return _combineTwoEphemerides(eph, EPH_MOON, EPH_EARTH, EPH_EARTH,
              EPH_SSB, 1.0, 1.0, date0, date1, xyz, vxyz);
    if (eph->have_moon_wrt_emb && eph->have_emb_wrt_ssb)
      return _combineTwoEphemerides(eph, EPH_MOON, EPH_EMB, EPH_EMB,
              EPH_SSB, 1.0, 1.0, date0, date1, xyz, vxyz);
    if (eph->have_moon_wrt_earth && eph->have_earth_wrt_emb && eph->have_emb_wrt_ssb)
      return _combineThreeEphemerides(eph, EPH_MOON, EPH_EARTH, EPH_EARTH,
              EPH_EMB, EPH_EMB, EPH_SSB, 1.0, 1.0, 1.0, date0, date1, xyz, vxyz);
  }
  else if (object == EPH_MOON && reference == EPH_EARTH)
  {
    if (eph->have_moon_wrt_emb && eph->have_earth_wrt_emb)
      return _combineTwoEphemerides(eph, EPH_MOON, EPH_EMB, EPH_EARTH, EPH_EMB,
                                    1.0, -1.0, date0, date1, xyz, vxyz);
  }
  else if (((reference >= EPH_MERCURY && reference <= EPH_SUN) && object == EPH_SSB) ||
           (((reference == EPH_MOON) || (reference == EPH_EARTH)) && object == EPH_SSB) ||
           (object == EPH_EARTH && reference == EPH_MOON))
  {
    // swap object and reference
    int rc = _calculateRectangular(eph, reference, object, date0, date1,
            xyz, vxyz);
    if (rc != 0)
      return rc;
    
    // change the sign of the obtained coordinates
    xyz[0] = -xyz[0];
    xyz[1] = -xyz[1];
    xyz[2] = -xyz[2];
    if (vxyz != NULL)
    {
      vxyz[0] = -vxyz[0];
      vxyz[1] = -vxyz[1];
      vxyz[2] = -vxyz[2];
    }
    SUCCESS(eph);
  }
  else if (reference != EPH_SSB && object != EPH_SSB)
  {
    // calculate both objects w.r.t. SSB and subtract
    double obj[6], ref[6];
    int i, rc1, rc2;
    
    rc1 = _calculateRectangular(eph, object, EPH_SSB, date0, date1, obj,
                                obj + 3);
    rc2 = _calculateRectangular(eph, reference, EPH_SSB, date0, date1, ref,
                                ref + 3);
    if (rc1 == 0 && rc2 == 0)
    {
      for (i = 0; i < 3; i++)
      {
        xyz[i] = obj[i] - ref[i];
        if (vxyz != 0)
          vxyz[i] = obj[i + 3] - ref[i + 3];
      }
      
      SUCCESS(eph);
    }
  }

  ERROR2(eph, EPH_ERROR_NOT_FOUND,
         "theory for object %d and reference %d not found", object, reference);
}

/* Calculate rectangular ephemerides and apply scaling factors */
CEXPORT int ephCalculateRectangular (EphAccess *eph, int object, int reference,
            double date0, double date1, double *xyz, double *vxyz)
{
  int rc = _calculateRectangular(eph, object, reference, date0, date1, xyz, vxyz);
  
  if (rc != 0)
    return rc;
  
  // apply time scaling factor to velocities
  // (discance scaling was already done in _calculateByTheory)
  if (vxyz != NULL)
  {
    vxyz[0] /= eph->time_scaling_factor;
    vxyz[1] /= eph->time_scaling_factor;
    vxyz[2] /= eph->time_scaling_factor;
  }
  return 0;
}

/* Calculate angular ephemerides and apply scaling factor */
CEXPORT int ephCalculateEulerAngles (EphAccess *eph, int frame,
            double date0, double date1, double *angles, double *rates)
{
  int i, rc;
  int single_frame = -1;
  int i_single_frame = -1;
  
  if (eph == NULL || angles == NULL)
    ERROR0(eph, EPH_ERROR_NULL_POINTER, "null pointer given");
  
  // Search for matching segment
  for (i = 0; i < eph->ntheories; i++)
  {
    // Check that given (JD0+JD1) minus frame's starting (JD0+JD1)
    // is within [0, length of segment]
    double diff = (date0 - eph->theories[i].jd0) + (date1 - eph->theories[i].jd1);
    if (diff < 0 || diff > eph->theories[i].nintervals * eph->theories[i].intlen)
      continue;
    
    if (eph->theories[i].object == frame)
      break;
    
    // check if we have a single PCK theory
    if (eph->theories[i].filetype == DAF_FILE_PCK)
    {
      if (single_frame == -1)
      {
        single_frame = eph->theories[i].object;
        i_single_frame = i;
      }
      else if (single_frame != eph->theories[i].object)
        single_frame = -2;
    }
  }

  if (i == eph->ntheories)
  {
    // check if we have frame = 0 and just one PCK theory
    if (frame == 0 && single_frame >= 0)
      i = i_single_frame;
    else
      ERROR1(eph, EPH_ERROR_NOT_FOUND, "theory for frame %d not found", frame);
  }
  
  rc = _calculateByTheory(eph, &eph->theories[i], date0, date1, angles, rates, 0);
  
  if (rc != 0)
    return rc;
  
  // apply scaling factor to rates
  if (rates != NULL)
  {
    rates[0] /= eph->time_scaling_factor;
    rates[1] /= eph->time_scaling_factor;
    rates[2] /= eph->time_scaling_factor;
  }
  return 0;
}

/* Calculate time ephemeride and apply scaling factor */
CEXPORT int ephCalculateTimeDiff (EphAccess *eph, int code, double date0, double date1,
                          double *diff)
{
  int i;
  double pos[3];
  int match = 0;
  
  if (eph == NULL || diff == NULL)
    ERROR0(eph, EPH_ERROR_NULL_POINTER, "null pointer given");

  // Search for matching segment
  for (i = 0; i < eph->ntheories; i++)
  {
    if (eph->theories[i].object == code)
    {
      match = 1;
      
      // Check that given (JD0+JD1) minus frame's starting (JD0+JD1)
      // is within [0, length of segment]
      double diff = (date0 - eph->theories[i].jd0) +
                    (date1 - eph->theories[i].jd1);
      if (diff < 0 ||
          diff > eph->theories[i].nintervals * eph->theories[i].intlen)
        continue;
    
      break;
    }
  }

  if (i == eph->ntheories)
  {
    if (match)
      ERROR1(eph, EPH_ERROR_DATE_OUT_OF_RANGE,
             "date is out of bounds for theory for time difference %d", code);
    else
      ERROR1(eph, EPH_ERROR_NOT_FOUND,
             "theory for time difference %d not found", code);
  }
  
  {
    int rc = _calculateByTheory(eph, &eph->theories[i], date0, date1, pos, NULL, 0);
    if (rc != 0)
      return rc;
  }
  
  // apply scaling factor to time difference
  // (TT-TDB is stored in seconds, so we divide by SEC_IN_DAY to scale to days,
  // and then apply our scaling factor)
  *diff = pos[0] * eph->time_scaling_factor / SEC_IN_DAY;
  SUCCESS(eph);
}

CEXPORT int ephSetDistanceUnits (EphAccess *eph, int units)
{
  if (units == EPH_KM)
    eph->distance_scaling_factor = 1.0; // the SPK files are already in km
  else if (units == EPH_AU)
    eph->distance_scaling_factor = 1.0 / 149597870.7;
  else
    ERROR1(eph, EPH_ERROR_BAD_INPUT, "unknown distance units: %d", units);
  eph->distance_units = units;
  SUCCESS(eph);
}

CEXPORT int ephSetTimeUnits (EphAccess *eph, int units)
{
  if (units == EPH_DAY)
    // internal units of SPK/PCK are days
    eph->time_scaling_factor = 1.0;
  else if (units == EPH_SEC)
    eph->time_scaling_factor = SEC_IN_DAY;
     // no nothing as the SPK files are already in sec
  else
    ERROR1(eph, EPH_ERROR_BAD_INPUT, "unknown time units: %d", units);
  SUCCESS(eph);
}

CEXPORT double ephGetLeftmostJD (EphAccess *eph)
{
  return eph->leftmost_jd;
}

CEXPORT double ephGetRightmostJD (EphAccess *eph)
{
  return eph->rightmost_jd;
}

CEXPORT void ephDestroy (EphAccess *eph)
{
  int i;
  
  if (eph == NULL)
    return;
  
  for (i = 0; i < eph->ndafs; i++)
    dafFreeData(&eph->dafs[i]);
  free(eph->theories);
  free(eph);
}
