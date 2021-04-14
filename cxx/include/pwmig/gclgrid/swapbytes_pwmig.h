#ifndef _GCLGRID_SWAPBYTES_H_
#define _GCLGRID_SWAPBYTES_H_
namespace pwmig::gclgrid
{
/* This is a local include file to define glp implementation of
   swapbyte routines in antelope.   Needed only for prototypes.*/
/*! Byte swap 32 bit integer vector of length n. */
void vectorswap4(int *in,int *out,int n);
/*! Byte swap a 32 bit integer in place */
void swap4(int *d);
/*! Byte swap a C float (32 bits) */
void swap4_float(float *d);
/*! Byte swap a double vector of length n and put result in out.
  Assumes in and out already allocated.  out can equal in */
void md2hd(double *in, double *out,int n);
/*! Returns true if this is a little endian machine */
bool IntelByteOrder();
/*! Byte swap a vector of doubles of length n in place */
void swapdvec(double *x,int nx);
/*! \brief Generic byte swap function for any simple type.

  This template is a generic method to byte swap standard numeric
  types in place.   It should only be used for integer and real
  numbers or results are unpredictable.   It reverses bytes in a
  range define by sizeof(T).

  Taken from:
http://stackoverflow.com/questions/105252/how-do-i-convert-between-big-endian-and-little-endian-values-in-c
*/
template <typename T>
void SwapEnd(T& var)
{
    T varSwapped;
    for(long i = 0; i < static_cast<long>(sizeof(var)); i++)
         ((char*)(&varSwapped))[sizeof(var) - 1 - i] = ((char*)(&var))[i];
    for(long i = 0; i < static_cast<long>(sizeof(var)); i++)
                          ((char*)(&var))[i] = ((char*)(&varSwapped))[i];
}
}
#endif
