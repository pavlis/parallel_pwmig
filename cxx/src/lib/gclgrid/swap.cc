#include <typeinfo>
#include <sys/types.h>
#include <iostream>
#include "pwmig/gclgrid/swapbytes_pwmig.h"


namespace pwmig::gclgrid
{
using namespace std;
using namespace pwmig::gclgrid;
/* This is the old antelope routine to swap bit endian 32 bit vectors  of
   length n*/
void vectorswap4(int *in,int *out,int n)
{
    int32_t so;  //so is short for swapped output
    int i;
    for(i=0;i<n;++i)
    {
        so=__builtin_bswap32(static_cast<int32_t>(in[i]));
        out[i]=static_cast<int>(so);
    }
}
void md2hd(double *in, double *out,int n)
{
    int i;
    //cout << "Byte swapped double result - length="<<n<<endl;
    for(i=0;i<n;++i)
    {
        //cout << in[i]<<" ";
        out[i]=in[i];
        SwapEnd(out[i]);
        //cout << out[i]<<endl;
    }
}
void swap4(int *d)
{
    int32_t so;
    so=__builtin_bswap32(static_cast<int32_t>(*d));
    *d=static_cast<int>(so);
}
void swap4_float(float *d)
{
    /*
    int32_t so;
    so=__builtin_bswap32(static_cast<int32_t>(*d));
    *d=static_cast<float>(so);
    */
    SwapEnd(*d);
}
/*! \brief Test for little endian condition.

To handle mixed processors it is essential to know if the
word structure of this machine you are on is little or big
endian.  Intel processors are dominate today and are little
endian while Sun machines, which are commonly used in geophysics,
are big endian.  Because it is common today to mix these platforms
on the same network a way to detect which byte order the current
machine is, is necessary.

\return true if this processor is little endian (Intel byte order).
	Conversely returns false if the processor is big endian.
\author  Gary L. Pavlis with the core algorithm stolen from the
	University of Texas supercomputer web site.
*/
bool IntelByteOrder()
{
        long i = 0x11223344; unsigned char* c = (unsigned char*) &i;
        if(*c != 0x44)
                return(false);
        else
                return(true);
}
/*! \brief Architecture indedependent procedure
to byte swap a vector of doubles.

In the seispp library most data are stored internally as doubles.
External data representations, however, are subject to byte order
issues.  This routine will take a vector of doubles and automatically
swap bytes using a method appropriate for the parent architecture.
It should always be preceded by logic to decide if byte swapping
is necessary as this will always swap bytes one way or the other.

\param x pointer to array of doubles to be byte swapped.
\param nx number of elements in x.  This is quietly assumed
	to be correct and not bounds checking is done by this procedure.
*/

void swapdvec(double *x,int nx)
{
    int i;
    for(i=0;i<nx;++i)
        SwapEnd(x[i]);
}

} // End mspass::gclgrid namespace declaration
