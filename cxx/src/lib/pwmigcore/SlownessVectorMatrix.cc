#include <vector>
#include "pwmig/pwmigcore/SlownessVectorMatrix.h"
#include "mspass/utility/MsPASSError.h"
namespace pwmig::pwmigcore {
using namespace pwmig::pwmigcore;
using namespace std;
using mspass::utility::MsPASSError;
using mspass::utility::ErrorSeverity;
using mspass::seismic::SlownessVector;

SlownessVectorMatrix::SlownessVectorMatrix()
{
    nrow=0;
    ncol=0;
}
SlownessVectorMatrix::SlownessVectorMatrix(const int nr, const int nc)
{
    nrow=nr;
    ncol=nc;
    int n=nr*nc;
    uarray.reserve(n);
    int i;
    for(i=0;i<n;++i)
        uarray.push_back(SlownessVector(0.0,0.0,0.0));
}
SlownessVectorMatrix::SlownessVectorMatrix(const SlownessVectorMatrix& parent)
                            : uarray(parent.uarray)
{
    nrow=parent.nrow;
    ncol=parent.ncol;
}
SlownessVector SlownessVectorMatrix::operator()(const int i, const int j) const
{
    if( (i<0) || (i>=nrow) || (j<0) || (j>=ncol) )
    {
        stringstream ss;
        ss << "SlownessVectorMatrix::operator(): row="<<i
            <<" column="<<j<<" is outside bounds of matrix size="
            << nrow <<"X"<<ncol<<endl;
        throw MsPASSError(ss.str(),ErrorSeverity::Fatal);
    }
    /* the original api returned a reference here.  That was dumb as a
    slowness vector is lightweight.  */
    return (this->uarray[i+nrow*j]);
}
SlownessVectorMatrix& SlownessVectorMatrix::operator=
                        (const SlownessVectorMatrix& parent)
{
    if(this!=&parent)
    {
        nrow=parent.nrow;
        ncol=parent.ncol;
        uarray=parent.uarray;
    }
    return *this;
}
/* Made a friend instead of part of object*/
ostream& operator<<(ostream& ostr, SlownessVectorMatrix& svm)
{
    ostr << "row index, column index, ux, uy"<<endl;
    int i,j;
    for(i=0;i<svm.nrow;++i)
        for(j=0;j<svm.ncol;++j)
        {
            SlownessVector uij=svm(i,j);
            ostr << i <<" "
                << j << " "
                << uij.ux << " "
                << uij.uy << endl;
        }
    return ostr;
}
}
