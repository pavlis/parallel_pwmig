#ifndef _SLOWNESSVECTORMATRIX_H_
#define _SLOWNESSVECTORMATRIX_H_
#include <iostream>
#include <sstream>
#include <vector>
#include "mspass/seismic/SlownessVector.h"
#include "mspass/utility/MsPASSError.h"
namespace pwmig::pwmigcore
{
class SlownessVectorMatrix
{
public:
    /*! Default constructor - 0x0 */
    SlownessVectorMatrix();
    /*! \brief Initializing constructor fills with 0 vectors.

      This is the main constructor for this object.   Because
      of the somewhat odd use it is totally nonstandard in its
      behavior.   It does not build a complete object, but initializes
      a matrix of slowness vectors of size nrowxncol to all zeros.
      The operator()method should be used to fill the matrix with
      actual data.   This was done because the entire purpose of this
      object is to divorce pwmigcore routines from global travel
      time calculators.   This was a maintenance issue choice
      created by the lack of a clean solution for C++ to the dsap
      version of the taup library and not common alternative.  It was
      carried forward with modification to the mspass version of pwmig
      with additions to allow python to fill in the data.

      \param nrow - number of rows grid (matrix) of SlownessVectors
      \param ncol - number of columns grid (matrix) of SlownessVectors
      */
    SlownessVectorMatrix(const int nrow, const int ncol);
    /*! Copy constructor */
    SlownessVectorMatrix(const SlownessVectorMatrix& parent);
    /*! \brief return number of rows in this matrix */
    int rows() const
    {
        return nrow;
    };
    /*! return number of columns in this matrix */
    int columns() const
    {
        return ncol;
    };
    /*! \brief overloaded subscript operator.

      \param i is row index (C order)
      \param j is col index (C order)

      \exception - Will throw a SeisppError if the index is outside
         row and column range.
      */
    mspass::seismic::SlownessVector operator()(const int i, const int j) const;
    /*! \brief Putter to set the slowness vector at specified grid index.

    This was added for the mspass version to allow the use of taup calculators
    in python (specifically obpsy's implementation) to set the slowness vectors
    by grid index.  It throws a MsPASSError exception if the index is out of
    range.

    \param u slowness vector to set for i and j index   (i.e. this is the value
      stored at [i[j]])
    \param i first array index of grid
    \param j second array index of grid
    */
    void set_slowness(const mspass::seismic::SlownessVector& u,const int i, const int j);
    SlownessVectorMatrix& operator=(const SlownessVectorMatrix& parent);
    friend std::ostream& operator<<(std::ostream& ostr, SlownessVectorMatrix& svm);
private:
    int nrow;
    int ncol;
    std::vector<mspass::seismic::SlownessVector> uarray;
};
} // End namespace encapsulation
#endif
