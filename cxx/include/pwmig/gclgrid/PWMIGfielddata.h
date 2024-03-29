#include "mspass/utility/ErrorLogger.h"
#include "pwmig/gclgrid/gclgrid.h"
#include "pwmig/pwmigcore/PWMIGmigrated_seismogram.h"

namespace pwmig::gclgrid
{
/*! \brief Special extension of GCLvectorfield for PwMIG output.

This is an extension of GCLvectorfield used only in the mspass version of
pwmig.  It is more specialized in that it always has 5-vectors as the field
grid data.  It extends its base in two ways.  First, it adds an accumulate
method.  That method does a lot more than just copy the contents of one "ray"
datum.  (see documentation of accumulate for what is done).  Second, it adds an
ErrorLogger component (symbol elog) that is the sum total of all error messages
to all inputs received through the accumulate function.
*/
class PWMIGfielddata : public pwmig::gclgrid::GCLvectorfield3d
{
public:
  mspass::utility::ErrorLogger elog;
  PWMIGfielddata();
  PWMIGfielddata(const pwmig::gclgrid::GCLgrid3d& g);
  PWMIGfielddata(const PWMIGfielddata& parent);
  PWMIGfielddata& operator=(const PWMIGfielddata& parent);
  /* this constructor probably does not need to be bound to python,  It is
  a tiny wrapper made to mesh with the pickle interface in pybind11.
  The args are the 3 grid dimensions and an ErrorLogger object from
  deserialization inside unpickling pybind11 function. */
  PWMIGfielddata(const int n1, const int n2, const int n3,
    const mspass::utility::ErrorLogger elog_to_clone)
      : GCLvectorfield3d(n1,n2,n3,5),elog(elog_to_clone){};

  /*! Accumulate the data encapsulated in one PWMIGmigrated_seismogram object.

  The inner loop of pwmig constructs a set of values it interpolates onto a
  raygrid that is assumed to be the same GCLgrid3d object used to create this
  object.   This method inserts that data into the grid to create the same
  5-vector structure as the original pwmig.  That means two things.  First,
  the data vector received in d is reversed so 0 becomes n-1, 1 becomes n-2,
  ...,n-1 becomes 0.  Second the migrated data values (0,1, and 2 components)
  are scaled by dweight*domega.  The latter is essential to get the relative
  amplitudes correct for the final plane wave sum.

  The other thing this function does is check the ErrorLogger contents of d
  and if it is not empty it is copied with a prefix id string and reposted
  to this object's elog.

  If the input d is marked dead only error log entries are copied and any
  data contained in d is ignored.
1
  The position the data are inserted (first two indices of the grid) are
  pulled from d.ix1. and d.ix2

  \exception - can throw a MsPASSError exception if d is marked live but has
    a size inconsistent with this object's internal n3 dimension.

  */
  void accumulate(const pwmig::pwmigcore::PWMIGmigrated_seismogram& d);
};

} // End namespace
