#include <vector>
#include "mspass/utility/dmatrix.h"
#include "mspass/utility/ErrorLogger.h"
#include "mspass/seismic/Seismogram.h"
#include <boost/archive/basic_archive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>

/*! Spectial class used in map operation with the function migrate_one_seismogram.

This class is the return of the migrate_one_seismogram function in the MsPASS
adaption of pwmig.  It is really a special struct designed specially for that
function.  It likely has no use outside of that one.
*/
class PWMIGmigrated_seismogram
{
public:
  /*! Grid index 1 - x1 generalized coordinate index - of parent pseudostation grid.*/
  int ix1;
  /*! Grid index 1 - x1 generalized coordinate index - of parent pseudostation grid.*/
  int ix2;
  /*! True if the data are valid.  Errors or a dead input can make this false. */
  bool live;
  /*! Holds migrated 3 data components with sample value=parent raygrid index 3 value.*/
  mspass::utility::dmatrix migrated_data;
  /*! Holds n3 length vector of delta omega values for each point in migrated data*/
  std::vector<double> domega;
  /*! Holds computed GRT weight for each sample of migrated_data.*/
  std::vector<double> dweight;
  /*! Any problems are posted to the error log held by this class.*/
  mspass::utility::ErrorLogger elog;

  /*! Default constructor.  Do not use. */
  PWMIGmigrated_seismogram():elog(){ix1=-1; ix2=-1;};
  /*! Memory allocating and simple initializing constructor.

  This is the main constructor used in pwmig.  It sets ix1=i and ix2=j and
  allocates space for dmatrix of length n.  reserve is called on vectors for
  to preallooc to length n. */
  PWMIGmigrated_seismogram(i,j,n);
  /*! Standard copy constructor. */
  PWMIGmigrated_seismogram(const PWMIGmigrated_seismogram& parent);
  /*! Standard assignmnt operator.*/
  PWMIGmigrated_seismogram& operator=(const PWMIGmigrated_seismogram& parent);
  /*! Copy error log from a parent Seismogram.

  In pwmig error messages in inputs can get list in the parallel environment
  if they aren't copied from th parent.  This method is called in migrate_one_seismogram
  to push any existing messages in a parent to this ErrorLogger copy in order
  to preserve them. That happens even if the input data are live. */
  int copy_elog(const Seismogram& d)
  {
    this->elog=d.elog;
  };

private:
  friend boost::serialization::access;
  template<class Archive>
     void serialize(Archive& ar,const unsigned int version)
  {
    ar & ix1;
    ar & ix2;
    ar & live;
    ar & migrated_data;
    ar & domega;
    ar & dweight;
    ar & elog;
  };
};
