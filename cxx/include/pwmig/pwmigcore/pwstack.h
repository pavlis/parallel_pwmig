#ifndef _PWSTACK_H_
#define _PWSTACK_H_
#include <string>
#include "mspass/utility/Metadata.h"
#include "mspass/algorithms/Taper.h"  // defines TopMute
#include "mspass/seismic/Ensemble.h"
#include "pwmig/pwmigcore/RectangularSlownessGrid.h"
#include "pwmig/pwmigcore/DepthDependentAperture.h"
#include "pwmig/dsap/coords.h"
namespace pwmig::pwmigcore {
/*
Main processing function for pwstack algorithm.  Takes an input
data ensemble and produce a complete suite of plane-wave stacks
defined by the RectangularSlownessGrid object.
This function uses an enhancement from the original Neal and Pavlis
and Poppeliers and Pavlis papers.  It allows the aperture to
be time variable.  This is generalized, but the normal expectation
is that the apeture would grow wider with time to compensate
somewhat for diffraction.  Most applications will want to use the
fresnel form of the aperture definition (defined in DepthDependentAperture)
but specialized applications may need to alter the definition.  An
example is reducing the aperture in the upper crustal section to reduce
horizontal smoothing to enhance the Moho to a nominal station spacing -
most appropriate for computing only a ccp section.

\param indata - input raw data ensemble (see object definition
the enscapsulates all this requires).  Note the contents are altered by
top mutes.  Also for efficiency the indata ensmble should only be data
within the range defined by the maximum cutoff distance.   Note the
input is not declared const for a reason - the data are altered.
That is consistent with mspass version that does a db query to assemble
independent input for each call to this function.
\param ugrid - object to define slowness grid for stacking
\param mute - mute applied to data before stacking
\param stackmute - mute applied to data after stacking
(This stack is aligned relative to latest
mute time of raw data in stack.)
\param lat0, lon0 - pseudostation grid point (IN RADIANS)
\param ux0, uy0 - slowness vector of input data (stacking
is relative to this vector, but the total slowness is posted to each output)
\param tstart and tend - define time period for output stack
The routine pretty much assumes relative timing
so this is normally time wrt the start of
each trace in the ensemble.
\param aperture - defines variable aperture stack weighting
(see above)
\param aperture_taper_length - defines a smoother length at the edge of the
 aperture (units of km)
\param centroid_cutoff is an important parameter for defining data holes.
 Practical experience showed that with irregularly sampled experimental
 data this parameter can reduce artifacts from bleeding the data into
 areas of poor coverage (holes)
\param mdlcopy - defines metadata to copy from the raw data
\param save_history - if true the mspass ProcessingHistory is enabled to
  save what inputs were stacked to produce each output.
\param algid - is the algorithm id string to assign to ProcessingHistory defining
  this instance of pwstack.   This parameter is only used if save_history is
  set true.

\return ensemble of seismograms produced from input for the range of
 specified plane wave relative slowness offsets.   An empty ensemble marked
 dead denotes a problem that should be handled.   That normally means just
 saving the error message posted to the ensemble definibng problem in MongoDB.

Changed July 1, 2008
Now returns count of fold for this grid point.  negative or 0
means did nothing.

Throws a MetadataError exception object if there are problems
parsing required metadata from any trace.  Current caller will
abort the program on this condition, but evolution might want
to produce a handler.
Change Jan 15,2004
Used to pass lat0, lon0 by arg list, now passed through the
ensemble metadata.
Change August,2009
Added two parameters to make variable depth aperture work correctly:
centroid_cutoff is used with a centroid calculation.  If the centroid
of the stations retained for a pseudostation fall outside this distance from
the pseudostation point, return -2 and do nothing.  The other is
aperture_taper_length.  A linear taper like a top mute of this length
(in seconds) is applied to weights whenever the first nonzero value is
not at zero lag.

spring 2021:  major changes to adapt the code to MsPASS.  Too many
changes to list.  Only the core algorithm was retained.
*/
mspass::seismic::LoggingEnsemble<mspass::seismic::Seismogram> pwstack_ensemble
(mspass::seismic::LoggingEnsemble<mspass::seismic::Seismogram>& indata,
  const pwmig::pwmigcore::RectangularSlownessGrid& ugrid,
   mspass::algorithms::TopMute& mute,
    mspass::algorithms::TopMute& stackmute,
     int stack_count_cutoff,
      const double tstart,
       const double tend,
        const pwmig::pwmigcore::DepthDependentAperture& aperture,
         const double aperture_taper_length,
          const double centroid_cutoff,
           const mspass::utility::MetadataList& mdlcopy,
            const bool save_history,
             const std::string algid);
} // end namespace pwmig::pwmigcore
#endif
