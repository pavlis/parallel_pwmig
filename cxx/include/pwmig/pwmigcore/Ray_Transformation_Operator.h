#ifndef _RAY_TRANSFORMATION_OPERATOR_H_
#define _RAY_TRANSFORMATION_OPERATOR_H_
#include <vector>
#include "mspass/utility/dmatrix.h"
#include "pwmig/gclgrid/gclgrid.h"
namespace pwmig::pwmigcore
{
class Ray_Transformation_Operator
{
public:
  int npoints;
  Ray_Transformation_Operator(const int np);
  // constructor for constant azimuth
  Ray_Transformation_Operator(pwmig::gclgrid::GCLgrid& g,
    mspass::utility::dmatrix& path, double azimuth);
  // constructor for depth dependent operator
  Ray_Transformation_Operator(pwmig::gclgrid::GCLgrid& g,
    mspass::utility::dmatrix& path,double azimuth, mspass::utility::dmatrix& nup);
  Ray_Transformation_Operator(const Ray_Transformation_Operator& pat);
  Ray_Transformation_Operator& operator=(const Ray_Transformation_Operator& parent);

  mspass::utility::dmatrix apply(mspass::utility::dmatrix& in);
private:
  vector<mspass::utility::dmatrix> U;
};

} // End namespace
#endif
