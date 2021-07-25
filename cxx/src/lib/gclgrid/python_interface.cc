/*This set of methods are extensions added in 2021 MsPASS conversion as a
less than elegant way to do setters with pybind11.  A more elegant method
would use subscripting on arrays, but I chose to not attempt that.
*/
#include "sstream"
#include "pwmig/gclgrid/GCLgridError.h"
#include "pwmig/gclgrid/gclgrid.h"
namespace pwmig::gclgrid
{
using namespace std;
std::string range_error_message(const string objname, const string objmethod,
  const int i, const int j, const int n1, const int n2)
{
  stringstream ss;
  ss<<objname<<"::"<<objmethod<<":  Received illegal indices:  i="<<i
       << " j=" << j<<endl
       <<"Grid n1="<<n1<<" and n2="<<n2;
  return ss.str();
}
std::string range_error_message(const string objname, const string objmethod,
  const int i, const int j, const int k, const int n1, const int n2, const int n3)
{
  stringstream ss;
  ss<<objname<<"::"<<objmethod<<":  Received illegal indices:  i="<<i
       << " j=" << j<<" k="<<k<<endl
       <<"Grid n1="<<n1<<", n2="<<n2<<", and n3="<<n3;
  return ss.str();
}
void GCLgrid::set_coordinates(const Cartesian_point& p, const int i, const int j)
{
  if( (i<0) || (i>=this->n1) || (j<0) || (j>this->n2))
	{
    string message=range_error_message("GCLscalarfield","set_coordinates(Cartesian)",
       i,j,this->n1,this->n2);
		throw GCLgridError(message);
	}
  this->x1[i][j]=p.x1;
  this->x2[i][j]=p.x2;
  this->x3[i][j]=p.x3;
}
void GCLgrid::set_coordinates(const Geographic_point& gp, const int i, const int j)
{
  if( (i<0) || (i>=this->n1) || (j<0) || (j>this->n2))
	{
    string message=range_error_message("GCLscalarfield","set_coordinates(Geographic)",
       i,j,this->n1,this->n2);
		throw GCLgridError(message);
	}
  Cartesian_point p=this->gtoc(gp);
  this->x1[i][j]=p.x1;
  this->x2[i][j]=p.x2;
  this->x3[i][j]=p.x3;
}
void GCLscalarfield::set_value(const double newvalue, const int i, const int j)
{
	if( (i<0) || (i>=this->n1) || (j<0) || (j>this->n2))
	{
    string message=range_error_message("GCLscalarfield","set_value",i,j,this->n1,this->n2);
		throw GCLgridError(message);
	}
	this->val[i][j]=newvalue;
}
void GCLvectorfield::set_value(const double *newvals, const int i, const int j)
{
	if( (i<0) || (i>=this->n1) || (j<0) || (j>this->n2))
	{
    string message=range_error_message("GCLvectorfield","set_value",
       i,j,this->n1,this->n2);
		throw GCLgridError(message);
	}
  for(auto iv=0;iv<this->nv;++iv) this->val[i][j][iv]=newvals[iv];
}
void GCLgrid3d::set_coordinates(const Cartesian_point& p,
  const int i, const int j, const int k)
{
  if( (i<0) || (i>=this->n1) || (j<0) || (j>this->n2)|| (k<0) || (k>=this->n3))
	{
    string message=range_error_message("GCLgrid3d","set_coordinates(Cartesian)",i,j,k,
        this->n1,this->n2,this->n3);
		throw GCLgridError(message);
	}
  this->x1[i][j][k]=p.x1;
  this->x2[i][j][k]=p.x2;
  this->x3[i][j][k]=p.x3;
}
void GCLgrid3d::set_coordinates(const Geographic_point& gp,
  const int i, const int j, const int k)
{
  if( (i<0) || (i>=this->n1) || (j<0) || (j>this->n2)|| (k<0) || (k>=this->n3))
	{
    string message=range_error_message("GCLgrid3d","set_coordinates(Geographic)",i,j,k,
        this->n1,this->n2,this->n3);
		throw GCLgridError(message);
	}
  Cartesian_point p=this->gtoc(gp);
  this->x1[i][j][k]=p.x1;
  this->x2[i][j][k]=p.x2;
  this->x3[i][j][k]=p.x3;
}
void GCLscalarfield3d::set_value(const double newvalue,
  const int i, const int j, const int k)
{
	if( (i<0) || (i>=this->n1) || (j<0) || (j>this->n2)|| (k<0) || (k>=this->n3))
	{
    string message=range_error_message("GCLscalarfield3d","set_value",i,j,k,
        this->n1,this->n2,this->n3);
		throw GCLgridError(message);
	}
	this->val[i][j][k]=newvalue;
}
void GCLvectorfield3d::set_value(const double *newvals,
  const int i, const int j, const int k)
{
	if( (i<0) || (i>=this->n1) || (j<0) || (j>this->n2) || (k<0) || (k>=this->n3))
	{
    string message=range_error_message("GCLvectorfield3d","set_value",i,j,k,
        this->n1,this->n2,this->n3);
    throw GCLgridError(message);
	}
	for(auto iv=0;iv<this->nv;++iv) this->val[i][j][k][iv]=newvals[iv];
}
/* getters for coordinates and field values follow.*/
Cartesian_point GCLgrid::get_coordinates(const int i, const int j) const
{
  if( (i<0) || (i>=this->n1) || (j<0) || (j>this->n2) )
	{
    string message=range_error_message("GCLgrid","get_coordinates",i,j,
        this->n1,this->n2);
    throw GCLgridError(message);
	}
  Cartesian_point result;
  result.x1=this->x1[i][j];
  result.x2=this->x2[i][j];
  result.x3=this->x3[i][j];
  return result;
}
Cartesian_point GCLgrid3d::get_coordinates(const int i, const int j, const int k) const
{
  if( (i<0) || (i>=this->n1) || (j<0) || (j>this->n2) || (k<0) || (k>=this->n3))
	{
    string message=range_error_message("GCLgrid3d","get_coordinates",i,j,k,
        this->n1,this->n2,this->n3);
    throw GCLgridError(message);
	}
  Cartesian_point result;
  result.x1=this->x1[i][j][k];
  result.x2=this->x2[i][j][k];
  result.x3=this->x3[i][j][k];
  return result;
}
/*! Get the field value at specified grid cell */
double GCLscalarfield::get_value(const int i, const int j) const
{
  if( (i<0) || (i>=this->n1) || (j<0) || (j>this->n2) )
	{
    string message=range_error_message("GCLscalarfield","get_value",i,j,
        this->n1,this->n2);
    throw GCLgridError(message);
	}
  return this->val[i][j];
}
std::vector<double> GCLvectorfield::get_value(const int i, const int j) const
{
  if( (i<0) || (i>=this->n1) || (j<0) || (j>this->n2) )
  {
    string message=range_error_message("GCLvectorfield","get_value",i,j,
        this->n1,this->n2);
    throw GCLgridError(message);
  }
  std::vector<double> retvector;
  retvector.reserve(this->nv);
  for(auto iv=0;iv<this->nv;++iv) retvector.push_back(this->val[i][j][iv]);
  return retvector;
}
double GCLscalarfield3d::get_value(const int i, const int j, const int k) const
{
  if( (i<0) || (i>=this->n1) || (j<0) || (j>this->n2) || (k<0) || (k>=this->n3))
	{
    string message=range_error_message("GCLscalarfield3d","get_value",i,j,k,
        this->n1,this->n2,this->n3);
    throw GCLgridError(message);
	}
  return this->val[i][j][k];
}
std::vector<double> GCLvectorfield3d::get_value(const int i, const int j, const int k) const
{
  if( (i<0) || (i>=this->n1) || (j<0) || (j>this->n2) || (k<0) || (k>=this->n3))
	{
    string message=range_error_message("GCLvectorfield3d","get_value",i,j,k,
        this->n1,this->n2,this->n3);
    throw GCLgridError(message);
	}
  std::vector<double> retvector;
  retvector.reserve(this->nv);
  for(auto iv=0;iv<this->nv;++iv) retvector.push_back(this->val[i][j][k][iv]);
  return retvector;
}
} //End namespace
