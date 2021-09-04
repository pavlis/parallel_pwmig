#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <string>
#include <math.h>
#include "pwmig/seispp/VelocityModel_1d.h"
#include "mspass/utility/MsPASSError.h"

namespace pwmig::seispp
{
using namespace std;
using namespace pwmig::seispp;
using namespace mspass::utility;

double VelocityModel_1d::getv(const double zin) const
{
	double dz;
        /* if above top, return top model value */
        if(zin<z[0]) return(v[0]);
	/* Assume at layer sigularity when within dz*this factor */
	const double dzeqfactor(0.01);
	for(int i=1;i<nlayers;++i)
	{
		if(zin<z[i])
		{
                        dz=zin-z[i-1];
			if(fabs((zin-z[i])/dz)<dzeqfactor)
				return(v[i]);
			else
				return(v[i-1]+dz*grad[i-1]);
		}
	}
	dz=zin-v[nlayers-1];
	return(v[nlayers-1]+dz*grad[nlayers-1]);
}
/* Read from a file constructor.  fname is the file name to be
read and form is a format name.  Currently supports two
names:  rbh - Herrmmann synthetics package format; and
plain - simple ascii depth, velocity pairs.
For plain the order is thickness,P,S.  i.e. three columns of
ascii numbers, free form with blank separators, layer thicknes
is in column 1, P velocity is in column 2, and S velocity
is in column 3.  (no gradients in this format allowed).
property must be "P" or "S".  IMPORTANT:  note the use of
layer THICKNESS not DEPTH. */
VelocityModel_1d::VelocityModel_1d(const std::string fname,
	const std::string form, const std::string property)
{
  const string base_error("VelocityModel_1d file constructor:  ");
	int i;

	ifstream input;

	input.open(fname.c_str(), ios::in);
	if(input.fail())
	  throw MsPASSError(base_error+"Cannot open file "+fname,ErrorSeverity::Invalid);
	char line[255];
	if(form=="rbh" || form=="plain")
	{
		// throw away the first few lines for rbh format
		if(form=="rbh")
		{
			for(i=0;i<12;++i) input.getline(line,255);
		}
		while(!input.eof())
		{
			double f1,f2,f3;
			double skipper;
			input >> f1;
			input >> f2;
			input >> f3;
      if(input.bad())
          throw MsPASSError(base_error+"read error.   Check data file="+fname,
					   ErrorSeverity::Invalid);
			if(form=="rbh")
				for(i=0;i<7;++i) input >> skipper;
			if(input.eof()) break;
			if(property=="P")
			{
				z.push_back(f1);
				v.push_back(f2);
			}
			else if(property == "S")
			{
				z.push_back(f1);
				v.push_back(f3);
			}
			else
			{
				input.close();
				throw MsPASSError(base_error+"Illegal property parameter = "+property,
					 ErrorSeverity::Invalid);
			}
			// all values are 0 for gradient here
			grad.push_back(0.0);
		}
		nlayers = z.size();
		// replace z values (currently intervals) with accumulated
		// depth to top of each layer
		for(i=1,z[0]=0.0;i<nlayers;++i)  z[i]+=z[i-1];
	}
  else if(form=="mod1d")
  {
		/* mod1d files can be obtained by subsetting or editing tables in
		an antelope database with the table name mod1d.
    We do a straight read of a mod1d table.  Require
    all rows have the same model name tag and are sorted
    by depth. */
            double zin,zlast,vin,gradin;
            string modname,mnlast,property_name;
            const string Ptag("Pvelocity");
            const string Stag("Svelocity");
            i=0;
            while(input.getline(line,255))
            {
                stringstream ss(line);
                ss>>modname;  ss>>property_name;
                ss>>zin;   ss>>vin;   ss>>gradin;
                if(input.bad())
								  throw MsPASSError(base_error+"read error.   Check data file="+fname);
                if(i==0)
                {
                    mnlast=modname;
                    zlast=zin;
                }
                if(mnlast!=modname)
                {
                    input.close();
										throw MsPASSError(base_error + "model name field changed from "
									     + mnlast + " to "+modname+"\nMust be constant for this format");
                }
                if(zin<zlast)
                {
                    input.close();
										throw MsPASSError(base_error + "Error in tabulated depths.  Data file="
									     + fname + " Has a depth that decreases or is constant\nCheck data file",
										   ErrorSeverity::Invalid);
                }
                if( ((property=="P") && (property_name==Ptag))
                    || ((property=="S") && (property_name==Stag)) )
                {
                    z.push_back(zin);
                    v.push_back(vin);
                    grad.push_back(gradin);
                    mnlast=modname;
                    zlast=zin;
                }
            }
            nlayers=z.size();
        }
	else
	{
		input.close();
		throw MsPASSError(base_error+"Unsupported format specified="+form,
			ErrorSeverity::Invalid);
	}
  input.close();
}
/* Standard copy constructor */
VelocityModel_1d::VelocityModel_1d(const VelocityModel_1d& old)
{
	nlayers=old.nlayers;
	grad=old.grad;
	v=old.v;
	z=old.z;
}
/* Standard assignment operator */
VelocityModel_1d& VelocityModel_1d::operator=(const VelocityModel_1d& old)
{
	if(this!=&old)
	{
		nlayers=old.nlayers;
		grad=old.grad;
		v=old.v;
		z=old.z;
	}
	return(*this);
}
} // Termination of namespace
