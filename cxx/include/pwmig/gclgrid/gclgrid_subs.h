#ifndef _GCLGRID_SUBS_H_
#define _GCLGRID_SUBS_H_
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <sstream>
#include "mspass/utility/Metadata.h"
#include "mspass/utility/SphericalCoordinate.h"
#include "pwmig/gclgrid/GCLgridError.h"
#include "pwmig/gclgrid/swapbytes_pwmig.h"
#include "pwmig/gclgrid/gclgrid.h"
namespace pwmig::gclgrid
{
bool byte_swap_is_needed(const mspass::utility::Metadata& md);

template <class T>
    void pfload_common_GCL_attributes(T& g,const mspass::utility::Metadata& par)
{
    try {
	    /* This template loads common attributes from BasicGCLgrid base class*/
			g.name=par.get_string("name");
			g.lat0=par.get_double("origin_latitude");
			g.lon0=par.get_double("origin_longitude");
			// Immediately convert these to radians
			g.lat0=mspass::utility::rad(g.lat0);
			g.lon0=mspass::utility::rad(g.lon0);
			g.r0=par.get_double("origin_radius");
      g.azimuth_y=par.get_double("azimuth_y");
      g.azimuth_y=mspass::utility::rad(g.azimuth_y);
      g.dx1_nom=par.get_double("dx1_nom");
      g.dx2_nom=par.get_double("dx2_nom");
			g.n1=par.get_int("n1");
			g.n2=par.get_int("n2");
			g.i0=par.get_int("i0");
			g.j0=par.get_int("j0");
      /* There perhaps should be a way to force these to be
           computed, but for now we assume they were set by
           a writer and we don't need to compute them. */
      g.x1low=par.get_double("x1low");
      g.x1high=par.get_double("x1high");
      g.x2low=par.get_double("x2low");
      g.x2high=par.get_double("x2high");
      g.x3low=par.get_double("x3low");
      g.x3high=par.get_double("x3high");
      /* This perhaps should be set by caller, but since all
           callers will need this do it here.*/
      g.set_transformation_matrix();
    } catch(...) {throw;};
}

template <class T>
    void pfload_3dgrid_attributes(T& g, const mspass::utility::Metadata& par)
{
    try {
        g.dx3_nom=par.get_double("dx3_nom");
        g.n3=par.get_int("n3");
        g.k0=par.get_int("k0");
    } catch(...) {throw;};
}

template <typename GCLtype> void read_GCL2d_coord_arrays(GCLtype& d, const mspass::utility::Metadata& md)
{
	try{
		std::string base_error("read_GCL2d_coord_arrays:  ");
		pfload_common_GCL_attributes(d,md);
		bool need_to_swap_bytes;
		need_to_swap_bytes=pwmig::gclgrid::byte_swap_is_needed(md);
		std::string fname,dir,dfile,dfileext;
		dir=md.get_string("dir");
		dfile=md.get_string("grid_data_file");
		fname=dir+"/"+dfile;
		if(md.is_defined("grid_data_file_extension"))
		{
			dfileext=md.get_string("grid_data_file_extension");
			fname+=".";
			fname+=dfileext;
		}
	  FILE *fp = fopen(fname.c_str(),"r");
	  if(fp == NULL)
		  throw GCLgridError(base_error
						+ "fopen failed on file "
						+ dfile);
    d.x1 = pwmig::gclgrid::create_2dgrid_contiguous(d.n1,d.n2);
	  d.x2 = pwmig::gclgrid::create_2dgrid_contiguous(d.n1,d.n2);
	  d.x3 = pwmig::gclgrid::create_2dgrid_contiguous(d.n1,d.n2);
	  /* Database stores data in geographic coordinates.
	  File stores Cartesian form.  A bit inconsistent,
	  but a design choice.  Storing geo coordinates
	  would be a future format choice.*/
	  int gridsize = d.n1*d.n2;
		if(md.is_defined("grid_data_foff"))
		{
			long foff;
			foff=md.get_long("grid_data_foff");
			if(fseek(fp,foff,SEEK_SET))
			{
				stringstream ss;
				ss << base_error
				  << "fseek to offset="<<foff<<" failed for file="<<dfile;
				throw GCLgridError(ss.str());
			}
		}
	  if(fread(d.x1[0],sizeof(double),gridsize,fp) != gridsize)
	  {
		  fclose(fp);
		  throw GCLgridError(base_error
						+ "fread failed on file reading x1 coordinate  array"
						+ dfile);
	  }
	  if(fread(d.x2[0],sizeof(double),gridsize,fp) != gridsize)
    {
	  	fclose(fp);
		  throw GCLgridError(base_error
						+ "fread failed on file reading x2 coordinate  array"
						+ dfile);
    }
    if(fread(d.x3[0],sizeof(double),gridsize,fp) != gridsize)
    {
		  fclose(fp);
		  throw GCLgridError(base_error
						+ "fread failed on file reading x3 coordinate array"
						+ dfile);
    }
	  fclose(fp);
	  if(need_to_swap_bytes)
	  {
		  swapdvec(d.x1[0],gridsize);
		  swapdvec(d.x2[0],gridsize);
		  swapdvec(d.x3[0],gridsize);
	  }
  }catch(...){throw;}
};
/* this is version for 3d grids - different because of dimensions of arrays */
template <typename GCLtype> void read_GCL3d_coord_arrays(GCLtype& d, const mspass::utility::Metadata& md)
{
	try{
		std::string base_error("read_GCL3d_coord_arrays:  ");
		pfload_common_GCL_attributes(d,md);
		bool need_to_swap_bytes;
		need_to_swap_bytes=pwmig::gclgrid::byte_swap_is_needed(md);
		std::string fname,dir,dfile,dfileext;
		dir=md.get_string("dir");
		dfile=md.get_string("grid_data_file");
		fname=dir+"/"+dfile;
		if(md.is_defined("grid_data_file_extension"))
		{
			dfileext=md.get_string("grid_data_file_extension");
			fname+=".";
			fname+=dfileext;
		}
	  FILE *fp = fopen(fname.c_str(),"r");
	  if(fp == NULL)
		  throw GCLgridError(base_error
						+ "fopen failed on file "
						+ dfile);
    d.x1 = pwmig::gclgrid::create_3dgrid_contiguous(d.n1,d.n2,d.n3);
	  d.x2 = pwmig::gclgrid::create_3dgrid_contiguous(d.n1,d.n2,d.n3);
	  d.x3 = pwmig::gclgrid::create_3dgrid_contiguous(d.n1,d.n2,d.n3);
	  /* Database stores data in geographic coordinates.
	  File stores Cartesian form.  A bit inconsistent,
	  but a design choice.  Storing geo coordinates
	  would be a future format choice.*/
	  int gridsize = d.n1*d.n2*d.n3;
		if(md.is_defined("grid_data_foff"))
		{
			long foff;
			foff=md.get_long("grid_data_foff");
			if(fseek(fp,foff,SEEK_SET))
			{
				stringstream ss;
				ss << base_error
				  << "fseek to offset="<<foff<<" failed for file="<<dfile;
				throw GCLgridError(ss.str());
			}
		}
	  if(fread(d.x1[0][0],sizeof(double),gridsize,fp) != gridsize)
	  {
		  fclose(fp);
		  throw GCLgridError(base_error
						+ "fread failed on file reading x1 coordinate  array"
						+ dfile);
	  }
	  if(fread(d.x2[0][0],sizeof(double),gridsize,fp) != gridsize)
    {
	  	fclose(fp);
		  throw GCLgridError(base_error
						+ "fread failed on file reading x2 coordinate  array"
						+ dfile);
    }
    if(fread(d.x3[0][0],sizeof(double),gridsize,fp) != gridsize)
    {
		  fclose(fp);
		  throw GCLgridError(base_error
						+ "fread failed on file reading x3 coordinate array"
						+ dfile);
    }
	  fclose(fp);
	  if(need_to_swap_bytes)
	  {
		  swapdvec(d.x1[0][0],gridsize);
		  swapdvec(d.x2[0][0],gridsize);
		  swapdvec(d.x3[0][0],gridsize);
	  }
  }catch(...){throw;}
};
template <typename GCLtype> void read_fielddata(const mspass::utility::Metadata& md,
	  double *buffer,size_t buffer_size)
{
	const std::string base_error("read_fieldata:  ");
	try{
		std::string fname,dir,dfile,dfileext;
		dir=md.get_string("dir");
		dfile=md.get_string("field_data_file");
		fname=dir+"/"+dfile;
		if(md.is_defined("field_data_file_extension"))
		{
			dfileext=md.get_string("field_data_file_extension");
			fname+=".";
			fname+=dfileext;
		}
	  FILE *fp = fopen(fname.c_str(),"r");
	  if(fp == NULL)
		  throw GCLgridError(base_error
						+ "fopen failed on file "
						+ fname);
		size_t foff;
		foff=md.get_long("field_data_foff");
		if(fseek(fp,foff,SEEK_SET))
    {
      fclose(fp);
			stringstream ss;
			ss<<base_error
			  << "fseek to start of data area at byte offset="<<foff
				<<" failed for file="<<fname<<endl;
      throw GCLgridError(ss.str());
    }
		if(fread(buffer,sizeof(double),buffer_size,fp) != buffer_size)
		{
			fclose(fp);
			stringstream ss;
			ss<<base_error
			  << "fread tried to read "<<buffer_size<<" bytes but failed"<<endl;
			throw GCLgridError(ss.str());
		}
		bool need_to_swap_bytes;
		need_to_swap_bytes=pwmig::gclgrid::byte_swap_is_needed(md);
		if(need_to_swap_bytes)
		{
			swapdvec(buffer,buffer_size);
		}
	}catch(...){throw;};
};

}  // End of namespace
#endif
