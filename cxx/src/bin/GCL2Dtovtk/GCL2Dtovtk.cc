#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <float.h>
#include "mspass/utility/AntelopePf.h"
#include "mspass/utility/MsPASSError.h"
#include "pwmig/gclgrid/gclgrid.h"
#include "pwmig/gclgrid/GCLMasked.h"

using namespace mspass::utility;
using namespace pwmig::gclgrid;
/* This is also in dsap but avoiding a large include for one line */
#define rad(d)    ((d) * M_PI/180.0)
int vtk_output_GCLgrid(GCLgrid& g, string ofile,string scalars_tag="Elevation");
int vtk_output_GCLgrid(GCLscalarfield& g, string ofile,string scalars_tag);



void usage()
{
  cerr << "GCL2Dtovtk infile outfile [-f ftype -r -pf pffile -V]"<<endl
    <<    "Read from infile and write to outfile (.vtk appended to this name)"<<endl
    <<    "  -f specify type of object infile is expected to contain (default is grid2d)"<<endl
    <<    "     (one of:  grid2d, scalar2d, MaskedScalar2d)"<<endl
    <<    "  -r remap using parameters form pffile"<<endl
    <<    "  -pf specify alternative pf to default GCL2Dtovtk.pf"<<endl
    <<    "  -V  echo this usage "<<endl;
  exit(-1);
}

int main(int argc, char **argv)
{
  int i,j;
  ios::sync_with_stdio();
  if(argc<3) usage();
  string infile(argv[1]);
  string outfile(argv[2]);
  string argstr;
  string pffile("GCL2Dtovtk.pf");
  string outfieldname;
  string fieldtype("grid2d");
  bool remap(false);
  for(i=3;i<argc;++i)
  {
    argstr=string(argv[i]);
    if(argstr=="-pf")
    {
      ++i;
      if(i>=argc)usage();
      pffile=string(argv[i]);
    }
    if(argstr=="-f")
    {
      ++i;
      if(i>=argc)usage();
      fieldtype=string(argv[i]);
    }
    else if(argstr=="-V")
    {
      usage();
    }
    else if(argstr=="-r")
    {
      remap=true;
    }
    else
    {
      usage();
    }
  }
  try {

    AntelopePf control=pfread(pffile);

    string scalars_tag=control.get_string("scalars_name_tag");
    /* Slightly odd logic here, but this allows remap off
    to be the default.  pf switch is ignored this way if
    the -r flag was used */
    if(!remap) remap=control.get_bool("remap_grid");
    BasicGCLgrid *rgptr;
    rgptr=NULL;
    if(remap)
    {
      cout << "Remapping enabled"<<endl;
      /* We create a small temporary GCLgrid to allow
                           use of remap_grid procedure.   */
      double lat0,lon0,r0,azm;
      lat0=control.get_double("latitude_origin");
      lon0=control.get_double("longitude_origin");
      r0=control.get_double("radius_origin");
      azm=control.get_double("azimuth_y_axis");
      cout << "Origin lat,lon,r="
                            <<lat0<<", "<<lon0<<", "<<r0<<endl;
      cout << "Coordinate system y axis rotation="
                            <<azm<<endl;
      lat0=rad(lat0);
      lon0=rad(lon0);
      azm=rad(azm);
      GCLgrid *gtmp=new GCLgrid(2,2,string("reference"),
                                lat0,lon0,r0,azm,1.0,1.0,0,0);
      rgptr=dynamic_cast<BasicGCLgrid*>(gtmp);

    }

    if(fieldtype=="grid2d")
    {
      int npoly;
      GCLgrid g;

      g=GCLgrid(infile);
      if(remap)
      {
        //if(g!=(*rgptr))
        remap_grid(g,*rgptr);
      }
      outfile=outfile+".vtk";
      npoly=vtk_output_GCLgrid(g,outfile,scalars_tag);
      cout << "Wrote "<<npoly<<" polygons to output file"<<endl;
    }
    else if(fieldtype=="scalar2d")
    {
      int npoly;
      GCLscalarfield field(infile);
      if(remap)
      {
          remap_grid(dynamic_cast<GCLgrid&>(field),
            *rgptr);
      }
      outfile=outfile+".vtk";
      npoly=vtk_output_GCLgrid(field,outfile,scalars_tag);
      cout << "Wrote "<<npoly<<" polygons to output file"<<endl;
    }
    else if( (fieldtype=="MaskedScalar2d")
                        || (fieldtype=="MaskedGrid2d"))
    {
      int npoly;
      GCLMaskedScalarField *field;
      if(fieldtype=="MaskedScalar2d")
          field=new GCLMaskedScalarField(infile);
      else
      {
        /* Logic assumes can only land here if fieldtype is
           masked 2d grid.  In that case we read the grid and
           add depth/elevation as an attribute */
        GCLMaskedGrid gtmp(infile);
        field=new GCLMaskedScalarField(gtmp);
        if(scalars_tag=="Elevation")
        {
          for(i=0;i<gtmp.n1;++i)
            for(j=0;j<gtmp.n2;++j)
            {
              if(field->point_is_valid(i,j))
                field->val[i][j]=(-gtmp.depth(i,j));
              else
                field->val[i][j]=0.0;
              }
        }
        else
        {
          for(i=0;i<gtmp.n1;++i)
            for(j=0;j<gtmp.n2;++j)
            {
              if(field->point_is_valid(i,j))
                field->val[i][j]=gtmp.depth(i,j);
              else
                field->val[i][j]=0.0;
            }
          }
      }
      if(remap)
      {
          remap_grid(dynamic_cast<GCLgrid&>(*field),
            *rgptr);
      }
                        outfile=outfile+".vtk";
      npoly=vtk_output_GCLgrid(*field,outfile,scalars_tag);
      cout << "Wrote "<<npoly<<" polygons to output file"<<endl;
    }
    else
    {
      cerr << "Unsupported fieldtype = "<<fieldtype<<endl
        << "Exiting with no output\n" << endl;
    }

  }
  catch (MsPASSError& merr)
  {
    cerr << "Caught MsPASSError object with this content:"<<endl;
    merr.log_error();
    exit(-1);
  }
  catch(GCLgridError& gerr)
  {
    cerr << "GCLgridError thrown"<<endl<<gerr.what();
    exit(-1);
  }
  catch(std::exception& sexcp)
  {
    cerr << sexcp.what();
    exit(-1);
  }
  catch (...)
  {
      cerr << "Something threw an unhandled exception"<<endl;
      exit(-1);
  }
}
