#include <fstream>
#include "pwmig/dsap/stock.h"
#include "pwmig/dsap/coords.h"
#include "mspass/utility/Metadata.h"
#include "mspass/utility/AntelopePf.h"
#include "mspass/utility/SphericalCoordinate.h"
#include "pwmig/gclgrid/gclgrid.h"
using namespace std;
using namespace mspass::utility;
using namespace pwmig::gclgrid;

void usage()
{
	cerr << "extract_section base_file_name "
		<< "[-pf pffile -o outfile -vectordata -v]"<<endl
	      << "(Default outfile is stdout)"<<endl;

	exit(-1);
}
double *extract_vertical(GCLscalarfield3d *g,Geographic_point p,double z0, double dz, int nz)
{
	double *result=new double[nz];
	Cartesian_point pc;
	/* p is not a reference so we alter it in the loop below */
	double r0=p.r;
	int i;
	for(i=0;i<nz;++i)
	{
		p.r=r0-z0-dz*static_cast<double>(i);
		pc=g->gtoc(p);
		if(g->lookup(pc.x1,pc.x2,pc.x3)==0)
		{
			result[i]=g->interpolate(pc.x1,pc.x2,pc.x3);
		}
		else
		{
			result[i]=0.0;
		}
	}
	return(result);
}

int main(int argc, char **argv)
{
	int i,j;
	ios::sync_with_stdio();
	if(argc<2) usage();
	string basename(argv[1]);
	string pfname("extract_section");
	bool vectordata(false);
	ofstream outstrm;
	bool out_to_other(false);
	string outfile;
	bool Verbose(false);

	for(i=2;i<argc;++i)
	{
		string argstr=string(argv[i]);
		if(argstr=="-pf")
		{
			++i;
			pfname=string(argv[i]);
		}
		else if(argstr=="-vectordata")
			vectordata=true;
		else if(argstr=="-o")
		{
			++i;
			out_to_other=true;
			outfile=string(argv[i]);
		}
		else if(argstr=="-v")
			Verbose=true;
		else
		{
			cerr << "Unknown argument = "<<argstr<<endl;
			usage();
		}
	}
	AntelopePf pf(pfname);
	/* format switch.  For now just ascii formats, but segy planned */
	enum OutputFormat {Dmatrix, GMT};
        /*
        char *pfkey;
        pfkey=strdup("output_format");
	char *formatname=pfget_string(pf,pfkey);
        free(pfkey);
	string stmp(formatname);
        */
        string stmp(pf.get_string("output_format"));
	OutputFormat odform;
	if(stmp=="GMT" || stmp=="gmt")
		odform=GMT;
	else
		odform=Dmatrix;
	if(Verbose)
	{
		cerr << "Using output format ";
		switch(odform)
		{
		case GMT:
			cerr << "GMT"<<endl;
			break;
		case Dmatrix:
			cerr << "dmatrix (matlab load ascii compatible)"<<endl;
		};
	}
	try {
		outstrm.open(outfile.c_str(),ios::out);
	} catch (ios::failure& var)
	{
		cerr << "Open failure on outfile="<<outfile
			<<endl
			<< "System message:  "
			<< var.what() <<endl;
		usage();
	}
	list<string> pointlist;
	pointlist=pf.get_tbl("section_points");
	if(pointlist.size()<=0)
	{
		cerr << "pf error on Tbl section_points"
			<<endl
			<<"Check parameter file"<<endl;
		exit(-1);
	}
	list<Geographic_point> section_points;
	for(auto lptr=pointlist.begin();lptr!=pointlist.end();++lptr)
	{
		double lat,lon;
		stringstream ss(*lptr);
		ss>>lat;
		ss>>lon;
		Geographic_point p;
		p.lat=rad(lat);
		p.lon=rad(lon);
		p.r=r0_ellipse(p.lat);
		section_points.push_back(p);
	}
	if(Verbose) cerr << "Read section path defined by "
		<< section_points.size()
		<< " control points"<<endl;
	try {
		double dx=pf.get_double("path_sample_interval");
		double dz=pf.get_double("depth_sample_interval");
		double zmax=pf.get_double("maximum_depth");
		double zmin=pf.get_double("minimum_depth");
		bool output_latlon=pf.get_bool("geographic_output");
		int nz=static_cast<int>((zmax-zmin)/dz) + 1;
		if(nz<=0)
		{
			cerr << "Illegal maximum_depth or minimum_depth "
				<< "parameter specified."<<endl
				<< "Check parameter file"<<endl;
			exit(-1);
		}
		bool use_points_directly=pf.get_bool("use_points_directly");
		if(use_points_directly && Verbose)
			cerr << "Using input points directly to define output data"<<endl;
		/* Now build the full path of control points as equally spaced
		as possible.  Using the latlon function from antelope. */
		Geographic_point p;
		double del,ddel,az,lat1,lon1,lat2,lon2;
		double delkm,r0,dr;
		int n;
		list<Geographic_point> path;
		list<Geographic_point>::iterator sptr;
		if(use_points_directly)
			path=section_points;
		else
		{
		    int ipath;
		    for(ipath=0,sptr=section_points.begin();
			ipath<(section_points.size() - 1);++ipath)
		    {
			if(sptr==section_points.begin())
			{
				lat1=sptr->lat;
				lon1=sptr->lon;
				r0=sptr->r;
			}
			else
			{
				lat1=p.lat;
				lon1=p.lon;
				r0=p.r;
			}
			++sptr;
			lat2=sptr->lat;
			lon2=sptr->lon;
			dr=(sptr->r)-r0;
			dist(lat1,lon1,lat2,lon2,&del,&az);
			delkm=del*r0;
			n=static_cast<int>(delkm/dx) + 1;
			n=(int)round(delkm/dx) + 1;
			ddel=dx/r0;
			dr=dr*dx/delkm;
			for(i=0;i<n;++i)
			{
				del=ddel*static_cast<double>(i);
				latlon(lat1,lon1,del,az,&lat2,&lon2);
				p.lat=lat2;
				p.lon=lon2;
				p.r=r0+dr*static_cast<double>(i);
				path.push_back(p);
			}
		    }
		}

		GCLscalarfield3d *g;
		if(vectordata)
		{
			int component=pf.get_int("component_number");
			GCLvectorfield3d *gvec=new GCLvectorfield3d(basename);
			g=extract_component(*gvec,component);
			delete gvec;
		}
		else
		{
			g=new GCLscalarfield3d(basename);
		}
		/* This will hold the output.  */
		dmatrix section(nz,path.size());
		section.zero();  /* this allows us to ignore edge effects*/
		double *seis;
		for(sptr=path.begin(),j=0;sptr!=path.end();++sptr,++j)
		{
			seis=extract_vertical(g,*sptr,zmin,dz,nz);
			for(i=0;i<nz;++i)
			{
				section(i,j)=seis[i];
			}
			delete [] seis;
		}
		switch (odform)
		{
		case GMT:
			for(j=0,sptr=path.begin();j<section.columns() && sptr!=path.end();++j,++sptr)
			{
				for(i=0;i<nz;++i)
				{
				    if(out_to_other)
				    {
					if(output_latlon)
					    outstrm << deg(sptr->lat)
						<< " "<<deg(sptr->lon)
						<< " "<<zmin+i*dz
						<< " "<<section(i,j)<<endl;
					else
					    outstrm << j*dx
						<< " "<<zmin+i*dz
						<< " "<<section(i,j)<<endl;
				    }
				    else
				    {
					if(output_latlon)
					    cout << deg(sptr->lat)
						<< " "<<deg(sptr->lon)
						<< " "<<zmin+i*dz
						<< " "<<section(i,j)<<endl;
					else
					    cout << j*dx
						<< " "<<zmin+i*dz
						<< " "<<section(i,j)<<endl;
				    }
				}
				if(out_to_other)
					outstrm<<">"<<endl;
				else
					cout << ">"<<endl;
			}
			break;
		case Dmatrix:
		default:
			if(out_to_other)
				outstrm << section;
			else
				cout << section;
		}
		delete g;
	}
	catch (MsPASSError& serr)
	{
		serr.log_error();
	}
        catch (exception& gerr)
        {
            cerr << "Error was throw.  This is the message:"<<endl
                << gerr.what()<<endl;
        }
}
