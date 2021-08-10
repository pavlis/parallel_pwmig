#include <fstream>
#include "stock.h"
#include "dbpp.h"
#include "Metadata.h"
#include "seispp.h"
#include "gclgrid.h"
#include "FixedFormatTrace.h"
using namespace std;
using namespace SEISPP;

void usage()
{
	cerr << "extract_section db gridname fieldname "
		<< "[-pf pffile -o outfile -vectordata -v]"<<endl
	      << "(Default outfile is stdout)"<<endl;

	exit(-1);
}
double *extract_vertical(GCLscalarfield3d *g,Geographic_point p,double z0, double dz, int nz)
{
	double *result=new double[nz];
	Cartesian_point pc;
	/* p is not a reference so we alter it in the loop below */
	double z;
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
void writesu(dmatrix& section, double zmin, double dz, list<Geographic_point>& path,
		ofstream& out)
{
    try {
	int i,j;
	int nz=section.rows();
	int nd=section.columns();
	list<Geographic_point>::iterator p0,p;
	p0=path.begin();
	double lat0,lon0;
	double lat,lon,azimuth,delta,distkm;
	lat0=p0->lat;
	lon0=p0->lon;
	FixedFormatTrace d(string("SEGYfloat"),nz);
	const double DEG2KM(111.320);
	vector<double> values;
	values.reserve(nz);
	for(j=0,p=path.begin();j<nd;++j,++p)
	{
		for(i=0;i<nz;++i) values.push_back(section(i,j));
		d.put(0,values);
		values.clear();
		/* These frozen names are a maintenance nightmare because
		these relate to an independent pf file that defines them
		(HeaderMap.pf) */
		d.put<int>("lineSeq",j+1);
		d.put<int>("reelSeq",1);
		d.put<int>("event_number",1);
		d.put<int>("cdpEns",1);
		d.put<int>("traceInEnsemble",j+1);
		d.put<int>("traceId",1);
		d.put<int>("coordUnits",10000);  //degree multiplier for lat lon
		lat=deg(p->lat);
		lon=deg(p->lon);
		dist(lat0,lon0,p->lat,p->lon,&delta,&azimuth);
		distkm=deg(delta)*DEG2KM;
		lat*=10000;
		lon*=10000;
		d.put<double>("recLongOrX",lon);
		d.put<double>("recLatOrY",lat);
		/* In SEGY/SU dt is microseconds.  This puts
		km scale depths to soomething closer to standard
		seismic.  */
		const double dz_multiplier(1000.0);
		d.put<double>("deltaSample",dz_multiplier*dz);
		d.put<int>("sampleLength",nz);
		d.put<double>("sourceToRecDist",distkm);
		/* This marks a trace live in SEGY */
		d.put<int>("dataUse",1);
		d.write(out);
	} 
	out.close();
    }catch(...){throw;};
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
	int i,j;
	ios::sync_with_stdio();
	if(argc<4) usage();
	string dbname(argv[1]);
	string gridname(argv[2]);
	string fieldname(argv[3]);
	string pfname("extract_section");
	bool vectordata(false);
	ofstream outstrm;
	bool out_to_other(false);
	string outfile;

	for(i=4;i<argc;++i)
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
			SEISPP_verbose=true;
		else
		{
			cerr << "Unknown argument = "<<argstr<<endl;
			usage();
		}
	}
	Pf *pf;
	if(pfread(const_cast<char *>(pfname.c_str()),&pf)) 
	{
		cerr << "pfread failed"<<endl;
		usage();
	}
	/* format switch.  For now just ascii formats, but segy planned */
	enum OutputFormat {Dmatrix, GMT, SEISMICUNIX};
        char *pfkey;
        pfkey=strdup("output_format");
	char *formatname=pfget_string(pf,pfkey);
        free(pfkey);
	string stmp(formatname);
	OutputFormat odform;
	if(formatname==NULL)
	{
		cerr << "Required parameter output_format is not in "
			<< "parameter file"<<endl;
		usage();
	}
	else if(stmp=="GMT" || stmp=="gmt")
		odform=GMT;
	else if(stmp=="SEISMICUNIX")
		odform=SEISMICUNIX;
	else
		odform=Dmatrix;
	if(SEISPP_verbose)
	{
		cerr << "Using output format ";
		switch(odform)
		{
		case GMT:
			cerr << "GMT"<<endl;
			break;
		case Dmatrix:
			cerr << "dmatrix (matlab load ascii compatible)"<<endl;
			break;
		case SEISMICUNIX:
			cerr << "SEISMICUNIX"<<endl;
		};
	}
	if(odform==SEISMICUNIX && !out_to_other)
	{
		cerr << "Illegal parameter combination.  "
			<< "You must use the -o argument when writing SEISMICUNIX data"
			<<endl;
		exit(-1);
	}
	try {
		if(odform==SEISMICUNIX)
			outstrm.open(outfile.c_str(),ios::out | ios::binary);
		else
			outstrm.open(outfile.c_str(),ios::out);
	} catch (ios::failure& var)
	{
		cerr << "Open failure on outfile="<<outfile
			<<endl
			<< "System message:  "
			<< var.what() <<endl;
		usage();
	}
		
	Tbl *pointlist;
        pfkey=strdup("section_points");
	pointlist=pfget_tbl(pf,pfkey);
        free(pfkey);
	if(pointlist==NULL)
	{
		cerr << "pf error on Tbl section_points"
			<<endl
			<<"Check parameter file"<<endl;
		exit(-1);
	}
	list<Geographic_point> section_points;
	for(i=0;i<maxtbl(pointlist);++i)
	{
		char *line;
		double lat,lon;
		line=(char *)gettbl(pointlist,i);
		sscanf(line,"%lf%lf",&lat,&lon);
		Geographic_point p;
		p.lat=rad(lat);
		p.lon=rad(lon);
		p.r=r0_ellipse(p.lat);
		section_points.push_back(p);
	}
	if(SEISPP_verbose) cerr << "Read section path defined by "
		<< section_points.size()
		<< " control points"<<endl;
	try {
		DatascopeHandle dbh(dbname,true);
		dbh.lookup(string("gclgdisk"));
		DatascopeHandle dbhg(dbh);
		dbhg.lookup(string("gclfield"));
		Dbptr db=dbh.db;
		Dbptr dbgrd=dbhg.db;
		
		Metadata control(pf);
		double dx=control.get_double("path_sample_interval");
		double dz=control.get_double("depth_sample_interval");
		double zmax=control.get_double("maximum_depth");
		double zmin=control.get_double("minimum_depth");
		bool output_latlon=control.get_bool("geographic_output");
		int nz=static_cast<int>((zmax-zmin)/dz) + 1;
		if(nz<=0)
		{
			cerr << "Illegal maximum_depth or minimum_depth "
				<< "parameter specified."<<endl
				<< "Check parameter file"<<endl;
			exit(-1);
		}
		bool use_points_directly=control.get_bool("use_points_directly");
		if(use_points_directly && SEISPP_verbose)
			cerr << "Using input points directly to define output data"<<endl;
		/* Now build the full path of control points as equally spaced 
		as possible.  Using the latlon function from antelope. */
		double del,ddel,az,lat1,lon1,lat2,lon2;
		double delkm,r0,dr;
		int n;
		Geographic_point p={0.0,0.0,0.0};
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
			n=SEISPP::nint(delkm/dx) + 1;
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
			int component=control.get_int("component_number");
                        int nvcomp=control.get_int("number_vector_components");
			GCLvectorfield3d *gvec=new GCLvectorfield3d(dbh,
                                gridname,fieldname,nvcomp);
			g=extract_component(*gvec,component);
			delete gvec;
		}
		else
		{
			g=new GCLscalarfield3d(dbh,gridname,fieldname);
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

		case SEISMICUNIX:
			writesu(section,zmin,dz,path,outstrm);
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
        catch (exception& gerr)
        {
            cerr << "Error was throw.  This is the message:"<<endl
                << gerr.what()<<endl;
        }
	catch (SeisppError& serr)
	{
		serr.log_error();
	}
}
