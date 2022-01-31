#include "fstream"
#include "mspass/utility/dmatrix.h"
#include "pwmig/gclgrid/GCLgridError.h"
#include "pwmig/gclgrid/gclgrid.h"
#include "pwmig/gclgrid/GCLMasked.h"
namespace pwmig::gclgrid
{
using namespace std;
using namespace mspass::utility;
using namespace pwmig::gclgrid;
/* First in this file - GCLMask code. */
GCLMask::GCLMask()
{
    n1_mask=0;
    n2_mask=0;
}
GCLMask::GCLMask(string base_name)
{
    string fname;
    fname=base_name+".mask";
    std::ifstream ifs(fname.c_str(),ios::in);
    boost::archive::text_iarchive ia(ifs);
    ia >> *this;
}
GCLMask::GCLMask(GCLgrid& parent)
{
    n1_mask=parent.n1;
    n2_mask=parent.n2;
    int ntotal=parent.n1*parent.n2;
    valid.reserve(ntotal);
    int k;
    for(k=0;k<ntotal;++k)valid.push_back(false);
}
GCLMask::GCLMask(const GCLMask& parent) : valid(parent.valid)
{
    n1_mask=parent.n1_mask;
    n2_mask=parent.n2_mask;
}
GCLMask& GCLMask::operator=(const GCLMask& parent)
{
    if(this != & parent)
    {
        this->valid=parent.valid;
        n1_mask=parent.n1_mask;
        n2_mask=parent.n2_mask;
    }
    return(*this);
}
bool GCLMask::mask_point(int i, int j)
{
    if(i<0 || (i>=n1_mask)) return false;
    if(j<0 || (j>=n2_mask)) return false;
    valid[this->voffset(i,j)] = false;
    return(true);
}
bool GCLMask::enable_point(int i, int j)
{
    if(i<0 || (i>=n1_mask)) return false;
    if(j<0 || (j>=n2_mask)) return false;
    valid[this->voffset(i,j)] = true;
    return(true);
}
bool GCLMask::point_is_valid(int i, int j)
{
    if(i<0 || (i>=n1_mask)) return false;
    if(j<0 || (j>=n2_mask)) return false;
    int k=voffset(i,j);
    return(valid[k]);
}



GCLMaskedGrid::GCLMaskedGrid(string fname) 
    : GCLgrid(fname),GCLMask(fname)
{
}
GCLMaskedGrid::GCLMaskedGrid(GCLgrid& parent) : GCLgrid(parent), GCLMask(parent)
{
}
GCLMaskedGrid::GCLMaskedGrid(GCLgrid& g,GCLMask& mask) : GCLgrid(g), GCLMask(mask)
{
}
GCLMaskedGrid::GCLMaskedGrid(GCLgrid& parent, dmatrix& mask, double maskmin)
                        : GCLgrid(parent)
{
    const string base_error("GCLMaskedGrid constructor with dmatrix mask:  ");
    if(mask.columns()!=parent.n2) throw GCLgridError(base_error
            + "number of columns in mask array does not match grid");
    if(mask.rows()!=parent.n1) throw GCLgridError(base_error
            + "number of columns in mask array does not match grid");
    int i,j;
    for(i=0;i<parent.n1;++i)
        for(j=0;j<parent.n2;++j)
        {
            if(mask(i,j)>maskmin)
                this->mask_point(i,j);
            else
                this->enable_point(i,j);
        }
}
GCLMaskedGrid::GCLMaskedGrid(GCLMaskedGrid& parent) 
    : GCLgrid(dynamic_cast<GCLgrid&>(parent)), GCLMask(dynamic_cast<GCLMask&>(parent))
{
}
GCLMaskedGrid& GCLMaskedGrid::operator=(const GCLMaskedGrid& parent)
{
    if(this != & parent)
    {
        this->GCLgrid::operator=(parent);
        this->GCLMask::operator=(parent);
    }
    return(*this);
}
/* This is the first of 3 save routines in this file.  They are painfully
   similar.  Considered messing with a template, but decided it was 
   an unnecessary complication and I'd deal with the repetitious code. */
void GCLMaskedGrid::save(string fname)
{
    try{
    /* Downcast to the parent GCLgrid and call it's save method to
       store base data.  Limitation here is I don't allow use
     of a directory as in GCLgrid library.  Alway current directory*/
        this->GCLgrid::save(fname,string("."));
        /* we have a boost serialization to save the mask data.
           We store that in a name derived from fname.  This 
           assume format with fname+".dat" as data and fname+".pf" 
           to store metadata. */
        string maskfile=fname+".mask";
        std::ofstream ofs(maskfile.c_str(),ios::out);
        boost::archive::text_oarchive oa(ofs);
        GCLMask *mask=dynamic_cast<GCLMask*>(this);
        oa << *mask;
        ofs.close();
    }catch(...){throw;};
}
GCLMaskedVectorField::GCLMaskedVectorField(GCLvectorfield& f, GCLMask& m) 
    : GCLvectorfield(f), GCLMask(m)
{
}
GCLMaskedVectorField::GCLMaskedVectorField(GCLMaskedGrid& g, int nvsize)
    : GCLvectorfield(dynamic_cast<GCLgrid&>(g),nvsize), GCLMask(dynamic_cast<GCLMask&>(g))
{
}
GCLMaskedVectorField::GCLMaskedVectorField(const GCLMaskedVectorField& parent)
    : GCLvectorfield(dynamic_cast<const GCLvectorfield&>(parent)),
        GCLMask(dynamic_cast<const GCLMask&>(parent))
{
}
GCLMaskedVectorField& GCLMaskedVectorField::operator=(const GCLMaskedVectorField& parent)
{
    if(this != & parent)
    {
       this->GCLvectorfield::operator=(parent);
       this->GCLMask::operator=(parent);
    }
    return(*this);
}
void GCLMaskedVectorField::save(string fname)
{
    try{
    /* Downcast to the parent GCLgrid and call it's save method to
       store base data.  Limitation here is I don't allow use
     of a directory as in GCLgrid library.  Alway current directory*/
        this->GCLvectorfield::save(fname,string("."));
        /* we have a boost serialization to save the mask data.
           We store that in a name derived from fname.  This 
           assume format with fname+".dat" as data and fname+".pf" 
           to store metadata. */
        string maskfile=fname+".mask";
        std::ofstream ofs(maskfile.c_str(),ios::out);
        boost::archive::text_oarchive oa(ofs);
        GCLMask *mask=dynamic_cast<GCLMask*>(this);
        oa << *mask;
        ofs.close();
    }catch(...){throw;};
}
GCLMaskedScalarField::GCLMaskedScalarField(string fname) 
    : GCLscalarfield(fname,default_output_format,false),GCLMask(fname)
{
}
GCLMaskedScalarField::GCLMaskedScalarField(GCLscalarfield& f, GCLMask& m) 
    : GCLscalarfield(f), GCLMask(m)
{
}
GCLMaskedScalarField::GCLMaskedScalarField(GCLMaskedGrid& g) 
    : GCLscalarfield(dynamic_cast<GCLgrid&>(g)), GCLMask(dynamic_cast<GCLMask&>(g))
{
}
GCLMaskedScalarField::GCLMaskedScalarField(const GCLMaskedScalarField& parent)
    : GCLscalarfield(dynamic_cast<const GCLscalarfield&>(parent)), 
        GCLMask(dynamic_cast<const GCLMask&>(parent))
{
}
GCLMaskedScalarField& GCLMaskedScalarField::operator=(const GCLMaskedScalarField& parent)
{
    if(this != & parent)
    {
       this->GCLscalarfield::operator=(parent);
       this->GCLMask::operator=(parent);
    }
    return(*this);
}
void GCLMaskedScalarField::save(string fname)
{
    try{
    /* Downcast to the parent GCLgrid and call it's save method to
       store base data.  Limitation here is I don't allow use
     of a directory as in GCLgrid library.  Alway current directory*/
        this->GCLscalarfield::save(fname,string("."));
        /* we have a boost serialization to save the mask data.
           We store that in a name derived from fname.  This 
           assume format with fname+".dat" as data and fname+".pf" 
           to store metadata. */
        string maskfile=fname+".mask";
        std::ofstream ofs(maskfile.c_str(),ios::out);
        boost::archive::text_oarchive oa(ofs);
        GCLMask *mask=dynamic_cast<GCLMask*>(this);
        oa << *mask;
        ofs.close();
    }catch(...){throw;};
}
}  // End namespace pwmig::gclgrid
