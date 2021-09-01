#include "pwmig/gclgrid/gclgrid.h"
#include "pwmig/gclgrid/PWMIGfielddata.h"
namespace pwmig::gclgrid {
using namespace pwmig::gclgrid;
using namespace mspass::utility;

PWMIGfielddata::PWMIGfielddata() : GCLvectorfield3d(),elog()
{};

PWMIGfielddata::PWMIGfielddata(const pwmig::gclgrid::GCLgrid3d& g)
  : GCLvectorfield3d(g,5),elog()
{};

PWMIGfielddata::PWMIGfielddata(const PWMIGfielddata& parent)
  : GCLvectorfield3d(dynamic_cast<const GCLvectorfield3d&>(parent)),elog(parent.elog)
{};

PWMIGfielddata& PWMIGfielddata::operator=(const PWMIGfielddata& parent)
{
  if(this != &parent)
  {
    this->GCLvectorfield3d::operator=(parent);
    elog=parent.elog;
  }
  return *this;
};

void PWMIGfielddata::accumulate(const pwmig::pwmigcore::PWMIGmigrated_seismogram& d)
{
  size_t i,j,n3(0);   // initializing to zero to assure 0 is returned for dead data
  if(d.live)
  {
    i=d.ix1;
    j=d.ix2;
    n3=this->n3;
    if( (n3!=d.dweight.size()) || (n3!=d.domega.size()) || (n3!=d.migrated_data.columns()))
    {
      stringstream ss;
      ss << "PWMIGfielddata::accumulate:  inconsistent array sizes"<<endl
         << "dweight vector length="<<d.dweight.size()
         << " domega vector length="<<d.domega.size()<<endl
         << "migrated data matrix size="<<d.migrated_data.rows()<<"X"<<d.migrated_data.columns()<<endl
         << "Internal grid n3 dimension="<<n3<<endl;
      throw MsPASSError(ss.str(),ErrorSeverity::Fatal);
    }
    size_t k,kk;
    for(k=0,kk=this->n3-1;k<this->n3;++k,--kk)
    {
      for(size_t l=0;l<3;++l)
      {
        this->val[i][j][k][l]=d.migrated_data(l,kk)
          *d.dweight[kk]*d.domega[kk];
      }
      this->val[i][j][k][3]=d.domega[kk];
      this->val[i][j][k][4]=d.dweight[kk];
    }
  }
  if(d.elog.size()>0)
  {
    list<LogData> logdata=d.elog.get_error_log();
    for(auto lptr=logdata.begin();lptr!=logdata.end();++lptr)
    {
      stringstream ss;
      ss << "Migrated data with ix1="<<i<<" and ix2="<<j<<" contained the following error message:"<<endl
          << lptr->message<<endl;
      this->elog.log_error(lptr->algorithm,string(ss.str()),lptr->badness);
    }
  }
}

} // End namespace
