#include <exception>
#include "pwmig/utility/coords.h"
#include "pwmig/utility/EventCatalog.h"
// This file should be completely ignored if antelope is not defined.
// This is here to allow simple compilation of all seispp with the NO_ANTELOPE option
#ifndef NO_ANTELOPE
using namespace std;
using namespace pwmig::utility;
namespace pwmig::utility
{
EventCatalog::EventCatalog(const MetadataList& mdl)
{
	mdloaded=mdl;
}


TimeWindow EventCatalog::range() const
{
	map<Hypocenter,Metadata,SpaceTimeCompare>::iterator hs,he;
	hs=catalog.begin();
	he=catalog.end();
	return(TimeWindow(hs->first.time,he->first.time));
}
bool EventCatalog::add(const Hypocenter& h,const Metadata& md)
{
	if(catalog.find(h)!=catalog.end())
	{
		return false;
	}
	else
	{
		catalog.insert(pair<Hypocenter,Metadata>(h,md));
		return  true;
	}
}
bool EventCatalog::add(const Metadata& md)
{
	try{
		Hypocenter h(md);
		return this->add(h,md);
	}catch(...){throw;};
}

bool EventCatalog::replace(const Hypocenter&h, const Metadata& md)
{
	bool result;
	result=false;
	if(catalog.find(h)==catalog.end()) result=true;
	catalog.insert(pair<Hypocenter,Metadata>(h,md));
	return result;
}

bool EventCatalog::replace(const Metadata& md)
{
	try{
		Hypocenter h(md);
		return this->replace(h,md);
	}catch(...){throw;};
}

bool EventCatalog::find(Hypocenter& h)
{
	bool result;
	map<Hypocenter,Metadata,SpaceTimeCompare>::iterator hptr;
	hptr=catalog.find(h);
	if(hptr==catalog.end())
	{
		result=false;
	}
	else
	{
		current_hypo=hptr;
		result=true;
	}
	return result;
}
Hypocenter EventCatalog::current()
{
	return(current_hypo->first);
}
Metadata EventCatalog::current_aux()
{
	return(current_hypo->second);
}
void EventCatalog::delete_current()
{
	catalog.erase(current_hypo);
	if(catalog.size()>0) current_hypo=catalog.begin();
}
void EventCatalog::rewind()
{
	current_hypo=catalog.begin();
}
void EventCatalog::operator++()
{
	if(current_hypo==catalog.end())
		throw MsPASSError(string("EventCatalog::opeator++:  ")
			+ "attempt to go outside valid range");
	++current_hypo;
}
int EventCatalog::size()
{
	return(catalog.size());
}
EventCatalog::EventCatalog(const EventCatalog& parent)
{
	catalog=parent.catalog;
	current_hypo=parent.current_hypo;
	mdloaded=parent.mdloaded;
}

EventCatalog& EventCatalog::operator=(const EventCatalog& parent)
{
    if(this!=&parent)
    {
	catalog=parent.catalog;
	current_hypo=parent.current_hypo;
	mdloaded=parent.mdloaded;
    }
    return(*this);
}

} /* End SEISPP namespace encapsulation*/
#endif
