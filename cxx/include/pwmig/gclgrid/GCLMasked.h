#ifndef _GCLMASK_H_
#define _GCLMASK_H_
#include <vector>
#include "mspass/utility/dmatrix.h"
#include "pwmig/gclgrid/gclgrid.h"
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
namespace pwmig::gclgrid
{
/*! Undefined values in a field are set to this magic value */
const double GCLFieldNullValue(-9.99999e-99);
/*! \brief  Used to apply a mask to a GCLgrid object.

  GCLgrid objects have a matrix structure that limits their natural
  extent to a distorted square in 3 space.  A mask is one way to
  build a more complex shape.   Grid cells are defined as either valid
  or invalid (true or false respectively).
  */
class GCLMask
{
    public:
        /*! Default constructor.  zero initializer */
        GCLMask();
        /*! \brief Construct from a file.

          This constructs the object from a file name.  This
          implementation uses boost's serialization to create
          the object from a file name convention.

          \param base_name is used to construct the file name for
            data as base_name + ".mask"
            */
        GCLMask(std::string base_name);
        /*! \brief Build an initial mask based on a pattern grid.

          A GCLMask object is used to turn components on or off of a GCLgrod
          or one of the field objects derived from it.  This constructor
          builds a pattern based on g and initializes all cells to
          invalid.   Use enable_point method to turn points back on. */
        GCLMask(pwmig::gclgrid::GCLgrid& g);
        /*! Standard copy constructor. */
        GCLMask(const GCLMask& parent);
        GCLMask& operator=(const GCLMask& parent);
        /*! \brief  Define point i,j of the grid as masked.

           When building a mask point by point this is a core method.
           \param i is x1 index of point to mask
           \param j is x2 index of point to mask
           \return true on success.  false if failed (outside range implied)
           */
        bool mask_point(int i, int j);
        /*! \brief  Define point i,j of the grid as valid.

           When building a mask point by point this is a core method.
           This is the inverse of mask_point as it makes a point valid.
           \param i is x1 index of point to mask
           \param j is x2 index of point to mask
           \return true on success.  false if failed (outside range implied)
           */
        bool enable_point(int i, int j);
        /*! \brief Query if a point is marked as visible.

           This is a core method.  Returns true if the point at grid
           position i,j is marked as visible(valid).  Note if the
           index is outside the grid silently returns false since
           by definition such a point is invisible.
           \param i is x1 index of point to query
           \param j is x2 index of point to query

           */
        bool point_is_valid(int i, int j);
        int nx1(){return n1_mask;};
        int nx2(){return n2_mask;};
    protected:
        /* Created as a single array, easily indexed
           because GCLgrid is n1xn2 */
        std::vector<bool> valid;
        /* These are the sizes used to create mask.  Constructors check these for
           consistency when applied to a parent GCLgrid */
        int n1_mask, n2_mask;
    private:
        /* These are used for input and output of this object. Uses base class definitions
           from gclgrid */
        friend class boost::serialization::access;
        template<class Archive>void serialize(Archive & ar, const unsigned int version)
        {
            ar & n1_mask;
            ar & n2_mask;
            ar & valid;
        };
        /* convenience routine */
        int voffset(int i, int j)
        {
            int result;
            result=j*n1_mask+i;
            return(result);
        }
};
class GCLMaskedGrid : public pwmig::gclgrid::GCLgrid, public GCLMask
{
    public:
        /*! \brief Default constructor.  */
        GCLMaskedGrid(){};
        /*! \brief Construct from a file.

          This constructor will read the output of the save
          method and recreate a clone of the file saved previously.
          This implementation uses boost serialization in a text file.

          \param fname is the file name where data is stored.
          */
        GCLMaskedGrid(std::string fname);
        /* \brief Construct with no original mask.

           Often we want to apply our own mask.  This constructor clones
           g and defines all points to be visible.
           \param g parent grid to clone
           */
        GCLMaskedGrid(pwmig::gclgrid::GCLgrid& g);
        /*! \brief Apply a mask m to g.

          Sometimes it can be useful, such as in build_masked_surface, to
          build a mask independently an then apply it to a given grid.
          This constructor does that.

          \param g - grid to which mask m is to be applied
          \param m - defining mask
          \exception GCLgridError will be thrown if g and m are not the same size.

          */
        GCLMaskedGrid(pwmig::gclgrid::GCLgrid& g,GCLMask& m);
        /* \brief Construct with a preassembled mask defined by a matrix of doubles.

           \param parent is a grid that will be cloned.
           \param mask is a double array used as a mask
           \param maskmin is a threshold used to avoid float ambiguity
               of dmatrix values.  Any value in dmatrix greater than maskmin
               is used to say mask that point.   Normally set values to some
               sensible positive number like 1.0
           \exception GCLgridError is thrown if mask size is not consistent
               with parent.
               */
        GCLMaskedGrid(pwmig::gclgrid::GCLgrid& parent,mspass::utility::dmatrix& mask,double maskmin=1.0e-13);
        /* \brief Standard copy constructor */
        GCLMaskedGrid(GCLMaskedGrid& parent);
        GCLMaskedGrid& operator=(const GCLMaskedGrid& parent);
        void save(std::string fname);
};
/* Variant of GCLscalarfield but with a smaller interface */
class GCLMaskedScalarField : public GCLscalarfield, public GCLMask
{
    public:
        /*! Construct from a file with root hame fname. */
        GCLMaskedScalarField(std::string fname);
        GCLMaskedScalarField(GCLscalarfield& g, GCLMask& m);
        GCLMaskedScalarField(GCLMaskedGrid& g);
        GCLMaskedScalarField(const GCLMaskedScalarField& parent);
        GCLMaskedScalarField& operator=(const GCLMaskedScalarField& parent);
        void save(std::string fname);
};
/* This is a variant of the GCLvectorfield but with a smaller
   interface.  */
class GCLMaskedVectorField : public GCLvectorfield, public GCLMask
{
    public:
        GCLMaskedVectorField(GCLvectorfield& g, GCLMask& m);
        GCLMaskedVectorField(GCLMaskedGrid& g, int nvsize);
        GCLMaskedVectorField(const GCLMaskedVectorField& parent);
        GCLMaskedVectorField& operator=(const GCLMaskedVectorField& parent);
        void save(std::string fname);
};
class TiePoint
{
    public:
        double radius;
        int i_tie, j_tie;
};
} // End namespace pwmig::gclgrid
#endif
