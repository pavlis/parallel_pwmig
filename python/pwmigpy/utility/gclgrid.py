import math
from pwmigpy.ccore.gclgrid import GCLscalarfield3d,GCLvectorfield3d
from mspasspy.ccore.utility import MsPASSError

def clip(f,lower_limit=0.0,upper_limit=10.0,clip_value=None,
                  components=None):
    """
    Clip field data values in a GCLvectorfield3d or GCLscalarfield3d and 
    (optionally) set the clipped data to a specific value.   To disable 
    upper or lower clip levels set the associated clip parameter to 
    a very small (lower) or very large (upper) number.  What defines 
    small and large is data dependent so what to use for that purpose is 
    the user's responsibility.
    
    :param f:   GCL field object to be processed.  Currently only 3d 
      versions are supported 
    :param lower_limit:   lower clip bound.  Any value smaller than this 
      is set to this value or the value defined by clip_value. 
    :param upper_limit:  upper clip limit.  Any value larger than this 
      value is set to upper_limit value or the value defined by clip_value. 
    :param clip_value:  optional clip value.  When set any clipped field 
      values are set to this value.   When not set (default when None) 
      the clip levels are used to set limits.  
    :param components:  optional list of vector components to which this 
      set of clip parameters is to be applied.  Vector fields can have 
      mixed data that would require different clip values for different 
      components of the field.   The default (defined with None) 
      applies the same values to all components.  When defined only 
      components in the list will be processed.  Note an exception 
      will be thrown immediately if any of the components is outside 
      the range of nv (number of vector components)
    :return:  number of cells to which the clip was applied.  
      Note user can always get the total number of cells as 
      f.n1*f.n2*f.n3.
    """
    if clip_value == None:
        use_clip_value=False
    else:
        use_clip_value=True
    number_clipped = 0
    if isinstance(f,GCLscalarfield3d):
        for i in range(f.n1):
            for j in range(f.n2):
                for k in range(f.n3):
                    v = f.get_value(i,j,k)
                    if v < lower_limit:
                        number_clipped += 1
                        if use_clip_value:
                            v=clip_value
                        else:
                            v=lower_limit
                        f.set_value(v,i,j,k)
                    if v > upper_limit:
                        number_clipped += 1
                        if use_clip_value:
                            v=clip_value
                        else:
                            v=upper_limit
                        f.set_value(v,i,j,k)
    elif isinstance(f,GCLvectorfield3d):
        components_to_clip=list()
        if components == None:
            for l in range(f.nv):
                components_to_clip.append(l)
        else:
            for testval in components:
                if testval<=0 or testval>=f.nv:
                    # TODO  this is a bit of a crytic error - should print more information
                    raise MsPASSError('clip received invalid components parameter','Invalid')
                else:
                    components_to_clip.append(testval)
                    
        for i in range(f.n1):
            for j in range(f.n2):
                for k in range(f.n3):
                    vec = f.get_value(i,j,k)
                    Changed = False
                    if l in components_to_clip:
                        v = vec[l]
                        if v < lower_limit:
                            Changeed=True
                            if use_clip_value:
                                v=clip_value
                            else:
                                v=lower_limit
                            vec[l] = v
                        if v > upper_limit:
                            Changed=True
                            if use_clip_value:
                                v=clip_value
                            else:
                                v=upper_limit
                            vec[l] = v
                        if Changed:
                            number_clipped += 1
                            f.set_value(vec,i,j,k)
    else:
        raise MsPASSError('clip:  invalid field data input - must be GCLvectorfield3d or GCLscalarfield3d object')
    return number_clipped

def find_small(f,threshold=1e-2,components=None):
    """
    Sometimes people use 0 as a definition of a null value.  That can 
    be very problematic for finding what data are valid since an approximate 
    0 is easy to get in computations.   This function uses a simple 
    threshold to detect zeros.  For scale data it uses abs of the value 
    while for vector data it uses the L2 norm of the components.   In 
    either case it results a GCLscalarfield3d filled with 1s and 0s.  
    0 means the data was small and 1 means the data were larger than 
    the threshold.   
    
    :param f:  is the field object to be processed. Must be either 
      GCLvectorfield3d or GCLscalarfield3d or the function will throw 
      an exception. 
    :param threshold:  postive number defining the value used to define 
      a hard zero.  Less than this value is used for the zero test. 
    :param components:  list of integer component numbers in a vector 
      field to be used to compute magnitude.  Default uses all component s
      of vector data.   A single number can be used to test one component 
      of a mixed vector of data (e.g. vp and vs from an EMC tomography model)
      while a [0,1,2] would be used for pwmig output where the first 3 
      components are vector samples of the image.
    """
    result = GCLscalarfield3d(f.n1,f.n2,f.n3)
    if isinstance(f,GCLscalarfield3d):
        for i in range(f.n1):
            for j in range(f.n2):
                for k in range(f.n3):
                    v = f.get_value(i,j,k)
                    coords=f.get_coordinates(i,j,k)
                    result.set_coordinates(coords,i,j,k)
                    if abs(v)<threshold:
                        result.set_value(0.0,i,j,k)
                    else:
                        result.set_value(1.0,i,j,k)
    elif isinstance(f,GCLvectorfield3d):
        components_to_test=list()
        if components == None:
            for l in range(f.nv):
                components_to_test.append(l)
        else:
            for testval in components:
                if testval<=0 or testval>=f.nv:
                    # TODO  this is a bit of a crytic error - should print more information
                    raise MsPASSError('find_small: received invalid components parameter','Invalid')
                else:
                    components_to_test.append(testval)
        for i in range(f.n1):
            for j in range(f.n2):
                for k in range(f.n3):
                    vec = f.get_value(i,j,k)
                    mag=0.0
                    for l in range(f.nv):
                        if l in components_to_test:
                            mag += vec[l]*vec[l]
                    mag=math.sqrt(mag)
                    if mag<threshold:
                        result.set_value(0.0,i,j,k)
                    else:
                        result.set_value(1.0,i,j,k)
                    
    result.compute_extents()
    return result
                    
                    