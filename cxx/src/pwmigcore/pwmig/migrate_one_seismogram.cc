bool weight_functions_set=false;
// Warning:  n3 must not change in the loop below
int n3=raygrid.n3;
vector<double> domega_ij(n3),dweight_ij(n3);
for(is=0;is<pwdata->member.size();++is)
{
  // convenient shorthand variables.  ns is data length
  // while n3 is the ray path onto which data are mapped
  int ns=pwdata->member[is].ns;
  double zmaxray,tmaxray;
  vector<double> Stime(n3), SPtime(n3);
  double t0;

  dt=pwdata->member[is].dt;
  // first decide in which grid cell to place this seismogram
  i = pwdata->member[is].get_int("ix1");
  j = pwdata->member[is].get_int("ix2");
  if(ApplyElevationStatics)
  {
    /* We assume elevation has been posted as elev.  This is done
    by PwmigFileHandle methods on loading, but if data handling method
    changes this needs to be watched.  Be aware this elevation is produced
    in pwstack as a weighted average of elevations of summed traces */
    double pselev=pwdata->member[is].get_double("elev");
    t0=pwdata->member[is].t0;
    ApplyGeometricStatic(dynamic_cast<BasicTimeSeries *>(&(pwdata->member[is])),
      static_velocity,pselev);
    if(SEISPP_verbose)
      cout << "Applying static time shift = "
        << pwdata->member[is].t0 - t0<<endl;
  }
  t0=pwdata->member[is].t0;
  if( (i<0) || (i>raygrid.n1) || (j<0) || (j>raygrid.n2) )
  {
    cerr << "grid index of data = ("
      << i << ","<< j
      << "is outside range of grid." <<endl
      << "Check grid definitions in database and pf"<<endl;
    exit(-1);
  }
  // We need a reference ray for the incident wavefield.
  // This ray is chopped up and transformed to GCLgrid coordinates
  // in the integration loop below.  We compute it here to
  // avoid constantly recomputing it
  zmaxray = raygrid.depth(i,j,0);
  zmaxray *= zpad;                      // Small upward adjustment to avoid rubber banding
  // problems in 3d media
  tmaxray=10000.0;                      // arbitrary large number.  depend on zmax

  // Now compute a series of vectorish quantities that
  // are parallel with the ray path associated with the
  // current seismogram (i,j).  First we compute the S wave
  // travel time down this ray path (i.e. ordered from surface down)
  Stime = compute_Stime(*Us3d,i,j,raygrid,use_3d_vmodel);
  // Now we compute the gradient in the S ray travel time
  // for each point on the ray.  Could have been done with
  // less memory use in the loop below, but the simplification
  // it provides seems useful to me
  gradTs = compute_gradS(raygrid,i,j,Vs1d);
  dmatrix gradTp(3,n3);
  dmatrix nup(3,n3);
  vector<double> zP(n3);
  // Now loop over each point on this ray computing the
  // p wave path and using it to compute the P time and
  // fill the vector gradP matrix
  // Note the double index, one running forward (k) and
  // a second (kk) running backward.  This complication is
  // caused by raygrid having upward oriented curves while
  // the S ray path (pathptr) is oriented downward.
  // Only need to compute P time to surface once for
  // each ray so we compute it before the loop below
  /* old form
                             double Tpr=compute_receiver_Ptime(raygrid,i,j,
                                      *TPptr,hypo);
                                      */
  double Tpr=TPptr->val[i+border_pad][j+border_pad][TPptr->n3-1];
  double tlag,Tpx,tdelta;
  bool needs_padding;
  int padmark=raygrid.n3-1;
  bool tcompute_problem;                // set true if tt calculation fails away from edges
  /* This is the position of the incident ray at the base of the model.
   * Common for all in next loop and needed for travel time calculation.
   * Note the ugly border pad construct necessary to account for padding
   * efficiently - avoids a lookup */
  Geographic_point rTP_gp=fetch_TP_point_at_base(*TPptr,i+border_pad,j+border_pad);
  /* It is a good idea to initialize these before each raty path */
  TPptr->reset_index();
  SPtime.clear();
  tcompute_problem=false;
  needs_padding=false;
  for(k=0,kk=raygrid.n3-1;k<raygrid.n3;++k,--kk)
  {
    vector<double>nu;
    double vp;
    SlownessVector u0=svm0(i,j);
    /* This is needed below to compute p*delta term */
    Geographic_point x_gp,rxTP_gp0;
    x_gp=raygrid.geo_coordinates(i,j,kk);

    nu = compute_unit_normal_P(*TPptr,raygrid.x1[i][j][kk],
      raygrid.x2[i][j][kk], raygrid.x3[i][j][kk]);
    for(l=0;l<3;++l) nup(l,kk)=nu[l];
    // zP and gradTp are computed in reverse order
    zP[k]=raygrid.depth(i,j,kk);
    vp=Vp1d.getv(zP[k]);
    for(l=0;l<3;++l) gradTp(l,k)=nu[l]/vp;
    /* This section used to be a procedure.  Inlined for
    speed during a cleanup June 2012 */
    int error_lookup;
    error_lookup=TPptr->lookup(raygrid.x1[i][j][kk],
      raygrid.x2[i][j][kk],raygrid.x3[i][j][kk]);
    switch(error_lookup)
    {
      case 0:
        Tpx=TPptr->interpolate(raygrid.x1[i][j][kk],
          raygrid.x2[i][j][kk],raygrid.x3[i][j][kk]);
        break;
      case 1:
        /* This means the point falls outside the P ray grid.
        Normally that means we should just break the loop.
        Issue a warning only if the kk index is small indicating
        the ray path is very short  */
        needs_padding=true;
        padmark=k;
        break;
      case -1:
      default:
        /* This is a convergence error in lookup and needs to be
        flagged as an error */
        tcompute_problem=true;
        needs_padding=true;
        padmark=k;
    }
    /* New procedure (Mar 2016) to fix error in approximation found for
     * no antelope version pwmig2.0 */
    try
    {
      rxTP_gp0=get_gp_base_TPx(*TPptr,x_gp);
    }catch(SeisppError& serr)
    {
      serr.log_error();
      tcompute_problem=true;
    }
    if(tcompute_problem || needs_padding) break;
    /* Original version had a different form here with 3 terms.
     * We have to add a new one for relative time calculation
     * here to correct for p*delta */
    tdelta=compute_delta_p_term(rxTP_gp0,rTP_gp,u0);
    tlag=Tpx+Stime[k]+tdelta-Tpr;
    //DEBUG
    /*cout << "k="<<k<<"Tpx="<<Tpx<<" Stime[k]="<<Stime[k]
        << " tdelta="<<tdelta<<" Tpr="<<Tpr<<" = tlag of "<<tlag<<endl;
        */
    //cout << "tdelta="<<tdelta <<" tlag="<<tlag<<endl;
    SPtime.push_back(tlag);
  }

  // skip this ray if there were travel time computation problems
  if(tcompute_problem)
  {
    cerr << "Warning:  slowness gridid "<< gridid
      << ", grid position index ("
      << i << ","
      << j << ").  GCLscalarfield3d lookup failed at ray index="
      << k <<endl
      << "Geographic location: "<<raygrid.lat(i,j,kk)
      <<", "<<deg(raygrid.lon(i,j,kk))
      <<", "<<raygrid.depth(i,j,kk)<<endl
      << "Data from failed ray index downward will be zeroed"<<endl
      << "This may leave ugly holes in the output image."
      <<endl;

  }
  if (needs_padding)
  {
    if(padmark==0)
    {
      cerr << "Warming:   gridid "<< gridid
        << ", grid position index ("
        << i << ","
        << j << ")." <<endl
        << "P time interpolation failed at earth's surface.  "<<endl
        << "This should not happen and may cause holes in the output."<<endl
        << "This trace will be projected as zeros"<<endl;
    }
    else
    {
      /* Arbitrarily increment the times by dt intervals.  Will
      zero data beyond pad mark below.  Not the most efficient way
      to do this, but don't expect this to be that common of an issue*/
      double tpadding0=Stime[padmark-1];
      for(k=padmark;k<raygrid.n3;++k)
        Stime.push_back(tpadding0+dt*static_cast<double>(k-padmark+1));
    }
  }
  // We now interpolate the data with tlag values to map
  // the data from time to an absolute location in space
  // Note carefully SPtime and work sense is inverted from
  // the raygrid.  i.e. they are oriented for a ray from the surface
  // to the bottom.  Below we have to reverse this
  dmatrix work(3,raygrid.n3);
  linear_vector_regular_to_irregular(t0,dt,pwdata->member[is].u,
    &(SPtime[0]),work);
  /* Cosine aper all data that was padded and zero the extension */
  if(needs_padding)
    cosine_taper_highend(work,padmark,taper_length);
  dmatrix raycoh(4,coh3cens->member[is].ns);
  int ncohcopy=coh3cens->member[i].ns;
  /* excessively paranoid, but murphy's law */
  if(cohens->member[is].ns != ncohcopy)
  {
    cerr << "Warning:  inconsistent coherence data size"<<endl
      << "Three-C data size="<<ncohcopy
      << " while all component sum size="
      << cohens->member[is].ns<<endl;
    ncohcopy=min(cohens->member[is].ns,ncohcopy);
    cerr << "Set to minimum="<<ncohcopy<<endl;
  }
  /* This component coherence into first 3 rows and sum in 4 */
  dcopy(ncohcopy,coh3cens->member[is].u.get_address(0,0),3,
    raycoh.get_address(0,0),4);
  dcopy(ncohcopy,coh3cens->member[is].u.get_address(1,0),3,
    raycoh.get_address(1,0),4);
  dcopy(ncohcopy,coh3cens->member[is].u.get_address(2,0),3,
    raycoh.get_address(2,0),4);
  dcopy(ncohcopy,&(cohens->member[is].s[0]),1,
    raycoh.get_address(3,0),4);
  int nout=SPtime.size();
  dmatrix raycohout(4,nout);
  /* we assume t0 and dt are consistent or scalar and vector
  coherence data */
  linear_vector_regular_to_irregular(cohens->member[is].t0,
    cohens->member[is].dt, raycoh,&(SPtime[0]),raycohout);
  /* This nontrivial effort may not be required, but best to
  have it available until proven to be a big issue. */
  raycohout=running_average(raycohout,cohsl);
  /* Now copy these results to the coherence raygrid.
  Note carefully the reverse order going from raycohout order
  to rayrid order and decimation by cohdecfac */
  for(k=0,kk=pwcohgrid.n3-1;(k<nout) && (kk>=0);
    k+=cohdecfac,--kk)
  {
    for(l=0;l<4;++l)
      pwcohgrid.val[i][j][kk][l]=raycohout(l,k);
  }

  //
  // Compute the transformation vector for each
  // ray point in space.  This comes in two flavors.
  // A fast, inexact version and a exact slow version
  //
  dmatrix *pathptr=extract_gridline(raygrid,i,j,raygrid.n3-1,3,true);
  if(use_depth_variable_transformation)
  {
    // the slow, exact method
    troptr=new Ray_Transformation_Operator(
      parent,
      *pathptr,
      ustack.azimuth(),
      nup);
  }
  else
  {
    // the fast approximation
    troptr=new Ray_Transformation_Operator(
      parent,
      *pathptr,ustack.azimuth());
  }
  Ray_Transformation_Operator& trans_operator=*troptr;
  work=trans_operator.apply(work);
  // done with these now
  delete pathptr;
  delete troptr;

  // This computes domega for a ray path in constant dz mode
  // We would have to interpolate anyway to mesh with
  // time data choose the faster, coarser dz method here
  //
  if(rcomp_wt || !weight_functions_set)
  {
    /* This is a gross inefficiency in stack_only but since I only expect it to be
    used for CCP stacking equivalent, this should not be a big deal.  Probably should
    do it right some day */
    if(stack_only)
    {
      domega_ij.clear();
      for(k=0;k<n3;++k) domega_ij.push_back(1.0);
    }
    else
    {
      /* This is the normal block for computing solid angles */
      domega_ij=compute_domega_for_path(ustack,dux,duy,
        Vs1d, zmaxray,dz,
        raygrid, i, j, gradTp,zP);
    }
    if(use_grt_weights && (!stack_only) )
      dweight_ij=compute_weight_for_path(gradTp,gradTs);
    else
      for(k=0;k<n3;++k)dweight_ij[k]=1.0;
    if(smooth_wt && (!stack_only))
    {
      domega_ij=running_average(domega_ij,nwtsmooth);
      // Unnecessary when not using grt weighting
      if(use_grt_weights)
        dweight_ij=running_average(dweight_ij,nwtsmooth);
    }
    weight_functions_set=true;
  }

  //
  // copy transformed data to vector field
  // copy weights and domega at same time
  //
  for(k=0,kk=raygrid.n3-1;k<raygrid.n3;++k,--kk)
  {
    for(l=0;l<3;++l)
    {
      pwdgrid.val[i][j][k][l]=work(l,kk)
        *dweight_ij[kk]*domega_ij[kk];
    }
    pwdgrid.val[i][j][k][3]=domega_ij[kk];
    pwdgrid.val[i][j][k][4]=dweight_ij[kk];
  }
  }
  delete pwdata;
