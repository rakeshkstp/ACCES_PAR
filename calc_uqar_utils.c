#include "uqar_header.h"

/********************************/
/* Get PAR bands for the sensor */
/********************************/
void get_bands(l1str *l1rec, int *bstart, int *bstop, int *blue, int *green, int *red, int *nirl, int *swir)
{   filehandle *l1file = l1rec->l1file;
    int nbands = l1file->nbands;
    int bstrt, bstp, r, indx, nir, gr, bl, sw;
    bstrt=-1; bstp=-1;
    for (indx = 0; indx < nbands; indx++)
        if (l1file->fwave[indx] > 400.0 && bstrt<0)
          {	bstrt=indx; continue;}
    for (indx = nbands-1; indx >= 0; indx--)
        if (l1file->fwave[indx] < 700.0 && bstp<0)
          {	bstp=indx;  continue;}
    if (bstrt>=0 && bstp>=0 && bstrt<=bstp)
      printf("\nCalculating PAR using %d bands starting from %d nm to %d nm.\n",bstp-bstrt+1,(int)l1file->fwave[bstrt],(int)l1file->fwave[bstp]);
    else
      {	printf("%s line %d: Can't find PAR bands for this sensor.\n",__FILE__,__LINE__);
        exit(1);
      }

    if ((bl = bindex_get(469)) < 0)
      if ((bl = bindex_get(490)) < 0)
        if ((bl = bindex_get(491)) < 0)
          if ((bl = bindex_get(482)) < 0)
            if ((bl = bindex_get(486)) < 0)
              if ((bl = bindex_get(488)) < 0)
              { printf("%s line %d: can't find blue band\n", __FILE__, __LINE__);
                exit(1);
              }
    if ((gr = bindex_get(555)) < 0)
      if ((gr = bindex_get(550)) < 0)
        if ((gr = bindex_get(560)) < 0)
          if ((gr = bindex_get(565)) < 0)
          { printf("%s line %d: can't find green band\n", __FILE__, __LINE__);
            exit(1);
          }
    if ((r = bindex_get(645)) < 0)
      if ((r = bindex_get(671)) < 0)
        if ((r = bindex_get(660)) < 0)
          if ((r = bindex_get(665)) < 0)
            if ((r = bindex_get(655)) < 0)
              if ((r = bindex_get(620)) < 0)
              {	printf("%s line %d: Can't find red band for thsi sensor.\n",__FILE__,__LINE__);
                exit(1);
              }
    if ((nir = bindex_get(859)) < 0)
      if ((nir = bindex_get(865)) < 0)
        if ((nir = bindex_get(867)) < 0)
          if ((nir = bindex_get(866)) < 0)
          {	printf("%s line %d: Can't find nirl band for thsi sensor.\n",__FILE__,__LINE__);
            exit(1);
          }
    if ((sw = bindex_get(2130)) < 0)
      if ((sw = bindex_get(2201)) < 0)
      {	printf("%s line %d: Can't find SWIR band for this sensor.\n",__FILE__,__LINE__);
        exit(1);
      }
    *bstart=bstrt;
    *bstop=bstp;
    *red=r;
    *nirl=nir;
    *green=gr;
    *blue=bl;
    *swir=sw;
}

/***************************************/
/* Calculating Ozone Optical thickness */
/***************************************/

float calc_O3(l1str *l1rec, int32_t ip, int16_t year, int16_t month, int16_t mday)
{	float O3;

	/*Reading O3 information from ancillary data, if available*/
		O3 = l1rec->oz[ip]*1000.0;

	/*If O3 information is not available from ancillary data,
	 *using climatology*/
	if (O3 < 100.0f)
		O3 = EstimateDobson(year, month, mday, l1rec->lat[ip]);
	/*DEfining boundaries such that O3 values does not exeed the LUT bounds*/
  if (O3 < 200.0)
    O3 = 200.01;
  if (O3 > 550.0)
    O3 = 549.99;

	return(O3);
}

/*******************************/
/* Detect Ice, Cloud and Water */
/*******************************/
int calc_uqar_icw(l1str *l1rec, int32_t ip, int nbands, int blue, int green, int nirl, int swir)
{ filehandle *l1file = l1rec->l1file;
  float *rhos = l1rec->rhos;
  float solz = l1rec->solz[ip];
  int icw;
  float ngb, nns;
  float rhosb, rhosg, rhosn, rhoss;
  float k;
  float intercept;

if (solz > 83.0)
	return(-1);

  rhosb = rhos[ip*nbands+blue];
  rhosg = rhos[ip*nbands+green];
  rhosn = rhos[ip*nbands+nirl];
  rhoss = rhos[ip*nbands+swir];

  intercept=(l1file->fwave[nirl]*rhoss - l1file->fwave[swir]*rhosn)/(l1file->fwave[nirl]-l1file->fwave[swir]);

  if (rhosg > rhosb)
    k = rhosg / rhosb;
  else
    k = 1.0;

  ngb = (rhosg - rhosb) / (rhosg + rhosb);

  nns = (rhosn - rhoss) / (rhosn + rhoss);

  icw =0;   //water
  if (intercept>0.1 && ngb < 0.1)
      icw=2;    //cloud
  //if (nns / k > 0.7 && rhosb > 0.06)
  if (nns / k > 0.6 && rhosb > 0.12)
      icw=1;    //ice

  return(icw);
}

/******************************/
/* Calculating Surface albedo */
/******************************/

float calc_salb(l2str *l2rec, int32_t ip, int nbands, int blue, int green,int red, int nirl, int swir, int ICW, int doy)
{	static float p0 = STDPR;
  static float pi = PI;
  l1str *l1rec = l2rec->l1rec;
  filehandle *l1file = l1rec->l1file;
  float It=0, Isen=0;
  float mu = l1rec->csenz[ip];
  float *Fo = l1rec->Fo;
  float *t_sol = l1rec->t_sol;
  float *t_sen = l1rec->t_sen;
  float *tg_sen = l1rec->tg_sen;
  float *tg_sol = l1rec->tg_sol;
  float *Lt = l1rec->Lt;
  float *Lr = l1rec->Lr;
  int32_t ipb;
  int ib;
  float salb=0;
  float mu0=l1rec->csolz[ip];
  float *taua = &l2rec->taua [ip * l1file->nbands];
  float taur[nbands];
  float tdbar, Tdbar, E0, sum_E0;
  int parbands[3]={blue,green,red};
  int i;
  float Aice;
  float SIC;

if(ICW == -1)
    return(BAD_FLT);

/* Get Sea ice fraction */
  SIC = ice_fraction(l1rec->lon[ip], l1rec->lat[ip]); //old or nsidc or oisst

 /* Coakley, J. A. (2003). Reflectance and albedo, Surface.
	* In J. R. Holton, J. A. Curry, & J. A. Pyle (Eds.),
	* Encyclopedia of Atmospheric Sciences (pp. 1914–1923).
	* https://doi.org/10.1016/B0-12-227090-8/00069-5 */

  /*Frouin, R.J., Franz, B.A., Werdell, P.J., 2003. The SeaWiFS PAR Product,
  in: Hooker, S.B., Firestone, E.R. (Eds.), SeaWiFS Postlaunch Technical Report Series, Volume 22,
  Algorithm Updates for the Fourth SeaWiFS Data Reprocessing.
  NASA Technical Memorandum 2003, Greenbelt, Maryland, pp. 46–50.*/

  if (ICW != 1)   //  cloud or water
  { tdbar=0;
    Tdbar=0;
    sum_E0=0;
    for(i=0; i<3; i++)
    { ib = parbands[i];
      taur[ib] = l1rec->pr[ip] / p0 * l1file->Tau_r[ib];
      E0=Fo[ib] * mu0;
      if (taua[ib] < 0.0 || taua[ib] > 1.0 || taua[ib] != taua[ib]) // taua is not valid
      { tdbar += exp(-(taur[ib])/mu0) * exp(0.52*taur[ib]/mu0) * E0;
        Tdbar += exp(-(taur[ib])/mu0) * E0;
        sum_E0+= E0;
      }
      else
      { tdbar += exp(-(taur[ib]+taua[ib])/mu0) * exp((0.52*taur[ib] + 0.83*taua[ib])/mu0) * E0;
        Tdbar += exp(-(taur[ib]+taua[ib])/mu0) * E0;
        sum_E0+= E0;
      }
    }

    tdbar/=sum_E0;
    Tdbar/=sum_E0;

    salb = (0.05/(1.1*pow(mu0,1.4)+0.15))*(Tdbar/tdbar) + 0.08*(1-(Tdbar/tdbar)); //Frouin et al., 2003

    if (ICW == 2)   //  cloud
    /* Perovich, D.K., Nghiem, S. V., Markus, T., Schweiger, A., 2007.
     * Seasonal evolution and interannual variability of the local solar
     * energy absorbed by the Arctic sea ice–ocean system. J. Geophys. Res.
     * 112, C03005. https://doi.org/10.1029/2006JC003558 */
    { Aice = 1.25/(1.25 + exp(-(pow((float)doy-208.0,2.0))/(1250.0)))-0.154;
      //Aice=0.864;
      salb=salb*(1.0 - SIC) + Aice * SIC;
    }
  }
  else // ice
  { sum_E0=0;
    salb=0;
    for(i=0; i<3; i++)
    { ib = parbands[i];
      ipb = ip * nbands + ib;
      E0=Fo[ib] * mu0;
      It = E0 * t_sol[ipb] * tg_sol[ipb] - (Lr[ipb] * pi * mu / t_sen[ipb] / tg_sen[ipb]);
      Isen = (Lt[ipb] - Lr[ipb]) * pi * mu / t_sen[ipb] / tg_sen[ipb];
      salb += (Isen/It) * E0;
      sum_E0+= E0;
    }
    salb/=sum_E0;
  }

  return(salb);
}


/***************************************/
/* Calculating Cloud Optical thickness */
/***************************************/
float calc_COT(l2str *l2rec, int32_t ip, int nbands, float *cldtr, int red, int blue, int green, int nirl, int swir, int ICW, float salb_in, int doy)
{ float t, COT;
  int32_t ipb;
  float I0, Ir;
  static float p0 = STDPR;
  l1str *l1rec = l2rec->l1rec;
  float mu0 = l1rec->csolz[ip];
  float mu = l1rec->csenz[ip];
  float *Fo = l1rec->Fo;
  //float fsol = l1rec->fsol;
  float *t_sol = l1rec->t_sol;
  float *t_sen = l1rec->t_sen;
  float *tg_sen = l1rec->tg_sen;
  float *tg_sol = l1rec->tg_sol;
  float *Lt = l1rec->Lt;
  float *Lr = l1rec->Lr;
  int cband = red;
  float salb=0;
  float *taua = &l2rec->taua [ip * l1rec->l1file->nbands];
  float taur;
  float td, Td, E0;
  float Aice;
  float SIC = ice_fraction(l1rec->lon[ip], l1rec->lat[ip]); //old or nsidc or oisst


E0=Fo[cband] * mu0;
ipb = ip * nbands + cband;

if(ICW == -1)
	return(BAD_FLT);

  if (ICW != 1)   //  cloud or water
  {   taur = l1rec->pr[ip] / p0 * l1rec->l1file->Tau_r[cband];
      if (taua[cband] < 0.0 || taua[cband] > 1.0 || taua[cband] != taua[cband]) // taua is not valid
      { td = exp(-(taur)/mu0) * exp(0.52*taur/mu0);
        Td = exp(-(taur)/mu0);
      }
      else
      { td = exp(-(taur+taua[cband])/mu0) * exp((0.52*taur + 0.83*taua[cband])/mu0);
        Td = exp(-(taur+taua[cband])/mu0);
      }

      salb = (0.05/(1.1*pow(mu0,1.4)+0.15))*(Td/td) + 0.08*(1.0-(Td/td)); //Frouin et al., 2003

      if (ICW == 2)   //  cloud
      /* Perovich, D.K., Nghiem, S. V., Markus, T., Schweiger, A., 2007.
        * Seasonal evolution and interannual variability of the local solar
        * energy absorbed by the Arctic sea ice–ocean system. J. Geophys. Res.
        * 112, C03005. https://doi.org/10.1029/2006JC003558 */
      { Aice = 1.25/(1.25 + exp(-(pow((float)doy-208.0,2.0))/(1250.0)))-0.154;
        I0 = E0 * t_sol[ipb] * tg_sol[ipb] - Lr[ipb] * PI * mu / t_sen[ipb] / tg_sen[ipb];
        Ir = (Lt[ipb] - Lr[ipb]) * PI * mu / t_sen[ipb] / tg_sen[ipb];
     
        if (salb*(1.0 - SIC) + Aice * SIC <= Ir/I0)
          salb=salb*(1.0 - SIC) + Aice * SIC;
      }
  }
  else // ice
  { I0 = E0 * t_sol[ipb] * tg_sol[ipb] - Lr[ipb] * PI * mu * t_sen[ipb] * tg_sen[ipb];
    Ir = (Lt[ipb] - Lr[ipb]) * PI * mu * t_sen[ipb] * tg_sen[ipb];
    salb = (Ir/I0);
  }

  if (ICW != 1)
  {
  /* Calculating cloud transmittance */
  I0 = E0 * t_sol[ipb] * tg_sol[ipb] - Lr[ipb] * PI * mu * t_sen[ipb] * tg_sen[ipb];
  Ir = Lt[ipb] * PI * mu * t_sen[ipb] * tg_sen[ipb];
  t=(E0-sqrt(E0*E0-4.0*I0*salb*(E0-Ir)))/(2.0*I0*salb);
  }
  else
    t=1.0;

  if (t>0.93)
    t=0.93;

*cldtr=t;

	/* Pandey, P., De Ridder, K., Gillotay, D., van Lipzig, N.P.M., 2012.
	 * Estimating cloud optical thickness and associated surface UV irradiance
	 * from SEVIRI by implementing a semi-analytical cloud retrieval algorithm.
	 * Atmos. Chem. Phys. 12, 7961–7975. https://doi.org/10.5194/acp-12-7961-2012
	 * COT = (1/0.75)*(1/t-1.07)/(1-g)*/

   /*	Yang, P., Baum, B.A., 2003. SATELLITE REMOTE SENSING | Cloud Properties, 
   	in: Holton, J.R., Curry, J.A., Pyle, J.A. (Eds.), Encyclopedia of Atmospheric Sciences. 
   	Elsevier, pp. 1956–1965. https://doi.org/10.1016/B0-12-227090-8/00348-1
   	0.65um is sensitive primarily to COT */
 COT = ((1.0/0.75)*(1.0/t-1.07))/(1.0-0.85);

  return(COT);
}

/********************************************************/
/* Calculating Kd                                       */
/* Lee, Z., Hu, C., Shang, S., Du, K., Lewis, M.,       */
/* Arnone, R., & Brewin, R. (2013). Penetration of      */
/* UV-visible solar radiation in the global oceans:     */
/* Insights from ocean color remote sensing. Journal    */
/* of Geophysical Research: Oceans, 118(9), 4241–4255.  */
/* https://doi.org/10.1002/jgrc.20308                   */
/********************************************************/
int calc_uqar_Kd(l2str *l2rec, int32_t ip, int nbands, int bstart, int bstop, int aer_b1, int aer_b2, int aer_g, int red, float parwave[], float Kd[])
{	l1str *l1rec = l2rec->l1rec;
  float rrs[nbands],u[nbands],a[nbands],bbp[nbands];
  float chi,etaw;
  int ib;
  float a_ref,bbp_ref;
  int ref;
  float eta;
  float solz = l1rec->solz[ip];
  float *Rrs=&l2rec->Rrs[ip * nbands];

  if (Rrs[bstart]<0)
    return(-1);

  for(ib = bstart; ib <= bstop; ib++)
  { rrs[ib] = Rrs[ib]/(0.52 + 1.7 * Rrs[ib]);
    u[ib] = (-0.089 + sqrt(0.089 * 0.089 + 4.0 * 0.125 * rrs[ib])) / (2.0 * 0.125);
  }

  //if (Rrs[red] < 0.0015)
  { ref = aer_g;
    chi = log10((rrs[aer_b1] + rrs[aer_b2]) / (rrs[aer_g] + 5.0 * (rrs[red] / rrs[aer_b2]) * rrs[red]));
    a_ref = l1rec->sw_a[ref] + pow(10.0, -1.146 - 1.366 * chi - 0.469 * chi * chi);
    bbp_ref = ((u[ref] * a_ref) / (1.0 - u[ref])) - l1rec->sw_bb[ref];
  }
  /*else
  { ref = red;
    a_ref = aw[ref] + 0.39 * pow(Rrs[red] / (Rrs[aer_b1] + Rrs[aer_b2]), 1.14);
    bbp_ref = ((u[ref] * a_ref) / (1.0 - u[ref])) - bbw[ref];
  }*/

  bbp_ref=fabs(bbp_ref);
  a_ref=fabs(a_ref);

  eta = 2.0 * (1.0 - 1.2 * exp(-0.9 * (rrs[aer_b1] / rrs[ref])));

  for(ib = bstart; ib <= bstop; ib++)
  { bbp[ib] = bbp_ref * pow(parwave[ref] / parwave[ib], eta);
    a[ib] = (1.0 - u[ib]) * (l1rec->sw_bb[ib] + bbp[ib])/u[ib];
    etaw = l1rec->sw_bb[ib] / (bbp[ib] + l1rec->sw_bb[ib]);
    Kd[ib] = (1.0 + 0.005 * solz) * a[ib] + (1.0 - 0.265 * etaw) * 4.259 * (1.0 - 0.52 * exp(-10.8 * a[ib])) * (bbp[ib] + l1rec->sw_bb[ib]);
    if (a[ib]<0)
      return(-1);
    }
	return(1);
}


/*******************************/
/* Calculate surface PAR       */
/*******************************/
float calc_par_surf(int32_t ip, l1str *l1rec, int16_t year, int16_t doy, float LUT[NWL][NTHETA][NO3][NCOT][NALB], float trise, float deltaT, int step, float O3, float COT, float salb, int choice)
{ float solz, sola;
  float par_temp[step+1];
  int intyear = year, intday = doy;
  int indx;
  float tmstep;
  float Ed[NWL];
  float par=0;

  /* Calculate PAR at each time step */
  for (indx = 0; indx <= step; indx++)
    {	tmstep = trise + indx * deltaT;
      sunangs_(&intyear, &intday, &tmstep, l1rec->lon+ip,l1rec->lat+ip, &solz, &sola); //Getting solar zenith angle for the time of the day
      if (solz >= 90.0 || solz < 0.0)
        { par_temp[indx] = 0.0;
          continue;
        }
      else
        { par_temp[indx] = calc_par_uqar_(&solz,&O3,&COT,&salb,LUT,Ed);	//The unit for PAR is muE m^-2 s^-1
          if (par_temp[indx]>0)
            par_temp[indx] *= 0.0036;		//The unit for PAR converted to E m^-2 hr^-1
          else
            par_temp[indx] = 0.0;

          if (!choice)
            par_temp[indx] *= (1.0 - salb);
        }
      }

   /* Integrate PAR through the day to get Daily PAR */
   for (indx = 0; indx < step; indx++)
      par += deltaT * (par_temp[indx] + par_temp[indx + 1]) / 2.0f;

   return(par);
}

/*******************************/
/* Calculate PAR at depth z    */
/*******************************/
float calc_par_z(l2str *l2rec, int32_t ip, int nbands, int bstart, int bstop, float par0m, float depth)
{   float parz;
    float kdpar;

    kdpar = calc_kdpar_uqar(l2rec, ip, nbands, bstart, bstop);
    parz = par0m * exp(-kdpar*depth);

    return(parz);
}


/*******************************/
/* Calculate surface iPAR      */
/*******************************/
/*	1. Wavelength = 290 : 700 : 5
	2. ThetaS = 0 : 90 : 5
	3. Ozone = 200 : 550 : 50
	4. Cloud optical Thickness = 0 to 64 = c(0,1,2,4,8,16,32,64)
	5. Surface Albedo = 0.05 : 0.95 : 0.15
*/

/* Laliberté, J., Bélanger, S., Frouin, R.J., 2016.
 * Evaluation of satellite-based algorithms to estimate
 * photosynthetically available radiation (PAR) reaching the
 * ocean surface at high northern latitudes.
 * Remote Sens. Environ. 184, 199–211.
 * https://doi.org/10.1016/j.rse.2016.06.014 */
float calc_ipar_surf(int32_t ip, l1str *l1rec, int16_t year, int16_t doy, double sec, float LUT[NWL][NTHETA][NO3][NCOT][NALB], float O3, float COT, float salb, int choice)
{ float solz, sola;
  int intyear = year, intday = doy;
  float Ed[NWL];
  float ipar=0.0;
  //float icefrac;
  float hour = sec/3600.0;

  /* Calculate iPAR */
  sunangs_(&intyear, &intday, &hour, l1rec->lon+ip,l1rec->lat+ip, &solz, &sola); //Getting solar zenith angle for the time of the day
  if (solz >= 90.0 || solz < 0.0)
    ipar = 0.0;
  else
  { ipar = calc_par_uqar_(&solz,&O3,&COT,&salb,LUT,Ed);	//The unit for PAR is muE m^-2 s^-1
    if (ipar>0)
      ipar /= 1000000.0;		//The unit for PAR converted to E m^-2 s^-1
    else
      ipar = 0.0;

    if (!choice)
      ipar *= (1.0 - salb);
  }

 return(ipar);
}


 /*********************************************************/
 /* Calculate Kd PAR                                      */
 /* Saulquin, B., Hamdi, A., Gohin, F., Populus, J.,      */
 /* Mangin, A., & D’Andon, O. F. (2013). Estimation of    */
 /* the diffuse attenuation coefficient KdPAR using MERIS */
 /* and application to seabed habitat mapping.            */
 /* Remote Sensing of Environment, 128, 224–233.          */
 /* https://doi.org/10.1016/j.rse.2012.10.002             */
 /*********************************************************/
 float calc_kdpar_uqar(l2str *l2rec, int32_t ip, int nbands, int bstart, int bstop)
 {  /********************************************************/
    /* Calculating Kd                                       */
    /* Lee, Z., Hu, C., Shang, S., Du, K., Lewis, M.,       */
    /* Arnone, R., & Brewin, R. (2013). Penetration of      */
    /* UV-visible solar radiation in the global oceans:     */
    /* Insights from ocean color remote sensing. Journal    */
    /* of Geophysical Research: Oceans, 118(9), 4241–4255.  */
    /* https://doi.org/10.1002/jgrc.20308                   */
    /********************************************************/
    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;
    float rrs[nbands],u[nbands],a[nbands],bbp[nbands];
    float chi,etaw;
    int ib;
    float a_ref,bbp_ref;
    int ref;
    float eta;
    float solz = l1rec->solz[ip];
    float *Rrs=&l2rec->Rrs[ip * nbands];
    float kdpar,kd490;
    static int firstCall = 0;
    static int b1, b2, g, r; //b1=443, b2=490, g=55x and r=66x

    if (!firstCall)
    {   firstCall=1;

        if ((b1 = bindex_get(443)) < 0)
          if ((b1 = bindex_get(442)) < 0)
            {	printf("%s line %d: Can't find 44x nm band for this sensor.\n",__FILE__,__LINE__);
              exit(1);
            }
        if ((b2 = bindex_get(490)) < 0)
          if ((b2 = bindex_get(491)) < 0)
            if ((b2 = bindex_get(482)) < 0)
              if ((b2 = bindex_get(486)) < 0)
                if ((b2 = bindex_get(488)) < 0)
                { printf("%s line %d: can't find 490 nm band ... Cannot calculate Kd490\n", __FILE__, __LINE__);
                  exit(1);
                }
        if ((g = bindex_get(547)) < 0)
          if ((g = bindex_get(555)) < 0)
          {	printf("%s line %d: Can't find 55x nm band for this sensor.\n",__FILE__,__LINE__);
            exit(1);
          }
        if ((r = bindex_get(667)) < 0)
          if ((r = bindex_get(671)) < 0)
            if ((r = bindex_get(660)) < 0)
              if ((r = bindex_get(665)) < 0)
                if ((r = bindex_get(655)) < 0)
                  if ((r = bindex_get(620)) < 0)
                  {	printf("%s line %d: Can't find 66x nm band for thsi sensor.\n",__FILE__,__LINE__);
                    exit(1);
                  }
    }

    if (Rrs[bstart]<0)
      return(BAD_FLT);

    for(ib = bstart; ib <= bstop; ib++)
    { rrs[ib] = Rrs[ib]/(0.52 + 1.7 * Rrs[ib]);
      u[ib] = (-0.089 + sqrt(0.089 * 0.089 + 4.0 * 0.125 * rrs[ib])) / (2.0 * 0.125);
    }

    ref = g;
    chi = log10((rrs[b1] + rrs[b2]) / (rrs[g] + 5.0 * (rrs[r] / rrs[b2]) * rrs[r]));
    a_ref = l1rec->sw_a[ref] + pow(10.0, -1.146 - 1.366 * chi - 0.469 * chi * chi);
    bbp_ref = ((u[ref] * a_ref) / (1.0 - u[ref])) - l1rec->sw_bb[ref];

    bbp_ref=fabs(bbp_ref);
    a_ref=fabs(a_ref);

    eta = 2.0 * (1.0 - 1.2 * exp(-0.9 * (rrs[b1] / rrs[ref])));

    bbp[b2] = bbp_ref * pow(l1file->fwave[ref] / l1file->fwave[b2], eta);
    a[b2] = (1.0 - u[b2]) * (l1rec->sw_bb[b2] + bbp[b2])/u[b2];
    etaw = l1rec->sw_bb[b2] / (bbp[b2] + l1rec->sw_bb[b2]);
    kd490 = (1.0 + 0.005 * solz) * a[b2] + (1.0 - 0.265 * etaw) * 4.259 * (1.0 - 0.52 * exp(-10.8 * a[b2])) * (bbp[b2] + l1rec->sw_bb[b2]);

    if (a[b2]<0 || kd490 <=0)
      return(BAD_FLT);

    if (kd490 > 0.115)
      kdpar=0.81*pow(kd490,0.8256);
    else
      kdpar=4.6051*kd490/(6.07*kd490+3.2);

  	return(kdpar);
 }


 /*******************************/
 /* Calculate isolume           */
 /*******************************/
 /* Letelier, R.M., Karl, D.M., Abbott, M.R., Bidigare, R.R., 2004. Light driven
  * seasonal patterns of chlorophyll and nitrate in the lower euphotic zone of the
  * North Pacific Subtropical Gyre. Limnol. Oceanogr. 49, 508–519.
  * https://doi.org/10.4319/lo.2004.49.2.0508
  */
 float calc_isolume_uqar(int32_t ip, l2str *l2rec, int16_t year, int16_t doy, float LUT[NWL][NTHETA][NO3][NCOT][NALB], float trise, float deltaT, int step, float O3, float COT, float salb, int nbands, int bstart, int bstop, float depth)
 { float isolume = BAD_FLT;
   float par0m, kdpar;
   l1str *l1rec = l2rec->l1rec;

   par0m=calc_par_surf(ip, l1rec, year, doy, LUT, trise, deltaT, step, O3, COT, salb, 0);
   kdpar=calc_kdpar_uqar(l2rec, ip, nbands, bstart, bstop);

   if (kdpar > 0 && depth > 0)
     isolume=(-1.0/kdpar)*log(0.415/par0m);
   else
    isolume = BAD_FLT;

  if (isolume>depth)
    isolume=depth;

     return(isolume);
 }

 /***************************************/
 /* Calculate depth for PAR 1% and 10%  */
 /***************************************/
 float calc_dPAR_uqar(l2str *l2rec, int32_t ip, int nbands, int bstart, int bstop, float depth, float pcnt)
 { float pdepth;
   float kdpar;

   kdpar=calc_kdpar_uqar(l2rec, ip, nbands, bstart, bstop);

   if (kdpar > 0 && depth > 0)
      pdepth=(-1.0/kdpar)*log(pcnt/100.0);
   else
      pdepth=BAD_FLT;

   if(pdepth>depth)
    pdepth=depth;

   return(pdepth);
 }
