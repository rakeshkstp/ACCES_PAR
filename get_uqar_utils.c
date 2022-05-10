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

#include "uqar_header.h"
/**
* @file get_PAR_UQAR.c
* @author Rakesh Kumar Singh
* @date 22 March 2021
* @copyright 2021 Rakesh Kumar Singh
* @brief Main file to calculate PAR using Laliberté et al. 2016
*/

void get_uqar_utils(l2str *l2rec, int prodnum, float prod[])
{   /**
      * @brief Main function to calculate PAR using Laliberté et al. 2016
      * @param [in] l2rec Level 2 record in SeaDAS
      * @return par, calculated PAR for the pixel
      */

    int32_t ip;
    int32_t npix = l2rec->l1rec->npix;
    l1str *l1rec = l2rec->l1rec;    /*!< L1 record of OCSSW */
    filehandle *l1file = l1rec->l1file; /*!< L1 file pointer */
    static int firstCall = 0;   /*!< Is it the first call to this function? Initial value FALSE */
    static int nbands; /*!< Number of bands available in the sensor */
    static int bstart, bstop;   /*!< Index for the first and last band to be used in PAR calculation. Bands correponding to 400 nm and 700 nm. */
    static int blue, green, red, nirl, swir;
    static float LUT0p[NWL][NTHETA][NO3][NCOT][NALB];
    static float LUT0m[NWL][NTHETA][NO3][NCOT][NALB];
    static int32_t mask = LAND;
    static int16_t year, month, mday, doy;
    static double sec;
    float trise, tset, deltaT;
    int step=10;
    int status = 1, choice =1;
    int ICW = -1;
    float COT, O3, salb, cldtr;
    float z, par0m;

    if (!firstCall)
      {   firstCall=1;
          /* Reading Ed0+ LUT (Laliberté et al., 2016)*/
          status = get_ed0_(LUT0p, &choice); //1 for PAR0+ and 0 for PAR0-
          if (status)
            exit(1);

          /* Reading Ed0- LUT (Laliberté et al., 2016)*/
          choice = 0;
          status = get_ed0_(LUT0m, &choice); //1 for PAR0+ and 0 for PAR0-
          if (status)
            exit(1);

          /* Extract Day, Month, Year from Scantime */
          unix2yds(l1rec->scantime, &year, &doy, &sec);
          unix2ymds(l1rec->scantime, &year, &month, &mday, &sec);

          /* Get all the PAR bands in the sensor (400 < Lambda < 700) */
          nbands = l1file->nbands;
          get_bands(l1rec, &bstart, &bstop, &blue, &green, &red, &nirl, &swir);
      }

    for (ip = 0; ip < npix; ip++)
    {   /* Get Depth */
        z = get_elev(l1rec->lat[ip], l1rec->lon[ip]);

        /* Skip pixel if masked */
        if (l1rec->Lt[ip * nbands] <= 0.0 || (l1rec->flags[ip] & mask) != 0
                || l1rec->solz[ip] > 83.0 || z >= 0)
        {   prod[ip] = BAD_FLT;
            l1rec->flags[ip] |= PRODFAIL;
            continue;
        }

        /* Get the rise and set time for this day */
        triseset(doy, l1rec->lon[ip], l1rec->lat[ip], &trise, &tset);

        /* Set up delta X for 10 intervals through the day in hours */
        deltaT = (tset - trise) / (float)step;

        /* Check if pixel has ice cloud or water */
        ICW = calc_uqar_icw(l1rec, ip, nbands, blue, green, nirl, swir);

        /* Calculate Ozone Optical Thickness (Unit: DU) */
        O3 = calc_O3(l1rec, ip, year, month, mday);

        /* Calculate surface albedo */
        salb = calc_salb(l2rec, ip, nbands, blue, green, red, nirl, swir, ICW, doy);

        /* Calculate Cloud Optical Thickness (Range: 0 - 64) */
        COT = calc_COT(l2rec, ip, nbands, &cldtr, blue, green, red, nirl, swir, ICW, salb, doy);

		  switch (prodnum)
      {  case CAT_uqar_par0p:
              prod[ip] = calc_par_surf(ip, l1rec, year, doy, LUT0p, trise, deltaT, step, O3, COT, salb, 1);
              break;
        case CAT_uqar_par0m:
              prod[ip] = calc_par_surf(ip, l1rec, year, doy, LUT0m, trise, deltaT, step, O3, COT, salb, 0);
              break;
        case CAT_uqar_parb:
              par0m=calc_par_surf(ip, l1rec, year, doy, LUT0m, trise, deltaT, step, O3, COT, salb, 0);
              prod[ip] = calc_par_z(l2rec, ip, nbands, bstart, bstop, par0m, -z);
              break;
        case CAT_uqar_kdpar:
              prod[ip] = calc_kdpar_uqar(l2rec, ip, nbands, bstart, bstop);
              break;
        case CAT_uqar_COT:
              prod[ip] = COT;
              break;
        case CAT_uqar_O3:
              prod[ip] = O3;
              break;
        case CAT_uqar_salb:
              prod[ip] = salb;
              break;
        case CAT_uqar_icw:
              prod[ip] = ICW;
              break;
        default:
              printf("Unknown product id: %d \n", prodnum);
              break;
      }
    }
}
