#include "l2_struc.h"
#include <timeutils.h>
#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"

#define NWL 84
#define NTHETA 20
#define NO3 11
#define NCOT 9
#define NALB 8

int get_ed0_(float LUT[NWL][NTHETA][NO3][NCOT][NALB], int *choice);

void get_bands(l1str *l1rec, int *bstart, int *bstop, int *blue, int *green, int *red, int *nirl, int *swir);

float calc_salb(l2str *l2rec, int32_t ip, int nbands, int blue, int green, int red, int nirl, int swir, int ICW, int doy);

float calc_O3(l1str *l1rec, int32_t ip, int16_t year, int16_t month, int16_t mday);

float calc_COT(l2str *l2rec, int32_t ip, int nbands, float *cldtr, int red, int blue, int green, int nirl, int swir, int ICW, float salb, int doy);

float calc_par_uqar_(float *solz, float *O3, float *COT, float *salb, float LUT[NWL][NTHETA][NO3][NCOT][NALB], float Ed[NWL]);

int calc_uqar_Kd(l2str *l2rec, int32_t ip, int nbands, int bstart, int bstop, int aer_b1, int aer_b2, int aer_g, int red, float parwave[], float Kd[]);

int calc_uqar_icw(l1str *l1rec, int32_t ip, int nbands, int blue, int green, int nirl, int swir);

float calc_par_surf(int32_t ip, l1str *l1rec, int16_t year, int16_t doy, float LUT[NWL][NTHETA][NO3][NCOT][NALB], float trise, float deltaT, int step, float O3, float COT, float salb, int choice);

float calc_par_z(l2str *l2rec, int32_t ip, int nbands, int bstart, int bstop, float par0m, float depth);

float calc_ipar_surf(int32_t ip, l1str *l1rec, int16_t year, int16_t doy, double sec, float LUT[NWL][NTHETA][NO3][NCOT][NALB], float O3, float COT, float salb, int choice);

float calc_kdpar_uqar(l2str *l2rec, int32_t ip, int nbands, int bstart, int bstop);

float calc_isolume_uqar(int32_t ip, l2str *l2rec, int16_t year, int16_t doy, float LUT[NWL][NTHETA][NO3][NCOT][NALB], float trise, float deltaT, int step, float O3, float COT, float salb, int nbands, int bstart, int bstop, float depth);

float calc_dPAR_uqar(l2str *l2rec, int32_t ip, int nbands, int bstart, int bstop, float depth, float pcnt);
