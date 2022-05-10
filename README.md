# ACCES_PAR
Scripts to estimate PAR reaching the coastal Arctic seafloor.

*For implementation of [SSP atmospheric correction](https://doi.org/10.1364/OE.27.0A1118) in [SeaDAS](https://seadas.gsfc.nasa.gov/), check [this repository](https://github.com/rakeshkstp/AtmosphericCorrection).*

* Add the following files to `$OCSSWROOT/ocssw_src/src/l2gen/CMakeLists.txt`
  * get_uqar_utils.c
  * calc_uqar_utils.c
  * interpol_ed0LUT_5nm_v2.f
  * calc_par_uqar.f
  * get_ed0_LUT_v2.f

## Author
**Rakesh Kumar Singh**

**[Simon Bélanger](https://github.com/belasi01)**
