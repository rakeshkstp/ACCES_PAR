# ACCES_PAR
Scripts to estimate PAR reaching the coastal Arctic seafloor.

*Note: These scripts are tested with SeaDAS v2020.1, and may need minor changes to work with newer versions.*

*For implementation of [SSP atmospheric correction](https://doi.org/10.1364/OE.27.0A1118) in [SeaDAS](https://seadas.gsfc.nasa.gov/), check [this repository](https://github.com/rakeshkstp/AtmosphericCorrection).*

## Where to place the scripts?
* Download and place contents of this repository in `$OCSSWROOT/ocssw_src/src/l2gen`.
* Add the following files to `$OCSSWROOT/ocssw_src/src/l2gen/CMakeLists.txt` in the `L2GEN_PRODUCT_FILES` section.
  * `get_uqar_utils.c`
  * `calc_uqar_utils.c`
  * `interpol_ed0LUT_5nm_v2.f`
  * `calc_par_uqar.f`
  * `get_ed0_LUT_v2.f`
* Add the following definitions to `l2prod.h`. You can choose your own product ID and name, if it is available. 
  ```
  #define CAT_uqar_par0p              348
  #define CAT_uqar_par0m              349
  #define CAT_uqar_parb               350
  #define CAT_uqar_icw                351
  #define CAT_uqar_COT                352
  #define CAT_uqar_salb               353
  #define CAT_uqar_O3                 354
  #define CAT_uqar_kdpar              360
  ```
* Add the definitions of these products in `$OCDATAROOT/common/product.xml`.
* Add the following line to `l12_proto.h`.
  ```
  void get_uqar_utils(l2str *l2rec, int prodnum, float prod[]);
  ```
* Add the following lines to `prodgen.c`
  ```
  case CAT_uqar_par0p:
  case CAT_uqar_par0m:
  case CAT_uqar_parb:
  case CAT_uqar_kdpar:
  case CAT_uqar_COT:
  case CAT_uqar_O3:
  case CAT_uqar_salb:
  case CAT_uqar_icw:
  get_uqar_utils(l2rec, p->cat_ix, fbuf);
  pbuf = (VOIDP) fbuf;
  break;       
  ```
* Compile the code as described [here](https://seadas.gsfc.nasa.gov/build_ocssw/#building-the-code).
* Now you can run `l2gen` and it will identify the product names described in `$OCDATAROOT/common/product.xml`.

## Author
**Rakesh Kumar Singh**

**[Simon Bélanger](https://github.com/belasi01)**
