        real*4 FUNCTION calc_par_uqar(rSZA,rO3,rTAUCL,rALB,ED_LUT,Ed)
          use iso_c_binding

!         Define parameters
          parameter (nwl=83, nthetas=19, nozone=10, ntaucl=8, nalb=7)
          real(c_float)    rSZA,   rO3,   rTAUCL, rALB
          real(c_float) deltax

!         Declare input Look-up-tables
          real(c_float) ED_LUT(nwl, nthetas, nozone, ntaucl, nalb)

!         Outputs:
          real(c_float) ED_inst(nwl)
          real(c_float), intent(out) :: Ed(nwl)

!         Interpolate the LUTs for the given set of SZA, O3 and TAUCL
          CALL interpol_ed0LUT_5nm_v2(ED_LUT,rSZA, rO3, rTAUCL,rALB, ED_inst)

!         Integrate over PAR and UV domains
          calc_par_uqar=0
          deltax = 5.0; 
          do i=23,nwl-1
              calc_par_uqar = calc_par_uqar + deltax*(ED_inst(i)+ED_inst(i+1))/2.0
          enddo
          do i=1,nwl
              Ed(i) = ED_inst(i)
          enddo
          return
        END
