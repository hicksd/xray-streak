PRO xstreak_savedata,imstr,filetosave


	OPENW,unit,filetosave,/GET_LUN

	 PRINTF,unit,"Saved on: "+SYSTIME()
	 PRINTF,unit,"Image name: "+imstr.file
   PRINTF,unit,"Target Name: "+imstr.targetdata.name
   
   PRINTF,unit,"Rotation (rads): "+STRCOMPRESS(STRING(imstr.lastrot),/REMOVE_ALL)
	 PRINTF,unit,"Shear: "+STRCOMPRESS(STRING(imstr.lastshear),/REMOVE_ALL)

   WIDGET_CONTROL,imstr.mag_id,GET_VALUE=micpx  ;mics/pixel
   WIDGET_CONTROL,imstr.mag0_id,GET_VALUE=d0    ;pixels
   WIDGET_CONTROL,imstr.fidpx_id,GET_VALUE=fidpx ;pixels
   WIDGET_CONTROL,imstr.fidtt_id,GET_VALUE=fidtt ;ps
   WIDGET_CONTROL,imstr.dx_id,GET_VALUE=dx      ;pixels
   WIDGET_CONTROL,imstr.dy_id,GET_VALUE=dy      ;pixels
   WIDGET_CONTROL,imstr.bkg_id,GET_VALUE=bkg    ;bkg level

   WIDGET_CONTROL,imstr.bpoly_id,GET_VALUE=Npoly;Polynomial order for background
   WIDGET_CONTROL,imstr.skew_id,GET_VALUE=skew ;skew for density
   WIDGET_CONTROL,imstr.mix_id,GET_VALUE=Mix ;(mics) mix thickness
   
;  WIDGET_CONTROL,imstr.rho0_id,GET_VALUE=rho0  ;g/cc
;  WIDGET_CONTROL,imstr.mu0_id,GET_VALUE=mu     ;mass attenuation coefficient, cm2/g
;  WIDGET_CONTROL,imstr.eta_id,GET_VALUE=eta    ;rho/rho0 of shell
;  WIDGET_CONTROL,imstr.R0_id,GET_VALUE=R0      ;(mics) outer shell edge position
   WIDGET_CONTROL,imstr.slit_id,GET_VALUE=slit  ;mic
   WIDGET_CONTROL,imstr.sloty0_id,GET_VALUE=sl_y0    ;(mics) hohl slot position (y0=0=centered on capsule)
   WIDGET_CONTROL,imstr.slotdy_id,GET_VALUE=sl_dy    ;(mics) hohl slot width

   WIDGET_CONTROL,imstr.R0_id,GET_VALUE=R0      ;(mics) outer shell edge position
   WIDGET_CONTROL,imstr.rho_id,GET_VALUE=rhom   ;rho_max of shell
   WIDGET_CONTROL,imstr.del_id,GET_VALUE=del    ;(mics) shell stdev thickness

   WIDGET_CONTROL,imstr.freetype_id,GET_VALUE=free ;Which parameters to vary [R0,rho,del] - 1=vary,0=fix

	 PRINTF,unit,"Polynomial coefficients: (a1 -> a"+STRCOMPRESS(STRING(FIX(N_ELEMENTS(imstr.tpoly)-1)),/REMOVE_ALL)+")"
	 PRINTF,unit,STRING(imstr.tpoly(1:N_ELEMENTS(imstr.tpoly)-1),/PRINT)
	 PRINTF,unit,"Mag (mic/px):"+micpx
	 PRINTF,unit,"x0 (pixel):"+d0
   PRINTF,unit,"Fidu (px):"+fidpx+" Fidu (ps):"+fidtt
	 PRINTF,unit,"dx(pix):"+dx+" dy(pix):"+dy
	 PRINTF,unit,"Background level:  "+bkg
;   PRINTF,unit,"Initial density (g/cc):"+rhom
;   PRINTF,unit,"Mass attenuation (cm2/g):"+mu
   PRINTF,unit,"Background polynomial order: "+Npoly
   PRINTF,unit,"Density function skew: "+skew
   PRINTF,unit,"Mix length (mic): "+Mix

   PRINTF,unit,"Slit size (mic):"+slit
   PRINTF,unit,"Masking slot position (mic):"+sl_y0
   PRINTF,unit,"Masking slot width (mic):"+sl_dy

   PRINTF,unit,"Analysis Method: ",imstr.dat.method
   PRINTF,unit,"Use average mu fit? (yes=1,no=0):"+STRCOMPRESS(STRING(imstr.avmufit))
   
   PRINTF,unit,"First guess - Position (mic):"+R0
   PRINTF,unit,"First guess - Density:"+rhom
   PRINTF,unit,"First guess - Thickness (mic):"+del
   
   PRINTF,unit,"Free parameter, Position: "+STRCOMPRESS(STRING(free(0)))
   PRINTF,unit,"Free parameter, Density: "+STRCOMPRESS(STRING(free(1)))
   PRINTF,unit,"Free parameter, Thickness: "+STRCOMPRESS(STRING(free(2)))
   
	 PRINTF,unit,"Plot ROI (xmin,xmax,ymin,ymax)"
	 PRINTF,unit,STRCOMPRESS(STRING(REFORM(TRANSPOSE(imstr.xyroi),4)))
	 PRINTF,unit,"Inversion ROI (xmin,xmax,ymin,ymax)"
	 PRINTF,unit,STRCOMPRESS(STRING(REFORM(TRANSPOSE(imstr.phroi),4)))
	 PRINTF,unit,"Intensity ROI (xmin,xmax,ymin,ymax)"
	 PRINTF,unit,STRCOMPRESS(STRING(REFORM(TRANSPOSE(imstr.inroi),4)))
   PRINTF,unit,"Block ROI (xmin,xmax,ymin,ymax)"
   PRINTF,unit,STRCOMPRESS(STRING(REFORM(TRANSPOSE(imstr.ghroi),4)))

;	 PRINTF,unit,"Warp region (xmin,xmax,ymin,ymax)"
;	 PRINTF,unit,STRCOMPRESS(STRING([imstr.warp.xmin,imstr.warp.xmax,imstr.warp.ymin,imstr.warp.ymax]))
;	 PRINTF,unit,"Warp matrix"
;	 PRINTF,unit,STRCOMPRESS(STRING(imstr.warp.c))

 	 PRINTF,unit,'***'
;***************************************
   dat=IMSTR.dat
   PRINTF,unit,STRING(['Time(ns)','rho(g/cc)','drho','R(mic)','dR','del(mic)','ddel',$
            'vel(km/s)','dvel','rhoR(g/cm2)','drhoR','mass(ug)','dmass',$
            'MassFrac','dMassFrac','Kapp(cm2/g)',$
            'a0','a1','a2','a3'],FORMAT='(20A14)')
	 PRINTF,unit,STRING([TRANSPOSE(dat.t),TRANSPOSE(dat.maxdensity),TRANSPOSE(dat.dmaxdensity),$
	          TRANSPOSE(dat.R0),TRANSPOSE(dat.dR0),TRANSPOSE(dat.del),TRANSPOSE(dat.ddel),$
	          TRANSPOSE(dat.v),TRANSPOSE(dat.dv),$
	          TRANSPOSE(dat.rhoR),TRANSPOSE(dat.drhoR),$
	          TRANSPOSE(dat.mass),TRANSPOSE(dat.dmass),$
	          TRANSPOSE(dat.mass)/MAX(dat.target.cmass),TRANSPOSE(dat.dmass)/MAX(dat.target.cmass),$
	          TRANSPOSE(dat.av_mu),$
	          TRANSPOSE(REFORM(dat.a[0,*])),TRANSPOSE(REFORM(dat.a[1,*])),$
	          TRANSPOSE(REFORM(dat.a[2,*])),TRANSPOSE(dat.dxmass)],$
	          FORMAT='(20G14.5)')
	CLOSE,unit
	FREE_LUN,unit
END