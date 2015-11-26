;This procedure is for use with xstreak.

FUNCTION XmuFromMass,tg,r,rho,av_mu,cmass
;Based on the provided rho(r) and mu(r) for the pre-shot configuration (from tg)
;it returns the mu(r) for the given configuration, assuming that mu vs cumulative mass
;is the same for both.
;r must be uniformly spaced and in ascending order.
;mu is in cm2/g

N=N_ELEMENTS(r)
;mass for r<0 is negative, for r>0 is positive
dr=ABS(r(1)-r(0)) 
mass=4.*!pi*rho*r^2*dr*1.E-6;*(r GE 0) - 4.*!pi*rho*r^2*dr*1.E-6*(r LT 0);ug
;cmass=TOTAL(mass,/CUMULATIVE)-0.5*(mass(0)+mass(N-1)) ;cumulative mass,ug

cmass=DBLARR(N)
w=WHERE(r LE 0) & IF w(0) NE -1 THEN cmass(w)=REVERSE(TOTAL(REVERSE(mass(w)),/CUMULATIVE))
w=WHERE(r GE 0) & IF w(0) NE -1 THEN cmass(w)=TOTAL(mass(w),/CUMULATIVE) ;cumulative mass,ug

;normalize cmass to be 0 at r=0 and positive for r<0 and r>0
;IF MIN(r) GT 0. OR MAX(r) LT 0. THEN BEGIN
;  Rmin=MIN(ABS(r)) & cmass0=cmass(WHERE(ABS(r) EQ Rmin))
;ENDIF ELSE cmass0=INTERPOL(cmass,r,0.) ;cumulative mass at r=0
;cmass=ABS(cmass-cmass0(0))

;make maximum inferred cmass be target max mass
w=WHERE(cmass GT MAX(tg.cmass)) & IF w(0) NE -1 THEN cmass(w) = MAX(tg.cmass) 

mu=INTERPOL(tg.mu,tg.cmass,cmass) ;mu of mass element at r
av_mu=INTERPOL(tg.av_mu,tg.cmass,MAX(cmass));average mu of total mass

RETURN,mu
END

FUNCTION GetMu,dat,sym,En
;Takes the structure 'dat', searches for the element 'sym' and finds mu at photon energy 'En'

exist=STRCMP(REPLICATE(sym,N_ELEMENTS(dat.symb)),dat.symb,/FOLD_CASE) ;finds position of that
w=MIN(WHERE(exist GT 0))                                              ;element in DAT
;w=WHERE(dat.symb EQ sym)
IF w(0) EQ -1 THEN RETURN,REPLICATE(0,N_ELEMENTS(En))
RETURN,INTERPOL(dat(w).mu,dat(w).E,En)
END

FUNCTION LoadTarget,XKAPPAFILE=XKAPPAFILE,TDAT=Tdat,En=En,Mix=Mix  ;returns the initial target details in a structure
;XKAPPAFILE=filename of xkappa.sav file containing Henke data.
;TDAT=target data structure output by procedure "target" into "xstreak_event".
;En=Photon energy in keV (overrides default); 
;Mix=scale length over which mix occurs, in microns (default=0)
;Scheme guesses rho(r) then fits the resulting mu*rho profile to the radiograph. 
;mu is calculated at each time step based on the distribution of elements in the 
;initial target and assuming 1-D ablation

IF KEYWORD_SET(XKAPPAFILE) THEN BEGIN
  RESTORE,xkappafile 
ENDIF ELSE BEGIN
  RESTORE,'c:\hicks\xstreak\xKappa.sav' ;loads the structure 'xray' containing Kappa info
ENDELSE

IF KEYWORD_SET(TDAT) THEN BEGIN
  IF NOT(KEYWORD_SET(En)) THEN En=Tdat.En ;keV. x-ray photon energy
  ;************Layer thicknesses and densities
  name=Tdat.name
  Rin=Tdat.Rin ;inner radius of capsule
  N_Region=Tdat.N_Region ;Number of regions
  d_Region=Tdat.d_region[0:N_region-1];thickness of each region 
  rho_Region=Tdat.rho_region[0:N_region-1];density of each region
  ;************Element details for each region
  ;f can be first listed in STOICHIOMETRIC form for convenience; automatically normalized afterwards
  ;  mat = tdat.mat(0:tdat.N_mat-1)
  mat =REPLICATE({sym:'' ,A:0.,mu:0.,f:FLTARR(N_Region)},tdat.N_mat)
  STRUCT_ASSIGN,tdat.mat,mat 
; print,name
; print,Rin
; print,N_Region
; print,d_Region
; print,rho_Region
; print,mat
  GOTO,jump
ENDIF 



;;OMEGA Be-Cu capsule; Shot 49780 or 49783
;IF NOT(KEYWORD_SET(En)) THEN En=5.2 ;keV. x-ray photon energy
;;************Layer thicknesses and densities
;name='Omega49780/49783'
;Rin=209.0 ;inner radius of capsule
;;d_Region=  [3.9,26.4,5.7] ;thickness of each region Shot49780
;d_Region= [4.5,26.9,21.] ;thickness of each region Shot49783
;rho_Region=[1.75,2.07,1.75];density of each region
;N_Region=N_ELEMENTS(d_Region) ;Number of regions
;;************Element details for each region
;;f can be first listed in STOICHIOMETRC form for convenience; automatically normalized afterwards
;mat =[$
;     {sym:'Be',A: 9.012,mu:0.,f:[1.0,0.97,1.0]},$
;     {sym:'Cu',A:63.546,mu:0.,f:[0.0,0.03,0.0]}$
;;     {sym:'Ge',A:  72.64,mu:0.,f:[0.0,0.00,0.00,0.0,0.00,0.0]}$
;     ]

;Scale 0.9; 155um wall
;IF NOT(KEYWORD_SET(En)) THEN En=8.35 ;keV. x-ray photon energy
;;************Layer thicknesses and densities
;Rin=789.2 ;inner radius of capsule
;d_Region=  [21.9,28.8,13.0,91.9];thickness of each region 
;rho_Region=[1.034,1.08,1.06,1.034];density of each region
;N_Region=N_ELEMENTS(d_Region) ;Number of regions
;;************Element details for each region
;;f can be first listed in STOICHIOMETRC form for convenience; automatically normalized afterwards
;mat =[$
;     {sym:'C' ,A:12.0107,mu:0.,f:REPLICATE(0.4235,N_Region)},$
;     {sym:'H' ,A: 1.0079,mu:0.,f:REPLICATE(0.5715,N_Region)},$
;     {sym:'O' ,A:15.9994,mu:0.,f:REPLICATE(0.0050,N_Region)},$
;     {sym:'Ge',A:  72.64,mu:0.,f:[0.,0.0057,0.0034,0.]}$
;     ]
     
;;Scale 0.9; 135um wall SHOT N091123
IF NOT(KEYWORD_SET(En)) THEN En=8.35 ;keV. x-ray photon energy
;************Layer thicknesses and densities
name='Scale 0.9; 135um wall SHOT N091123'
Rin=784.15 ;inner radius of capsule
d_Region=  [23.1,29.4,12.8,68.3];thickness of each region 
rho_Region=[1.02,1.07,1.05,1.02];density of each region
N_Region=N_ELEMENTS(d_Region) ;Number of regions
;************Element details for each region
;f can be first listed in STOICHIOMETRC form for convenience; automatically normalized afterwards
mat =[$
     {sym:'C' ,A:12.0107,mu:0.,f:REPLICATE(0.4235,N_Region)},$
     {sym:'H' ,A: 1.0079,mu:0.,f:REPLICATE(0.5715,N_Region)},$
     {sym:'O' ,A:15.9994,mu:0.,f:REPLICATE(0.0050,N_Region)},$
     {sym:'Ge',A:  72.64,mu:0.,f:[0.,0.0062,0.0039,0.]}$
;     {sym:'Ge',A:  72.64,mu:0.,f:[0.0,0.00,0.00,0.0,0.00,0.0]}$
     ]

;;Rev5 ConvAbl, 3/29/10; changed to 8.95 keV on 4/13/10
;IF NOT(KEYWORD_SET(En)) THEN En=8.95 ;keV. x-ray photon energy
;;************Layer thicknesses and densities
;Rin=902.7 ;inner radius of capsule
;d_Region=  [20.30,5.0000,34.000,13.000,133.0];thickness of each region 
;rho_Region=[1.069,1.1077,1.1464,1.1077,1.069];density of each region
;N_Region=N_ELEMENTS(d_Region) ;Number of regions
;;************Element details for each region
;;f can be first listed in STOICHIOMETRC form for convenience; automatically normalized afterwards
;mat =[$
;     {sym:'C' ,A:12.0107,mu:0.,f:REPLICATE(0.42337,N_Region)},$
;     {sym:'H' ,A: 1.0079,mu:0.,f:REPLICATE(0.57155,N_Region)},$
;     {sym:'O' ,A:15.9994,mu:0.,f:REPLICATE(0.00508,N_Region)},$
;     {sym:'Ge',A:  72.64,mu:0.,f:[0.0,0.00503,0.0101,0.00503,0.0]}$
;;     {sym:'Ge',A:  72.64,mu:0.,f:[0.0,0.00,0.00,0.0,0.00,0.0]}$
;     ]

;Rev5 THD ConvAbl, 4/22/10; 
;IF NOT(KEYWORD_SET(En)) THEN En=8.95 ;keV. x-ray photon energy
;;************Layer thicknesses and densities
;Rin=923.8 ;inner radius of capsule
;d_Region=  [5.300,5.3000,33.900,13.100,129.6];thickness of each region 
;rho_Region=[1.069,1.1077,1.1464,1.1077,1.069];density of each region
;N_Region=N_ELEMENTS(d_Region) ;Number of regions
;;************Element details for each region
;;f can be first listed in STOICHIOMETRIC form for convenience; automatically normalized afterwards
;mat =[$
;     {sym:'C' ,A:12.0107,mu:0.,f:REPLICATE(0.42337,N_Region)},$
;     {sym:'H' ,A: 1.0079,mu:0.,f:REPLICATE(0.57155,N_Region)},$
;     {sym:'O' ,A:15.9994,mu:0.,f:REPLICATE(0.00508,N_Region)},$
;     {sym:'Ge',A:  72.64,mu:0.,f:[0.0,0.00502,0.0098,0.00502,0.0]}$
;;     {sym:'Ge',A:  72.64,mu:0.,f:[0.0,0.00,0.00,0.0,0.00,0.0]}$
;     ]

;Hammel Scale 1.07 test
;IF NOT(KEYWORD_SET(En)) THEN  En=10.69;14.43 ;keV. x-ray photon energy
;;************Layer thicknesses and densities
;Rin=902.0;844.0 ;inner radius of capsule IGNORE DTH LAYER
;d_Region=  [  5.0,  5.0, 38.0, 14.0,128.0];58.0;thickness of each region
;rho_Region=[1.082,1.096,1.116,1.096,1.082];0.255;density of each region
;N_Region=N_ELEMENTS(d_Region) ;Number of regions
;;************Element details for each region
;;f can be first listed in STOICHIOMETRIC form for convenience; automatically normalized afterwards
;mat =[$
;     {sym:'C' ,A:12.0107,mu:0.,f:[0.42337,0.42252,0.42125,0.42252,0.42337]},$
;     {sym:'H' ,A: 1.0079,mu:0.,f:[0.57155,0.57041,0.56869,0.57041,0.57155]},$
;     {sym:'O' ,A:15.9994,mu:0.,f:[0.00508,0.00507,0.00507,0.00506,0.00508]},$
;     {sym:'Ge',A:  72.64,mu:0.,f:[0.00000,0.00200,0.00200,0.00500,0.00000]}$
;     ]

;IF NOT(KEYWORD_SET(En)) THEN En=8.35 ;keV. x-ray photon energy
;;************Layer thicknesses and densities
;Rin=900. ;inner radius of capsule
;d_Region=  [ 10.0, 10.0,  10.0,  30.0,  15.0, 70.0];thickness of each region 
;rho_Region=[1.082,1.082,1.0956, 1.116,1.0956,1.082];density of each region
;N_Region=N_ELEMENTS(d_Region) ;Number of regions
;;************Element details for each region
;;f can be first listed in STOICHIOMETRC form for convenience; automatically normalized afterwards
;mat =[$
;     {sym:'C' ,A:12.0107,mu:0.,f:REPLICATE(1.,N_Region)},$
;     {sym:'H' ,A: 1.0079,mu:0.,f:REPLICATE(1.35,N_Region)},$
;     {sym:'O' ,A:15.9994,mu:0.,f:REPLICATE(0.012,N_Region)},$
;     {sym:'Ge',A:  72.64,mu:0.,f:[0.0,0.0,0.0048,0.012,0.0048,0.0]}$
;;     {sym:'Ge',A:  72.64,mu:0.,f:[0.0,0.00,0.00,0.0,0.00,0.0]}$
;     ]

;IF NOT(KEYWORD_SET(En)) THEN En=7.0 ;keV. x-ray photon energy
;;************Layer thicknesses and densities
;Rin=1027.4 ;inner radius of capsule
;d_Region=  [11.0, 6.0, 6.0, 34.1,11.3,100.2];thickness of each region 
;rho_Region=[1.80,1.80,1.85, 1.90,1.85,1.80];density of each region
;N_Region=N_ELEMENTS(d_Region) ;Number of regions
;;************Element details for each region
;;f can be first listed in STOICHIOMETRC form for convenience; automatically normalized afterwards
;mat =[$
;     {sym:'Be',A:9.01218,mu:0.,f:REPLICATE(1.000,N_Region)},$
;     {sym:'Ar',A: 39.948,mu:0.,f:REPLICATE(0.002,N_Region)},$
;     {sym:'O' ,A:15.9994,mu:0.,f:REPLICATE(0.004,N_Region)},$
;     {sym:'Cu',A: 63.546,mu:0.,f:[0.0,0.0,0.00506,0.0101,0.00506,0.0]}$
;     ]

jump:     
;************Load mu values if not provided
FOR i=0,N_ELEMENTS(mat)-1 DO BEGIN
  IF mat(i).mu EQ 0 THEN mat(i).mu=GetMu(XRAY,mat(i).sym,En)
ENDFOR

;*************Normalize stoichiometry to give atomic fraction for each element
FOR i=0,N_Region-1 DO mat(*).f[i]=mat(*).f[i]/TOTAL(mat(*).f[i])
;print,mat.f(0:N_Region-1)*100.

;************Edge and midpoint radii of each region
R_Mesh=Rin+[0,TOTAL(d_Region,/CUMULATIVE)] ;edge points of each layer
R_Region=FLTARR(N_Region) ;midpoint radius of each layer 
FOR i=0,N_Region-1 DO R_Region(i)=(R_Mesh(i)+R_Mesh(i+1))/2.

;*******rho(r) and mu(r) along finely-resolved radius, r
N=FIX(TOTAL(d_Region))*2. ;Number of radial elements ;1 element per half micron
;1st element=inner radius; last element=outer radius
r=FINDGEN(N)/(N-1)*(MAX(R_Mesh)-MIN(R_Mesh))+MIN(R_Mesh) 
dr=r(1)-r(0)
;Create rho(r) and mu(r)
rho=FLTARR(N) & mu=FLTARR(N)
FOR i=0,N_Region-1 DO BEGIN 
  w=WHERE(r GE R_Mesh(i) AND r LE R_Mesh(i+1))
  rho(w)=rho_Region(i)
  newmu=TOTAL(mat(*).mu*mat(*).A*mat(*).f[i])/TOTAL(mat(*).A*mat(*).f[i])
  mu(w)=newmu
  print,"Region "+STRING(i,FORMAT='(I2)')+" kappa= "+STRING(newmu,FORMAT='(F8.5)')
ENDFOR
;***********Add Mix, if any
IF KEYWORD_SET(mix) THEN BEGIN
;  npix=2.*mix/dr ;double "mix" to form boxcar 
;  rho=SMOOTH(rho,npix) & mu=SMOOTH(mu,npix,/EDGE_TRUNCATE)
  fconvol=FLTARR(N)
  sig_hi=mix/2. & w=WHERE(r GE MEAN(r)) & fconvol(w)=EXP(-(r(w)-MEAN(r))^2/(2*sig_hi^2)) ;
  sig_lo=mix/2.  & w=WHERE(r LE MEAN(r)) & fconvol(w)=EXP(-(r(w)-MEAN(r))^2/(2*sig_lo^2))
  rho=CONVOL(rho,fconvol,TOTAL(fconvol),EDGE_TRUNCATE=1,EDGE_WRAP=0) ;conserving blurring function
  mu =CONVOL(mu ,fconvol,TOTAL(fconvol),EDGE_TRUNCATE=1,EDGE_WRAP=0)
ENDIF ELSE BEGIN
  rho=CONVOL(rho,[1,1,1],3,/EDGE_TRUNCATE) ;smooth by 1 pixel even when no mix is specified
  mu =CONVOL(mu ,[1,1,1],3,/EDGE_TRUNCATE) ;this tries to prevent discontinuity problems
ENDELSE

;***********Cumulative mass
mass=4.*!pi*rho*r^2*dr*1.E-6 ;ug
cmass=TOTAL(mass,/CUMULATIVE)-0.5*(mass(0)+mass(N-1)) ;cumulative mass,ug
w=WHERE(cmass LE 0) & IF w(0) NE -1 THEN cmass(w)=0.1 ;avoid having cmass LE 0; set minimum to 0.1 ug
;***********Average mu
mumass=TOTAL(mu*mass,/CUMULATIVE)-0.5*(mu(0)*mass(0)+mu(N-1)*mass(N-1)) ;cumulative mu*mass,ug*cm2/g
w=WHERE(mumass LE 0) & IF w(0) NE -1 THEN mumass(w)=0.0 ;avoid having mumass LE 0; set minimum to 0.0

av_mu=mumass/cmass ;average mu of remaining mass, cm2/g
av_mu(0)=av_mu(1) ;avoid anomalies with first element by setting it equal to the second
;***********Deliver structure containing mu vs cmass
tg={name:name,N:N,r:r,rho:rho,mu:mu,cmass:cmass,mumass:mumass,av_mu:av_mu,mat:mat} ;target parameters

RETURN,tg
END



FUNCTION read_xdat,ncols,file
;Reads x-ray attentuation length data from LBL site
header=''
line=FLTARR(ncols)
OPENR,unit,file,/GET_LUN
WHILE NOT(EOF(unit)) DO BEGIN
    ON_IOERROR,skip
    READF,unit,header & READS,header,line
    IF SIZE(data,/TYPE) NE 0 THEN data=[[data],[line]] ELSE data=line
    skip:
ENDWHILE
CLOSE,unit & FREE_LUN,unit
;E in keV, mu in cm2/g
dat={E:REFORM(data(0,*))/1.E3,mu:1.E4/REFORM(data(1,*))}
RETURN,dat
END

PRO CreateXKappa
;This creates a *.sav file containing cold x-ray opacities
;which can be quickly read by xstreak

dir='c:\hicks\eos\x-ray\'
el=['H','Be','C','O','Al','Si','Ar','Mn','Fe','Ni','Cu','Ge'];,'C22H10N2O5']
;el=['H','Be','C']

dat=read_xdat(2,dir+el(0)+'-length.dat') 
xray=REPLICATE({symb:el(0),E:dat.E,mu:dat.mu},N_ELEMENTS(el))
FOR i=1,N_ELEMENTS(el)-1 DO BEGIN
  dat=read_xdat(2,dir+el(i)+'-length.dat') &  xray(i)={symb:el(i),E:dat.E,mu:dat.mu}
ENDFOR
SAVE,FILENAME='c:\hicks\xstreak\'+'xKappa.sav',xray
;E=FINDGEN(100)/99*10.
;PLOT,E,INTERPOL(xray(2).dat.mu,xray(2).dat.E,E),/YLOG
 
RETURN
END

PRO comparemix,eps=eps
;Examine how the average opacity changes with mix
;Compares Planck spectrum vs absorption of e.g. Ge
RED=255 & GREEN=255*256L & BLUE=255*256L*256L & YELLOW=RED+GREEN & PINK=RED+BLUE & LBLUE=GREEN+BLUE
ORANGE=255+126*256L & PURPLE=126+BLUE & PLAIN=RED+GREEN+BLUE & BACKG=0 & BROWN=174+71*256L
charsz=1.2 & thk=1.5 & symsz=1.2

!P.FONT=-1
IF KEYWORD_SET(eps) OR KEYWORD_SET(ps) THEN BEGIN
!P.FONT=1
    SET_PLOT,'ps'
    !X.MARGIN=[6,2] & !Y.MARGIN=[3.2,2]    
;    DEVICE,ENCAPSULATED=0,LANDSCAPE=1,XSIZE=9,YSIZE=6,/INCHES,/COLOR,/HELVETICA,FILENAME="C:\Documents and Settings\hicks13\DESKTOP\IDL.EPS"
;    DEVICE,ENCAPSULATED=1,LANDSCAPE=0,XSIZE=15,YSIZE=12,/CM,/COLOR,/HELVETICA,FILENAME="C:\Documents and Settings\hicks13\DESKTOP\IDL.EPS"
    DEVICE,/ENCAPSULATED,LANDSCAPE=0,XSIZE=15,YSIZE=7,/COLOR,/HELVETICA,FILENAME="C:\Documents and Settings\hicks13\DESKTOP\IDL.EPS"
;    DEVICE,/ENCAPSULATED,LANDSCAPE=0,XSIZE=8,YSIZE=14,/COLOR,/HELVETICA,FILENAME="C:\Documents and Settings\hicks13\DESKTOP\IDL.EPS"
    TVLCT, [0,255,0,0,255,255,0,255,126,255,174], [0,0,255,0,255,0,255,126,0,255,71], [0,0,0,255,0,255,255,0,255,255,0]
    PLAIN=0 & RED=1 & GREEN=2 & BLUE=3 & YELLOW=4 & PINK=5 & LBLUE=6 & ORANGE=7 & PURPLE=8 & BACKG=9 & BROWN=10
    FILL=1 & charsz=1.1 & thk=3.0*thk & symsz=0.7*symsz
    !P.MULTI=[0,2,2]
ENDIF ELSE !P.MULTI=[0,2,2]
clr=[PLAIN,RED,GREEN,BLUE,BROWN,PINK,LBLUE,ORANGE,PURPLE,YELLOW]
b = FINDGEN(17)*(!PI*2/16.) & USERSYM,1*COS(b),1*SIN(b),THICK=thk,FILL=0
mic=STRING(["265B])+'m'

;tg={r:r,rho:rho,mu:mu,cmass:cmass,mumass:mumass,av_mu:av_mu,mat:mat} ;target parameters
En=8.35 & En1=11.2 ;En=11.2 for above Ge-edge
tg=loadtarget(En=En) &  tgm=loadtarget(En=En,Mix=20.)
tg1=loadtarget(En=En1) &  tgm1=loadtarget(En=En1,Mix=20.)


!P.MULTI=[0,2,1]

PLOT,tg.r-MIN(tg.r),tg.mu,THICK=thk,CHARSIZE=charsz,CHARTHICK=thk,XTHICK=thk,YTHICK=thk,$
    XTITLE='Initial distance from inner radius, '+mic,YTITLE='Opacity, cm!u2!n/g';$
;    ,TITLE=STRING(En,FORMAT='(F5.2)')+' keV'
OPLOT,tgm.r-MIN(tgm.r),tgm.mu,THICK=thk,COLOR=RED,LINESTYLE=2
;
;PLOT,tg.cmass/MAX(tg.cmass),tg.av_mu,THICK=thk,CHARSIZE=charsz,CHARTHICK=thk,XTHICK=thk,YTHICK=thk,XRANGE=[0,1],$
;    XTITLE='Remaining mass fraction',YTITLE='Mass-averaged opacity (cm!u2!n/g)',$
;    TITLE=STRING(En,FORMAT='(F5.2)')+' keV'
;OPLOT,tgm.cmass/MAX(tgm.cmass),tgm.av_mu,THICK=thk,COLOR=RED

PLOT,tg.cmass/MAX(tg.cmass),tg.av_mu/tg.av_mu,THICK=thk,CHARSIZE=charsz,$
    CHARTHICK=thk,XTHICK=thk,YTHICK=thk,XRANGE=[0,1],YRANGE=[0.8,1.4],/XSTYLE,/YSTYLE,$
    XTITLE='Remaining mass fraction',YTITLE='Inferred mass/Real Mass'
;OPLOT,tg1.cmass/MAX(tg1.cmass),tgm1.av_mu/tg1.av_mu,THICK=thk,COLOR=RED
OPLOT,tg.cmass/MAX(tg.cmass),tgm.av_mu/tg.av_mu,THICK=thk,COLOR=RED,LINESTYLE=2


;PLOT,tg.cmass/MAX(tg.cmass),tgm.av_mu/tg.av_mu,THICK=thk,CHARSIZE=charsz,$
;    CHARTHICK=thk,XTHICK=thk,YTHICK=thk,XRANGE=[0,1],YRANGE=[0.7,2.7],/XSTYLE,/YSTYLE,$
;    XTITLE='Remaining mass fraction',YTITLE='Inferred mass/Real Mass'
;OPLOT,tg1.cmass/MAX(tg1.cmass),tgm1.av_mu/tg1.av_mu,THICK=thk,COLOR=RED
;OPLOT,tg.cmass/MAX(tg.cmass),tg.av_mu/tg.av_mu,THICK=thk,LINESTYLE=2
;
;PLOT,tgm.av_mu/tg.av_mu*tg.cmass/MAX(tg.cmass),$ 
;  tgm1.av_mu/tg1.av_mu*tg1.cmass/MAX(tg1.cmass),THICK=thk,CHARSIZE=charsz,$
;    CHARTHICK=thk,XTHICK=thk,YTHICK=thk,/XSTYLE,/YSTYLE,XRANGE=[0,1],YRANGE=[0,1],$
;    XTITLE='Inferred mass fraction (8.35 keV)',$
;    YTITLE='Inferred mass fraction (12.3 keV)'
;OPLOT,tg.cmass/MAX(tg.cmass),tg1.cmass/MAX(tg1.cmass),THICK=thk,COLOR=PLAIN,LINESTYLE=2

!P.MULTI=0
IF KEYWORD_SET(eps) OR KEYWORD_SET(ps)  THEN BEGIN
    DEVICE,/CLOSE
    SET_PLOT,'WIN'
    !P.FONT=-1
ENDIF
RETURN
END

PRO examineCu,eps=eps
;examines how much trace Cu confuses measurements using a Zn He-a backlighter.
RED=255 & GREEN=255*256L & BLUE=255*256L*256L & YELLOW=RED+GREEN & PINK=RED+BLUE & LBLUE=GREEN+BLUE
ORANGE=255+126*256L & PURPLE=126+BLUE & PLAIN=RED+GREEN+BLUE & BACKG=0 & BROWN=174+71*256L
charsz=1.2 & thk=1.5 & symsz=1.2

!P.FONT=-1
IF KEYWORD_SET(eps) OR KEYWORD_SET(ps) THEN BEGIN
!P.FONT=1
    SET_PLOT,'ps'
    !X.MARGIN=[6,2] & !Y.MARGIN=[3.2,2]    
;    DEVICE,ENCAPSULATED=0,LANDSCAPE=1,XSIZE=9,YSIZE=6,/INCHES,/COLOR,/HELVETICA,FILENAME="C:\Documents and Settings\hicks13\DESKTOP\IDL.EPS"
;    DEVICE,ENCAPSULATED=1,LANDSCAPE=0,XSIZE=15,YSIZE=12,/CM,/COLOR,/HELVETICA,FILENAME="C:\Documents and Settings\hicks13\DESKTOP\IDL.EPS"
    DEVICE,/ENCAPSULATED,LANDSCAPE=0,XSIZE=15,YSIZE=7,/COLOR,/HELVETICA,FILENAME="C:\Documents and Settings\hicks13\DESKTOP\IDL.EPS"
;    DEVICE,/ENCAPSULATED,LANDSCAPE=0,XSIZE=8,YSIZE=14,/COLOR,/HELVETICA,FILENAME="C:\Documents and Settings\hicks13\DESKTOP\IDL.EPS"
    TVLCT, [0,255,0,0,255,255,0,255,126,255,174], [0,0,255,0,255,0,255,126,0,255,71], [0,0,0,255,0,255,255,0,255,255,0]
    PLAIN=0 & RED=1 & GREEN=2 & BLUE=3 & YELLOW=4 & PINK=5 & LBLUE=6 & ORANGE=7 & PURPLE=8 & BACKG=9 & BROWN=10
    FILL=1 & charsz=1.1 & thk=3.0*thk & symsz=0.7*symsz
    !P.MULTI=[0,2,1]
ENDIF ELSE !P.MULTI=[0,2,1]
clr=[PLAIN,RED,GREEN,BLUE,BROWN,PINK,LBLUE,ORANGE,PURPLE,YELLOW]
b = FINDGEN(17)*(!PI*2/16.) & USERSYM,1*COS(b),1*SIN(b),THICK=thk,FILL=0

;Tri-doped CH: Si + Ge + Cu, Aug 15 2011
;En=8.980 ;8.976, 8.980 ;keV. x-ray photon energy
;************Layer thicknesses and densities
name='CH(SiGeCu)'
Rin=931.0 ;inner radius of capsule
d_Region=  [6.00000,6.00000,34.0000,10.0000,139.000];thickness of each region 
rho_Region=[1.07100,1.08800,1.09900,1.07600,1.06500];density of each region
N_Region=N_ELEMENTS(d_Region) ;Number of regions
N_Mat=6 ;Number of different materials (=elements)
;************Element details for each region
;f can be first listed in STOICHIOMETRIC form for convenience; automatically normalized afterwards
mat =[$
     {sym:'C' ,A:12.0107,mu:0.,f:REPLICATE(0.42337,N_Region)},$
     {sym:'H' ,A: 1.0079,mu:0.,f:REPLICATE(0.57155,N_Region)},$
     {sym:'O' ,A:15.9994,mu:0.,f:REPLICATE(0.00508,N_Region)},$
     {sym:'Si',A:28.0855,mu:0.,f:[0.000000,0.007000,0.017000,0.010000,0.000000]},$
     {sym:'Ge',A:72.6400,mu:0.,f:[0.000000,0.001500,0.001500,0.000000,0.000000]},$
     {sym:'Cu',A:63.5460,mu:0.,f:[0.001000,0.000000,0.000000,0.000000,0.000000]}$     
     ]
Tdat={name:name,En:0.,Rin:Rin,N_Region:N_Region,N_Mat:N_Mat,d_Region:d_region,rho_Region:rho_Region,mat:mat}

En0=8.976 & En1=8.980 ;below and above Cu K-edge
q0=loadtarget(En=En0,Tdat=Tdat)
q1=loadtarget(En=En1,Tdat=Tdat)
HELP,/struct,q0

PLOT,q0.R,q0.mu,XRANGE=[920,1000],$
  CHARSIZE=charsz,XTHICK=thk,YTHICK=thk,THICK=thk,CHARTHICK=thk,$ 
  XTITLE="Radius (!Mmm)",YTITLE="Opacity (cm!u2!n/g)"
OPLOT,q1.R,q1.mu,THICK=thk,LINESTYLE=2,COLOR=RED

fmass=q0.cmass/MAX(q0.cmass)
PLOT,[0,1],[1,1],XTHICK=thk,YTHICK=thk,THICK=thk,CHARTHICK=thk,$
  /XSTYLE,/YSTYLE,XRANGE=[0,0.3],YRANGE=[0.9,2.0],$
  XTITLE="Remaining mass fraction",YTITLE="Inferred Mass/Real Mass"
OPLOT,fmass,q1.mumass/q0.mumass,THICK=thk,LINESTYLE=2,COLOR=RED

!P.MULTI=0
IF KEYWORD_SET(eps) OR KEYWORD_SET(ps)  THEN BEGIN
    DEVICE,/CLOSE
    SET_PLOT,'WIN'
    !P.FONT=-1
ENDIF
RETURN
END