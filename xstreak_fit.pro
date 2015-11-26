FUNCTION shelldy,x,P,y0,dy,N
;takes a calculated projected depth P(x), which is assumed to be
;taken along a slice through the center of the sphere, and calculates the new projected depth
;assuming the slice has a finite width (dy) and is actually centered at a position (y0)
;which may not be the center of the sphere. 'N' is the number of slices, each taken at a 
;different position within the width 'dy', over which P(x,y) is averaged.

IF N LE 1 THEN RETURN,P

y=DINDGEN(N)/(N-1)*(dy-dy/N)-dy/2.+dy/N/2.+y0 ;range of offset positions, centered on y0

PP=DINDGEN(N_ELEMENTS(x),N) ;array of new P
FOR i=0,N-1 DO BEGIN
  xx=SQRT(x^2+y(i)^2)
  PP(*,i)=INTERPOL(P,x,xx)
ENDFOR

Pnew=REFORM(REBIN(PP,N_ELEMENTS(x),1)) ;Average over all the various P slices

RETURN,Pnew
END

FUNCTION shell,N,a,xm,rho,mu,M,func,$
  INIT=INIT,TARGET=TARGET,MAXRHO=MAXRHO,MEANR=MEANR,TWOSIGMA=TWOSIGMA,MASS=MASS,RHOR=RHOR,FITMURHO=FITMURHO,SKEW=SKEW
;Calculates the optical depth of a shell with parameters given by a.
;Returns the optical depth, as well as the density profile, rho 
;a=[R0,eta,del,x0,slit,micpx,sl_y0,sl_dy]

R0=ABS(a(0))  & rhomax=ABS(a(1)) &  del=ABS(a(2))  
x0=a(3)  & pin=a(4) & micpx=a(5)    & sl_y0=a(6) & sl_dy=a(7)
;rhomaxT=ABS(a[8]) ;max of tail density

del=ABS(R0*(1-EXP(-del/R0)))+1.E-4;make sure that del is less than R0 and not zero
;rhomaxb=(rhomax*1.)*(1-EXP(-ABS(a(8))/(rhomax*1.))) ;ablation wave max density (make sure it's less than rhomax)

;IF SIZE(func,/type) EQ 0 THEN func='Gaussian !4q!3'
;CONSTANT DENSITY
;rho=maxrho*DOUBLE(ABS(ABS(xm)-ABS(R0)) LE ABS(del)/2.) & func='Constant !4q!3'   ;Constant density shell
;GAUSSIAN DENSITY
;rho=maxrho*EXP(-(ABS(xm)-R0)^2/2/(del/2)^2)   & func='Gaussian !4q!3'        ;Gaussian density shell
;TRIANGULAR DENSITY
;rho=maxrho*(ABS(xm)-ABS(R0))/ABS(del)*DOUBLE((ABS(xm)-ABS(R0)) GE 0. AND (ABS(xm)-ABS(R0)) LE ABS(del)) & func='Triangle'
;
;rho=maxrho*EXP(-((xm^2-R0^2)/(R0*del))^2)   & func='d!4q!3/dr(0)=0'
;rho=maxrho*EXP(-((xm^2-R0^2)/(2.*del)^2)^2) & func='2d!4q!3/dr(0)=0'
;SKEWED GAUSSIAN
IF NOT(KEYWORD_SET(SKEW)) THEN skew=0. ;0 is Gaussian; -7 is skewed heavily towards larger radius.
alpha=skew & s=(ABS(xm)-R0)/(del/2.) & rho=rhomax*2.*SQRT(2.*!PI)*(EXP(-s^2/2.)/SQRT(2.*!PI))*(0.5*(1+IMSL_ERF(alpha*s/SQRT(2)))) & func='Skewed Gaussian' 
;;;;alpha=skew & s=(ABS(xm)-R0)/(del/2.) & rho=rhomax*2.*SQRT(2.*!PI)*(EXP(-s^2/2.)/SQRT(2.*!PI))*(0.5*(1+ERF(alpha*s/SQRT(2)))) & func='Skewed Gaussian' 

;skewT=0. & R0T=R0+del & delT=del*2. ;& rhomaxT=1.;rhomax
;alpha=skewT & s=(ABS(xm)-R0T)/(delT/2.) & tail=rhomaxT*2.*SQRT(2.*!PI)*(EXP(-s^2/2.)/SQRT(2.*!PI))*(0.5*(1+IMSL_ERF(alpha*s/SQRT(2)))) & func='Skewed Gaussian' 
;rho=SQRT(rho^2+tail^2)

;*******Physically motivated density distribution inside ablation front
;g=5./3. ;adiabatic exponent
;nu=-0. ;acceleration index (du/dt=r^nu)
;w=WHERE(xm^2 GE (R0-del)^2 AND xm^2 LE R0^2)
;rho=DBLARR(N)
;;IF w(0) NE -1 THEN rho(w)=rhomax*((ABS(xm(w))^(nu+1)-(R0-del)^(nu+1))/(R0^(nu+1)-(R0-del)^(nu+1)))^(1./(g-1))
;IF w(0) NE -1 THEN rho(w)=g/(g-1)*rhomax/del*(1+(ABS(xm(w))-R0)/del)^(1./(g-1)) ;here rhomax is actually rho*R
;**********

;prevent distribution from spilling outside ROI by using Fermi function 
;to smoothly drop edges to zero.
rho=rho/(1+EXP((ABS(xm)-0.9*MAX(ABS(xm)))/(0.015*MAX(ABS(xm))))) 
;Stop distribution from spilling into center
core=20. & rho=rho/(1+EXP(-(ABS(xm)-core)/(0.1*core)))



;THIS IS WHERE MU IS CALCULATED; mu=1 if fitting to mu*rho
IF KEYWORD_SET(fitmurho) THEN mu=REPLICATE(1.,N) ELSE BEGIN
  IF KEYWORD_SET(target) THEN mu=XmuFromMass(target,xm,rho) ELSE mu=REPLICATE(1.,N)  
ENDELSE

;IF KEYWORD_SET(TARGET) AND NOT(KEYWORD_SET(fitmurho)) THEN BEGIN
;  mu=XmuFromMass(target,xm,rho,av_mu)  
  ;mu=REPLICATE(1.,N);*av_mu ;use average mu over entire profile to stabilize fitting routine.Not as accurate but faster.10/17/09
;ENDIF ELSE BEGIN
;  mu=REPLICATE(1.,N) ;defaults to mu=1 when no target info is specified
;  av_mu=1.           ;defaults to mu=1 when no target is specified
;ENDELSE 

;IF func EQ 'Constant !4q!3' THEN BEGIN ;DON'T USE CONSTANT RHO UNTIL THIS IS FIXED!!!!!!!
;    OD=DBLARR(N) & del=ABS(del)
;    w=WHERE(ABS(xm) LE R0+del/2.)
;    IF w(0) NE -1 THEN OD(w)=2*maxrho*SQRT((R0+del/2)^2-xm(w)^2) ;THIS NEEDS TO BE FIXED TO CALCULATE MU*RHO!!!!!!!!!
;    w2=WHERE(ABS(xm) LE R0-del/2)
;    IF w2(0) NE -1 THEN OD(w2)=2*maxrho*(SQRT((R0+del/2.)^2-xm(w2)^2)-SQRT((R0-del/2.)^2-xm(w2)^2) )
;ENDIF ELSE OD=abeltransform(xm,mu*rho,M,init=init)/1.E4 ;Need to make unitless: um*cm2/g*g/cm3=um/cm

;Smoothing before abeltransform to improve inversion of profiles with sharp edges (e.g. constant top hat profile)
;rho=SMOOTH(rho,pin/micpx)    ;Average over imaging slit width (x-steps)
;mu=SMOOTH(mu,pin/micpx)


posx=WHERE(xm GE 0) & negx=WHERE(xm LT 0)
IF N_ELEMENTS(posx) GT N_ELEMENTS(negx) THEN w=WHERE(xm GE 0) $
  ELSE w=WHERE(xm LE 0) ;this doesn't work for calculating quantities below when all xm LE 0.

IF KEYWORD_SET(Maxrho) THEN Maxrho=MAX(rho(w)) ;g/cc
IF KEYWORD_SET(MeanR) THEN MeanR=INT_TABULATED(xm(w),rho(w)*ABS(xm(w))*xm(w)^2)/INT_TABULATED(xm(w),rho(w)*xm(w)^2);center of mass
;IF KEYWORD_SET(MeanR) THEN MeanR=INT_TABULATED(xm(w),rho(w)*ABS(xm(w)))/INT_TABULATED(xm(w),rho(w)) ;um
IF KEYWORD_SET(TwoSigma) THEN TwoSigma=2.*SQRT(INT_TABULATED(xm(w),rho(w)*(ABS(xm(w))-MeanR)^2)/INT_TABULATED(xm(w),rho(w))) ;um
IF KEYWORD_SET(RHOR) THEN rhoR=INT_TABULATED(xm(w),rho(w))/1.E4 ;g/cm2
IF KEYWORD_SET(MASS) THEN mass=4.*!pi*INT_TABULATED(xm(w),rho(w)*xm(w)^2)/1.E12 ;g
;PRINT,A(0:2)

;****ablated wave
;;rhomaxb=rhomax*EXP(-ABS(a(8))) ;max density
;rhomaxb=0.
;delb=180. ; scale length
;nub=0.7 ;superexponential power
;w=WHERE(xm^2 GT R0^2)
;IF w(0) NE -1 THEN rho(w)=rhomaxb*EXP(-((ABS(xm(w))-R0)/delb)^nub)
;print,rhomaxb/rhomax
;****************
;rho=SMOOTH(rho,3) ;remove infinitely sharp ablation front

;************CONVOLUTION
OD=abeltransform(xm,mu*rho,M,init=init)/1.E4 ;Need to make unitless: um*cm2/g*g/cm3=um/cm
OD=-ALOG(SMOOTH(EXP(-OD),pin/micpx)) ;Smooth the intensity profile (not the OD).

;GXD1 convolution function 
;yy=(FINDGEN(N)-N/2.+0.5*(N MOD 2))*micpx ;CCD pixels THIS CANNOT BE DEFINED WITHOUT dy. NEED TO FIX!!!!!
;cc=[0.78,0.2,0.02]& kk=[10.,300.,1400.]
;kernal=FLTARR(N) & FOR j=0,N_ELEMENTS(cc)-1 DO kernal=kernal+SQRT(!pi)*cc(j)/kk(j)*EXP(-!pi^2/kk(j)^2*(yy)^2)
;kernal=kernal/TOTAL(kernal)
;OD=-ALOG(IMSL_CONVOL1D(SHIFT(kernal,-N/2),EXP(-OD),/PERIODIC)) 
;***************
;Average over slot width
num=FIX(sl_dy/pin) ;number of y-steps over which to average OD
IF num GE 1 THEN OD=shelldy(xm,OD,sl_y0,sl_dy,num) ;Average over finite hohlraum slit width
;*********************

RETURN,OD
END

FUNCTION splinesmooth,x,xx,y,sig=sig
;returns a spline smoothed version of y(x), given at the new values xx
;sig is roughly equal to error^2 * N, where N is number of points in x.
Ys=IMSL_CSSMOOTH(x,y,SMPAR=sig) & yy = IMSL_SPVALUE(xx, Ys)
RETURN,yy
END

FUNCTION velocity,t,R,tt,smth=smth

IF KEYWORD_SET(smth) THEN BEGIN
    ;RR=SMOOTH(R,3)
    RR=splinesmooth(t,tt,R)
    v=DERIV(tt,-RR)
ENDIF ELSE BEGIN
    IF N_ELEMENTS(R) GT 3 THEN BEGIN
      RR=SMOOTH(R,1) 
      v=DERIV(t,-RR)
    ENDIF ELSE BEGIN
      v=REPLICATE(0.,N_ELEMENTS(R))
    ENDELSE
ENDELSE

RETURN,v
END

FUNCTION fitbkg,xm,y,OD,Npoly=Npoly,Iguess=Iguess,I0guess=I0guess,block=block
;fits a function to the inferred bkg, where y is the measured intensity (or ln(I))
;OD is the optical depth.
;'block' specifies which portion of the profile has no real data.

LogY=y
w=WHERE(y LE 0,COMPLEMENT=w0) & ymin=MIN(y(w0))
IF w(0) NE -1 THEN LogY(w)=ALOG(ymin) 
LogY(w0)=ALOG(y(w0)) 
;LogY=ALOG(y)

Iguess=OD
I0guess=LogY+Iguess

IF NOT(KEYWORD_SET(Npoly)) THEN Npoly=4

IF KEYWORD_SET(block) THEN BEGIN
  w=WHERE(xm GE MIN(block.y) AND xm LE MAX(block.y))
  err=xm/xm 
  IF w(0) NE -1 THEN err(w)=1.E4
  cff=POLY_FIT(xm,I0guess,Npoly,YFIT=I0fit,MEASURE_ERRORS=err)
ENDIF ELSE BEGIN
  err=xm/xm ;&  err(0)=1.E-4 & err(N_ELEMENTS(err)-1)=1.E-4
  cff=POLY_FIT(xm,I0guess,Npoly,YFIT=I0fit,MEASURE_ERRORS=err)
ENDELSE

RETURN,I0fit
END

PRO curvefit_convabl,y,a,f
;Function to be used with CURVEFIT
COMMON xparams,M,xm,blk,tg,ffitmurho,NNpoly,sskew,zz7,zz8,zz9,zz10,zz11,zz12

;a(3)=ABS(a(3))
;eta=a(0) & R0=a(1) & del=a(2) & x0=a(3) & pin=a(4) & micpx=a(5) ;& mur0=a(6)
N=N_ELEMENTS(y)

;a=[R0,maxrho,del,x0,slit,micpx,sl_y0,sl_dy,fbkg]
eta=a(0) & R0=a(1) & del=a(2) & x0=a(3) & pin=a(4) & micpx=a(5) & sl_y0=a(6) & sl_dy=a(7) ;& fbkg=a(8)

OD=shell(N,a,xm,rho,mu,M,TARGET=tg,FITMURHO=ffitmurho,SKEW=sskew)

yb=y;-fbkg ; substract guess of background

IF SIZE(blk,/TYPE) EQ 8 THEN I0fit=fitbkg(xm,yb,OD,I0guess=I0guess,block=blk,Npoly=NNpoly) $
ELSE I0fit=fitbkg(xm,yb,OD,I0guess=I0guess,Npoly=NNpoly)

f=(I0fit-I0guess) ;& print,TOTAL(bkguess*ALOG(bkguess))
;f=f+REPLICATE(del,N)^2*3.E-7;1E-5;add a regularization functional that biases towards thinner shells

RETURN
END

;FUNCTION lm_convabl,y,a,f
;;LEVENBERG-MARQUARDT DOES NOT SEEM TO WORK. 
;;Function to be used with LMFIT
;COMMON xparams,M,xm,zz2,zz3,zz4,zz5,zz6,zz7,zz8,zz9,zz10,zz11,zz12
;
;Na=N_ELEMENTS(a)
;a(3)=ABS(a(3))
;eta=a(0) & R0=a(1) & del=a(2) & x0=a(3) & pin=a(4) & micpx=a(5) & mur0=a(6)
;N=N_ELEMENTS(y)
;
;aa=a 
;P=shell(N,aa,xm,M) & bfit=fitbkg(xm,y,P,mur0,bkguess=bkguess) 
;f=TOTAL((bfit-bkguess)^2)
;
;da=0.001 & aa=[(1.+da)*a(0),a(1:Na-1)] 
;P=shell(N,aa,xm,M) & bfit=fitbkg(xm,y,P,mur0,bkguess=bkguess) 
;dfda0=(TOTAL((bfit-bkguess)^2)-f)/(da*a(0))
;
;da=0.001 & aa=[a(0),(1.+da)*a(1),a(2:Na-1)] 
;P=shell(N,aa,xm,M) & bfit=fitbkg(xm,y,P,mur0,bkguess=bkguess) 
;dfda1=(TOTAL((bfit-bkguess)^2)-f)/(da*a(1))
;
;da=0.001 & aa=[a(0),a(1),(1.+da)*a(2),a(3:Na-1)] 
;P=shell(N,aa,xm,M) & bfit=fitbkg(xm,y,P,mur0,bkguess=bkguess) 
;dfda2=(TOTAL((bfit-bkguess)^2)-f)/(da*a(2))
;
;print,f,dfda0,dfda1,dfda2
;RETURN,[f,dfda0,dfda1,dfda2,REPLICATE(0.,Na-3)]
;END


FUNCTION xstreak_fit,imstr,sub,tsub,psub,a,NPOLY=NPOLY,SKEW=SKEW,MIX=MIX,BLOCK=BLOCK,INIT=INIT,FITMURHO=FITMURHO
;sub(time,position) is the array containing the re-binned streaked image
;tsub is a vector of times at which to perform the inversions,
;psub are the corresponding pixel positions.
;a=[eta,R0,del,x0,slit,micpx,mur0]
;FITMURHO fits the product mu*rho (without having to determine mu during the fit). The average mu
;and hence rho, mass, rho*r etc. are determined afterwards using TARGET information. This procedure is 
;faster, smoother, but (presumably) slightly less accurate.
COMMON xparams

;**********
;IF KEYWORD_SET(TARGET) THEN tg=TARGET ELSE tg=0;
IF KEYWORD_SET(MIX) THEN mix=mix ELSE mix=0. ;amount of atomic-scale mix in target
IF KEYWORD_SET(FITMURHO) THEN ffitmurho=1. ELSE ffitmurho=0.;if=1 then fit mu*rho only, get average(mu) after fitting.
IF KEYWORD_SET(NPOLY) THEN NNpoly=Npoly ELSE NNpoly=4 ;I0 polynomial order
IF KEYWORD_SET(SKEW) THEN sskew=skew ELSE sskew=0.
WIDGET_CONTROL,imstr.freetype_id,GET_VALUE=free ;Which parameters to vary [R0,rho,del] - 1=vary,0=fix

tg=imstr.targetdata
;tg=LoadTarget(mix=mix) ;OLD version of xstreak uses this hardwired target input

;Also add:
;1) Gaussian skew
;2) Which parameters to freeze in fit

;********************Convert to a log image for analysis
;Imin=MIN(ABS(sub)) ;this is the minimum absolute value in the image. All values less than Imin are set to Imin.
;w=WHERE(sub LE 0,COMPLEMENT=w0) 
;Imin=MIN(sub(w0))
;IF w(0) NE -1 THEN sub(w)=Imin 
;sub=ALOG(sub) 
;********************

;a=[R0,maxrho,del,x0,slit,micpx,sl_y0,sl_dy,fbkg]
x0=a(3) & slit=a(4) & micpx=a(5) & sl_y0=a(6) & sl_dy=a(7) ;& fbkg=a(8)

sz=SIZE(sub) & snx=sz(1) & sny=sz(2)

aa=DBLARR(N_ELEMENTS(a),snx) ;create array to store fitted values of 'a'
xm=(DINDGEN(sny)-DOUBLE(x0))*micpx ;create spatial vector
OD=shell(sny,a,xm,rho,mu,M,func,TARGET=tg,init=init,FITMURHO=ffitmurho,SKEW=sskew) ;create the transformation matrix M

;********************
z=DBLARR(snx) & z2=DBLARR(snx,sny)  ;initializing array for structure
dat={target:tg,method:'Forward',nx:snx,ny:sny,t:z,pix:z,chisq:z,$
      mix:mix,fitmurho:ffitmurho,Npoly:NNpoly,skew:sskew,$
      a:aa,xm:xm,im:REBIN(z2,snx,sny,3),$ ;'im' contains data for the running profile plots
      maxdensity:z,R0:z,del:z,x0:z,slit:z,micpx:z,sl_y0:z,sl_dy:z,$
      v:z,rhoR:z,mass:z,av_mu:z,$
      xrhoR:z,xmass:z,dxrhoR:z,dxmass:z,$ ;rhoR/mur0 and mass/mur0 with errors
      dmaxdensity:z,dR0:z,ddel:z,dv:z,drhoR:z,dmass:z,$ ;errors
      inten:z2,bkg:z2,OD:z2,density:z2,mu:z2}

FOR i=0,snx-1 DO BEGIN
    y=REFORM(sub(i,*)) 
    ;w=WHERE(ABS(xm) GE ABS(a(1))-ABS(a(2)) AND ABS(xm) LE ABS(a(1))+ABS(a(2)))
    ;wghts=REPLICATE(1.,sny) & IF w(0) NE -1 THEN wghts(w)=wghts(w)+1.
    weights=REPLICATE(1.,sny);1./y^4;;wghts;1./y;REPLICATE(1.,sny) ;xm^2;
    ;weights=1./ABS(y)
    IF KEYWORD_SET(block) THEN BEGIN
    ;only send block roi info if block roi is within current portion of invert roi
      IF psub(i) GE MIN(block.p) AND psub(i) LE MAX(block.p) AND $
         MAX(xm) GE MIN(block.y) AND MIN(xm) LE MAX(block.y) THEN blk=block ELSE blk=0 
      w=WHERE(xm GE MIN(block.y) AND xm LE MAX(block.y))
      IF w(0) NE -1 THEN weights(w)=0.
    ENDIF ELSE blk=0
    ;SETTING blk=0 does two things:
    ;(i) CURVEFIT weights the blocked roi region with zero.
    ;(ii)POLY_FIT for the background weights the blocked roi region with large errors.
    ;y=y+y*randomn(seed,n_elements(y))*0.05 ;introduce some noise
                                               ;free=whether [R0,Rho,del] is allowed to vary
    result=CURVEFIT(y,DBLARR(sny),weights,a,siga,FITA=[free,0,0,0,0,0,0],CHISQ=chisqr,$
        /NODERIVATIVE,STATUS=stat,FUNCTION_NAME='curvefit_convabl',TOL=1.E-6,ITMAX=2000,/DOUBLE)
;    result=CURVEFIT(y,DBLARR(sny),weights,a,siga,FITA=[1,0,0,0,0,0,0,0,1],CHISQ=chisqr,$
;        /NODERIVATIVE,STATUS=stat,FUNCTION_NAME='curvefit_convabl',TOL=1.E-6,ITMAX=2000,/DOUBLE)
;    print,stat,a(8)
;    result=CURVEFIT(y,DBLARR(sny),weights,a,siga,FITA=[free,0,0,0,0,0,0],CHISQ=chisqr,$
;        /NODERIVATIVE,STATUS=stat,FUNCTION_NAME='curvefit_convabl',TOL=1.E-6,ITMAX=2000,/DOUBLE)
    IF stat NE 0 THEN BEGIN
      PRINT,"Failed to fit at t = "+STRING(tsub(i),FORMAT='(F15.2)')+" ps";+"  a(0:2)="+STRING(a(0:2))
    ENDIF
    ;print,"Tail density, g/cc = "+STRING(a[8])
    print,"Chi-Squared = "+STRING(chisqr*1000.)
    ;PRINT,a(0),a(1),a(2)
;    result=CURVEFIT(y,FLTARR(sny),weights,a,siga,FITA=[1,1,1,0,0,0,0,0],CHISQ=chisqr,$
;        /NODERIVATIVE,STATUS=stat,FUNCTION_NAME='curvefit_convabl',TOL=1.E-6,ITMAX=2000,/DOUBLE)
;    IF stat NE 0 THEN PRINT,"Failed to fit at t = "+STRING(tsub(i),FORMAT='(F15.2)')+" ps"

;LEVENBERG-MARQUARDT DOES NOT SEEM TO WORK. 
;   result=LMFIT(y,FLTARR(sny),a,CHISQ=chisqr,CONVERGENCE=stat,/DOUBLE,$
;      FITA=[1,1,1,0,0,0,0,0,0,0],FUNCTION_NAME='lm_convabl',ITER=Niter,$
;      ITMAX=2000,ITMIN=10,MEASURE_ERRORS=weights,SIGMA=siga,TOL=1.E-6) 
;    IF stat NE 1 THEN PRINT,"Failed to fit at t = "+STRING(tsub(i),FORMAT='(F15.2)')+" ps"


;a=[R0,eta,del,x0,slit,micpx,sl_y0,sl_dy]

    a(0)=ABS(a(0)) & a(1)=ABS(a(1)) & a(2)=ABS(a(2))
    dat.a(*,i)=a
    dat.t(i)=tsub(i)  & dat.pix(i)=psub(i)
    dat.x0(i)=a(3)   & dat.slit(i)=a(4) & dat.micpx=a(5) 
    dat.sl_y0=a(6)   & dat.sl_dy=a(7) 
    dat.chisq(i)=chisqr
    ;****Calculate projected thickness (P) and background profile for this best fit
    maxrho=1. & meanR=1. & twosigma=1. & mass=1. & rhoR=1. & av_mu=1.;initialize parameters to trigger KEYWORD_SET
    OD  =shell(sny,a,xm,rho,mu,M,TARGET=tg,FITMURHO=ffitmurho,SKEW=sskew,$
              MAXRHO=maxrho,MEANR=meanR,TWOSIGMA=twosigma,MASS=mass,RHOR=rhoR)
    I0fit=fitbkg(xm,y,OD,Iguess=Iguess,I0guess=I0guess,block=blk,Npoly=NNpoly)
    
    IF KEYWORD_SET(fitmurho) THEN BEGIN ;calculate average mu after the fit
;      IF KEYWORD_SET(target) THEN BEGIN ;only do it if target data is available
        
        ;If max mu*mass is greater than value at the maximum mass then 
        ;set av_mu equal to that at the maximum mass to avoid extrapoloation problems.
        IF mass*1.E6 GT MAX(tg.mumass) THEN av_mu=tg.av_mu(tg.N-1) $
        ELSE av_mu=INTERPOL(tg.av_mu,tg.mumass,mass*1.E6) ;get av_mu based on fitted mu*mass
        
        maxrho=maxrho/av_mu & rhor=rhor/av_mu & mass=mass/av_mu
        rho=rho/av_mu       & mu=REPLICATE(av_mu,N_ELEMENTS(mu))
 ;     ENDIF
    ENDIF ELSE BEGIN
      ;IF KEYWORD_SET(target) THEN 
      mmu=XmuFromMass(tg,xm,rho,av_mu) ;run this just to extract av_mu (not mu itself)
    ENDELSE
    
    ;dat.maxdensity(i)=a(0)*dat.rho0(i)& dat.R0(i)=a(1) & dat.del(i)=a(2)
    dat.R0(i)=ABS(meanR)
    dat.maxdensity(i)=maxrho 
    dat.del(i)=twosigma ;This allows a(0),a(1) etc. to be used for parameters other than maxeta, R0, etc.

;mm=INT_TABULATED(xm*1.E-4,P*ABS(xm)*1.E-8)*1*!pi*dat.rho0(i)*1.E6 ;this is identical to the density-integrated mass
;w=WHERE(ABS(xm) LE dat.R0(i)+1.5*dat.del(i)) & ;proj=(I0fit(w)-y(w))/(dat.mur0(i)*1.E4)*dat.rho0(i)
;mm=INT_TABULATED(xm(w)*1.E-4, proj(w)*ABS(xm(w))*1.E-4 )*1.*!pi*1.E6 ;this is noisier than the density-integrated mass
    ;****Store results for this time step
    dat.rhoR(i)=rhoR*1.E3 ;mg/cm^2 
    dat.mass(i)=mass*1.E6 ;ug
    dat.av_mu(i)=av_mu ;cm2/g

;    dat.xrhoR(i)=rhoR*dat.mur0(i)*1.E4 ;no units ("radiographed rhoR")
;    dat.xmass(i)=mass*dat.mur0(i)*1.E4 ;cm^2 ("radiographed mass")
    dat.inten(i,*)= EXP(-OD);fitted x-ray intensity (normalized)
    dat.bkg(i,*)= EXP(I0guess)            ;fitted background intensity
    dat.OD(i,*)= OD                            ;fitted optical depth
    dat.density(i,*) = rho    ;density
    dat.mu(i,*) = mu    ;attenuation coefficient, cm2/g
    ;****Calculate and store errors
    ;dat.dR0(i)=siga(0)
    ;dat.dmaxdensity(i)=siga(1)
    ;dat.ddel(i)=siga(2)
    ;**********
    OUTPLOT,xm,(y),imstr,/newplot,$
      XTITL="Distance (mic)",YTITL="Intensity",$
      TITL="t = "+STRCOMPRESS(STRING(tsub(i),FORMAT='(F8.1)'),/REMOVE_ALL)+" ns"
    ;OUTPLOT,xm,EXP(I0guess),imstr,/overplot,color=[255,0,255]
    OUTPLOT,xm,EXP(I0fit),imstr,/overplot,color=[0,0,255]
    OUTPLOT,xm,EXP(I0fit-OD),imstr,/overplot,color=[255,0,0]
 
    yfit=EXP(I0fit-OD)
    ;OUTPLOT,[-dat.R0[i],dat.R0[i]],INTERPOL(yfit,xm,[-dat.R0[i],dat.R0[i]]),$
    ;  imstr,/overplot,COLOR=[255,0,0],PSY=24
    ;OUTPLOT,xm,rho*0.05,imstr,/overplot,COLOR=[0,255,0]
    
    dat.im[i,*,0]=(y)
    dat.im[i,*,1]=EXP(I0fit)
    dat.im[i,*,2]=EXP(I0fit-OD)
    
;    PLOT,xm,y,/YSTYLE,/XSTYLE,CHARSIZE=charsz,$
;        XTITLE="!4l!3m",TITLE="t = "+STRING(tns(i),FORMAT='(F5.2)')+" ns"
;    OPLOT,xm,I0fit-Iguess,COLOR=RED
;    OPLOT,xm,I0fit,COLOR=RED,LINESTYLE=2
ENDFOR
   
dat.v=velocity(dat.t,dat.R0,dat.t,smth=0)*1.E0 ;velocity from t vs R0 fit, convert to mic/ns

RETURN,dat
END


;******************************
;******************************
;******************************
;******************************
;Here begins a region to test new inversion techniques

FUNCTION readdat,file,ncols
;Reads data in columns
;'file' is the filename
;'ncols' is the number of columns in the data
;This is an efficient code that assumes anything which is not formatted 
;as floating point with 'ncols' columns is just a comment
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

RETURN,data
END

PRO xconvablfunc2,y,a,f
COMMON xparams,M,xm,zz2,zz3,zz4,zz5,zz6,zz7,zz8,zz9,zz10,zz11,zz12

a(3)=ABS(a(3))
eta=a(0) & R0=a(1) & del=a(2) & x0=a(3) & pin=a(4) & micpx=a(5) & mur0=a(6)
N=N_ELEMENTS(y)

P=shell(N,a,xm,M)
bkguess=(y+mur0*P)
;f=(bkguess-MEAN(bkguess))^2 + REPLICATE(del,N)^2*1.E-8
f=ABS(DERIV(bkguess))+REPLICATE(del,N)^2*1.E-8

;bfit=fitbkg(xm,y,P,mur0,bkguess=bkguess)
;bkguess=y+mur0*P
;q=1. & dx=250.
;wlo=WHERE(xm GE -R0-q*dx AND xm LE -R0+q*dx)
;whi=WHERE(xm GE R0-q*dx AND xm LE R0+q*dx)
;f=FLTARR(N)
;IF wlo(0) NE -1 THEN BEGIN
;  bfitlo=fitbkg(xm(wlo),y(wlo),P(wlo),mur0,bkguess=bkguesslo) 
;  f(wlo)=(bfitlo-bkguesslo)^2
;ENDIF 
;IF whi(0) NE -1 THEN BEGIN
;  bfithi=fitbkg(xm(whi),y(whi),P(whi),mur0,bkguess=bkguesshi)
;  f(whi)=(bfithi-bkguesshi)^2
;ENDIF


RETURN
END


PRO testxmuxstreak
;Tests iteration method that calculates mu for each element of density distribution based on 
;cumulative mass from inside out

N=501
xmin=-1000. & xmax=800.
xm=DINDGEN(N)/(N-1)*(xmax-xmin)+xmin
 
r0=800. & del=40. & maxrho=7. & x0=(N-1)/2. & slit=5. & micpx=xm(1)-xm(0) & sl_y0=0. & sl_dy=0.

alpha=-0. & s=(ABS(xm)-R0)/(del/2.) 
rho=maxrho*2.*SQRT(2.*!PI)*(EXP(-s^2/2.)/SQRT(2.*!PI))*(0.5*(1+IMSL_ERF(alpha*s/SQRT(2))))

;rho=maxrho*DOUBLE(ABS(ABS(xm)-ABS(R0)) LE ABS(del)/2.)
a=[R0,maxrho,del,x0,slit,micpx,sl_y0,sl_dy]

target=loadtarget()
mu=xmufrommass(target,xm,rho,av_mu,cmass)
maxrho=1. & meanr=1. & twosigma=1. & mass=1. & rhor=1.
OD=shell(N,a,xm,rho,mu,M,func,INIT=1,TARGET=TARGET,MAXRHO=MAXRHO,MEANR=MEANR,TWOSIGMA=TWOSIGMA,MASS=MASS,RHOR=RHOR)

!P.MULTI=[0,2,2]

PLOT,xm,EXP(-od),XTITLE="Position (um)",YTITLE="Transmission"
PLOT,xm,mu,XTITLE="Position (um)",YTITLE="mu (cm2/g)"
OPLOT,xm,av_mu,LINESTYLE=2,COLOR=255
OPLOT,xm,rho,LINESTYLE=3
PLOT,xm,cmass

!P.MULTI=0
RETURN
END

PRO CreateStreak
;This creates a radiograph with known input density profiles. It can then be 
;analyzed with xstreak to see if the inputs can be recovered.

Nr=1024 & Nt=512

maxrho=9. & del=40. 
Rmax=600. 
R0max=500. & R0min=100.

pin=25. ;blurring/resolution

tmax=1000.

r  =DINDGEN(Nr/2.)/(Nr/2.-1)*Rmax
R0 =REVERSE(DINDGEN(Nt)/(Nt-1)*(R0max-R0min)+R0min)
t  =DINDGEN(Nt)/(Nt-1)*tmax

dr=r(1)-r(0)
dt=t(1)-t(0)

rho=DBLARR(Nt,Nr/2.)
mu=DBLARR(Nt,Nr/2.)
OD1=DBLARR(Nt,Nr/2.)
OD=DBLARR(Nt,Nr-1)

vel=DERIV(t,R0)
rhor=DBLARR(Nt)
mass=DBLARR(Nt)

tg=loadtarget()

i=0.
rho(i,*)=maxrho*EXP(-(ABS(r)-R0(i))^2/2/(del/2)^2)
mu(i,*) =XmuFromMass(tg,r,REFORM(rho(i,*)),av_mu)
OD1(i,*)=abeltransform(r,REFORM(mu(i,*))*REFORM(rho(i,*)),M,/init)/1.E4 ;Make unitless: um*cm2/g*g/cm3=um/cm
;OD(i,*) =[REVERSE(REFORM(OD1(i,1:Nr/2.-1))),REFORM(OD1(i,*))]
;OD(i,*)=SMOOTH(OD(i,*),pin/dr)

FOR i=0,Nt-1 DO BEGIN
  rho(i,*)=maxrho*EXP(-(ABS(r)-R0(i))^2/2/(del/2)^2)
  mu(i,*) =XmuFromMass(tg,r,REFORM(rho(i,*)),av_mu)
  OD1(i,*)=abeltransform(r,REFORM(mu(i,*))*REFORM(rho(i,*)),M)/1.E4 ;Make unitless: um*cm2/g*g/cm3=um/cm
  OD(i,*) =[REVERSE(REFORM(OD1(i,1:Nr/2.-1))),REFORM(OD1(i,*))]
  OD(i,*)=SMOOTH(OD(i,*),pin/dr)
  rhor(i)=INT_TABULATED(r,rho(i,*))/1.E1 ;mg/cm2
  mass(i) =4.*!pi*INT_TABULATED(r,rho(i,*)*r^2)/1.E12 *1.E6;ug
ENDFOR
image=EXP(-OD)

PRINT,"mic/pixel= ",dr
PRINT,"ps/pixel= ",dt
PRINT,"slit (um)= ",pin
PRINT,"vel (um/ns)= ",-1.E3*mean(vel)
SAVE,image,FILENAME='c:\hicks\nic\nif-radiography\simulations\simcam\testimage.sav'

;PLOT,t,-1.E3*vel
;PLOT,t,rhor
PLOT,t,mass
;PLOT,r,OD1(Nt/2,*)
;TVSCL,mu
;PM,[[t],[mass]]
;print,(mass)
RETURN
END

PRO testconvol,eps=eps
;tests convolution function

RED=255 & GREEN=255*256L & BLUE=255*256L*256L & YELLOW=RED+GREEN & PINK=RED+BLUE & LBLUE=GREEN+BLUE
ORANGE=255+126*256L & PURPLE=126+BLUE & PLAIN=RED+GREEN+BLUE & BACKG=0 & BROWN=174+71*256L
charsz=1.2 & thk=1.1 & symsz=1.3

!P.FONT=-1
IF KEYWORD_SET(eps) OR KEYWORD_SET(ps) THEN BEGIN
!P.FONT=1
    SET_PLOT,'ps'
    !X.MARGIN=[6.5,1.5] & !Y.MARGIN=[3.2,2]    
    DEVICE,ENCAPSULATED=1,LANDSCAPE=0,XSIZE=15,YSIZE=14,/COLOR,/HELVETICA,FILENAME="C:\Documents and Settings\hicks13\DESKTOP\IDL.EPS"
    TVLCT, [0,255,0,0,255,255,0,255,126,255,174], [0,0,255,0,255,0,255,126,0,255,71], [0,0,0,255,0,255,255,0,255,255,0]
    PLAIN=0 & RED=1 & GREEN=2 & BLUE=3 & YELLOW=4 & PINK=5 & LBLUE=6 & ORANGE=7 & PURPLE=8 & BACKG=9 & BROWN=10
    FILL=1 & charsz=1.2 & thk=3.0*thk & symsz=0.7*symsz
ENDIF
!P.MULTI=[0,2,2]
clr=[PLAIN,RED,GREEN,BLUE,BROWN,PINK,LBLUE,ORANGE,PURPLE,YELLOW]
b = FINDGEN(17)*(!PI*2/16.) & USERSYM,1*COS(b),1*SIN(b),THICK=thk,FILL=1
mic=STRING(["265B])+'m'

N=1024
;ymin=-10. & ymax=10.
y=FINDGEN(N)-N/2.+0.5*(N MOD 2);/(N-1)*(ymax-ymin)+ymin
;I=0.2*FLOAT(ABS(y) GT 2. AND ABS(y) LT 50.)
I0=EXP(-y^2/(2*(300.)^2))

;a=[R0,eta,del,x0,slit,micpx,sl_y0,sl_dy]
R0=300.  & rhomax=20. &  del=70.  & x0=N/2. & pin=25. & micpx=1.    & sl_y0=0. & sl_dy=0.
a=[R0,rhomax,del,x0,pin,micpx,sl_y0,sl_dy]
OD=shell(N,a,y,rho,mu,M,func,$
  /INIT,TARGET=TARGET,MAXRHO=MAXRHO,MEANR=MEANR,TWOSIGMA=TWOSIGMA,MASS=MASS,RHOR=RHOR,FITMURHO=0,SKEW=0)

Rmax=25. & murho=1./4.15 ;3.47um for 8.35keV ;4.15um for 8.95 keV; 25um for 30 keV
w=WHERE(ABS(y) LT 25.) & OD(w)=2.*murho*Rmax*SQRT(1-(y(w)/Rmax)^2) ;add central wire

I=I0*EXP(-OD)
;I=SMOOTH(I,35,/EDGE_TRUNCATE)

c=[0.78,0.2,0.02]& k=[10.,300.,1400.]
kernal=FLTARR(N) & FOR j=0,N_ELEMENTS(c)-1 DO kernal=kernal+SQRT(!pi)*c(j)/k(j)*EXP(-!pi^2/k(j)^2*y^2)
;kernal=EXP(-yk^2/2./sig^2) 
kernal=kernal/TOTAL(kernal)
;kernal=SHIFT(kernal,-N/2)
;conv=fftconvolve(I,kernal)
;conv=BLK_CON(kernal(1:N-2),I)
;conv=CONVOL(I,SHIFT(kernal,-N/2),/EDGE_ZERO)
;conv=IMSL_CONVOL1D(SHIFT(kernal,-N/2),I,/PERIODIC) 
;conv=CONVOLVE(I,SHIFT(kernal,-N/2),FT_PSF=psf_ft,/CORREL, NO_FT=1 )
;conv=I0*conv & I=I0*I

Niter=20
FOR ii=1,Niter DO Max_Entropy,conv,kernal,deconv,multipliers
reconv=IMSL_CONVOL1D(SHIFT(kernal,-N/2),deconv,/PERIODIC)
;deconv=deconv_SOP(conv,kernal,reconv) ;ry is the reconvolved data
;deconv=REAL_PART(deconv)

PLOT,y,kernal,YLOG=1,yrange=[1.E-6,MAX(I)],XTHICK=thk,YTHICK=thk,CHARTHICK=thk,THICK=thk,$
  XTITLE='Pixels @ CCD',YTITLE=''
;OPLOT,y,I

PLOT,y,I,YRANGE=[-0.05,0.8],/XSTYLE,/YSTYLE,XTHICK=thk,YTHICK=thk,CHARTHICK=thk,THICK=thk,$
  XTITLE='Pixels @ CCD',YTITLE=''
OPLOT,y,conv,THICK=thk,COLOR=RED
OPLOT,y,SMOOTH(I,35,/EDGE_TRUNCATE),THICK=thk,COLOR=BLUE,LINESTYLE=1
;OPLOT,y,SMOOTH(I,80,/EDGE_TRUNCATE),COLOR=255
OPLOT,y,deconv,COLOR=GREEN,LINESTYLE=1
OPLOT,y,reconv,COLOR=GREEN,LINESTYLE=2

limbI=I(where(y GT 200 and y lt 400))
limbC=conv(where(y GT 200 and y lt 400))
print,MIN(limbI)/MAX(limbI),MIN(limbC)/MAX(limbC)

!P.MULTI=0
IF KEYWORD_SET(eps) OR KEYWORD_SET(ps)  THEN BEGIN
    DEVICE,/CLOSE
    SET_PLOT,'WIN'
ENDIF
!P.FONT=-1
END

FUNCTION create2Dpsf,saveimage=saveimage,noshow=noshow
;creates a 2D psf image

filename='c:\hicks\nif\convabl\results\GXD1_psf_v1.tif'
Nx=4200 & Ny=4200
x=REBIN(FINDGEN(Nx)-Nx/2.+0.5*(Nx MOD 2),Nx,Ny)
y=TRANSPOSE(REBIN(FINDGEN(Ny)-Ny/2.+0.5*(Ny MOD 2),Ny,Nx))

c=[0.78,0.2,0.02]& k=[10.,300.,1400.]
kernal=FLTARR(Nx,Ny) 
FOR j=0,N_ELEMENTS(c)-1 DO kernal=kernal+!pi*c(j)/k(j)^2*EXP(-!pi^2/k(j)^2*(x^2+y^2))
kernal=kernal/TOTAL(kernal)

IF KEYWORD_SET(saveimage) THEN WRITE_TIFF,filename,kernal,/FLOAT

;PRINT,MAX(kernal),MIN(kernal),TOTAL(kernal)
IF NOT(KEYWORD_SET(noshow)) THEN TVSCL,CONGRID(ALOG(kernal),500,500)

RETURN,kernal
END

PRO test2Dpsf,eps=eps
;explores convolution/deconvolution of 2D psf
RED=255 & GREEN=255*256L & BLUE=255*256L*256L & YELLOW=RED+GREEN & PINK=RED+BLUE & LBLUE=GREEN+BLUE
ORANGE=255+126*256L & PURPLE=126+BLUE & PLAIN=RED+GREEN+BLUE & BACKG=0 & BROWN=174+71*256L
charsz=1.0 & thk=1.1 & symsz=1.3

!P.FONT=-1
IF KEYWORD_SET(eps) OR KEYWORD_SET(ps) THEN BEGIN
!P.FONT=1
    SET_PLOT,'ps'
    !X.MARGIN=[6.5,1.5] & !Y.MARGIN=[3.2,2]    
    DEVICE,ENCAPSULATED=1,LANDSCAPE=0,XSIZE=8,YSIZE=7,/COLOR,/HELVETICA,FILENAME="C:\Documents and Settings\hicks13\DESKTOP\IDL.EPS"
    TVLCT, [0,255,0,0,255,255,0,255,126,255,174], [0,0,255,0,255,0,255,126,0,255,71], [0,0,0,255,0,255,255,0,255,255,0]
    PLAIN=0 & RED=1 & GREEN=2 & BLUE=3 & YELLOW=4 & PINK=5 & LBLUE=6 & ORANGE=7 & PURPLE=8 & BACKG=9 & BROWN=10
    FILL=1 & charsz=1.0 & thk=3.0*thk & symsz=0.7*symsz
ENDIF
;!P.MULTI=[0,2,2]
clr=[PLAIN,RED,GREEN,BLUE,BROWN,PINK,LBLUE,ORANGE,PURPLE,YELLOW]

;**************Create 2D PSF
Nx0=4200 & Ny0=4200 ;original PSF size
Nx=1024 & Ny=1024 ;PSF is scaled to this size
x=REBIN(FINDGEN(Nx)-Nx/2.+0.5*(Nx MOD 2),Nx,Ny)*Nx0/Nx
y=TRANSPOSE(REBIN(FINDGEN(Ny)-Ny/2.+0.5*(Ny MOD 2),Ny,Nx))*Ny0/Ny
c=[0.78,0.2,0.02]& k=[10.,300.,1400.]
psf=FLTARR(Nx,Ny) 
FOR j=0,N_ELEMENTS(c)-1 DO psf=psf+!pi*c(j)/k(j)^2*EXP(-!pi^2/k(j)^2*(x^2+y^2))
psf=psf/TOTAL(psf)
;***************Create second PSF (for imaging slot)
sigx=35./2.35 & sigy=2.
psf2=EXP(-x^2/(2.*sigx^2)-y^2/(2.*sigy^2)) & psf2=psf2/TOTAL(psf2)
;***************Create a single PSF by convolving the two PSF's
psf=convolve(psf2,psf) & psf=psf/TOTAL(psf)
;***************Create 2D image
N=Nx
R0=400.  & rhomax=90. &  del=15.  & x0=N/2. & pin=2. & micpx=2.    & sl_y0=0. & sl_dy=0.
y=(FINDGEN(N)-N/2.+0.5*(N MOD 2))*micpx
I0=EXP(-y^2/(2*250.^2))
a=[R0,rhomax,del,x0,pin,micpx,sl_y0,sl_dy]
OD=shell(N,a,y,rho,mu,M,func,$
  /INIT,TARGET=TARGET,MAXRHO=MAXRHO,MEANR=MEANR,TWOSIGMA=TWOSIGMA,MASS=MASS,RHOR=RHOR,FITMURHO=0,SKEW=0)
Rmax=25. & murho=1./4.15 ;3.47um for 8.35keV ;4.15um for 8.95 keV; 25um for 30 keV
w=WHERE(ABS(y) LT 25.) & OD(w)=2.*murho*Rmax*SQRT(1-(y(w)/Rmax)^2) ;add central wire
Ix=I0*EXP(-OD)
;Ix=SMOOTH(Ix,20) 
image=REBIN(Ix,Nx,Ny) ;turn into a 2D image
;***************Convolve image and PSF
imconv=convolve(image,psf)
imconvn=imconv+RANDOMN(SYSTIME(1),Nx,Ny)*0.0
;***************rebin
Nx=1024 & Ny=1024
y=REBIN(y,Nx)
image=REBIN(image,Nx,Ny)
psf=REBIN(psf,Nx,Ny) & psf=psf/TOTAL(psf)
imconv=REBIN(imconv,Nx,Ny)
imconvn=REBIN(imconvn,Nx,Ny)

;***************Deconvolve image
Niter=15 & FOR ii=1,Niter DO Max_Entropy,imconvn,psf,imdeconv,multipliers,FT_PSF=psf_ft
;imdeconv=image_deconvolve(imconvn,psf,0.01*imconv)

;***************
PLOT,y,image(*,Ny/2.),XRANGE=[-600,600],THICK=thk,XTHICK=thk,YTHICK=thk,CHARTHICK=thk,CHARSIZE=charsz,$
  XTITLE="Radius, !Mmm",YTITLE="Intensity (arb.)"
;OPLOT,y,imconv(*,Ny/2.),THICK=thk,COLOR=RED
OPLOT,y,imconvn(*,Ny/2.),THICK=thk,COLOR=RED
OPLOT,y,imdeconv(*,Ny/2.),COLOR=GREEN,THICK=thk,LINESTYLE=2

;***************

!P.MULTI=0
IF KEYWORD_SET(eps) OR KEYWORD_SET(ps)  THEN BEGIN
    DEVICE,/CLOSE
    SET_PLOT,'WIN'
ENDIF
!P.FONT=-1

RETURN
END

PRO test2Dpsf2
;creates a test Gaussian image, and blurs it with psf and a slit function.
;Saves image to be used as a test on the xstreak deconvolution process

;**************Create 2D PSF
Nx0=4200 & Ny0=4200 ;original PSF size
Nx=4200 & Ny=4200 ;PSF is scaled to this size
x=REBIN(FINDGEN(Nx)-Nx/2.+0.5*(Nx MOD 2),Nx,Ny)*Nx0/Nx
y=TRANSPOSE(REBIN(FINDGEN(Ny)-Ny/2.+0.5*(Ny MOD 2),Ny,Nx))*Ny0/Ny
psf=FLTARR(Nx,Ny) 
;*****Create image
sigx=8. & sigy=500.
im=EXP(-x^2/(2.*sigx^2)-y^2/(2.*sigy^2))
;****GXD1
;c=[0.78,0.2,0.02]& k=[10.,300.,1400.]
;psf=FLTARR(Nx,Ny) 
;FOR j=0,N_ELEMENTS(c)-1 DO psf=psf+!pi*c(j)/k(j)^2*EXP(-!pi^2/k(j)^2*(x^2+y^2))
;psf=psf/TOTAL(psf)
;****GXD3
c=[0.719,0.173,0.0267,0.0813]& k=[9.7,319.,1176.,670.]
FOR j=0,2 DO psf=psf+!pi*c(j)/k(j)^2*EXP(-!pi^2/k(j)^2*(x^2+y^2))
psf=psf+c(3)*0.548/(40000+(SQRT(x^2+y^2))^2.4)
psf=psf/TOTAL(psf)
;***************Create second PSF (for imaging slot)
slit=25. & micpx=1.
;sigx=slit/micpx/2.7 & sigy=4.;The 2.7 was calibrated to properly deconvolve a boxcar using a Gaussian
;psf2=EXP(-x^2/(2.*sigx^2)-y^2/(2.*sigy^2)) & psf2=psf2/TOTAL(psf2)

xblur=slit/micpx & yblur=3.
psf2=FLTARR(Nx,Ny) & psf2[WHERE(ABS(x) LE xblur/2. AND ABS(y) LE yblur/2.)]=1. & psf2=psf2/TOTAL(psf2)

psf=convolve(psf,psf2) & psf=psf/TOTAL(psf) & psf2=0
;psf=psf2 & psf=psf/TOTAL(psf) & psf2=0

;slit=25.
;ims=SMOOTH(im,[slit,1])
;sigx=35./2.35 & sigy=2.
;psf2=EXP(-x^2/(2.*sigx^2)-y^2/(2.*sigy^2)) & psf2=psf2/TOTAL(psf2)
;***************Create a single PSF by convolving the two PSF's
imconv=convolve(im,psf)

image=imconv
SAVE,image,FILENAME="c:\Hicks\NIF\ConvAbl\SimCam\Hyades\8umGaussian+GXD3+25umBoxcar.sav"
image=im
SAVE,image,FILENAME="c:\Hicks\NIF\ConvAbl\SimCam\Hyades\8umGaussian.sav"


RETURN
END

PRO testTail
;examines an adapted density profile to take into account high opacity blow-off
RED=255 & GREEN=255*256L & BLUE=255*256L*256L & YELLOW=RED+GREEN & PINK=RED+BLUE & LBLUE=GREEN+BLUE
ORANGE=255+126*256L & PURPLE=126+BLUE & PLAIN=RED+GREEN+BLUE & BACKG=0 & BROWN=174+71*256L
charsz=1.0 & thk=1.1 & symsz=1.3



N=1000 & xmax=500. & xmin=-500. & xm=FINDGEN(N)/(N-1)*(xmax-xmin)+xmin

skew=0. & R0=350. & del=50. & rhomax=6.
alpha=skew & s=(ABS(xm)-R0)/(del/2.) & rho=rhomax*2.*SQRT(2.*!PI)*(EXP(-s^2/2.)/SQRT(2.*!PI))*(0.5*(1+IMSL_ERF(alpha*s/SQRT(2)))) & func='Skewed Gaussian' 

skew=0. & R0=R0+del & del=del*2. & rhomax=1.;rhomax
alpha=skew & s=(ABS(xm)-R0)/(del/2.) & tail=rhomax*2.*SQRT(2.*!PI)*(EXP(-s^2/2.)/SQRT(2.*!PI))*(0.5*(1+IMSL_ERF(alpha*s/SQRT(2)))) & func='Skewed Gaussian' 

;tail=1/(1+EXP(-(ABS(xm)-R0-del)/(0.5*del)))
PLOT,xm,sqrt(rho^2+tail^2),XRANGE=[200,500],THICK=thk
OPLOT,xm,rho,LINESTYLE=1
OPLOT,xm,tail,LINESTYLE=2

RETURN
END

;
;PRO testinversion
;;tests a new inversion scheme
;COMMON xparams
;
;RED=255 & GREEN=255*256L & BLUE=255*256L*256L & YELLOW=RED+GREEN & PINK=RED+BLUE & LBLUE=GREEN+BLUE
;ORANGE=255+126*256L & PURPLE=126+BLUE & PLAIN=RED+GREEN+BLUE & BACKG=0 & BROWN=174+71*256L
;charsz=1.2 & thk=1.5 & symsz=0.4
;
;file='c:\hicks\omega\data\0712\results\49780_IvsY_03.txt'
;;dat=readxvisdat(file,x,y,ysig)
;dat=readdat(file,4)
;xm=REFORM(dat(0,*)) & y=REFORM(dat(1,*))
;y=-ALOG(y/MAX(y))
;N=N_ELEMENTS(y)
;
;!P.MULTI=[0,1,2]
;
;;************Plot measured intensity distribution and FFT
;PLOT,xm,y,/YSTYLE,xrange=[0,300]
;pp0=!p & px0=!x & py0=!y; capture system variables for this plot
;PLOT,FINDGEN(N)-N/2,SHIFT(ABS(FFT(y)),N/2),YLOG=1,XRANGE=[0,100]
;pp1=!p & px1=!x & py1=!y; capture system variables for this plot
;
;;Setup array of shell parameters with initial guesses
;eta=1.8/2.015 & R0=215. & del=25. & x0=1032 & slit=10.
;micpx=0.875 & mur0=28.9350*1.E-4 & rho0=2.015
;sl_y0=0. & sl_dy=0.
;a=[eta,R0,del,x0,slit,micpx,mur0,rho0,sl_y0,sl_dy]
;P=shell(N,a,xm,M,func,INIT=1,MASS=MASS,RHOR=RHOR)
;bkguess=y+mur0*P
;
;;************Plot intensity ditribution after subtracting initial guess
;!p=pp0 & !x=px0 & !y=py0
;OPLOT,xm,bkguess,COLOR=GREEN
;!p=pp1 & !x=px1 & !y=py1
;OPLOT,FINDGEN(N)-N/2,SHIFT(ABS(FFT(bkguess)),N/2),COLOR=GREEN
;OPLOT,FINDGEN(N)-N/2,SHIFT(ABS(FFT(mur0*P)),N/2),COLOR=GREEN,LINESTYLE=2
;
;;************Calculate best fit
;weights=REPLICATE(1.,N) ;1/y;xm^2;
;result=CURVEFIT(y,FLTARR(N),weights,a,siga,FITA=[1,1,1,0,0,0,0,0,0,0],CHISQ=chisqr,$
;        /NODERIVATIVE,STATUS=stat,FUNCTION_NAME='xconvablfunc2',TOL=1.E-5,/DOUBLE)
;
;maxeta=1. & meanR=1. & twosigma=1. & mass=1. & rhoR=1. & compression=1 ;initialize parameters to trigger KEYWORD_SET
;P=shell(sny,a,xm,M,MAXETA=maxeta,MEANR=meanR,TWOSIGMA=twosigma,MASS=mass,RHOR=rhoR,compression=compression)
;bkguess=y+mur0*P
;print,maxeta,meanR,twosigma
;
;;*************Plot final intensity distribution
;!p=pp0 & !x=px0 & !y=py0
;OPLOT,xm,bkguess,COLOR=RED
;!p=pp1 & !x=px1 & !y=py1
;OPLOT,FINDGEN(N)-N/2,SHIFT(ABS(FFT((bkguess))),N/2),COLOR=RED
;OPLOT,FINDGEN(N)-N/2,SHIFT(ABS(FFT(mur0*P)),N/2),COLOR=RED,LINESTYLE=2
;
;!P.MULTI=0
;RETURN
;END
