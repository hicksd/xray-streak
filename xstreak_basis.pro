@xmufrommass
@xstreak_fit
@xinvert

FUNCTION get_H,N,order
;Constructs a matrix H for linear regularization
;Order parameter, 1:u'=constant, 2:u"=constant etc.

IF N LE 6 THEN BEGIN
  PRINT,'H must be greater than 6x6 array'
  RETURN,0
END

B=DBLARR(N-order,N)
CASE order OF
1:BEGIN
  B(0,0)=-1 & B(0,1)=1
  FOR i=1,N-1-order DO B(i,*)=SHIFT(B(0,*),i)
END
2:BEGIN
  B(0,0)=-1 & B(0,1)=2 & B(0,2)=-1
  FOR i=1,N-1-order DO B(i,*)=SHIFT(B(0,*),i) 
END
3:BEGIN
  B(0,0)=-1 & B(0,1)=3 & B(0,2)=-3 & B(0,3)=1
  FOR i=1,N-1-order DO B(i,*)=SHIFT(B(0,*),i) 
END
ELSE:
ENDCASE

RETURN,MATRIX_MULTIPLY(B,B,/ATRANSPOSE)
END


FUNCTION funcbasis,k,x,orig=orig
COMMON bparams,del,zz3,zz4,zz5
;Basis function to be used in IMSL_FCNLSQ to return projected function, P_k
;If ORIG=1 then return value of original (density) function, f_k
IF SIZE(del,/TYPE) EQ 0 THEN del=5. ;If not set, then make 5 pixels the width of each basis
R0=(k-0.5)*del ;k=1 is at del/2.
;print,k,x ;Gaussian basis function seem to give less smooth result than square function
;and are MUCH slower (using the numerical method applied here).
;f_k=EXP(-(ABS(r)-R0)^2/2/(del/2)^2)   ;Gaussian density function
;IF k(0) EQ 1 AND x(0) EQ 0 THEN P_k=abeltransform(r,f_k,M,/init) ELSE P_k=abeltransform(r,f_k,M)
;P_k=INTERPOL(P_k,r,ABS(x))

IF NOT(KEYWORD_SET(orig)) THEN BEGIN
  P_k=DBLARR(N_ELEMENTS(x))
  w=WHERE(ABS(x) LE R0+del/2.) & $
    IF w(0) NE -1 THEN P_k(w)=2*SQRT((R0+del/2.)^2-x(w)^2)
  w2=WHERE(ABS(x) LE R0-del/2) & $
    IF w2(0) NE -1 THEN P_k(w2)=2*(SQRT((R0+del/2.)^2-x(w2)^2)-SQRT((R0-del/2.)^2-x(w2)^2) )
  RETURN,P_k
ENDIF ELSE BEGIN
  ;make sure not
  ;to double count by carefully setting GE and LT limits
  RETURN,DOUBLE((ABS(x)-ABS(R0)) GE -ABS(del)/2. AND (ABS(x)-ABS(R0)) LT ABS(del)/2.) ;f_k
ENDELSE

END



FUNCTION AbelBasis,x,y,ddel,gfit,ccor,eedge
COMMON bparams
;Uses the basis function approach to extract the best fit Abel inverse of y vs x
;x must be positive and start at zero.

dx=ABS(x(1)-x(0))
Nx=N_ELEMENTS(x) ;

del=ddel ;width of basis function
est=gfit  ;first guess for Gaussian
cor=ccor ;minimum radius for use in Gaussfit
edge=eedge ;fraction of data outside which is considered background (usually around 0.95)

y_end=MEAN(y(edge*Nx-1:Nx-1)) & y=y/y_end ;force transmission @ edge to be 1 (krho=0 there)
P=-ALOG(ABS(y))/1.E-4 ;Take natural log, multiply by -1; convert to cm from um.

nbasis=ROUND(MAX(x)/del)-1 

x=DOUBLE(x)
c_k = IMSL_FCNLSQ('funcbasis',nbasis,x,P,/DOUBLE,SSE=sig) 
tau =DBLARR(Nx)  & FOR k=1,nbasis DO tau =tau +c_k(k-1)*funcbasis(k,x) ;builds fitted tau
krho=DBLARR(Nx)  & FOR k=1,nbasis DO krho=krho+c_k(k-1)*funcbasis(k,x,/orig) ;builds fitted krho

;;interpolate between basis coefficients to recover radial density profile (for square basis functions only)
;;INTERPOL(c_k,FINDGEN(nbasis)*del,x) 
err=REPLICATE(1,Nx) 
core=WHERE(ABS(x) LT cor) & err(core)=1000. & err(edge*Nx:Nx-1)=1000
cfit=GAUSSFIT(x,krho,a,NTERMS=6,MEASURE_ERRORS=err,ESTIMATES=est)

krhob=krho ;krho before fitted background is subtracted
cfitb=cfit ;cfit before fitted background is subtracted
bkgfit=a(3)+a(4)*x+a(5)*x^2 ;isolates fitted background

krho=krho-bkgfit ;subtract background
cfit=cfit-bkgfit ;subtract background

krho(core)=0. & krho(edge*Nx-1:Nx-1)=0.

q={x:x,y:y,del:del,est:est,cor:cor,nbasis:nbasis,tau:tau,krho:krho,krhob:krhob,cfit:cfit,cfitb:cfitb,bkgfit:bkgfit,a:a}

RETURN,q
END

FUNCTION AbelBasis2,xx,yy,delp,filter=filter,est=est,Rmin=Rmin,edge=edge,q0=q0,q1=q1
;Takes both +ve and -ve x-data, does Abel inverse on each independently, recombines
;filter is width, in um, of filtering window
;delp is width, in PIXELS, of square basis function
;est is initial estimate for Gaussfit
;Rmin is minimum radius to use in Gaussfit

!EXCEPT=0 ;don't print floating point underflow errors

IF NOT(KEYWORD_SET(est)) THEN est=[10.,300.,25.,0,0,0]   ;first guess for Gaussian
IF NOT(KEYWORD_SET(Rmin)) THEN Rmin=100.

dx=ABS(xx(1)-xx(0)) ;get x-spacing
del=delp*dx ;width of basis function in um

IF filter GT 0 THEN BEGIN
  filt=ROUND((filter/dx-1)/2.) 
  IF filt GT 2 AND filt LT N_ELEMENTS(yy) THEN ys = LEEFILT(yy,filt,/DOUBLE) ELSE ys=yy
ENDIF ELSE ys=yy

;********** +ve limb
w0=WHERE(xx GE 0) 
q0=AbelBasis(xx(w0),ys(w0),del,est,Rmin,edge)
;********** -ve limb
w1=WHERE(xx LE 0) 
q1=AbelBasis(-REVERSE(xx(w1)),REVERSE(ys(w1)),del,est,Rmin,edge)
;********** Find average
N0=N_ELEMENTS(q0.x) & N1=N_ELEMENTS(q1.x) & N=MAX([N0,N1]) 
x=(y=(tau=(krho=FLTARR(N))))
IF N0 GT N1 THEN BEGIN
  x=q0.x
  y[0:N1-1]   =0.5*(q0.y[0:N1-1]   +q1.y[0:N1-1])    & y[N1:N0-1]   =q0.y[N1:N0-1]    ;doesn't make sense
  tau[0:N1-1] =0.5*(q0.tau[0:N1-1] +q1.tau[0:N1-1])  & tau[N1:N0-1] =q0.tau[N1:N0-1]  ;doesn't make sense
  krho[0:N1-1]=0.5*(q0.krho[0:N1-1]+q1.krho[0:N1-1]) & krho[N1:N0-1]=q0.krho[N1:N0-1]
END
IF N1 GT N0 THEN BEGIN
  x=q1.x
  y[0:N0-1]   =0.5*(q0.y[0:N0-1]   +q1.y[0:N0-1])    & y[N0:N1-1]   =q1.y[N0:N1-1] ;doesn't make sense
  tau[0:N0-1] =0.5*(q0.tau[0:N0-1] +q1.tau[0:N0-1])  & tau[N0:N1-1] =q1.tau[N0:N1-1] ;doesn't make sense
  krho[0:N0-1]=0.5*(q0.krho[0:N0-1]+q1.krho[0:N0-1]) & krho[N0:N1-1]=q1.krho[N0:N1-1] 
ENDIF
IF N1 EQ N0 THEN BEGIN
  x=q0.x
  y=0.5*(q0.y+q1.y) ;doesn't make sense
  tau=0.5*(q0.tau+q1.tau) ;doesn't make sens
  krho=0.5*(q0.krho+q1.krho)
ENDIF
nbasis=N
;x=0.5*(q0.x[0:N-1]+q1.x[0:N-1]) ;These should be the same but average them anyway
;y=0.5*(q0.y[0:N-1]+q1.y[0:N-1]) ;This step probably doesn't make sense
;nbasis=MIN([q0.nbasis,q1.nbasis])
;tau=0.5*(q0.tau[0:N-1]+q1.tau[0:N-1]) ;this step doesn't make sense because of background subtraction differences between limbs
;krho=0.5*(q0.krho[0:N-1]+q1.krho[0:N-1])

;est=0.5*(q0.a+q1.a) & est(3)=0. & est(4)=0. & est(5)=0.

;Redo the Gaussian fit on the average krho
err=REPLICATE(1,N) & core=WHERE(ABS(x) LT Rmin) & err(core)=1000. ;& err(edge*N:N-1)=1000
cfit=GAUSSFIT(x,krho,a,NTERMS=6,MEASURE_ERRORS=err,ESTIMATES=est)

krho=krho-(a(3)+a(4)*x+a(5)*x^2) ;subtract background
cfit=cfit-(a(3)+a(4)*x+a(5)*x^2) ;subtract background

a(0:2)=ABS(a(0:2)); make sure Gaussian is +ve
q={xx:xx,yy:yy,ys:ys,x:x,y:y,del:del,est:est,Rmin:Rmin,nbasis:nbasis,tau:tau,krho:krho,cfit:cfit,a:a,q0:q0,q1:q1}


RETURN,q
END

FUNCTION xstreak_basis,imstr,sub,tsub,psub,a,BLOCK=block,MIX=mix,Npoly=Npoly,SKEW=skew,FITMURHO=FITMURHO

IF KEYWORD_SET(MIX) THEN mix=mix ELSE mix=0. ;amount of atomic-scale mix in target
IF KEYWORD_SET(FITMURHO) THEN ffitmurho=1. ELSE ffitmurho=0.;if=1 then fit mu*rho only, get average(mu) after fitting.
IF KEYWORD_SET(NPOLY) THEN NNpoly=Npoly ELSE NNpoly=4 ;I0 polynomial order
IF KEYWORD_SET(SKEW) THEN sskew=skew ELSE sskew=0.
WIDGET_CONTROL,imstr.freetype_id,GET_VALUE=free ;Which parameters to vary [R0,rho,del] - 1=vary,0=fix
;WIDGET_CONTROL,imstr.slit_id,GET_VALUE=slit & slit=FLOAT(slit[0])

tg=imstr.targetdata

WIDGET_CONTROL,imstr.dpix_id,GET_VALUE=dpix & dpix=FLOAT(dpix[0])
WIDGET_CONTROL,imstr.Rmin_id,GET_VALUE=Rmin & Rmin=FLOAT(Rmin[0])
WIDGET_CONTROL,imstr.edge_id,GET_VALUE=edge & edge=FLOAT(edge[0])

;;filter=25. ; um of smoothing filter
;dpix=5. ;# of PIXELS for basis function
;Rmin=70. ;(um) innermost radius inside of which is considered polluted
;edge=0.95 ;fraction of data outside which is considered background

;a=[R0,maxrho,del,x0,slit,micpx,sl_y0,sl_dy,fbkg]
x0=a(3) & slit=ABS(a(4)) & micpx=ABS(a(5)) & sl_y0=a(6) & sl_dy=a(7) ;& fbkg=b(8)

sz=SIZE(sub) & snx=sz(1) & sny=sz(2)
filter=slit/2. ; um of smoothing filter

aa=DBLARR(N_ELEMENTS(a),snx) ;create array to store fitted values of 'a'
xm=(DINDGEN(sny)-DOUBLE(x0))*micpx ;create spatial vector
;OD=shell(sny,a,xm,rho,mu,M,func,TARGET=tg,init=init,FITMURHO=ffitmurho,SKEW=sskew) ;create the transformation matrix M

;********************
z=DBLARR(snx) & z2=DBLARR(snx,sny)  ;initializing array for structure
dat={target:tg,method:'Inverse',nx:snx,ny:sny,t:z,pix:z,chisq:z,$
      mix:mix,fitmurho:ffitmurho,Npoly:NNpoly,skew:sskew,$
      a:aa,xm:xm,im:REBIN(z2,snx,sny,3),$ ;'im' contains data for the running profile plots
      maxdensity:z,R0:z,del:z,x0:z,slit:z,micpx:z,sl_y0:z,sl_dy:z,$
      v:z,rhoR:z,mass:z,av_mu:z,$
      xrhoR:z,xmass:z,dxrhoR:z,dxmass:z,$ ;rhoR/mur0 and mass/mur0 with errors
      dmaxdensity:z,dR0:z,ddel:z,dv:z,drhoR:z,dmass:z,$ ;errors
      inten:z2,bkg:z2,OD:z2,density:z2,mu:z2}

!p.multi=[0,4,4]

est=[ABS(a(1)),ABS(a(0)),ABS(a(2))/2.,0,0,0] ;first guess of fit parameters

FOR i=0,snx-1 DO BEGIN
    x=xm
    y=REFORM(sub(i,*)) 
    
      ;***********Deconvolve
   IF slit/micpx GT 1 THEN BEGIN 
    Nx=N_ELEMENTS(x) & xpsf=FINDGEN(Nx)-Nx/2.+0.5*(Nx MOD 2) ;construct psf
;    ;    sigx=motion(j)/2.35/umperpix ; 
;    ;    psf=EXP(-xpsf^2/(2.*sigx^2)) & psf=psf/TOTAL(psf);/umperpix
    psf=FLTARR(Nx) & psf[WHERE(ABS(xpsf) LE slit/micpx/2. )]=1. & psf=psf/TOTAL(psf)
;    ;***************Max Entropy deconvolve image
    multipliers=FLTARR(Nx) 
    Niter=10 & FOR ii=1,Niter DO Max_Entropy,y,psf,yyd,multipliers,FT_PSF=psf_ft,/ONED
    y=yyd
  ENDIF
  
    q=AbelBasis2(x,y,dpix,filter=filter,est=est,Rmin=Rmin,edge=edge,q0=q0,q1=q1) 
    est=ABS(q.a) ;use best fit parameters as guess for next round

    ;a=[R0,eta,del,x0,slit,micpx,sl_y0,sl_dy]
    a(0)=ABS(est(1)) & a(1)=ABS(est(0)) & a(2)=ABS(est(2))*2.
    dat.a[*,i]=a
    dat.t[i]=tsub(i)  & dat.pix[i]=psub(i)
    dat.x0[i]=a(3)   & dat.slit[i]=a(4) & dat.micpx=a(5) 
    dat.sl_y0=a(6)   & dat.sl_dy=a(7) 
    dat.chisq[i]=0.;chisqr
    
    xx=q.x
    krho=q.cfit ;Should replace this with q.krho to reduce dependencies on fit function
    w=WHERE(xx LT MIN(xx) OR xx GT MAX(xx)) ;This is a placeholder for when limits are to be placed on the integral
    IF w(0) NE -1 THEN krho(w)=0.
    
    ;********************
    kmass=4.*!pi*INT_TABULATED(xx,krho*xx^2)/1.E12 ;g
    ;If max mu*mass is greater than value at the maximum mass then 
    ;set av_mu equal to that at the maximum mass to avoid extrapoloation problems.
    IF kmass*1.E6 GT MAX(tg.mumass) THEN av_mu=tg.av_mu[tg.N-1] $
    ELSE av_mu=INTERPOL(tg.av_mu,tg.mumass,kmass*1.E6) ;get av_mu based on observed mu*mass
    mass=kmass/av_mu
    ;********************
    
    IF KEYWORD_SET(fitmurho) THEN BEGIN ;calculate average mu after the fit
      mu=REPLICATE(av_mu,N_ELEMENTS(xx))
      rho=krho/av_mu ;here rho preserves the shape of krho (i.e. full mix case)
    ENDIF ELSE BEGIN
      mumass=TOTAL(krho,/CUMULATIVE)
      mu=INTERPOL(tg.mu,tg.mumass,mumass)
      rho=krho/mu
    ENDELSE
       
    Maxrho=MAX(rho) ;g/cc
    MeanR=INT_TABULATED(xx,rho*ABS(xx)*xx^2)/INT_TABULATED(xx,rho*xx^2);center of mass
    ;MeanR=INT_TABULATED(xx,rho*ABS(xx))/INT_TABULATED(xx,rho) ;density centroid, um
    TwoSigma=2.*SQRT(INT_TABULATED(xx,rho*(ABS(xx)-MeanR)^2)/INT_TABULATED(xx,rho)) ;um
    rhoR=INT_TABULATED(xx,rho)/1.E4 ;g/cm2
    
    dat.R0[i]=ABS(meanR)
    dat.maxdensity[i]=maxrho 
    dat.del[i]  =twosigma ;This allows a(0),a(1) etc. to be used for parameters other than maxeta, R0, etc.
    dat.rhoR[i] =rhoR*1.E3 ;mg/cm^2 
    dat.mass[i] =mass*1.E6 ;ug
    dat.av_mu[i]=av_mu ;cm2/g


;    dat.inten(i,*)= EXP(-OD);fitted x-ray intensity (normalized)
;    dat.bkg(i,*)= EXP(I0guess)            ;fitted background intensity
;    dat.OD(i,*)= OD                            ;fitted optical depth
;    dat.density(i,*) = rho    ;density
;    dat.mu(i,*) = mu    ;attenuation coefficient, cm2/g

    OUTPLOT,q.x,q.krho,imstr,/newplot,$
      XTITL="Radius (mic)",YTITL="kappa*rho (/cm)",$
      TITL="t = "+STRCOMPRESS(STRING(tsub(i),FORMAT='(F8.1)'),/REMOVE_ALL)+" ps"
    ;OUTPLOT,xm,EXP(I0guess),imstr,/overplot,color=[255,0,255]
    OUTPLOT,q.x,q.cfit,imstr,/overplot,color=[0,0,255]
    OUTPLOT,[MIN(q.x),MAX(q.x)],[0,0],imstr,/overplot,color=[0,0,0]

    IF SIZE(qq,/TYPE) EQ 0 THEN qq=q ELSE qq=[qq,q]
   
;    plot,q.x,q.krho
;    oplot,[min(q.x),max(q.x)],[0,0]
;    
    
;    yfit=EXP(I0fit-OD)
;    OUTPLOT,[-dat.R0[i],dat.R0[i]],INTERPOL(yfit,xm,[-dat.R0[i],dat.R0[i]]),$
;      imstr,/overplot,COLOR=[255,0,0],PSY=24
;    OUTPLOT,xm,rho/MAX(rho)*0.2,imstr,/overplot,COLOR=[0,255,0]
    
ENDFOR    

dat.v=velocity(dat.t,dat.R0,dat.t,smth=0)*1.E0 ;velocity from t vs R0 fit, convert to mic/ns
  
dat=imagetostructure(dat,'QQ',qq)

!p.multi=0
RETURN,dat
END
        
PRO testbasis,s,eps=eps,qq=qq,avK=avK
;tests a new inversion scheme using decomposition into basis functions
;COMMON bparams

RED=255 & GREEN=255*256L & BLUE=255*256L*256L & YELLOW=RED+GREEN & PINK=RED+BLUE & LBLUE=GREEN+BLUE
ORANGE=255+126*256L & PURPLE=126+BLUE & PLAIN=RED+GREEN+BLUE & BACKG=0 & BROWN=174+71*256L
charsz=1.2 & thk=1.5 & symsz=1.3

!P.FONT=-1
IF KEYWORD_SET(eps) OR KEYWORD_SET(ps) THEN BEGIN
!P.FONT=1
    SET_PLOT,'ps'
;    !X.MARGIN=[5.5,1.5] & !Y.MARGIN=[3.,1.]    
    !X.MARGIN=[5.5,1.5] & !Y.MARGIN=[3.,1.]    
    DEVICE,ENCAPSULATED=1,LANDSCAPE=0,XSIZE=17,YSIZE=12,/COLOR,/HELVETICA,FILENAME="C:\Documents and Settings\hicks13\DESKTOP\IDL.EPS"
    TVLCT, [0,255,0,0,255,255,0,255,126,255,174], [0,0,255,0,255,0,255,126,0,255,71], [0,0,0,255,0,255,255,0,255,255,0]
    PLAIN=0 & RED=1 & GREEN=2 & BLUE=3 & YELLOW=4 & PINK=5 & LBLUE=6 & ORANGE=7 & PURPLE=8 & BACKG=9 & BROWN=10
    FILL=1 & charsz=1.3 & thk=3.0*thk & symsz=0.7*symsz
ENDIF
clr=[PLAIN,RED,GREEN,BLUE,BROWN,ORANGE,PINK,LBLUE,PURPLE,YELLOW]

;***********Load streak profile
;dir='c:\hicks\nif\convabl\results\110630\invert\' & g=[500.,490.,460.,430,410,390,315,290,270,220,200,190]
;files=dir+['strip1a.txt','strip1b.txt','strip1c.txt','strip2a.txt','strip2b.txt','strip2c.txt',$
;  'strip3a.txt','strip3b.txt','strip3c.txt','strip4a.txt','strip4b.txt','strip4c.txt'] 
;dir='c:\hicks\nif\convabl\results\111011\invert\'
;files=dir+['strip1b.txt','strip2b.txt','strip3b.txt','strip4b.txt'] & g=[430.,390.,322.,230]
;dir='c:\hicks\nif\convabl\results\111009\invert\'
;files=dir+['strip1c.txt','strip2b.txt','strip3b.txt','strip4b.txt'] & g=[430.,410.,330.,220]

;dir='c:\hicks\nif\convabl\results\110630\' & g=[500.,490.,460.,420,400,390,360,300,260,220,200,190]
;file=dir+['N110630_gstreak3.sav'] & umperpix=1.027 & d0=579 & tbin=500 & xmin=0 & xmax=1150
;t=[20.8860,20.9600,21.0340,21.2110,21.2850,21.3590,21.6390,21.7130,21.7870,21.9550,22.0290,22.1030]
;TgName='N110630-CHSi'

dir='c:\hicks\nif\convabl\results\111007\' & g=[865.,813.,736.,623.]
file=dir+'N111007_gstreak2.sav' & umperpix=1.0321 & d0=1195 & tbin=500 & xmin=50 & xmax=2349 
t=[18.54,19.115,19.719,20.327]
Rlim=[500.,1000] & rholim=[-2,6] 
TgName='N111007-CHSi'

;dir='c:\hicks\nif\convabl\results\111009\' & g=[488.,477.,460.,420,410,390,350,330,310,270,250,220]
;file=dir+'N111009_gstreak1.sav' & umperpix=1.0296 & d0=579 & tbin=500 & xmin=10 & xmax=960
;t=[20.871,20.95,21.029,21.156,21.235,21.314,21.455,21.534,21.613,21.759,21.838,21.917]
;TgName='N111009-CHSi'

;dir='c:\hicks\nif\convabl\results\111011\' & g=[500.,490.,470.,430,410,390,350,330,310,270,250,220]
;file=dir+'N111011_gstreak1.sav' & umperpix=1.0237 & d0=615 & tbin=400 & xmin=85 & xmax=1168
;t=[20.88,20.96,21.04,21.173,21.252,21.332,21.473,21.552,21.632,21.776,21.855,21.935]
;TgName='N111011-CHSi'

filter=0. ; um of smoothing filter
del=2. ;# of PIXELS for basis function
Rmin=120. ;(um) innermost radius inside which is considered polluted
edge=0.95 ;fraction of data outside which is considered background


;***********Load target info, if available
IF SIZE(TgName[0],/TYPE) NE 0 THEN BEGIN
  TgStrFile='c:\hicks\xstreak\xstreak_target.sav' & xrayfile='c:\hicks\xstreak\xKappa.sav'
  RESTORE,file=TgStrfile
  w=WHERE(TgStr.name EQ TgName)
  tg=loadtarget(XKAPPAFILE=xrayfile,TDAT=TgStr[w])
ENDIF 
;************

RESTORE,file & im=image & sz=SIZE(im) & nx=sz(2) & nt=sz(1)/tbin

z=DBLARR(nt)
s={t:t,maxdensity:z,R0:z,Rmass:z,Rin:z,Rabl:z,del:z,mass:z,fmass:z}


;!P.MULTI=[0,3,4]
!P.MULTI=[0,2,2]
;!P.MULTI=0

;FOR j=0,N_ELEMENTS(files)-1 DO BEGIN
  ;dat=readdat(files(j),4) ;dat=readxvisdat(file,x,y,ysig)
  ;xx=REFORM(dat(0,*)) & yy=REFORM(dat(1,*)) & N=N_ELEMENTS(yy)

FOR j=0,nt-1 DO BEGIN
  xx=(FINDGEN(nx)-d0)*umperpix &   xx=xx(xmin:xmax)
  yy=REFORM(REBIN(im(j*tbin:(j+1)*tbin-1,xmin:xmax),1,xmax-xmin+1)) ;both limbs

  ;***********Deconvolve
;  Nyy=N_ELEMENTS(yy) & xpsf=FINDGEN(Nyy)-Nyy/2.+0.5*(Nyy MOD 2) ;construct psf
;  sigx=29./2.35/umperpix ; 
;  psf=EXP(-xpsf^2/(2.*sigx^2)) & psf=psf/TOTAL(psf);/umperpix
;      ;***************Max Entropy deconvolve image
;  multipliers=FLTARR(Nyy) 
;  Niter=5 & FOR ii=1,Niter DO Max_Entropy,yy,psf,yyd,multipliers,FT_PSF=psf_ft,/ONED
;  yy=yyd
  ;*******************
  
  est=[10,g(j),25.,0,0,0]
  q=AbelBasis2(xx,yy,del,filter=filter,est=est,Rmin=Rmin,edge=edge,q0=q0,q1=q1)
  IF j EQ 0 THEN qq=q ELSE qq=[qq,q]
  ;*******************
  IF SIZE(TgName[0],/TYPE) NE 0 THEN BEGIN
    kmass=4.*!pi*INT_TABULATED(q.x,q.cfit*q.x^2)/1.E12 ;g
    ;If max mu*mass is greater than value at the maximum mass then 
    ;set av_mu equal to that at the maximum mass to avoid extrapoloation problems.
    IF kmass*1.E6 GT MAX(tg.mumass) THEN av_mu=tg.av_mu[tg.N-1] $
    ELSE av_mu=INTERPOL(tg.av_mu,tg.mumass,kmass*1.E6) ;get av_mu based on observed mu*mass
    mass=kmass/av_mu
    ;********************
    IF j EQ 0 THEN avK=av_mu ELSE avK=[avK,av_mu]
    ;********************
    q.krho=q.krho/av_mu
    q.a[0]=q.a[0]/av_mu
    q.cfit=q.cfit/av_mu
    q0.krho=q0.krho/av_mu
    q1.krho=q1.krho/av_mu
    
    ;ymin=(nt-1-j)*drho & ymax=(nt-j)*drho
    ;ymin=rholim[0] & ymax=rholim[1]
    
    PLOT,[Rlim[0],Rlim[1]],[rholim[0],rholim[1]],/NODATA,/YSTYLE,$
    CHARSIZE=charsz,THICK=thk,XTHICK=thk,YTHICK=thk,CHARTHICK=thk,$;YRANGE=[-1,MAX(eta)],
    XTITLE="Radius (!Mmm)",YTITLE="Density (g/cc)" 
    pp1=!p & px1=!x & py1=!y; capture system variables for this plot
    
  ENDIF ELSE BEGIN
    PLOT,[Rlim[0],Rlim[1]],[rholim[0],rholim[1]],/NODATA,$
    CHARSIZE=charsz,THICK=thk,XTHICK=thk,YTHICK=thk,CHARTHICK=thk,$;YRANGE=[-1,MAX(eta)],
    XTITLE="Radius (!Mmm)",YTITLE="!Mk!d!Mn!n !Mr (cm!u-1!n)"
    pp1=!p & px1=!x & py1=!y; capture system variables for this plot  
  ENDELSE
  XYOUTS,0.02*(Rlim[1]-Rlim[0])+Rlim[0],rholim[0]+0.85*(rholim(1)-rholim(0)),$
    STRCOMPRESS(STRING(t(j),FORMAT='(F7.3)'))+' ns',CHARSIZE=charsz*0.5,CHARTHICK=thk
  
  s.Rmass[j]=INT_TABULATED(q.x,q.cfit*ABS(q.x)*q.x^2)/INT_TABULATED(q.x,q.cfit*q.x^2);center of mass
  s.maxdensity[j]=q.a[0]
  s.R0[j]=q.a[1]
  s.del[j]=2.*q.a[2]
  s.Rin[j]=q.a[1]-q.a[2]*2.
  s.Rabl[j]=q.a[1]+q.a[2]*2.

  OPLOT,q.x,q.krho,THICK=thk
  OPLOT,[0,MAX(q.x)],[0,0],THICK=thk,LINESTYLE=1
  OPLOT,q0.x,q0.krho,THICK=thk,COLOR=RED,LINESTYLE=1
  OPLOT,q1.x,q1.krho,THICK=thk,COLOR=BLUE,LINESTYLE=1
  OPLOT,q.x,q.cfit,THICK=thk,LINESTYLE=2
  
  OPLOT,[s.Rmass[j],s.Rmass[j]],[rholim[0],rholim[1]],THICK=thk,LINESTYLE=0
  OPLOT,[s.Rin[j],s.Rin[j]],[rholim[0],rholim[1]],THICK=thk,LINESTYLE=1
  OPLOT,[s.Rabl[j],s.Rabl[j]],[rholim[0],rholim[1]],THICK=thk,LINESTYLE=1
    
jump:
ENDFOR

pm,[[s.R0],[s.del],[s.Rin],[s.Rabl]]

;PLOT,xx,yy,CHARSIZE=charsz,THICK=thk,XTHICK=thk,YTHICK=thk,CHARTHICK=thk,$
;    XTITLE="Position (!Mmm)",YTITLE="Intensity (arb units)"
;OPLOT,q.xx,q.ys,THICK=thk,LINESTYLE=1

;w=WHERE(xm GE 0) & x=xm(w) & y=ys(w)
;;w=WHERE(xm LE 0) & x=-REVERSE(xm(w)) & y=REVERSE(ys(w))
;r=x ;define radial positions for funcbasis().
;Nx=N_ELEMENTS(x) ;
;
;y_end=MEAN(y(Nx-0.05*Nx-1:Nx-1)) & y=y/y_end
;
;P=-ALOG(y)/mur0 ;Take natural log, multiply by -1.
;
;del=5.*dx 
;nbasis=ROUND(MAX(x)/del)-1 & c_k = IMSL_FCNLSQ('funcbasis',nbasis,x,P,/DOUBLE,SSE=sig) 
;func=DBLARR(Nx) & FOR k=1,nbasis DO func=func+c_k(k-1)*funcbasis(k,x) ;builds fitted P
;krho=DBLARR(Nx)  & FOR k=1,nbasis DO krho=krho+c_k(k-1)*funcbasis(k,x,/orig) ;builds fitted eta
;
;;;interpolate between basis coefficients to recover radial density profile (for square basis functions only)
;;;INTERPOL(c_k,FINDGEN(nbasis)*del,x) 
;
;;m=MAX(eta,wmax) & print,wmax
;est=[10,g,25.,0,0,0]
;err=REPLICATE(1,N_ELEMENTS(x)) & core=WHERE(ABS(x) LT 100.) & err(core)=1000.
;cfit=GAUSSFIT(x,krho,a,NTERMS=6,MEASURE_ERRORS=err,ESTIMATES=est)
;print,a
;;dx=ABS(x(1)-x(0)) & print,dx
;;nf=2.*ABS(a(2))/dx ;# points on left and right to fit in Savitzky-Golay filter
;;nf=ROUND(nf)
;;savgolFilter = SAVGOL(nf, nf, 0, 4) 
;;etaFilt=CONVOL(eta, savgolFilter, /EDGE_TRUNCATE)
;
;
;;PLOT,x,func,CHARSIZE=charsz,THICK=thk,XTHICK=thk,YTHICK=thk,CHARTHICK=thk,$
;;    XTITLE="Position (!Mmm)",YTITLE="Projected depth (!Mmm)"
;;pp0=!p & px0=!x & py0=!y; capture system variables for this plot
;;OPLOT,x,P,COLOR=RED,THICK=thk,LINESTYLE=2


;w=WHERE(xm GE 0) & x=xm(w) & y=ys(w)
;krho=AbelBasis(x,y,ddel=del,gfit=est,ccor=Rmin,q=q)


;PLOT,q.x,q.krho,XRANGE=[0,MAX(xx)],$;YRANGE=[-40,40],$
;  CHARSIZE=charsz,THICK=thk,XTHICK=thk,YTHICK=thk,CHARTHICK=thk,$;YRANGE=[-1,MAX(eta)],
;  XTITLE="Radius (!Mmm)",YTITLE="!Mk!d!Mn!n !Mr (cm!u-1!n)"
;OPLOT,[0,1000],[0,0],THICK=thk
;pp1=!p & px1=!x & py1=!y; capture system variables for this plot
;OPLOT,q0.x,q0.krho,THICK=thk,COLOR=RED,LINESTYLE=1
;OPLOT,q1.x,q1.krho,THICK=thk,COLOR=BLUE,LINESTYLE=1

;w=WHERE(xm LE 0) & x=-REVERSE(xm(w)) & y=REVERSE(ys(w))
;krho=AbelBasis(x,y,ddel=del,gfit=est,ccor=Rmin,q=q)
;OPLOT,x,krho,THICK=thk,COLOR=BLUE
;OPLOT,x,q.cfit,THICK=thk,COLOR=RED,LINESTYLE=2



;;w=WHERE(xm GE 0) & x=xm(w) & y=ys(w)
;w=WHERE(xm LE 0) & x=-REVERSE(xm(w)) & y=REVERSE(ys(w))
;r=x ;define radial positions for funcbasis().
;
;Nx=N_ELEMENTS(x) ;
;y_end=MEAN(y(Nx-0.05*Nx-1:Nx-1)) & y=y/y_end
;P=-ALOG(y)/mur0 ;Take natural log, multiply by -1.
;
;
;del=5*dx & nbasis=ROUND(MAX(x)/del)-1 & c_k = IMSL_FCNLSQ('funcbasis',nbasis,x,P,/DOUBLE,SSE=sig) 
;func=DBLARR(Nx) & FOR k=1,nbasis DO func=func+c_k(k-1)*funcbasis(k,x) ;builds fitted P
;eta=DBLARR(Nx)  & FOR k=1,nbasis DO eta=eta+c_k(k-1)*funcbasis(k,x,/orig) ;builds fitted eta
;
;;;interpolate between basis coefficients to recover radial density profile (for square basis functions only)
;;;INTERPOL(c_k,FINDGEN(nbasis)*del,x) 
;
;;m=MAX(eta,wmax) & print,wmax
;;est=[2.1585516,210.07641,9.1109823,-0.33608420,0.00011869420,6.3929830e-006]
;;est=[10,x(wmax),25.,0,0,0]
;est=[10,g,25.,0,0,0]
;err=REPLICATE(1,N_ELEMENTS(x)) & core=WHERE(ABS(x) LT 100.) & err(core)=1000.
;cfit=GAUSSFIT(x,eta,a,NTERMS=6,MEASURE_ERRORS=err,ESTIMATES=est)
;print,a
;
;
;x=-REVERSE(x) & eta=REVERSE(eta) & cfit=REVERSE(cfit)
;OPLOT,x,eta,THICK=thk
;OPLOT,x,cfit,THICK=thk,COLOR=RED,LINESTYLE=2
;
;
;
;
;;w=WHERE(x LT 70.)
;;print,mean(eta(w)),stddev(eta(w))
;;OPLOT,x,f,COLOR=RED,LINESTYLE=2
;;
;;;Try other methods for Abel inversion
;;;Normal equations 
;;junk=abeltransform(x,P,M,/init) ;Used to generate M matrix only
;;;;;MM=2.*dx*M & f=LA_LEAST_SQUARES(TRANSPOSE(MM)#MM, TRANSPOSE(MM)#SMOOTH(P,1),/DOUBLE) 
;;;;; Compute the LU decomposition:
;;lambda=1.E6
;;H=get_H(Nx,2)
;;MM=2.*dx*M
;;Aludc=MATRIX_MULTIPLY(MM,MM,/ATRANSPOSE)+lambda*H 
;;print,TRACE(MATRIX_MULTIPLY(MM,MM,/ATRANSPOSE))/TRACE(H)
;;LA_LUDC,Aludc,index 
;;b=MATRIX_MULTIPLY(MM,P,/ATRANSPOSE) & f=LA_LUSOL(Aludc, index, b,/DOUBLE) 
;;OPLOT,x,f,COLOR=PINK,THICK=thk
;;!p=pp0 & !x=px0 & !y=py0
;;OPLOT,x,abeltransform(x,f,M),COLOR=PINK,THICK=thk,LINESTYLE=1
;;!p=pp1 & !x=px1 & !y=py1
;;
;;;Filtered backprojection
;;P=INTERPOL(P,x,x+dx/2.) ;extract values of P from x+dx/2. in prep. for creating a profile evenly symmetric about x=0
;;x=INTERPOL(x,x,x+dx/2.) ;make sure x starts at dx/2. in prep. for creating a profile evenly symmetric about x=0
;;x=[-REVERSE(x),x] & P=[REVERSE(P),P]
;;
;;q=fbp(x,P)
;;OPLOT,x,q,THICK=thk,COLOR=BLUE,LINESTYLE=1


!P.MULTI=0
IF KEYWORD_SET(eps) OR KEYWORD_SET(ps)  THEN BEGIN
    DEVICE,/CLOSE
    SET_PLOT,'WIN'
ENDIF
!P.FONT=-1
!X.MARGIN=[10,3] & !Y.MARGIN=[4,2] 


RETURN
END


PRO testbasis2,s,dat,eps=eps,qq=qq,avK=avK,density=density,peakonly=peakonly
;tests a new inversion scheme using decomposition into basis functions
;Different from testbasis in that it overplots results on one graph.
;"Peakonly" forces a plot only of the peak region
;COMMON bparams

RED=255 & GREEN=255*256L & BLUE=255*256L*256L & YELLOW=RED+GREEN & PINK=RED+BLUE & LBLUE=GREEN+BLUE
ORANGE=255+126*256L & PURPLE=126+BLUE & PLAIN=RED+GREEN+BLUE & BACKG=0 & BROWN=174+71*256L
charsz=1.2 & thk=1.5 & symsz=1.3

!P.FONT=-1
IF KEYWORD_SET(eps) OR KEYWORD_SET(ps) THEN BEGIN
!P.FONT=1
    SET_PLOT,'ps'
;    !X.MARGIN=[5.5,1.5] & !Y.MARGIN=[3.,1.]    
    !X.MARGIN=[5.5,1.5] & !Y.MARGIN=[3.,2.]    
    DEVICE,ENCAPSULATED=1,LANDSCAPE=0,XSIZE=11,YSIZE=9,/COLOR,/HELVETICA,FILENAME="C:\Documents and Settings\hicks13\DESKTOP\IDL.EPS"
    TVLCT, [0,255,0,0,255,255,0,255,126,255,174], [0,0,255,0,255,0,255,126,0,255,71], [0,0,0,255,0,255,255,0,255,255,0]
    PLAIN=0 & RED=1 & GREEN=2 & BLUE=3 & YELLOW=4 & PINK=5 & LBLUE=6 & ORANGE=7 & PURPLE=8 & BACKG=9 & BROWN=10
    FILL=1 & charsz=1.2 & thk=3.0*thk & symsz=0.7*symsz
ENDIF
clr=[PLAIN,RED,GREEN,BLUE,BROWN,ORANGE,PINK,LBLUE,PURPLE,YELLOW]

s=0 & dat=0
;***********Load streak profile
;dir='c:\hicks\nif\convabl\results\110630\invert\' & g=[500.,490.,460.,430,410,390,315,290,270,220,200,190]
;files=dir+['strip1a.txt','strip1b.txt','strip1c.txt','strip2a.txt','strip2b.txt','strip2c.txt',$
;  'strip3a.txt','strip3b.txt','strip3c.txt','strip4a.txt','strip4b.txt','strip4c.txt'] 
;dir='c:\hicks\nif\convabl\results\111011\invert\'
;files=dir+['strip1b.txt','strip2b.txt','strip3b.txt','strip4b.txt'] & g=[430.,390.,322.,230]
;dir='c:\hicks\nif\convabl\results\111009\invert\'
;files=dir+['strip1c.txt','strip2b.txt','strip3b.txt','strip4b.txt'] & g=[430.,410.,330.,220]

;dir='c:\hicks\nif\convabl\results\110630\' & g=[500.,490.,460.,420,400,390,360,300,260,220,200,190]
;file=dir+['N110630_gstreak3.sav'] & umperpix=1.027 & d0=579 & tbin=500 & xmin=0 & xmax=1150
;t=[20.8860,20.9600,21.0340,21.2110,21.2850,21.3590,21.6390,21.7130,21.7870,21.9550,22.0290,22.1030]
;TgName='N110630-CHSi'

;dir='c:\hicks\nif\convabl\results\111009\' & g=[488.,477.,460.,420,410,390,350,330,310,270,250,220]
;file=dir+'N111009_gstreak1.sav' & umperpix=1.0296 & d0=579 & tbin=500 & xmin=10 & xmax=960
;t=[20.871,20.95,21.029,21.156,21.235,21.314,21.455,21.534,21.613,21.759,21.838,21.917]
;TgName='N111009-CHSi'

;dir='c:\hicks\nif\convabl\results\111011\' & g=[500.,490.,470.,430,410,390,350,330,310,270,250,220]
;file=dir+'N111011_gstreak1.sav' & umperpix=1.0237 & d0=615 & tbin=400 & xmin=85 & xmax=1168
;t=[20.88,20.96,21.04,21.173,21.252,21.332,21.473,21.552,21.632,21.776,21.855,21.935]
;TgName='N111011-CHSi'

;dir='c:\hicks\nif\convabl\results\111007\' & g=[865.,813.,736.,623.]
;file=dir+'N111007_gstreak2.sav' & umperpix=1.0321 & d0=1199 & tbin=500 & xmin=50 & xmax=2349 
;t=[18.54,19.115,19.719,20.327]
;Rlim=[500.,1000] & rholim=[-5,30] 
;TgName='N111007-CHSi' & Shot='N111007'

;dir='c:\hicks\nif\convabl\results\111220-2\' & g=[970.,960.,950.,940.]
;file=dir+'N111220-2_gstreak1.sav' & umperpix=1.0467 & d0=1201 & tbin=500 & xmin=0 & xmax=2399 
;t=13.78+[0.,0.8,1.6,2.4]
;Rlim=[850.,1050] & rholim=[-5,20] 
;motion=[2.,2,2,4] ;um smoothing filter
;TgName='N111220-02-CHSi' & Shot='N111220-2'

;dir='c:\hicks\nif\convabl\results\111218\' & g=[933.,912.,875.,809.]
;file=dir+'N111218_gstreak1.sav' & umperpix=1.0485 & d0=1201 & tbin=500 & xmin=50 & xmax=2349 
;t=16.24+[0.,0.8,1.6,2.4]
;Rlim=[700.,1000] & rholim=[-5,35] 
;motion=[2.,4,7,11] ;um smoothing filter
;TgName='N111218-CHSi' & Shot='N111218'

;dir='c:\hicks\nif\convabl\results\111219\' & g=[820.,730.,580.,380.]
;file=dir+'N111219_gstreak1.sav' & umperpix=1.0355 & d0=1199 & tbin=500 & xmin=50 & xmax=2349 
;t=18.61+[0.,0.8,1.6,2.4]
;Rlim=[200.,1000] & rholim=[-5,35]
;motion=[10.,19,26.,31.]
;TgName='N111219-CHSi' & Shot='N111219'

dir='c:\hicks\nif\convabl\results\111220-1\' & g=[820.,730.,580.,380.]
file=dir+'N111220-1_gstreak1.sav' & umperpix=1.0056 & d0=1197 & tbin=500 & xmin=50 & xmax=2349 
t=18.63+[0.,0.8,1.6,2.4]
Rlim=[200.,1000] & rholim=[-5,35]
motion=[12.,18,24.,32.]
TgName='N111220-01-CHSi' & Shot='N111220-1'


;dir='c:\hicks\nif\convabl\results\120306\' & g=[475.835,405.200,334.445,250.765]
;file=dir+'N120306_gstreak1.sav' & umperpix=0.978 & d0=604 & tbin=500 & xmin=10 & xmax=1190
;t=[21.4100,21.6870,21.9680,22.2740]
;Rlim=[50.,600] & rholim=[-50,100]
;motion=[29.,29,29.,29.]
;TgName='N120306-CHSi' & Shot='N120306'

;dir='c:\hicks\nif\convabl\results\120324\'
;g=[485.365,463.135,444.535,424.250,398.580,378.255,357.650,328.690,307.870,286.080,263.145,236.885,212.255]
;file=dir+'N120324_Deconv72x63y_ff&wirebkg.raw' & umperpix=0.6114 & d0=2212 & tbin=140 & xmin=1110 & xmax=3125 & tmin=990 & tmax=2879
;t=[21.3140,21.4070,21.4990,21.5905,21.6820,21.7730,21.8660,21.9595,22.0540,22.1490,22.2450,22.3420,22.4380]
;Rlim=[50.,600] & rholim=[-50,100]
;motion=REPLICATE(25.,N_ELEMENTS(t))
;TgName='N120324-CHSi' & Shot='N120324'

;dir='c:\hicks\nif\convabl\results\120329\' & g=[356.,286.,211.,137.]
;file=dir+'N120329_gstreak1.sav' & umperpix=0.6063 & d0=959 & tbin=500 & xmin=90 & xmax=1800
;t=[21.7460,22.0230,22.3040,22.6100]
;Rlim=[50.,500] & rholim=[-50,100]
;motion=[29.,29,29.,29.]
;TgName='N120329-CHSiTHD' & Shot='N120329'

;dir='c:\hicks\nif\convabl\results\120408\' 
;g=[503.870,484.065,462.495,439.615,416.250,393.505,371.725,347.595,321.365,297.965,274.695,249.915,227.625,201.840]
;file=dir+'N120408_Deconv72x38y_wirebkg.raw' & umperpix=0.5955 & d0=2023 & tbin=140 & xmin=1065 & xmax=2980 & tmin=935 & tmax=2964
;t=[21.2775,21.3705,21.4630,21.5550,21.6460,21.7370,21.8295,21.9225,22.0165,22.1120,22.2080,22.3040,22.4000,22.4970]
;Rlim=[50.,600] & rholim=[-50,100]
;motion=REPLICATE(25.,N_ELEMENTS(t))
;TgName='N120408-CHSi' & Shot='N120408'

;dir='c:\hicks\nif\convabl\results\120409\' 
;g=[532.600,506.975,486.265,459.645,439.900,414.345,386.125,360.445,341.580,310.635,283.850,258.175,234.425,205.570]
;file=dir+'N120409_Deconv72x37y_ff_bkg1000.raw' & umperpix=0.5969 & d0=2190 & tbin=140 & xmin=1195 & xmax=3188 & tmin=565 & tmax=2524
;t=[21.0310,21.1240,21.2175,21.3105,21.4035,21.4960,21.5875,21.6785,21.7700,21.8625,21.9555,22.0505,22.1460,22.2420]
;Rlim=[50.,600] & rholim=[-50,100]
;motion=REPLICATE(25.,N_ELEMENTS(t))
;TgName='N120409-CHSi' & Shot='N120409'

;dir='c:\hicks\nif\convabl\results\120421\' 
;g=[478.175,459.865,438.330,413.585,395.530,373.410,348.115,329.550,305.510,280.825,253.515,231.325,210.600,191.245,169.725]
;file=dir+'N120421_Deconv70x22y_wirebkg.raw' & umperpix=0.5923 & d0=2175 & tbin=140 & xmin=1010 & xmax=3149 & tmin=400 & tmax=2569
;t=[20.9210,21.0145,21.1075,21.2005,21.2940,21.3870,21.4800,21.5710,21.6620,21.7540,21.8460,21.9395,22.0335,22.1290,22.2250]
;Rlim=[50.,600] & rholim=[-50,200]
;motion=REPLICATE(25.,N_ELEMENTS(t))
;TgName='N120421-CHSi-Uni' & Shot='N120421'

filter=5. ; um of smoothing filter
del=2. ;# of PIXELS for basis function
Rmin=90. ;(um) innermost radius inside which is considered polluted
edge=0.95 ;fraction of data outside which is considered background


;***********Load target info, if available
IF SIZE(TgName,/TYPE) NE 0 THEN BEGIN
  TgStrFile='c:\hicks\xstreak\xstreak_target.sav' & xrayfile='c:\hicks\xstreak\xKappa.sav'
  RESTORE,file=TgStrfile
  w=WHERE(TgStr.name EQ TgName)
  tg=loadtarget(XKAPPAFILE=xrayfile,TDAT=TgStr[w])
ENDIF 
;************

IF STRPOS(file,'.sav') NE -1 THEN RESTORE,file
IF STRPOS(file,'.raw') NE -1 THEN image=READ_RAW(file)

im=image & sz=SIZE(im) & nx=sz(2) & nt=sz(1)/tbin

IF SIZE(tmin,/TYPE) EQ 0 THEN tmin=0.
IF SIZE(tmax,/TYPE) EQ 0 THEN nt=sz(1)/tbin ELSE nt=(tmax-tmin)/tbin

z=DBLARR(nt)
s={t:t,maxdensity:z,R0:z,Rmass:z,Rin:z,Rabl:z,del:z,mass:z,fmass:z}

;!P.MULTI=[0,3,4]
;!P.MULTI=[0,2,2]
;!P.MULTI=0

;FOR j=0,N_ELEMENTS(files)-1 DO BEGIN
  ;dat=readdat(files(j),4) ;dat=readxvisdat(file,x,y,ysig)
  ;xx=REFORM(dat(0,*)) & yy=REFORM(dat(1,*)) & N=N_ELEMENTS(yy)

IF SIZE(TgName,/TYPE) NE 0 AND KEYWORD_SET(density) THEN BEGIN
    PLOT,[Rlim[0],Rlim[1]],[rholim[0],rholim[1]],/NODATA,/YSTYLE,$
    CHARSIZE=charsz,THICK=thk,XTHICK=thk,YTHICK=thk,CHARTHICK=thk,$;YRANGE=[-1,MAX(eta)],
    XTITLE="Radius (!Mmm)",YTITLE="Density (g/cc)",TITLE=Shot
ENDIF ELSE BEGIN
    PLOT,[Rlim[0],Rlim[1]],[rholim[0],rholim[1]],/NODATA,$
    CHARSIZE=charsz,THICK=thk,XTHICK=thk,YTHICK=thk,CHARTHICK=thk,$;YRANGE=[-1,MAX(eta)],
    XTITLE="Radius (!Mmm)",YTITLE="!Mk!d!Mn!n !Mr (cm!u-1!n)",TITLE=Shot
ENDELSE

FOR j=0,nt-1 DO BEGIN
;FOR j=2,2 DO BEGIN
  xx=(FINDGEN(nx)-d0)*umperpix &   xx=xx(xmin:xmax)
  yy=REFORM(REBIN(im(tmin+j*tbin:tmin+(j+1)*tbin-1,xmin:xmax),1,xmax-xmin+1)) ;both limbs

  ;***********Deconvolve
  IF motion(j) GT 0 THEN BEGIN
    Nyy=N_ELEMENTS(yy) & xpsf=FINDGEN(Nyy)-Nyy/2.+0.5*(Nyy MOD 2) ;construct psf
;    sigx=motion(j)/2.35/umperpix ; 
;    psf=EXP(-xpsf^2/(2.*sigx^2)) & psf=psf/TOTAL(psf);/umperpix
    psf=FLTARR(Nyy) & psf[WHERE(ABS(xpsf) LE motion(j)/umperpix/2. )]=1. & psf=psf/TOTAL(psf)
    ;***************Max Entropy deconvolve image
    multipliers=FLTARR(Nyy) 
    Niter=10 & FOR ii=1,Niter DO Max_Entropy,yy,psf,yyd,multipliers,FT_PSF=psf_ft,/ONED
    yy=yyd
  ENDIF
  ;*******************
  
  est=[10,g(j),25.,0,0,0]
  q=AbelBasis2(xx,yy,del,filter=filter,est=est,Rmin=Rmin,edge=edge,q0=q0,q1=q1)
  IF N_ELEMENTS(qq) EQ 0 THEN qq=q ELSE qq=[qq,q]
  ;*******************
  IF SIZE(TgName,/TYPE) NE 0 AND KEYWORD_SET(density) THEN BEGIN
    kmass=4.*!pi*INT_TABULATED(q.x,q.cfit*q.x^2)/1.E12 ;g
    ;If max mu*mass is greater than value at the maximum mass then 
    ;set av_mu equal to that at the maximum mass to avoid extrapoloation problems.
    IF kmass*1.E6 GT MAX(tg.mumass) THEN av_mu=tg.av_mu[tg.N-1] $
    ELSE av_mu=INTERPOL(tg.av_mu,tg.mumass,kmass*1.E6) ;get av_mu based on observed mu*mass
    mass=kmass/av_mu
    ;********************
    IF j EQ 0 THEN avK=av_mu ELSE avK=[avK,av_mu]
    ;********************
    q.krho=q.krho/av_mu
    q.a[0]=q.a[0]/av_mu
    q.cfit=q.cfit/av_mu
    q0.krho=q0.krho/av_mu
    q1.krho=q1.krho/av_mu
  
  ENDIF 
;  XYOUTS,0.02*(Rlim[1]-Rlim[0])+Rlim[0],rholim[0]+0.85*(rholim(1)-rholim(0)),$
;    STRCOMPRESS(STRING(t(j),FORMAT='(F7.3)'))+' ns',CHARSIZE=charsz*0.5,CHARTHICK=thk
  
  s.Rmass[j]=INT_TABULATED(q.x,q.cfit*ABS(q.x)*q.x^2)/INT_TABULATED(q.x,q.cfit*q.x^2);center of mass
  s.maxdensity[j]=q.a[0]
  s.R0[j]=q.a[1]
  s.del[j]=2.*q.a[2]
  s.Rin[j]=q.a[1]-q.a[2]*3.
  s.Rabl[j]=q.a[1]+q.a[2]*3.

  ;w=WHERE(q.x GE s.Rin[j] AND q.x LE s.Rabl[j])
  w=INDGEN(N_ELEMENTS(q.x))
  
  OPLOT,q.x[w],q.krho[w],THICK=thk,COLOR=clr[j MOD N_ELEMENTS(clr)]
  OPLOT,[0,MAX(q.x)],[0,0],THICK=thk,LINESTYLE=1
  OPLOT,q0.x[w],q0.krho[w],THICK=thk/2.,LINESTYLE=1,COLOR=clr[j MOD N_ELEMENTS(clr)]
  OPLOT,q1.x[w],q1.krho[w],THICK=thk/2.,LINESTYLE=2,COLOR=clr[j MOD N_ELEMENTS(clr)]

  ;OPLOT,q0.x,q0.krhob,COLOR=clr[j],THICK=thk*2
  ;OPLOT,q0.x,q0.bkgfit,COLOR=clr[j],THICK=thk*2
  
  IF SIZE(dat,/TYPE) NE 8 THEN BEGIN
    a=REPLICATE({t:t[j],x:q.x,krho:q.krho,x0:q0.x,krho0:q0.krho,x1:q1.x,krho1:q1.krho},nt)
    dat={shot:shot,TgName:TgName,a:a}
  ENDIF ELSE BEGIN
    dat.a[j]={t:t[j],x:q.x,krho:q.krho,x0:q0.x,krho0:q0.krho,x1:q1.x,krho1:q1.krho}
  ENDELSE
  
;  OPLOT,q.x,q.cfit,THICK=thk,LINESTYLE=2
  
;  OPLOT,[s.Rmass[j],s.Rmass[j]],[rholim[0],rholim[1]],THICK=thk,LINESTYLE=0
;  OPLOT,[s.Rin[j],s.Rin[j]],[rholim[0],rholim[1]],THICK=thk,LINESTYLE=1
;  OPLOT,[s.Rabl[j],s.Rabl[j]],[rholim[0],rholim[1]],THICK=thk,LINESTYLE=1
    
  
jump:
ENDFOR

IF SIZE(TgName,/TYPE) NE 0 AND NOT(KEYWORD_SET(density)) THEN BEGIN
 OPLOT,tg.r,tg.rho*tg.mu,THICK=thk,COLOR=PURPLE
ENDIF

pm,[[s.R0],[s.del],[s.Rin],[s.Rabl]]

!P.MULTI=0
IF KEYWORD_SET(eps) OR KEYWORD_SET(ps)  THEN BEGIN
    DEVICE,/CLOSE
    SET_PLOT,'WIN'
ENDIF
!P.FONT=-1
!X.MARGIN=[10,3] & !Y.MARGIN=[4,2] 


RETURN
END

PRO writekrho,dat
;writes dat output to a file


nt=N_ELEMENTS(dat.a.t)


  filetosave="C:\Documents and Settings\hicks13\DESKTOP\Shot"+dat.shot+"krho.txt"
  OPENW,unit,filetosave,/GET_LUN
  PRINTF,unit,"Output on: "+SYSTIME()
  PRINTF,unit,''
  
  f='('+STRING(nt+1)+'A12)'
  g='('+STRING(nt+1)+'G12.5)'
  
  PRINTF,unit,STRING(['Time(ns):',STRING(dat.a.t)],FORMAT=f)
  PRINTF,unit,STRING(['R(um)',REPLICATE('krho(/cm)',nt)],FORMAT=f)

  FOR i=0,N_ELEMENTS(dat.a[0].x)-1 DO BEGIN
   PRINTF,unit,STRING([dat.a[0].x[i],dat.a.krho[i]],FORMAT=g)
  ENDFOR
  
;  PRINTF,unit,'Time (ns)'
;  PRINTF,unit,STRING(['',dat.a[0].t,dat.a[1].t,dat.a[2].t,dat.a[3].t],FORMAT='(5A12)')
;  PRINTF,unit,''
;  PRINTF,unit,STRING(['R(um)','krho(/cm)','krho(/cm)','krho(/cm)','krho(/cm)'],FORMAT='(5A12)')
;  PRINTF,unit,STRING([$
;            TRANSPOSE(dat[0].x),TRANSPOSE(dat[0].krho),TRANSPOSE(dat[1].krho),$
;            TRANSPOSE(dat[2].krho),TRANSPOSE(dat[3].krho)],FORMAT='(5G12.5)') ;
  CLOSE,unit
  FREE_LUN,unit

RETURN
END

FUNCTION PieceDeriv,x,y,ipiece,erry=erry,errdydxmin=errdydxmin
;Extracts piecewise linear derivatives
;errdydxmin is the miminum relative error in dydx (e.g. 0.05 for velocity)

Ndat=N_ELEMENTS(x)
z=FLTARR(Ndat) & zp=FLTARR(N_ELEMENTS(ipiece)/2)
xp=zp & dxp=zp & yp=zp & dyp=zp & dydxp=zp & ddydxp=zp

FOR j=0,N_ELEMENTS(zp)-1 DO BEGIN
  imin=ipiece(2*j) & imax=ipiece(2*j+1)
  xp(j)=MEAN(x(imin:imax))

  IF KEYWORD_SET(erry) THEN BEGIN
    c_y=POLY_FIT(x(imin:imax),y(imin:imax),1,COVAR=cov,/DOUBLE,MEASURE_ERRORS=erry(imin:imax)) 
  ENDIF ELSE BEGIN
    c_y=POLY_FIT(x(imin:imax),y(imin:imax),1,COVAR=cov,/DOUBLE);,MEASURE_ERRORS=erry(imin:imax))   
  ENDELSE
  
  yp(j)=POLY(xp(j),c_y) 
  dyp(j)=SQRT(cov(0,0)+xp(j)^2*cov(1,1)+2.*xp(j)*cov(0,1))

  c_dydx=c_y(1)
  dydxp(j)=POLY(xp(j),c_dydx)
  ddydxp(j)=SQRT(cov(1,1))
ENDFOR

IF KEYWORD_SET(errdydxmin) THEN ddydxp=SQRT(ddydxp^2+(dydxp*errdydxmin)^2)

RETURN,{xx:x,yy:y,x:xp,dx:dxp,y:yp,dy:dyp,dydx:dydxp,ddydx:ddydxp}
END

PRO compare_velocities,s,eps=eps
;Compares velocities for CoM, inside edge, ablation surface etc.

RED=255 & GREEN=255*256L & BLUE=255*256L*256L & YELLOW=RED+GREEN & PINK=RED+BLUE & LBLUE=GREEN+BLUE
ORANGE=255+126*256L & PURPLE=126+BLUE & PLAIN=RED+GREEN+BLUE & BACKG=0 & BROWN=174+71*256L
charsz=1.2 & thk=1.5 & symsz=1.3

!P.FONT=-1
IF KEYWORD_SET(eps) OR KEYWORD_SET(ps) THEN BEGIN
!P.FONT=1
    SET_PLOT,'ps'
    !X.MARGIN=[5.5,1.5] & !Y.MARGIN=[3.,1.]    
    DEVICE,ENCAPSULATED=1,LANDSCAPE=0,XSIZE=17,YSIZE=8,/COLOR,/HELVETICA,FILENAME="C:\Documents and Settings\hicks13\DESKTOP\IDL.EPS"
    TVLCT, [0,255,0,0,255,255,0,255,126,255,174], [0,0,255,0,255,0,255,126,0,255,71], [0,0,0,255,0,255,255,0,255,255,0]
    PLAIN=0 & RED=1 & GREEN=2 & BLUE=3 & YELLOW=4 & PINK=5 & LBLUE=6 & ORANGE=7 & PURPLE=8 & BACKG=9 & BROWN=10
    FILL=1 & charsz=1.3 & thk=3.0*thk & symsz=0.7*symsz
ENDIF
clr=[PLAIN,RED,GREEN,BLUE,BROWN,ORANGE,PINK,LBLUE,PURPLE,YELLOW]
b = FINDGEN(17)*(!PI*2/16.) & USERSYM,1*COS(b),1*SIN(b),THICK=thk,FILL=1

tlim=[20.8,22.2] & Rlim=[0,600] & Vlim=[200,370]


;********************

!p.multi=[0,2,1]

PLOT,[tlim[0],tlim[1]],[Rlim[0],Rlim[1]],/NODATA,/XSTYLE,/YSTYLE,XTHICK=thk,YTHICK=thk,CHARTHICK=thk,CHARSIZE=charsz,$
  XTITLE='Time (ns)',YTITLE='Radius (!Mmm)'
OPLOT,s.t,s.Rin,THICK=thk,SYMSIZE=symsz,PSYM=2,COLOR=RED
OPLOT,s.t,s.Rabl,THICK=thk,SYMSIZE=symsz,PSYM=5,COLOR=BLUE
OPLOT,s.t,s.Rmass,THICK=hk,SYMSIZE=symsz,PSYM=8,COLOR=PLAIN

PLOT,[tlim[0],tlim[1]],[Vlim[0],Vlim[1]],/NODATA,/XSTYLE,/YSTYLE,XTHICK=thk,YTHICK=thk,THICK=hk,CHARTHICK=thk,CHARSIZE=charsz,$
  XTITLE='Time (ns)',YTITLE='Velocity (!Mmm/ns)'
p=piecederiv(s.t,s.Rin,[0,5,3,8,6,10],errdydxmin=0.05)
XYERRPLOT,p.x,-p.dydx,p.dx,p.ddydx,THICK=thk,SYMSIZE=symsz,PSYM=2,COLOR=RED
p=piecederiv(s.t,s.Rabl,[0,5,3,8,6,10],errdydxmin=0.05)
XYERRPLOT,p.x,-p.dydx,p.dx,p.ddydx,THICK=thk,SYMSIZE=symsz,PSYM=5,COLOR=BLUE
p=piecederiv(s.t,s.Rmass,[0,5,3,8,6,10],errdydxmin=0.05)
XYERRPLOT,p.x,-p.dydx,p.dx,p.ddydx,THICK=thk,SYMSIZE=symsz,PSYM=8,COLOR=PLAIN

!P.MULTI=0
IF KEYWORD_SET(eps) OR KEYWORD_SET(ps)  THEN BEGIN
    DEVICE,/CLOSE
    SET_PLOT,'WIN'
ENDIF
!P.FONT=-1
!X.MARGIN=[10,3] & !Y.MARGIN=[4,2] 

RETURN
END


PRO testbasis_omega,eps=eps
;tests a new inversion scheme using decomposition into basis functions
COMMON bparams

RED=255 & GREEN=255*256L & BLUE=255*256L*256L & YELLOW=RED+GREEN & PINK=RED+BLUE & LBLUE=GREEN+BLUE
ORANGE=255+126*256L & PURPLE=126+BLUE & PLAIN=RED+GREEN+BLUE & BACKG=0 & BROWN=174+71*256L
charsz=1.2 & thk=1.5 & symsz=0.4
clr=[PLAIN,RED,GREEN,BLUE,BROWN,PINK,LBLUE,ORANGE,PURPLE,YELLOW]

IF KEYWORD_SET(eps) THEN BEGIN
    SET_PLOT,'ps'
    DEVICE,/ENCAPSULATED,LANDSCAPE=0,XSIZE=9.5,YSIZE=7.5,/INCHES,/COLOR,FILENAME="C:\Documents and Settings\hicks13\DESKTOP\IDL.EPS"
    TVLCT, [0,255,0,0,255,255,0,255,126,255,174], [0,0,255,0,255,0,255,126,0,255,71], [0,0,0,255,0,255,255,0,255,255,0]
    PLAIN=0 & RED=1 & GREEN=2 & BLUE=3 & YELLOW=4 & PINK=5 & LBLUE=6 & ORANGE=7 & PURPLE=8 & BACKG=9 & BROWN=10
    clr=[PLAIN,RED,GREEN,BLUE,BROWN,PINK,LBLUE,ORANGE,PURPLE,YELLOW]
    FILL=1 & charsz=1.0 & thk=3.0*thk & symsz=0.7*symsz
ENDIF



;***********Load streak profile
mur0=28.9350*1.E-4
file='c:\hicks\omega\data\0712\results\Original Analysis\49780_IvsY_02.txt'
dat=readdat(file,4) ;dat=readxvisdat(file,x,y,ysig)
xm=REFORM(dat(0,*)) & y=REFORM(dat(1,*)) & N=N_ELEMENTS(y)
w=WHERE(xm GE 0) & sign=1.
;w=WHERE(xm LE 0) & sign=-1.
x=xm(w) 
r=x ;define radial positions for funcbasis().
dx=ABS(x(1)-x(0)) ;get x-spacing

;*********Noise filtering
;P=SMOOTH(P,4/0.875)
filt=FIX((10./dx-1)/2.) & y = LEEFILT(y,filt,/DOUBLE) 
;;Result = WV_TOOL_DENOISE(P); [, X] [, Y] [, GROUP_LEADER=widget_id] [, TITLE=string] [, UNITS=string] [, XTITLE=string] [, XUNITS=string] [, YTITLE=string] [, YUNITS=string]) 
;Pw=WV_DENOISE(P,"Daubechies",2,DWT_FILT=dwt_filt,PERCENT=98,THRESHOLD=1,WPS_FILT=wps_filt)
;P=Pw(0:Nx-1)
;**********************

P=-ALOG(y(w))/mur0 ;Take natural log, multiply by -1.
x=sign*xm(w) & Nx=N_ELEMENTS(x) ;

;*******Create test profile using Gaussian density profile
;del0=20. & k0=200.
;Nx=500 & xmax=250.
;x=DINDGEN(Nx)/(Nx-1)*xmax
;;f=DOUBLE(ABS(ABS(x)-ABS(k0)) LE ABS(del0)/2.)
;f=4.01*EXP(-(ABS(x)-k0)^2/2/(del0/2)^2)
;bkg=50+0.007*x^2-2.E-5*x^3
;noise=RANDOMN(5,Nx)*100.
;P=abeltransform(x,f,M,/init) + bkg + noise

;***Fix outside of projection to zero to avoid violent uptick in eta at end. This is critical
P_end=MEAN(P(Nx-0.05*Nx-1:Nx-1)) & P=P-P_end 

!P.MULTI=[0,2,2]
;PLOT,xm,y,/YSTYLE

del=2 & nbasis=ROUND(MAX(x)/del)-1 & c_k = IMSL_FCNLSQ('funcbasis',nbasis,x,P,/DOUBLE,SSE=sig) 
;PRINT,SQRT(sig/Nx)

func=DBLARR(Nx) & FOR k=1,nbasis DO func=func+c_k(k-1)*funcbasis(k,x)
eta=DBLARR(Nx)  & FOR k=1,nbasis DO eta=eta+c_k(k-1)*funcbasis(k,x,/orig)
;interpolate between basis coefficients to recover radial density profile (for square basis functions only)
;INTERPOL(c_k,FINDGEN(nbasis)*del,x) 

est=[2.1585516,210.07641,9.1109823,-0.33608420,0.00011869420,6.3929830e-006]
cfit=GAUSSFIT(x,eta,a,NTERMS=6,ESTIMATES=est)
print,a
dx=ABS(x(1)-x(0)) & print,dx
nf=2.*ABS(a(2))/dx ;# points on left and right to fit in Savitzky-Golay filter
nf=ROUND(nf)
savgolFilter = SAVGOL(nf, nf, 0, 4) 
etaFilt=CONVOL(eta, savgolFilter, /EDGE_TRUNCATE)

PLOT,x,y(w),YRANGE=[0,35000],CHARSIZE=charsz,THICK=thk,XTHICK=thk,YTHICK=thk,CHARTHICK=thk,$
    XTITLE="Position (!4l!3m)",YTITLE="Pixel brightness (arb units)"

PLOT,x,func,CHARSIZE=charsz,THICK=thk,XTHICK=thk,YTHICK=thk,CHARTHICK=thk,$
    XTITLE="Position (!4l!3m)",YTITLE="Projected depth (!4l!3m)"
pp0=!p & px0=!x & py0=!y; capture system variables for this plot
OPLOT,x,P,COLOR=RED,THICK=thk

PLOT,x,eta,YRANGE=[-1,MAX(eta)],CHARSIZE=charsz,THICK=thk,XTHICK=thk,YTHICK=thk,CHARTHICK=thk,$
    XTITLE="Radius (!4l!3m)",YTITLE="!4q/q!d0!4"
pp1=!p & px1=!x & py1=!y; capture system variables for this plot
OPLOT,x,cfit,THICK=thk,COLOR=green,LINESTYLE=1
OPLOT,x,etaFilt,THICK=thk,COLOR=RED,LINESTYLE=2
;w=WHERE(x LT 70.)
;print,mean(eta(w)),stddev(eta(w))
;OPLOT,x,f,COLOR=RED,LINESTYLE=2
;
;Try other methods for Abel inversion
;Normal equations 
junk=abeltransform(x,P,M,/init) ;Used to generate M matrix only
;;;MM=2.*dx*M & f=LA_LEAST_SQUARES(TRANSPOSE(MM)#MM, TRANSPOSE(MM)#SMOOTH(P,1),/DOUBLE) 
;;; Compute the LU decomposition:
lambda=1.E6
H=get_H(Nx,2)
MM=2.*dx*M
Aludc=MATRIX_MULTIPLY(MM,MM,/ATRANSPOSE)+lambda*H 
print,TRACE(MATRIX_MULTIPLY(MM,MM,/ATRANSPOSE))/TRACE(H)
LA_LUDC,Aludc,index 
b=MATRIX_MULTIPLY(MM,P,/ATRANSPOSE) & f=LA_LUSOL(Aludc, index, b,/DOUBLE) 
OPLOT,x,f,COLOR=PINK,THICK=thk
!p=pp0 & !x=px0 & !y=py0
OPLOT,x,abeltransform(x,f,M),COLOR=PINK,THICK=thk,LINESTYLE=1
!p=pp1 & !x=px1 & !y=py1

;Filtered backprojection
P=INTERPOL(P,x,x+dx/2.) ;extract values of P from x+dx/2. in prep. for creating a profile evenly symmetric about x=0
x=INTERPOL(x,x,x+dx/2.) ;make sure x starts at dx/2. in prep. for creating a profile evenly symmetric about x=0
x=[-REVERSE(x),x] & P=[REVERSE(P),P]

q=fbp(x,P)
OPLOT,x,q,THICK=thk,COLOR=BLUE,LINESTYLE=1

!P.MULTI=0

IF KEYWORD_SET(eps)  THEN BEGIN
    DEVICE,/CLOSE
    SET_PLOT,'WIN'
ENDIF

RETURN
END
PRO testdeconv
;Attempts to deconvolve - creates higher amplitude noise. Doesn't work for now
RED=255 & GREEN=255*256L & BLUE=255*256L*256L & YELLOW=RED+GREEN & PINK=RED+BLUE & LBLUE=GREEN+BLUE
ORANGE=255+126*256L & PURPLE=126+BLUE & PLAIN=RED+GREEN+BLUE & BACKG=0 & BROWN=174+71*256L
charsz=1.2 & thk=1.5 & symsz=0.4

mur0=28.9350*1.E-4
file='c:\hicks\omega\data\0712\results\49780_IvsY_02.txt'
dat=readdat(file,4) ;dat=readxvisdat(file,x,y,ysig)
xm=REFORM(dat(0,*)) & y=REFORM(dat(1,*)) & N=N_ELEMENTS(y)

sig=15./2.35
psf=EXP(-xm^2/(2.*sig^2)) & psf=psf/TOTAL(psf)

!P.MULTI=[0,4,4]
y=MAX(y)-y
Niter=16 ;# of iterations of maximum entropy procedure to recover original function
FOR i=1,Niter DO BEGIN
    Max_Entropy,y,psf,deconv,multipliers
    PLOT,xm,y
    OPLOT,xm,deconv,COLOR=RED
ENDFOR

;PLOT,xm,y
;OPLOT,xm,deconv,COLOR=RED
RETURN
END