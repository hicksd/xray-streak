PRO abel_integral,ii,jj,I0,I1,I2
;Calculates the integrals appearing in the Abel transformation matrix

i=FLOAT(ii) & j=FLOAT(jj)

IF i EQ 0 AND j EQ 0 THEN BEGIN
    I0=0. & I1=0 & GOTO,jump
ENDIF

IF j GT i THEN BEGIN
    a=SQRT((j+0.5)^2-i^2)
    b=SQRT((j-0.5)^2-i^2)
    I0=ALOG(((j+0.5)+a)/((j-0.5)+b))
    I1=a-b-j*I0
ENDIF ELSE BEGIN
    a=SQRT((j+0.5)^2-i^2)
    b=SQRT(j^2-i^2)
    I0=ALOG(((j+0.5)+a)/(j+b))
    I1=a-b-j*I0
ENDELSE
jump:

RETURN
END

FUNCTION abeltransform,rr,ff,M,INIT=INIT
;Forward Abel transform (i.e. projection)
;Integration performed using 1st order Taylor expansion of integrand (3 pt integration)
;rr must contain one and only one zero and be linearly increasing.
;M=Transformation matrix
;INIT= requests that M be initialized.
rr=REFORM(rr) & ff=REFORM(ff)
NN=N_ELEMENTS(rr) & dr=rr(1)-rr(0)

;************Shift rr and ff to include r=0 as one of the elements
Roffset=MIN(rr(WHERE(ABS(rr) EQ MIN(ABS(rr)))));this is the value of rr closest to zero
;keep a record of original rr and ff
rr_orig=rr  & ff_orig=ff
;Shift rr and ff to include r=0
rr=rr_orig-Roffset  & ff=INTERPOL(ff_orig,rr_orig,rr)

;****Check to make sure that rr contains one zero and is uniformly spaced
;****THIS CHECK ISN'T REALLY NEEDED ANYMORE BECAUSE OF Roffset CORRECTION ABOVE
;asc=rr-SHIFT(rr,1) & lin=TOTAL(ABS(asc(1:NN-1)))/(NN-1) ;lin NE dr implies non-linear
w0=WHERE(rr EQ 0) & srt=SORT(rr) & mono=TOTAL(ABS(rr(srt)-rr))
IF w0(0) EQ -1 OR N_ELEMENTS(w0) GT 1 OR mono NE 0 THEN BEGIN
    PRINT,"rr must contain one and only one zero and be linearly increasing"
    RETURN,FLTARR(NN)
ENDIF

;****Choose which part of rr is longer, and thus to be used for calculating transformation matrix
wp=WHERE(rr GE 0) & wm=WHERE(rr LE 0)
IF N_ELEMENTS(wp) EQ 1 THEN Np=0 ELSE Np=N_ELEMENTS(wp)
IF N_ELEMENTS(wm) EQ 1 THEN Nm=0 ELSE Nm=N_ELEMENTS(wm)
IF Np GE Nm THEN N=Np ELSE N=Nm

;*******Calculate transformation matrix (for r, starting at zero)
;The only input to M is the number of radial elements, N. (does not need f, r, or dr)
IF N_ELEMENTS(M) EQ 0 OR KEYWORD_SET(init) THEN BEGIN
M=FLTARR(N,N) ;Transformation matrix
   FOR i=0,N-1 DO BEGIN
    FOR j=i,N-1 DO BEGIN
        IF j EQ -1 THEN GOTO,jump

        IF i EQ 0 AND j EQ 0 THEN BEGIN
         M(i,j)=0.5
         GOTO,jump
        ENDIF

;*****original matrix: includes j=i-1 (must set in loop above)
;*****by taking central diff approx everywhere
;       abel_integral,i,j+1,I0p,I1p
;       IF j EQ i-1 THEN BEGIN
;         M(i,j)=-((j+1)/2*I1p)
;         GOTO,jump
;       ENDIF
;       abel_integral,i,j,I0,I1
;       IF j EQ i THEN BEGIN
;         M(i,j)=(j*I0+I1)-((j+1)/2*I1p)
;         GOTO,jump
;       ENDIF
;       abel_integral,i,j-1,I0m,I1m
;       M(i,j)=((j-1)/2*I1m)+(j*I0+I1)-((j+1)/2*I1p)
;*****New matrix: avoids j=i-1 by taking forward difference approx @ j=i, central diff elsewhere
       IF j EQ i-1 THEN GOTO,jump
       abel_integral,i,j+1,I0p,I1p
       abel_integral,i,j,I0,I1
       IF j EQ i THEN BEGIN
         M(i,j)=(j*I0+(1-j)*I1)-((j+1)/2*I1p) ;avoids using j=i-1
            GOTO,jump
       ENDIF
       abel_integral,i,j-1,I0m,I1m
       IF j EQ i+1 THEN BEGIN
       M(i,j)=((j-1)*I1m)+(j*I0+I1)-((j+1)/2*I1p)
       GOTO,jump
       ENDIF

       M(i,j)=((j-1)/2*I1m)+(j*I0+I1)-((j+1)/2*I1p)


       jump:
    ENDFOR
   ENDFOR
ENDIF

;****Calculate the projection vector in radius
IF Nm EQ 0 THEN RETURN,2*dr*M#ff
IF Np EQ 0 THEN RETURN,REVERSE(2*dr*M#REVERSE(ff))
;****Then reflect to create remaining half of final vector
PP=FLTARR(NN)
IF Np GE Nm THEN BEGIN
    Pp=2*dr*M#ff(wp) & Pm=Pp(1:Nm-1)
ENDIF ELSE BEGIN
    Pm=2*dr*M#REVERSE(ff(wm)) & Pp=Pm(1:Np-1)
ENDELSE
PP=[REVERSE(Pm),Pp]

;*********Convert back to original rr by adding Rpffset offset
PP=INTERPOL(PP,rr,rr_orig) & rr=rr_orig & ff=ff_orig

RETURN,PP
END


PRO testabel
;tests abel transform with known analytic solutions

N=397
x=FINDGEN(N)/(N-1)*20
x=-REVERSE(x)
;x=x-x(N/2)

;******Solid sphere
f=FLTARR(N_ELEMENTS(x))
R0=5. & w=WHERE(ABS(x) LE R0)
f(w)=1.
P=abeltransform(x,f,M) ;forward Abel transform
PLOT,x,P
Ptheory=FLTARR(N)
w=WHERE(ABS(x) LE R0) & IF w(0) NE -1 THEN Ptheory(w)=2*SQRT(R0^2-x(w)^2)
OPLOT,x,Ptheory,COLOR=255

;*******Annulus
f=FLTARR(N_ELEMENTS(x))
R0=5. & del=1. & w=WHERE(ABS(x) LE R0 AND ABS(x) GE R0-del)
f(w)=1.
P=abeltransform(x,f,M) ;forward Abel transform
OPLOT,x,P
Ptheory=FLTARR(N)
w=WHERE(ABS(x) LE R0 AND ABS(x) GE R0-del) & IF w(0) NE -1 THEN Ptheory(w)=2*SQRT(R0^2-x(w)^2)
w2=WHERE(ABS(x) LE R0-del) & IF w2(0) NE -1 THEN Ptheory(w2)=2*(SQRT(R0^2-x(w2)^2)-SQRT((R0-del)^2-x(w2)^2) )
OPLOT,x,Ptheory,COLOR=255

;*******Gaussian shell
f=FLTARR(N_ELEMENTS(x))
R0=5. & s=1.
f=EXP(-(ABS(x)-R0)^2/2/s^2)
P=abeltransform(x,f,M)
OPLOT,x,P

RETURN
END