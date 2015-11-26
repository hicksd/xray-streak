


FUNCTION flatfield,data,flat
;divides out the flatfield response
    data=FLOAT(data)/FLOAT(flat)
RETURN,data
END


FUNCTION darkimage,data,dark
;subtract the dark image
    data=data-dark
RETURN,data
END

PRO xrfc

path='c:\hicks\radiography\OMEGAdata\0410\'
fim1=path+'xrfccd_xrfc4_t2_37626.spe'
fdk1=path+'xrfccd_xrfc4_t2_37626_bkg.spe'
fim2=path+'xrfccd_xrfc4_t2_37627.spe'
fdk2=path+'xrfccd_xrfc4_t2_37627_bkg.spe'

image=read_spe(fim1) & dk=read_spe(fdk1) & image=image-dk
;bkg=read_spe(fim2) & dk=read_spe(fdk2) & bkg=bkg-dk

;bkg(WHERE(bkg EQ 0))=1
;image=FLOAT(image)/FLOAT(bkg)

dk=0 & bkg=0
SAVE,image,FILENAME=path+'37626.sav'

RETURN
END

PRO testbkg,dat,bkg,sz,im,init=init

bin=3.

IF NOT(KEYWORD_SET(init)) THEN GOTO,start

path='c:\hicks\radiography\omegadata\'

sub=500. & xm=[543,1942] & ym=[204,1603]
file='xrf4t2_36563.hdf'
dat=omegahdf(path+file) & dat=MEDIAN(dat,5)
dat=dat(xm(0):xm(1),ym(0):ym(1))
dat=IfromOD(dat,sub)
;sz=SIZE(dat) ;& IF bin GT 1 THEN dat=CONGRID(dat,sz(1)/bin,sz(2)/bin)

sub=550. & xm=543+[0,xm(1)-xm(0)] & ym=204+[0,ym(1)-ym(0)]
file='serpentine_xrf4t2_36547.hdf'
bkg=omegahdf(path+file) & bkg=MEDIAN(bkg,5)
bkg=bkg(xm(0):xm(1),ym(0):ym(1))
bkg=IfromOD(bkg,sub)
;sz=SIZE(bkg) ;& IF bin GT 1 THEN bkg=CONGRID(bkg,sz(1)/bin,sz(2)/bin)

start:

sz=SIZE(dat)

dsub=dat(100:200,400:500)
bsub=bkg(100:200,400:500)
tvscl,dsub
tvscl,bsub,111,0
;TVSCL,CONGRID(((dat)),sz(1)/bin,sz(2)/bin)
;TVSCL,CONGRID(((bkg)),sz(1)/bin,sz(2)/bin),sz(1)/bin,0

;fdat=FFT(dat) & fdat(670:729,670:729)=0
;fbkg=FFT(bkg) & fbkg(670:729,670:729)=0
;
;TVSCL,CONGRID(FFT(fdat,/INVERSE),sz(1)/bin,sz(2)/bin)
;TVSCL,CONGRID(FFT(fbkg,/INVERSE),sz(1)/bin,sz(2)/bin),sz(1)/bin,0

RETURN
END

PRO testfft

y=FLTARR(512)
y(0:200)=3 & y(201:400)=1
;y=y+randomn(0,512)/3.

!p.multi=[0,2,2]
plot,y
fy=FFT(y)
plot,abs(fy),/ylog
fy(251:255)=0 & fy(256:260)=0
plot,fft(fy,/inverse)

!p.multi=0
RETURN
END