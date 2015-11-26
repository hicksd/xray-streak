PRO wshow_images,imstr,imagename,origin=origin

WIDGET_CONTROL,/HOURGLASS

iname=WHERE(TAG_NAMES(imstr) EQ imagename)
im=imstr.(iname(0))

imstr.oimwindow -> IDLgrWindow::GetProperty,DIMENSIONS=D ;D=[xsize/factor,ysize/factor]
imstr.oimview -> IDLgrView::SetProperty,VIEWPLANE_RECT=[0,0,D(0),D(1)],COLOR=[0,0,0]


TVLCT,rval,gval,bval,/GET
imstr.oimpalette -> SetProperty,RED=rval,GREEN=gval,BLUE=bval
imstr.oimwindow -> SetProperty, PALETTE=imstr.oimpalette

IF imstr.roiclr EQ 1 THEN BEGIN
    xmin=imstr.inroi(0,0) & ymin=imstr.inroi(1,0)
    xmax=imstr.inroi(0,1) & ymax=imstr.inroi(1,1)
    IF imagename EQ 'IMAGE' AND $
       xmax GT 0 AND xmin GT 0 AND xmax GT xmin AND $
       ymax GT 0 AND ymin GT 0 AND ymax GT ymin AND $
       xmax LT imstr.xsize*imstr.factor AND $
       ymax LT imstr.ysize*imstr.factor THEN BEGIN
    cmin=MIN(im(xmin:xmax,ymin:ymax))
    cmax=MAX(im(xmin:xmax,ymin:ymax))
    ENDIF ELSE BEGIN
       PRINT,"ROI is incompatible for color typing"
       cmin=MIN(im) & cmax=MAX(im)
    ENDELSE
ENDIF ELSE BEGIN
    cmin=MIN(im) & cmax=MAX(im)
ENDELSE

ctable=!D.TABLE_SIZE

sz=SIZE(im)
nx=sz(1)/imstr.factor & ny=sz(2)/imstr.factor

CASE imstr.logplot OF
0: BEGIN
    dat=(CONGRID(im, nx, ny)-cmin)/(cmax-cmin)
    lo=WHERE(dat LE 0)
    hi=WHERE(dat GE 1)
    IF lo(0) NE -1 THEN dat(lo)=0
    IF hi(0) NE -1 THEN dat(hi)=1
    imstr.oimimage -> SetProperty, $
      DATA=(ctable-1)*dat;(CONGRID(im, nx, ny)-cmin)/(cmax-cmin)
   END
1: BEGIN
    imstr.oimimage -> SetProperty, $
      DATA=ctable*(ALOG(CONGRID(im, nx, ny)-cmin+1))/ALOG(cmax-cmin+2)
   END
ELSE:
ENDCASE

IF KEYWORD_SET(origin) THEN BEGIN
    imstr.oimimage -> IDLgrImage::SetProperty,LOCATION=origin
ENDIF ELSE imstr.oimimage -> IDLgrImage::SetProperty,LOCATION=[0,0]


imstr.oimwindow -> Draw, imstr.oimview


END