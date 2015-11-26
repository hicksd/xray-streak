@XmuFromMass
@max_entropy

PRO xstreak_event, event


WIDGET_CONTROL, event.id, GET_UVALUE = eventval     ;find the user value
                   ;of the widget where
                   ;the event occured
IF N_ELEMENTS(eventval) EQ 0 THEN RETURN
CASE eventval OF
'THEMENU': BEGIN
    CASE event.value OF

    "Open":BEGIN
       xvisopen,event
       END
    "Copy Image":BEGIN;copy image to clipboard
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
     IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No image loaded")
         GOTO,break_copyimage
       ENDIF

       oimclipb=OBJ_NEW('IDLgrClipboard')
       oimclipb -> SetProperty, PALETTE=imstr.oimpalette

       oimclipb -> Draw, imstr.oimview, VECTOR=0

       break_copyimage:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
       END
    "Copy Plot":BEGIN;copy plot to clipboard
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
     IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No image loaded")
         GOTO,break_copyplot
       ENDIF

       oimclipb=OBJ_NEW('IDLgrClipboard')

       imstr.oplwindow -> GetProperty, DIMENSIONS=dims ;IDLgrWindow object
       ;oimclipb -> GetProperty, DIMENSIONS=odims
       oimclipb -> SetProperty, DIMENSIONS=dims

       oimclipb -> Draw, imstr.oplview, VECTOR=0;,POSTSCRIPT=0,VECT_TEXT_RENDER_METHOD=1

       break_copyplot:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
       END
    "TIF (no bin)":BEGIN
      savetif,event,bin=0
    END
    "TIF (binned)":BEGIN
      savetif,event,bin=1
    END
    "Raw binary (no bin)":BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
     IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No image loaded")
         GOTO,break_save_binary
       ENDIF

       file=DIALOG_PICKFILE(PATH=imstr.pathsave,FILE="*.raw",$
         FILTER="*.raw",/WRITE,GET_PATH=path,GROUP=event.handler)

       IF file EQ '' THEN GOTO,break_save_binary

       imstr.pathsave=path
     OPENU,unit,imstr.pathfile,/GET_LUN
         PRINTF,unit,imstr.pathopen
         PRINTF,unit,imstr.pathsave
     CLOSE,unit
     FREE_LUN,unit

;print,file
;print,imstr.factor,imstr.xsize,imstr.ysize
       OPENW,unit,file,/GET_LUN
           isz=SIZE(imstr.image)
print,isz
           WRITEU,unit,UINT(isz(1)),UINT(isz(2))
           WRITEU,unit,imstr.image
       CLOSE,unit
       FREE_LUN,unit

       break_save_binary:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
        "Raw ROI":BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
     IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No image loaded")
         GOTO,break_save_rawroi
       ENDIF

       file=DIALOG_PICKFILE(PATH=imstr.pathsave,FILE="*.raw",$
         FILTER="*.raw",/WRITE,GET_PATH=path,GROUP=event.handler)

       IF file EQ '' THEN GOTO,break_save_rawroi

       imstr.pathsave=path
     OPENU,unit,imstr.pathfile,/GET_LUN
         PRINTF,unit,imstr.pathopen
         PRINTF,unit,imstr.pathsave
     CLOSE,unit
     FREE_LUN,unit

       sz=SIZE(imstr.image)
       xmin=imstr.xyroi(0,0) & xmax=imstr.xyroi(0,1)
       ymin=imstr.xyroi(1,0) & ymax=imstr.xyroi(1,1)

       IF xmax-xmin LE 0 OR xmin LT 0 OR xmax LT 0 THEN BEGIN
         xmin=0 & xmax=sz(1)-1
       END
       IF ymax-ymin LE 0 OR ymin LT 0 OR ymax LT 0 THEN BEGIN
         ymin=0 & ymax=sz(2)-1
       END

       OPENW,unit,file,/GET_LUN
           isz=SIZE(imstr.image(xmin:xmax,ymin:ymax))
print,isz
           WRITEU,unit,UINT(isz(1)),UINT(isz(2))
           WRITEU,unit,imstr.image(xmin:xmax,ymin:ymax)
       CLOSE,unit
       FREE_LUN,unit

       break_save_rawroi:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    "Save Scaled Image":BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
     IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No image loaded")
         GOTO,break_save_scaled_image
       ENDIF

       file=DIALOG_PICKFILE(PATH=imstr.pathsave,FILE="*.eps",$
         FILTER="*.eps",/WRITE,GET_PATH=path,GROUP=event.handler)

       IF file EQ '' THEN GOTO,break_save_scaled_image

       imstr.pathsave=path
     OPENU,unit,imstr.pathfile,/GET_LUN
         PRINTF,unit,imstr.pathopen
         PRINTF,unit,imstr.pathsave
     CLOSE,unit
     FREE_LUN,unit

       sz=SIZE(imstr.image)

       xmin=imstr.xyroi(0,0) & xmax=imstr.xyroi(0,1)
       ymin=imstr.xyroi(1,0) & ymax=imstr.xyroi(1,1)

       IF xmax-xmin LE 10 OR xmin LT 0 OR xmax LT 0 THEN BEGIN
         xmin=0 & xmax=sz(1)-1
       END
       IF ymax-ymin LE 10 OR ymin LT 0 OR ymax LT 0 THEN BEGIN
         ymin=0 & ymax=sz(2)-1
       END

       xsz=xmax-xmin+1 & ysz=ymax-ymin+1
;print,xmin,xmax,ymin,ymax
       im=imstr.image(xmin:xmax,ymin:ymax)
        px=FINDGEN(xsz)+xmin & py=FINDGEN(ysz)+ymin
        tt=time(imstr,px,py,dd)

       IF imstr.factor GT 1. THEN BEGIN
            tt=CONGRID(tt,xsz/imstr.factor)
            dd=CONGRID(dd,ysz/imstr.factor)
            im=CONGRID(im,xsz/imstr.factor,ysz/imstr.factor)
       ENDIF

       cmin=MIN(im) & cmax=MAX(im)
       IF imstr.logplot EQ 1 THEN BEGIN
         im=(ALOG(im-cmin+1))/ALOG(cmax-cmin+2)
       ENDIF ELSE BEGIN
         im=(im-cmin)/(cmax-cmin)
       ENDELSE

       savscaleimage,file,im,tt,dd

       break_save_scaled_image:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    "Save Plot as Text":BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
       
       IF imstr.nplots LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("Insufficient plots available" )
         GOTO,break_saveplottotext
       ENDIF

       filetosave=DIALOG_PICKFILE(PATH=imstr.pathsave,GET_PATH=path,GROUP=event.handler)

       IF filetosave NE '' THEN BEGIN
         preexist=FINDFILE(filetosave)

         IF preexist(0) NE '' THEN BEGIN
          answer=DIALOG_MESSAGE("File already exists. Overwrite?",/QUESTION)
          IF answer EQ 'No' THEN GOTO,break_saveplottotext
         ENDIF
         plottotext,imstr,filetosave,/streakdump

         imstr.pathsave=path
           OPENU,unit,imstr.pathfile,/GET_LUN
            PRINTF,unit,imstr.pathopen
            PRINTF,unit,imstr.pathsave
         CLOSE,unit
         FREE_LUN,unit
       ENDIF

       break_saveplottotext:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    "Save all data":BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
       IF SIZE(imstr.dat,/TYPE) NE 8 THEN BEGIN
        junk=DIALOG_MESSAGE("No Inversion data exists")
        GOTO,break_savealldata
       ENDIF

       filetosave=DIALOG_PICKFILE(PATH=imstr.pathsave,GET_PATH=path,GROUP=event.handler)

       IF filetosave NE '' THEN BEGIN
         preexist=FINDFILE(filetosave)

         IF preexist(0) NE '' THEN BEGIN
          answer=DIALOG_MESSAGE("File already exists. Overwrite?",/QUESTION)
          IF answer EQ 'No' THEN GOTO,break_savealldata
         ENDIF
         xstreak_savedata,imstr,filetosave
         
         answer=DIALOG_MESSAGE("Do you want to save the entire data structure as an IDL save file? (may be large)",/QUESTION)
         IF answer EQ 'Yes' THEN BEGIN
          alldat=imstr.dat
          SAVE,alldat,FILENAME=filetosave+'.sav'
         ENDIF
         
         imstr.pathsave=path
           OPENU,unit,imstr.pathfile,/GET_LUN
            PRINTF,unit,imstr.pathopen
            PRINTF,unit,imstr.pathsave
         CLOSE,unit
         FREE_LUN,unit
       ENDIF

       break_savealldata:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    "Print Image":BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
     IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No image loaded")
         GOTO,break_printimage
       ENDIF

       oprint = OBJ_NEW('IDLgrPrinter',COLOR_MODEL=1)
       oprint -> SetProperty, PALETTE=imstr.oimpalette, LANDSCAPE=1
       result=DIALOG_PRINTERSETUP(oprint)
       ;result=DIALOG_PRINTJOB(oprint)

       IF result NE 0 THEN BEGIN
         oprint -> Draw, imstr.oimview, VECTOR=0
         oprint -> NewPage
         oprint -> NewDocument
       END

       break_printimage:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    "Print Plot":BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
     IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No image loaded")
         GOTO,break_printplot
       ENDIF

       oprint = OBJ_NEW('IDLgrPrinter',COLOR_MODEL=1)
       result=DIALOG_PRINTERSETUP(oprint)
       ;result=DIALOG_PRINTJOB(oprint)

       IF result NE 0 THEN BEGIN
         oprint -> Draw, imstr.oplview, VECTOR=1
         ;oprint -> NewPage
         oprint -> NewDocument
       END

       break_printplot:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    "Exit": BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
       OBJ_DESTROY,imstr.oimwindow
       OBJ_DESTROY,imstr.oimview
       OBJ_DESTROY,imstr.oplwindow
       OBJ_DESTROY,imstr.oplview
     WIDGET_CONTROL, event.top, /DESTROY
       END
    "Show ROI":BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
     IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No image loaded")
         GOTO,break_showroi
       ENDIF

       imstr.oimpolyxyroi -> SetProperty, HIDE=0
       imstr.oimwindow -> Draw, imstr.oimview

       break_showroi:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    "Log/Linear":BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
     IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No image loaded")
         GOTO,break_log
       ENDIF

       IF imstr.logplot EQ 0 THEN imstr.logplot=1 ELSE imstr.logplot=0

       CASE imstr.whichimage OF
       'IMAGE':wshow_images,imstr,'IMAGE'
       'PHASE':wshow_images,imstr,'PHASE',origin=[imstr.phroi(0,0),imstr.phroi(1,0)]/imstr.factor
       'CUMULPHASE':wshow_images,imstr,'CUMULPHASE',origin=[imstr.phroi(0,0),imstr.phroi(1,0)]/imstr.factor
       'AMPL':wshow_images,imstr,'AMPL',origin=[imstr.phroi(0,0),imstr.phroi(1,0)]/imstr.factor
       'NORMALAMPL':BEGIN
          logplot=imstr.logplot & imstr.logplot=0 ;make sure that logplot is not on (since density has lots of zeros)
          wshow_images,imstr,'NORMALAMPL',origin=[imstr.phroi(0,0),imstr.phroi(1,0)]/imstr.factor
          imstr.logplot=logplot ;restore logplot setting
        END
       ENDCASE
       
       break_log:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    "Refresh":BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
     IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No image loaded")
         GOTO,break_ref
       ENDIF

       wshow_images,imstr,'IMAGE'
       imstr.whichimage='IMAGE'

       break_ref:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    "Color Table": BEGIN
     WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
     IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No image loaded")
         GOTO,break_ct
       ENDIF
       XLOADCT, GROUP = event.handler,/MODAL
       wshow_images,imstr,'IMAGE'
       help,imstr.image

       break_ct:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY

       END
    "Color Palette": BEGIN
     WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
       IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No image loaded")
         GOTO,break_pal
       ENDIF
       XPALETTE, GROUP = event.top;,/MODAL
       wshow_images,imstr,'IMAGE'

       break_pal:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY

       END
    "Linear Fit":BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
       oplot0=imstr.oplplot(0) & oplot0->GetProperty, DATA=dat
       IF N_ELEMENTS(dat) EQ 0 THEN BEGIN
         junk=DIALOG_MESSAGE("A plot must exist before a fit can be performed")
         GOTO,break_lfit
       ENDIF

       x=REFORM(dat(0,*)) & y=REFORM(dat(1,*))

       IF imstr.nplots GT 1 THEN BEGIN
         oplot1=imstr.oplplot(1) & oplot1-> GetProperty, DATA=dat
         yerr=ABS(REFORM(dat(1,*))-y)
       ENDIF ELSE yerr=SQRT(y)
       ;IF N_ELEMENTS(dat) EQ 0 THEN BEGIN
       ; junk=DIALOG_MESSAGE("Data must contain errors to perform a fit")
       ; GOTO,break_lfit
       ;ENDIF

       c=POLY_FIT(x,y,1,MEASURE_ERRORS=yerr,SIGMA=sig)
       ;PRINT,"Y = ("+STRCOMPRESS(STRING(c(1)))+"+/-"+STRCOMPRESS(STRING(sig(1)))+")*X + "+STRCOMPRESS(STRING(c(0)))+"+/-"+STRCOMPRESS(STRING(sig(0)))
       PRINT,"Coeffs = ",REFORM(c)
       PRINT,"Errors = ",REFORM(sig)
       imstr.oplplot(imstr.nplots) = OBJ_NEW('IDLgrPlot')
       imstr.oplmodel -> Add, imstr.oplplot(imstr.nplots)

       daty=FLTARR(N_ELEMENTS(x))
       FOR i=0,N_ELEMENTS(c)-1 DO daty=daty+c(i)*x^i

       oplot0=imstr.oplplot(0) & oplot0->GetProperty, XRANGE=xr0,YRANGE=yr0
       oplotn=imstr.oplplot(imstr.nplots) & oplotn -> SetProperty, DATAX=x,DATAY=daty,$;c(0)+c(1)*x, $
              MIN_VALUE=yr0(0),MAX_VALUE=yr0(1),LINESTYLE=2 & imstr.oplplot(imstr.nplots)=oplotn

       imstr.nplots=imstr.nplots+1

       imstr.oplwindow -> Draw, imstr.oplview

       break_lfit:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
       END
    "Quadratic Fit":BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
       oplot0=imstr.oplplot(0) & oplot0-> GetProperty, DATA=dat
       IF N_ELEMENTS(dat) EQ 0 THEN BEGIN
         junk=DIALOG_MESSAGE("A plot must exist before a fit can be performed")
         GOTO,break_qfit
       ENDIF

       x=REFORM(dat(0,*)) & y=REFORM(dat(1,*))

       IF imstr.nplots GT 1 THEN BEGIN
         oplot1=imstr.oplplot(1) & oplot1->GetProperty, DATA=dat
         yerr=ABS(REFORM(dat(1,*))-y)
       ENDIF ELSE yerr=SQRT(y)
       ;IF N_ELEMENTS(dat) EQ 0 THEN BEGIN
       ; junk=DIALOG_MESSAGE("Data must contain errors to perform a fit")
       ; GOTO,break_lfit
       ;ENDIF

       c=POLY_FIT(x,y,2,SIGMA=sig,MEASURE_ERRORS=yerr)
       ;PRINT,"Y = ("+STRCOMPRESS(STRING(c(1)))+"+/-"+STRCOMPRESS(STRING(sig(1)))+")*X + "+STRCOMPRESS(STRING(c(0)))+"+/-"+STRCOMPRESS(STRING(sig(0)))
       PRINT,"Coeffs = ",REFORM(c)

       imstr.oplplot(imstr.nplots) = OBJ_NEW('IDLgrPlot')
       imstr.oplmodel -> Add, imstr.oplplot(imstr.nplots)

       daty=FLTARR(N_ELEMENTS(x))
       FOR i=0,N_ELEMENTS(c)-1 DO daty=daty+c(i)*x^i
;daty=splinesmooth(x,x,y);,sig=SQRT(N_ELEMENTS(x)))
;daty=SMOOTH(y,N_ELEMENTS(x)/10,/EDGE_TRUNCATE)
;daty=TS_SMOOTH(y,N_ELEMENTS(x)/10)
       oplot0=imstr.oplplot(0) & oplot0->GetProperty, XRANGE=xr0,YRANGE=yr0
       oplotn=imstr.oplplot(imstr.nplots) & oplotn->SetProperty, DATAX=x,DATAY=daty,$;c(0)+c(1)*x, $
              MIN_VALUE=yr0(0),MAX_VALUE=yr0(1),LINESTYLE=2 & imstr.oplplot(imstr.nplots)=oplotn

       imstr.nplots=imstr.nplots+1

       imstr.oplwindow -> Draw, imstr.oplview

       break_qfit:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
       END
    "Gaussian Fit":BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
       oplot0=imstr.oplplot(0) & oplot0-> GetProperty, DATA=dat
       IF N_ELEMENTS(dat) EQ 0 THEN BEGIN
         junk=DIALOG_MESSAGE("A plot must exist before a fit can be performed")
         GOTO,break_gfit
       ENDIF

       x=REFORM(dat(0,*)) & y=REFORM(dat(1,*))

       ;IF imstr.nplots GT 1 THEN BEGIN
       ; imstr.oplplot(1) -> GetProperty, DATA=dat
       ; yerr=ABS(REFORM(dat(1,*))-y)
       ;ENDIF ELSE yerr=SQRT(y)
       ;IF N_ELEMENTS(dat) EQ 0 THEN BEGIN
       ; junk=DIALOG_MESSAGE("Data must contain errors to perform a fit")
       ; GOTO,break_lfit
       ;ENDIF

       yfit = GAUSSFIT(x,y,a,NTERMS=6)
       PRINT,"Y = "+STRCOMPRESS(STRING(a(0)))+"exp((x-"+$
         STRCOMPRESS(STRING(a(1)))+")^2/2*"$
         +STRCOMPRESS(STRING(a(2)))+"^2)"

       print,a
       imstr.oplplot(imstr.nplots) = OBJ_NEW('IDLgrPlot')
       imstr.oplmodel -> Add, imstr.oplplot(imstr.nplots)

       oplot0=imstr.oplplot(0) & oplot0-> GetProperty, XRANGE=xr0,YRANGE=yr0
       oplotn=imstr.oplplot(imstr.nplots) & oplotn-> SetProperty, DATAX=x,DATAY=yfit, $
              MIN_VALUE=yr0(0),MAX_VALUE=yr0(1),LINESTYLE=2 & imstr.oplplot(imstr.nplots)=oplotn

       imstr.nplots=imstr.nplots+1

       imstr.oplwindow -> Draw, imstr.oplview

       break_gfit:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
       END
;    "Find Warp Matrix":BEGIN
;       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
;       WIDGET_CONTROL,/HOURGLASS
;
;       phxmin=imstr.phroi(0,0) & phxmax=imstr.phroi(0,1)
;       phymin=imstr.phroi(1,0) & phymax=imstr.phroi(1,1)
;
;       xmin=imstr.xyroi(0,0) & xmax=imstr.xyroi(0,1)
;       ymin=imstr.xyroi(1,0) & ymax=imstr.xyroi(1,1)
;
;       sz=SIZE(imstr.image)
;       IF xmin LT 0 OR xmax GT sz(1) OR ymin LT 0 OR ymax GT sz(2) OR $
;          phxmin LT 0 OR phxmax GT sz(1) OR phymin LT 0 OR phymax GT sz(2)   THEN BEGIN
;         junk=DIALOG_MESSAGE("Plot ROI or Phase ROI exceeds bounds of image")
;         GOTO,break_findwarp
;       ENDIF
;       IF N_ELEMENTS(imstr.cumulphase) LE 1 THEN BEGIN
;         junk=DIALOG_MESSAGE("A phase image must be calculated before warping can be performed")
;         GOTO,break_findwarp
;       END
;
;       IF xmin GE phxmin AND xmax LE phxmax AND $
;         ymin GE phymin AND ymax LE phymax THEN BEGIN
;
;         cmlphim=imstr.cumulphase(xmin-phxmin:xmax-phxmin,ymin-phymin:ymax-phymin)
;
;         sz=SIZE(cmlphim) & npx=sz(1) & npy=sz(2)
;         px=INDGEN(npx) & py=INDGEN(npy)
;
;         ;Generate transformation matrix elements, c(npolyX,npolyY)
;         ;First, fit each column of data (y-direction) to a polynomial
;         ;of order npolyY
;         npolyY=imstr.warp.nPolyY
;         cY=FLTARR(npx,npolyY+1)
;         FOR ii=0,npx-1 DO cY(ii,*)=POLY_FIT(py,cmlphim(ii,*),npolyY)
;         ;Now fit each coeffcient from the above fit to a polynomial
;         ;of order npolyX.
;         npolyX=imstr.warp.nPolyX
;         c=FLTARR(npolyX+1,npolyY+1)
;         FOR ii=0,npolyY DO c(*,ii)=POLY_FIT(px,cY(*,ii),npolyX)
;
;         ;Now estimate the conversion from fringe shift to y-pixel shift
;         ;by taking the very first column only and finding the average phase
;         ;difference between each pixel. This is potentially the weakest point
;         ;in this analysis.
;         dphi=cumulative_phase(imstr.phase(xmin-phxmin,ymin-phymin:ymax-phymin),/RETURN_DPHI)
;         dphi=REFORM(dphi)
;         ;Find the average phase change per pixel at xmin in y-direction
;         av=MEAN(dphi) & sig=STDDEV(dphi)
;         norm=WHERE(ABS(dphi - av) LT 2.*sig)
;         IF norm(0) NE -1 THEN dphi0=MEAN(dphi(norm)) ELSE dphi0=MEAN(dphi0)
;
;         dF0=dphi0/2./!pi
;         PRINT,"Average Pixels/Fringe, dY/dF = "+STRCOMPRESS(STRING(1/dF0))
;
;
;       ENDIF ELSE BEGIN
;         junk=DIALOG_MESSAGE("Plot ROI must be contained within Phase ROI")
;         WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
;         RETURN
;       ENDELSE
;
;       ;Convert the transformation matrix c (presently describing the fringe
;       ;shift error dF) to a pixel shift error dY. dF0=fringes per pixel.
;       c=c/dF0
;
;       imstr.warp.c=c
;       imstr.warp.xmin=xmin & imstr.warp.xmax=xmax
;       imstr.warp.ymin=ymin & imstr.warp.ymax=ymax
;       imstr.warp.exist=1
;print,c
;print,xmin,xmax,ymin,ymax
;       break_findwarp:
;       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
;       RETURN
;
;    END
;    "Apply Warp":BEGIN
;       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
;       IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
;         junk=DIALOG_MESSAGE("No image loaded")
;         GOTO,break_applywarp
;       ENDIF
;       IF imstr.warp.exist EQ 0 THEN BEGIN
;         junk=DIALOG_MESSAGE("No warp matrix exists")
;         GOTO,break_applywarp
;       ENDIF
;
;       WIDGET_CONTROL,/HOURGLASS
;
;       xmin=imstr.warp.xmin & xmax=imstr.warp.xmax
;       ymin=imstr.warp.ymin & ymax=imstr.warp.ymax
;
;       nPolyX=imstr.warp.nPolyX & nPolyY=imstr.warp.nPolyY
;
;       c=imstr.warp.c
;
;       ;Apply the dY shift transformation to each column of data
;       sz=SIZE(imstr.image)
;
;       IF xmin LT 0 OR xmax GT sz(1) OR ymin LT 0 OR ymax GT sz(2) THEN BEGIN
;         junk=DIALOG_MESSAGE("Edges of warp region exceed current image size")
;         GOTO,break_applywarp
;       ENDIF
;
;       oldim=imstr.image(xmin:xmax,ymin:ymax)
;
;       y=FINDGEN(ymax-ymin+1)
;       dy=FINDGEN(ymax-ymin+1)
;
;       FOR ii=0,xmax-xmin DO BEGIN
;         IF nPolyX EQ 1 AND nPolyY EQ 1 THEN BEGIN
;          dy=(c(0,0)+c(1,0)*ii)+(c(0,1)+c(1,1)*ii)*y
;         ENDIF ELSE BEGIN
;          IF nPolyX EQ 2 AND nPolyY EQ 2 THEN BEGIN
;            dy=(c(0,0)+c(1,0)*ii+c(2,0)*ii^2)+(c(0,1)+c(1,1)*ii+c(2,1)*ii^2)*y+$
;            (c(0,2)+c(1,2)*ii+c(2,2)*ii^2)*y^2
;          ENDIF ELSE BEGIN
;              FOR jj=0,ymax-ymin DO BEGIN
;                 dy(jj)=TOTAL(c##TRANSPOSE(REPLICATE(ii,npolyX+1)^INDGEN(npolyX+1))$
;                   *y(jj)^INDGEN(npolyY+1) )
;              ENDFOR
;          ENDELSE
;         ENDELSE
;         imstr.image(xmin+ii,ymin:ymax)=INTERPOL(oldim(ii,*),y,y+dY)
;       ENDFOR
;
;       break_applywarp:
;       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
;    END
    "Camera Settings":BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
       result = CAMERA(imstr,event)
       IF result(0) EQ -1 THEN GOTO,break_cam

       WIDGET_CONTROL,imstr.rot_id,SET_VALUE=result(0)
       WIDGET_CONTROL,imstr.shr_id,SET_VALUE=result(1)
       WIDGET_CONTROL,imstr.mag_id,SET_VALUE=result(2)
       ;WIDGET_CONTROL,imstr.delay_id,SET_VALUE=result(3)
;       WIDGET_CONTROL,imstr.a_id,SET_VALUE=result(4)
;       WIDGET_CONTROL,imstr.b_id,SET_VALUE=result(5)
;       WIDGET_CONTROL,imstr.c_id,SET_VALUE=result(6)
;       WIDGET_CONTROL,imstr.fidpx_id,SET_VALUE=result(7)
;       WIDGET_CONTROL,imstr.fidtt_id,SET_VALUE=result(8)
      imstr.MinPix=result(3)
      imstr.MaxPix=result(4)
      FOR iii=1,N_ELEMENTS(imstr.tpoly)-1 DO imstr.tpoly(iii)=result(iii+5-1)
      FOR jjj=1,N_ELEMENTS(imstr.tpoly)-1 DO imstr.spoly(jjj-1)=imstr.tpoly(jjj)*jjj

       break_cam:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
   "Target Description":BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY

       Tgdata = TARGET(imstr,event) 
       IF SIZE(Tgdata,/TYPE) EQ 8 THEN BEGIN ;create tg structure if supplied by TARGET.
        ;help,/struct,Tgdata
        tg=loadtarget(XKAPPAFILE=imstr.xrayfile,TDAT=Tgdata)
        ;help,/struct,tg
        imstr=imagetostructure(imstr,'TARGETDATA',tg)
        imstr.loadedtarget=1 ;flip switch to indicate target has been loaded (stays loaded forever)
        ;print,imstr.targetdata.mat(0)
       ENDIF ELSE junk=DIALOG_MESSAGE("Warning! No new target information was loaded!")
       break_target:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    "Rotation":BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
       IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No image loaded")
         GOTO,break_rot
       ENDIF
       WIDGET_CONTROL,imstr.rot_id,GET_VALUE=rads
       WIDGET_CONTROL,/HOURGLASS

       drad=FLOAT(rads(0))-imstr.lastrot
       imstr.image=ROT(imstr.image,FLOAT(drad)/!pi*180.,/INTERP)
       imstr.lastrot=FLOAT(rads(0))

       wshow_images,imstr,'IMAGE'
       break_rot:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    "Right 90":BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
       IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No image loaded")
         GOTO,break_rotR
       ENDIF
       WIDGET_CONTROL,/HOURGLASS

       newimstr=RotateFlip(imstr,'Right 90')

       WIDGET_CONTROL,event.handler, SET_UVALUE=newimstr, /NO_COPY
       break_rotR:
    END
    "Left 90":BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
       IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No image loaded")
         GOTO,break_rotL
       ENDIF
       WIDGET_CONTROL,/HOURGLASS

       newimstr=RotateFlip(imstr,'Left 90')

       WIDGET_CONTROL,event.handler, SET_UVALUE=newimstr, /NO_COPY
       break_rotL:
    END
    "Flip Horizontal":BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
       IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No image loaded")
         GOTO,break_flipH
       ENDIF
       WIDGET_CONTROL,/HOURGLASS

       newimstr=RotateFlip(imstr,'Flip Horizontal')

       WIDGET_CONTROL,event.handler, SET_UVALUE=newimstr, /NO_COPY
       break_flipH:
    END
    "Flip Vertical":BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
       IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No image loaded")
         GOTO,break_flipV
       ENDIF
       WIDGET_CONTROL,/HOURGLASS

       newimstr=RotateFlip(imstr,'Flip Vertical')

       WIDGET_CONTROL,event.handler, SET_UVALUE=newimstr, /NO_COPY
       break_flipV:
    END
    "Shear":BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
       IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No image loaded")
         GOTO,break_shear
       ENDIF
       WIDGET_CONTROL,imstr.shr_id,GET_VALUE=shearr

       shr=FLOAT(shearr(0))-imstr.lastshear
       imstr.lastshear=FLOAT(shearr(0))

       nrows=N_ELEMENTS(imstr.image(0,*)-1)
       ncols=N_ELEMENTS(imstr.image(*,0)-1)

       WIDGET_CONTROL,/HOURGLASS
       FOR ii=0,nrows-1 DO BEGIN
         row=imstr.image(*,ii)
         dx=-shr*ii
         n=FIX(dx) & f=dx MOD 1
         imstr.image(*,ii)=(1-f)*SHIFT(row,n)+f*SHIFT(row,n+1)
       ENDFOR

       wshow_images,imstr,'IMAGE'
       break_shear:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    'Median Filter':BEGIN
    WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
      IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
       junk=DIALOG_MESSAGE("No image loaded")
       GOTO,break_medfil
      ENDIF
    WIDGET_CONTROL,imstr.medfilval_id,GET_VALUE=med

    WIDGET_CONTROL,/HOURGLASS

    imstr.image=MEDIAN(imstr.image,FLOAT(med(0)))
    ;imstr.image=SMOOTH(imstr.image,FLOAT(med(0)))
    ;imstr.image=wvfilter(imstr.image,FLOAT(med(0)),FLOAT(med(0)))
    ;imstr.image=LEEFILT(imstr.image,FLOAT(med(0)))
    ;imstr.image = WV_DENOISE(imstr.image,"Daubechies",2,PERCENT=97,$;COEFFICIENTS=1024,$
    ;    DWT_FILT=dwt_filt,THRESHOLD=0,WPS_FILT=wps_filt)

    wshow_images,imstr,'IMAGE'
    break_medfil:
    WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    'Subtract Image':BEGIN
        WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
      IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
       junk=DIALOG_MESSAGE("No image loaded")
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
       GOTO,break_subbkg
      ENDIF ELSE WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY

    xvisopen,event,/BKG
    break_subbkg:
    END
    'Subtract DC offset':BEGIN
    WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
      IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
       junk=DIALOG_MESSAGE("No image loaded")
       GOTO,break_subDC
      ENDIF
    WIDGET_CONTROL,imstr.bkg_id,GET_VALUE=sbkg

    WIDGET_CONTROL,/HOURGLASS

    b=FIX(sbkg(0))
    lb=WHERE(imstr.image LT b)
    gb=WHERE(imstr.image GE b)
    IF lb(0) NE -1 THEN imstr.image(lb) = 0
    IF gb(0) GT -1 THEN imstr.image(gb) = imstr.image(gb) - b

    wshow_images,imstr,'IMAGE'
    break_subDC:
    WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    'Subtract plot ROI':BEGIN
    WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
      IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
       junk=DIALOG_MESSAGE("No image loaded")
       GOTO,break_subplotroi
      ENDIF
    
    WIDGET_CONTROL,/HOURGLASS

    xmin=imstr.xyroi[0,0] & xmax=imstr.xyroi[0,1]
    ymin=imstr.xyroi[1,0] & ymax=imstr.xyroi[1,1]
    
    im=TEMPORARY(imstr.image)
    sz=SIZE(im) & nx=sz(1) & ny=sz(2)
    
    bkg=REFORM(REBIN(im(xmin:xmax,ymin:ymax),xmax-xmin+1,1))
    bkg=SMOOTH(bkg,25)
    bkg=REBIN(bkg,xmax-xmin+1,ny)
    im(xmin:xmax,0:ny-1)=im(xmin:xmax,0:ny-1)-bkg
    
    ;w=WHERE(im LT 0) & IF w(0) NE -1 THEN im(w)=0
    imstr.image=TEMPORARY(im)
        
    wshow_images,imstr,'IMAGE'
    break_subplotroi:
    WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    'Subtract Int ROI (wire)':BEGIN
    WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
      IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
       junk=DIALOG_MESSAGE("No image loaded")
       GOTO,break_subintroi
      ENDIF
    
    ask=DIALOG_MESSAGE("Do you want to subtract the bkg after looking at the bkg plots?",/DEFAULT_NO,/QUESTION)
    
    WIDGET_CONTROL,/HOURGLASS

    im=TEMPORARY(imstr.image)
    sz=SIZE(im) & nx=sz(1) & ny=sz(2)

    ixmin=imstr.inroi[0,0] & ixmax=imstr.inroi[0,1]
    pxmin=imstr.xyroi[0,0] & pxmax=imstr.xyroi[0,1]

    dw=100. ;# pixels on either side of ixmin, ixmid, ixmax over which to search for the wire minimum
    axmin=pxmin; 150. ;min extent of useful data region
    axmax=pxmax; 4050. ;max extent of useful data region
   

    IF ixmin LT dw OR ixmax GE nx-dw-1 THEN BEGIN
      PRINT,"Automatically find rough positions of wires"
      ixmin=700 & ixmax=3500
    ENDIF
    
    ixmid=MEAN([ixmin,ixmax]) ;middle wire position
    
    spc=FIX(ABS(ixmax-ixmin)/4.) ;spacing between wire and inter-image minima
    ;wmins=[ixmin-spc,ixmin,ixmin+spc,ixmid,ixmid+spc,ixmax,ixmax+spc] ;ConA: use wire and inter-image and edge of image
    wmins=[ixmin-spc,ixmin,ixmid,ixmax,ixmax+spc] ;ConA: Use wire and edge of image only
    ;wmins=[ixmin,ixmid,ixmax] ;ConAw
    
    
;    mask=read_raw(imstr.xpath+'GXD1_mask_v1.raw')
    mask=read_raw(imstr.xpath+'GXD3_mask_v1.raw')
    ylim=FLTARR(2,4) ;y limits for each strip
    x=FINDGEN(nx)
    ;cc=FLTARR(3,4) ;array of quadratic polynomial coefficientss for bkgs on each strip
    
    ;ERASE,255+255*256L+255*256L*256L
    !p.multi=[0,2,2]
    FOR s=1,4 DO BEGIN ;loop through each strip
      w=WHERE(mask EQ s) & a=ARRAY_INDICES(mask,w) & ylim(0,s-1)=MIN(a(1,*)) & ylim(1,s-1)=MAX(a(1,*)) ;y-max and y-min for this strip
      sp=MEDIAN(im[*,ylim[0,s-1]:ylim[1,s-1]],DIMENSION=2) ;strip, w/ y-dimension collapsed
      w=FLTARR(N_ELEMENTS(wmins))-1 ;array to hold final wire or inter-image positions
      FOR i=0,N_ELEMENTS(w)-1 DO BEGIN ;Loop thru rough positions of each minima to determine exact minima positions
        lo=wmins(i)-dw & hi=wmins(i)+dw
        IF lo LT axmin THEN lo=axmin & IF lo GT axmax THEN lo=axmax
        IF hi LT axmin THEN hi=axmin & IF hi GT axmax THEN hi=axmax
        IF lo GT hi THEN lo=hi
        w(i)=MEAN(WHERE(sp EQ MIN(sp[lo:hi]) AND x GE lo AND x LE hi))
      ENDFOR
      IF MIN(w) EQ -1 THEN GOTO,skip
      
      yb=FLTARR(nx)
;      FOR j=0,2 DO BEGIN ;Loop through each subimage
;        a=2.*j & b=2.*(j+1)
;        cf=POLY_FIT(w[a:b],sp[w[a:b]],2) ;quadratic fit to wire & each side of sub-image
;        ;cf=POLY_FIT([w[a],w[b]],[sp[w[a]],sp[w[b]]],1) ;linear fit to each side of sub-image
;        yb[w[a]:w[b]]=POLY(x[w[a]:w[b]],cf)
;      ENDFOR
      order=2

      cf=POLY_FIT(w,sp(w),order) & yb=POLY(x,cf)
      PLOT,x,sp,BACKGROUND=255+255*256L+255*256L*256L,COLOR=0 
      OPLOT,x,yb,COLOR=255

      ;cc(*,s-1)=POLY_FIT(w,sp(w),order) & yb=POLY(x,cc(*,s-1))
;      order=N_ELEMENTS(w)-1
;      cff=POLY_FIT(w,sp(w),order) & yb=POLY(x,cff)
;      PLOT,sp & oplot,yb
      IF ask EQ 'Yes' THEN BEGIN
        bkg=REBIN(yb,nx,ylim(1,s-1)-ylim(0,s-1)+1)
        ;bkg=min(sp(where(x GT 1800 AND x LT 2200))) ;N120329 only
        im(*,ylim(0,s-1):ylim(1,s-1))=im(*,ylim(0,s-1):ylim(1,s-1))-bkg
      ENDIF
      skip:
    ENDFOR
    !p.multi=0
    
;    FOR s=1,4 DO BEGIN
;      bkg=REBIN(POLY(x,cc(*,s-1)),nx,ylim(1,s-1)-ylim(0,s-1)+1)
;      im(*,ylim(0,s-1):ylim(1,s-1))=im(*,ylim(0,s-1):ylim(1,s-1))-bkg
;    ENDFOR

;    w=WHERE(im LT 1) & IF w(0) NE -1 THEN im(w)=1
    imstr.image=TEMPORARY(im)
        
    wshow_images,imstr,'IMAGE'
    break_subintroi:
    WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    'Subtract DISC':BEGIN
      WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
      IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
       junk=DIALOG_MESSAGE("No image loaded")
       GOTO,break_subdisc
      ENDIF

    ask=DIALOG_MESSAGE("Do you want to subtract the bkg after looking at the bkg profiles?",/DEFAULT_NO,/QUESTION)
    
    WIDGET_CONTROL,/HOURGLASS

    im=(imstr.image)
    sz=SIZE(im) & nx=sz(1) & ny=sz(2)

;    ixmin=imstr.inroi[0,0] & ixmax=imstr.inroi[0,1]
    pxmin=imstr.xyroi[0,0] & pxmax=imstr.xyroi[0,1]
    gxmin=imstr.ghroi[0,0] & gxmax=imstr.ghroi[0,1] 
    pymin=imstr.xyroi[1,0] & pymax=imstr.xyroi[1,1]
    gymin=imstr.ghroi[1,0] & gymax=imstr.ghroi[1,1]
  
    dy=10 ;# of rows over which to average for top and bottom
    Npoly=10  ;order of polynomial to fit
    
    btm=REBIN(im[*,pymin:pymin+dy-1],nx,1) ;average dy rows above pymin
    top=REBIN(im[*,pymax-dy+1:pymax],nx,1) ;average dy rows below pymax
    ctr=REBIN(im[*,gymin:gymax],nx,1) ;average over ghroi region
    x=FINDGEN(nx) & y=FINDGEN(ny)
    
    ;**** fits in x-direction to the bottom, top, and center lineouts used to estimate bkg
    bcffx=POLY_FIT(x,btm,Npoly,YFIT=bfitx) 
    tcffx=POLY_FIT(x,top,Npoly,YFIT=tfitx)
    ccffx=POLY_FIT(x,ctr,Npoly,YFIT=cfitx)
    
;    PLOT,x,btm,BACKGROUND=255+255*256L+255*256L*256L,COLOR=0,YRANGE=[0,MAX([btm,top,ctr])] & OPLOT,x,bfitx,COLOR=0
;    OPLOT,x,top,COLOR=255      & OPLOT,x,tfitx,COLOR=255
;    OPLOT,x,ctr,COLOR=255*256L & OPLOT,x,cfitx,COLOR=255*256L
  
    OUTPLOT,x,ctr,imstr,/newplot,XTITL="Time (pixels)",YTITL="Intensity",TITL="Blk=ctr; Red=btm; Blue=top"
    OUTPLOT,x,cfitx,imstr,/overplot
    OUTPLOT,x,top,imstr,/overplot,COLOR=[0,0,255]
    OUTPLOT,x,tfitx,imstr,/overplot,COLOR=[0,0,255]
    OUTPLOT,x,btm,imstr,/overplot,COLOR=[255,0,0]
    OUTPLOT,x,bfitx,imstr,/overplot,COLOR=[255,0,0]
    
    bkg=FLTARR(nx,ny) ;create bkg array
    FOR i=0,nx-1 DO BEGIN
      cffy=POLY_FIT([pymin,MEAN([gymin,gymax]),pymax],[bfitx[i],cfitx[i],tfitx[i]],2)
      bkg[i,*]=POLY(y,cffy)
    ENDFOR
 ;   TVSCL,CONGRID(bkg,210,210)
    
    IF ask EQ 'Yes' THEN BEGIN
      im=im-bkg
      w=WHERE(im LT 1) & IF w[0] NE -1 THEN im[w]=1
    ENDIF
      
    imstr.image=TEMPORARY(im)
        
    wshow_images,imstr,'IMAGE'
    break_subdisc:
    WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY

    END
;    'Remove Ghost Fringes':BEGIN
;      WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
;
;      IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
;       junk=DIALOG_MESSAGE("No image loaded")
;       GOTO,break_ghost
;      ENDIF
;      IF imstr.ghroi(0,0) EQ -1 THEN BEGIN
;        junk=DIALOG_MESSAGE("Define Ghost ROI")
;        GOTO,break_ghost
;      ENDIF
;
;    WIDGET_CONTROL,/HOURGLASS
;
;    ghxmin=imstr.ghroi(0,0) & ghxmax=imstr.ghroi(0,1)
;    ghymin=imstr.ghroi(1,0) & ghymax=imstr.ghroi(1,1)
;
;    im=imstr.image(ghxmin:ghxmax,ghymin:ghymax)
;    sz=SIZE(im) & nx=sz(1) & ny=sz(2)
;
;Fim=FFT(im)
;;    F0=Fim(0,0)
;;    Fim(0,*)=0.5*(Fim(1,*)+Fim(nx-1,*)) ;get rid of Kx=0 (or set to adjacent Kx=1 values)
;;    Fim(0,0)=F0
;;Fim(0,*)=COMPLEX(0,1)*Fim(0,*)
;
;;fmin=12
;;fcol=Fim(0,*)
;;infcol=COMPLEXARR(ny)
;;infcol(0:fmin-1)=fcol(0:fmin-1)
;;infcol(ny-1-(fmin-1):ny-1)=fcol(ny-1-(fmin-1):ny-1)
;;Fim(0,*)=infcol
;
;;eps=MinGhost(Fim)
;;
;F0=Fim(0,0) & med0=MEDIAN(im)
;theta=!pi/2.;+eps
;Fim(0,*)=Fim(0,*)*COMPLEX(COS(theta),SIN(theta))
;Fim(0,0)=F0 ;set overall DC level to same as before ghost subtraction.
;;Fim(0,*)=0. ;& Fim(1,*)=0. & Fim(nx-1,*)=0.
;
;imf=(FFT(Fim,/INVERSE))
;
;    imstr.image(ghxmin:ghxmax,ghymin:ghymax)=imf
;
;    wshow_images,imstr,'IMAGE'
;
;    break_ghost:
;    WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
;    END
'Correct Gain Droop':BEGIN
      WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
      IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
       junk=DIALOG_MESSAGE("No image loaded")
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
       GOTO,break_gaindroop
      ENDIF
    WIDGET_CONTROL,/HOURGLASS

;    gaincorrect=read_raw(imstr.xpath+'GXD1_droop_corr_v1.raw')
;    mask=read_raw(imstr.xpath+'GXD1_mask_v1.raw')
    gaincorrect=read_raw(imstr.xpath+'GXD3_droop_corr_v1.raw')
    mask=read_raw(imstr.xpath+'GXD3_mask_v1.raw')
    
;;   w=WHERE(mask EQ 3.) &  IF w(0) NE -1 THEN gaincorrect(w)=gaincorrect(w)/1.65 ;fix unusually bright 3rd strip
;    w=WHERE(mask EQ 1.) & IF w(0) NE -1 THEN gaincorrect(w)=gaincorrect(w)/1. 
;    w=WHERE(mask EQ 2.) & IF w(0) NE -1 THEN gaincorrect(w)=gaincorrect(w)*1.0 
;    w=WHERE(mask EQ 3.) & IF w(0) NE -1 THEN gaincorrect(w)=gaincorrect(w)*2. 
;    w=WHERE(mask EQ 4.) & IF w(0) NE -1 THEN gaincorrect(w)=gaincorrect(w)*2. 

;    nx=imstr.xsize*imstr.factor
;    ny=imstr.ysize*imstr.factor
;    gaincorrect=REBIN(REVERSE(FINDGEN(nx)/(nx-1)*2.+1),nx,ny)
    imstr.image=gaincorrect*imstr.image
    
    
    wshow_images,imstr,'IMAGE'    
    break_gaindroop:
    WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
'Flatfield DISC':BEGIN
      WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
      IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
       junk=DIALOG_MESSAGE("No image loaded")
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
       GOTO,break_discflatfield
      ENDIF
    WIDGET_CONTROL,/HOURGLASS

    gaincorrect=read_raw(imstr.xpath+'DISC3 FF 111011_Manson_normalized_dewarped.raw')
    gaincorrect=SMOOTH(gaincorrect,[50,50])
    w=WHERE(gaincorrect LT 0.1) & IF w[0] NE -1 THEN gaincorrect[w]=0.1 ;avoid tiny values
    imstr.image=imstr.image/gaincorrect
    
    wshow_images,imstr,'IMAGE'    
    break_discflatfield:
    WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END    
'Deconvolve':BEGIN
      WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
      IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
       junk=DIALOG_MESSAGE("No image loaded")
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
       GOTO,break_deconvolve
      ENDIF
    WIDGET_CONTROL,/HOURGLASS

    WIDGET_CONTROL,imstr.deconx_id,GET_VALUE=deconx  & deconx=FLOAT(deconx[0]) ;pixels to deconvolve in x
    WIDGET_CONTROL,imstr.decony_id,GET_VALUE=decony  & decony=FLOAT(decony[0]) ;pixels to deconvolve in y
    
    IF deconx LE 1 OR decony LE 0 THEN BEGIN
      junk=DIALOG_MESSAGE("Deconv x and Deconv y must be >1")
      GOTO,break_deconvolve
    ENDIF
    
    im=TEMPORARY(imstr.image)
    sz=SIZE(im) & Nx0=sz(1) & Ny0=sz(2) ;original image size

    ;**************Create 2D PSF
;   Nx=2048. & Ny=2048.
   Nx=1024. & Ny=1024.
;    Nx=Nx0 & Ny=Ny0
    x=REBIN(FINDGEN(Nx)-Nx/2.+0.5*(Nx MOD 2),Nx,Ny)*FLOAT(Nx0/Nx)             ;units of pixels in original image
    y=TRANSPOSE(REBIN(FINDGEN(Ny)-Ny/2.+0.5*(Ny MOD 2),Ny,Nx))*FLOAT(Ny0/Ny)  ;units of pixels in original image
    psf=FLTARR(Nx,Ny) 
    ;****GXD1
    ;c=[0.78,0.2,0.02]& k=[10.,300.,1400.]
    ;FOR j=0,N_ELEMENTS(c)-1 DO psf=psf+!pi*c(j)/k(j)^2*EXP(-!pi^2/k(j)^2*(x^2+y^2))
    ;****GXD3
    c=[0.719,0.173,0.0267,0.0813]& k=[9.7,319.,1176.,670.]
    FOR j=0,2 DO psf=psf+!pi*c(j)/k(j)^2*EXP(-!pi^2/k(j)^2*(x^2+y^2))
    psf=psf+c(3)*0.548/(40000+(SQRT(x^2+y^2))^2.4)
    ;****
    psf=psf/TOTAL(psf)
    ;***************Create second PSF (for imaging slot)
    WIDGET_CONTROL,imstr.mag_id,GET_VALUE=micpx  & micpx=FLOAT(micpx(0)) ;mics/pixel
    WIDGET_CONTROL,imstr.slit_id,GET_VALUE=slit  &  slit=FLOAT(slit(0)) ;mic
;    sigx=slit/micpx/2.7 & sigy=4.;The 2.7 was calibrated to properly deconvolve a boxcar using a Gaussian
;    psf2=EXP(-x^2/(2.*sigx^2)-y^2/(2.*sigy^2)) & psf2=psf2/TOTAL(psf2)
;    xblur=slit/micpx & yblur=3.
    ;xblur=63. & yblur=68. ;DISC 120324: 30 um slit
    ;xblur=37. & yblur=69. ;DISC 120408/120409/120418: 17 um slit
    ;xblur=22. & yblur=69. ;DISC 120421: 10 um slit
    ;xblur=72. & yblur=63. ;DISC 120324; 30 um slit, 2.8 ns sweep, 500 um pc slit
;    xblur=72. & yblur=38. ;DISC 120408; 17 um slit, 2.8 ns sweep, 500 um pc slit
    ;xblur=72. & yblur=37. ;DISC 120409; 17 um slit, 2.8 ns sweep, 500 um pc slit
    ;xblur=72. & yblur=41. ;DISC 120418; 17 um slit, 2.8 ns sweep, 500 um pc slit
    ;xblur=70. & yblur=22. ;DISC 120421 ;10 um slit, 2.8 ns sweep, 500 um pc slit
    xblur=deconx & yblur=decony
    print,xblur,yblur
    im=SMOOTH(TEMPORARY(im),[xblur/2.,yblur/2.]) ;smooth over a half a deconvolution width
    psf2=FLTARR(Nx,Ny) & psf2[WHERE(ABS(x) LE xblur/2. AND ABS(y) LE yblur/2.)]=1. & psf2=psf2/TOTAL(psf2)
    ;***************Create a single PSF by convolving the two PSF's
;    psf=convolve(psf,psf2) & psf=psf/TOTAL(psf) & psf2=0.
    psf=psf2 & psf=psf/TOTAL(psf) & psf2=0 ;use this to skip deconvolution of the GXD psf portion
    ;;psf=SMOOTH(psf,[8,0]) & psf=psf/TOTAL(psf)
    ;***************Re-grid image to size of psf
    im=CONGRID(TEMPORARY(im),Nx,Ny,CUBIC=-0.5)
    ;***************Max Entropy deconvolve image
    Niter=10 & FOR ii=1,Niter DO Max_Entropy,im,psf,imd,multipliers,FT_PSF=psf_ft
    ;***************Make deconvolved image the original image size and insert into structure
    imstr.image=CONGRID(imd,Nx0,Ny0,CUBIC=-0.5)
    ;***************
    
    wshow_images,imstr,'IMAGE'    
    break_deconvolve:
    WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    'Remove Nonlinearity':BEGIN
    WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
      IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
       junk=DIALOG_MESSAGE("No image loaded")
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
       GOTO,break_nonlinear
      ENDIF
    WIDGET_CONTROL,/HOURGLASS
    
    ;3rd order polynomial fit to nonlinearity curve from G. Kyrala for GXD2, Oct 10, 2010.
    c=[5.76941,0.970545,6.92451e-006,2.95540e-010] & Imin=450.
    im=TEMPORARY(imstr.image)
    w=WHERE(im LE Imin)
    ;only apply correction to data above Imin
    im(w)=POLY(im(w),c)
    imstr.image=im
    
    wshow_images,imstr,'IMAGE'    
    break_nonlinear:
    WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END

    'Flatfield':BEGIN
        WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
      IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
       junk=DIALOG_MESSAGE("No image loaded")
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
       GOTO,break_flatfield
      ENDIF ELSE WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY

    xvisopen,event,/flatfield
    break_flatfield:
    END
    'Find Sweep':BEGIN
        WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
      IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
       junk=DIALOG_MESSAGE("No image loaded")
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
       GOTO,break_findsweep
      ENDIF ELSE WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY

    junk=DIALOG_MESSAGE("Find sweep procedure not yet implemented...suckah!")
    break_findsweep:
    END
    'Linearize Scale':BEGIN
      WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
      oplot0=imstr.oplplot(0) & oplot0-> GetProperty, DATA=dat
      IF N_ELEMENTS(dat) EQ 0 THEN BEGIN
        junk=DIALOG_MESSAGE("A plot must exist before a fit can be performed")
        GOTO,break_linearizescale
      ENDIF
      ;Find polynomial fit to x vs y and plot
       x=REFORM(dat(0,*)) & y=REFORM(dat(1,*))
       IF imstr.nplots GT 1 THEN BEGIN
         oplot1=imstr.oplplot(1) & oplot1->GetProperty, DATA=dat
         yerr=ABS(REFORM(dat(1,*))-y)
       ENDIF ELSE yerr=SQRT(y)

       c=POLY_FIT(y,x,4,SIGMA=sig,STATUS=stat) & PRINT,"Coeffs = ",REFORM(c)
       IF stat NE 0 THEN GOTO,break_linearizescale
       
       imstr.oplplot(imstr.nplots) = OBJ_NEW('IDLgrPlot')
       imstr.oplmodel -> Add, imstr.oplplot(imstr.nplots)

       datx=POLY(y,c)
       oplot0=imstr.oplplot(0) & oplot0->GetProperty,XRANGE=xr0,YRANGE=yr0
       oplotn=imstr.oplplot(imstr.nplots) & oplotn->SetProperty,DATAX=datx,DATAY=y,MIN_VALUE=yr0(0),MAX_VALUE=yr0(1),LINESTYLE=2
       imstr.nplots=imstr.nplots+1 & imstr.oplwindow -> Draw, imstr.oplview
      ;Map image to this new scale
       linearimage=POLY(imstr.image,c)
       imstr.image=linearimage
    wshow_images,imstr,'IMAGE'
    break_linearizescale:
    WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    'Film Scale':BEGIN
    WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
      IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
       junk=DIALOG_MESSAGE("No image loaded")
       GOTO,break_filmscale
      ENDIF

    WIDGET_CONTROL,/HOURGLASS
    imstr.image=10.^(0.0016*imstr.image-2.66)

    wshow_images,imstr,'IMAGE'
    break_filmscale:
    WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    'Image':BEGIN ;Display image
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY

   IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No image loaded")
         WIDGET_CONTROL,event.handler,SET_UVALUE=imstr,/NO_COPY
         RETURN
       ENDIF
       wshow_images,imstr,'IMAGE'

       imstr.whichimage = 'IMAGE'

       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    'Normalized intensity':BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
        IF N_ELEMENTS(imstr.phase) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No inversion image exists")
         WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
              RETURN
        ENDIF

       wshow_images,imstr,'PHASE',origin=[imstr.phroi(0,0),imstr.phroi(1,0)]/imstr.factor

       imstr.whichimage = 'PHASE'

       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    'Background':BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
        IF N_ELEMENTS(imstr.cumulphase) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No inversion image exists")
         WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
              RETURN
        ENDIF
       wshow_images,imstr,'CUMULPHASE',origin=[imstr.phroi(0,0),imstr.phroi(1,0)]/imstr.factor

       imstr.whichimage = 'CUMULPHASE'

       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    'Optical Depth':BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
        IF N_ELEMENTS(imstr.ampl) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No inversion image exists")
         WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
              RETURN
        ENDIF
       wshow_images,imstr,'AMPL',origin=[imstr.phroi(0,0),imstr.phroi(1,0)]/imstr.factor

       imstr.whichimage = 'AMPL'

       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    'Density':BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
        IF N_ELEMENTS(imstr.normalampl) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No inversion image exists")
         WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
              RETURN
        ENDIF
       logplot=imstr.logplot & imstr.logplot=0 ;make sure that logplot is not on (since density has lots of zeros)
       wshow_images,imstr,'NORMALAMPL',origin=[imstr.phroi(0,0),imstr.phroi(1,0)]/imstr.factor
       imstr.logplot=logplot ;restore logplot setting
       
       imstr.whichimage = 'NORMALAMPL'

       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END
    'About':BEGIN
       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
       junk=DIALOG_MESSAGE(TITLE="  ",$
       [imstr.title,"D. G. Hicks, May 2012","hicks13@llnl.gov",$
                    "Lawrence Livermore National Laboratory","DO NOT DISTRIBUTE WITHOUT PERMISSION"],/INFORMATION)
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
    END

    ELSE: BEGIN
     junk=DIALOG_MESSAGE("Program Error: Event User Value Not Found")
    ENDELSE
ENDCASE
END
'Roitype':BEGIN
    WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
    break_roitype:
    WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
END
'XYTABLE':BEGIN
    WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
    IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
       junk=DIALOG_MESSAGE('Image must be loaded first')
       GOTO,break_xytable
    ENDIF

    IF event.type EQ 0 THEN BEGIN ;For end-of-line insertion events
       xcell=event.x & ycell=event.y
       oldvalue=imstr.xyvals(xcell,ycell)

       WIDGET_CONTROL,imstr.xytable,GET_VALUE=a
       imstr.xyvals=FLOAT(a)

       IF xcell NE 4 THEN BEGIN
         CASE xcell OF
         0:BEGIN
          imstr.xyroi(0,0)=imstr.xyvals(0,0) & imstr.xyroi(0,1)=imstr.xyvals(0,1)
          imstr.xyroi(1,0)=imstr.xyvals(0,2) & imstr.xyroi(1,1)=imstr.xyvals(0,3)
          roi=imstr.xyroi
           END
         1:BEGIN
          IF imstr.whichimage NE 'IMAGE' AND imstr.roilist(xcell) EQ 'Set ion ROI' THEN BEGIN
              junk=DIALOG_MESSAGE("Changing of Inversion ROI can be performed in image view only")
              break
          ENDIF

          IF N_ELEMENTS(imstr.phase) GT 1 AND imstr.roilist(xcell) EQ 'Set Inversion ROI' THEN BEGIN
              answer=DIALOG_MESSAGE("Changing Inversion ROI will delete all existing inversion data. Do you wish to continue?",/QUESTION)
              IF answer EQ 'Yes' THEN BEGIN
                 WIDGET_CONTROL,/HOURGLASS
                 imstr1  =imagetostructure(imstr,'PHASE',[-1])
                 imstr2  =imagetostructure(imstr1,'CUMULPHASE',[-1])
                 imstr3  =imagetostructure(imstr2,'NORMALAMPL',[-1])
                 newimstr=imagetostructure(imstr3,'AMPL',[-1])
                 WIDGET_CONTROL,event.handler, SET_UVALUE=newimstr, /NO_COPY
                 WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
              ENDIF ELSE BEGIN
                 imstr.xyvals(xcell,ycell)=oldvalue
                 WIDGET_CONTROL,imstr.xytable,SET_VALUE=imstr.xyvals
                 GOTO,break_xytable
              ENDELSE
          ENDIF

          imstr.phroi(0,0)=imstr.xyvals(1,0) & imstr.phroi(0,1)=imstr.xyvals(1,1)
          imstr.phroi(1,0)=imstr.xyvals(1,2) & imstr.phroi(1,1)=imstr.xyvals(1,3)
          roi=imstr.phroi
           END
         2:BEGIN
          imstr.inroi(0,0)=imstr.xyvals(2,0) & imstr.inroi(0,1)=imstr.xyvals(2,1)
          imstr.inroi(1,0)=imstr.xyvals(2,2) & imstr.inroi(1,1)=imstr.xyvals(2,3)
          roi=imstr.inroi
           END
         3:BEGIN
          imstr.ghroi(0,0)=imstr.xyvals(3,0) & imstr.ghroi(0,1)=imstr.xyvals(3,1)
          imstr.ghroi(1,0)=imstr.xyvals(3,2) & imstr.ghroi(1,1)=imstr.xyvals(3,3)
          roi=imstr.ghroi
           END
           
         ENDCASE
         set_roi,imstr,roi,xcell
       ENDIF ELSE BEGIN
         imstr.limindex=1
         WIDGET_CONTROL,imstr.limtype,SET_DROPLIST_SELECT=imstr.limindex
       ENDELSE

    END
    break_xytable:
    WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
END
'VImage Window':BEGIN
    WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY

    IF N_ELEMENTS(imstr.image) LE 1 THEN GOTO,break_fw

    IF (event.x) LT 0 OR (event.x GE imstr.xsize) THEN GOTO,break_fw
    IF (event.y) LT 0 OR (event.y GE imstr.ysize) THEN GOTO,break_fw

    WIDGET_CONTROL, imstr.vimageID, GET_DRAW_VIEW=drawview

    ;Find and display the cursor position & value in pixels and time & distance
       px=event.x*imstr.factor
       py=event.y*imstr.factor

       IF imstr.whichimage EQ 'IMAGE' THEN BEGIN
        zvalue=STRING(imstr.image(px,py))
       ENDIF ELSE BEGIN
        IF px GE imstr.phroi(0,0) AND px LE imstr.phroi(0,1) AND $
         py GE imstr.phroi(1,0) AND py LE imstr.phroi(1,1) THEN BEGIN
          CASE imstr.whichimage OF
                'PHASE':zvalue=STRING(imstr.phase(px-imstr.phroi(0,0),py-imstr.phroi(1,0)))
                'CUMULPHASE':zvalue=STRING(imstr.cumulphase(px-imstr.phroi(0,0),py-imstr.phroi(1,0)))
                'AMPL':zvalue=STRING(imstr.ampl(px-imstr.phroi(0,0),py-imstr.phroi(1,0)))
                'NORMALAMPL':zvalue=STRING(imstr.normalampl(px-imstr.phroi(0,0),py-imstr.phroi(1,0)))
          ENDCASE
        ENDIF ELSE zvalue=''
       ENDELSE
       xyz=STRING(px,FORMAT='(I)')+STRING(py,FORMAT='(I)')+ zvalue

       WIDGET_CONTROL,imstr.pixID,SET_VALUE='          [x,y,z]='+xyz

       tt=STRCOMPRESS(STRING(time(imstr,px,py,dst)))
       dst=STRCOMPRESS(STRING(dst))

       WIDGET_CONTROL,imstr.dtID,$
       SET_VALUE="                    d (mic) ="+dst(0)+"   t (ns) ="+tt(0)

       x=FIX(event.x*imstr.factor)
       y=FIX(event.y*imstr.factor)
       ;z=imstr.image(x,y)

       WIDGET_CONTROL, imstr.vimageID, SET_DRAW_VIEW=drawview

    ;If right button is pressed then print (x,y,z)
       IF event.press EQ 4 THEN BEGIN
         PRINT,x,y,imstr.image(x,y)
       ENDIF

    ;If left button is pressed then define ROI
       IF event.press EQ 1 THEN BEGIN ;start drag on button press
         imstr.drag = 1
         imstr.roi = [[x,y],[0,0]]
         imstr.oldphroi=imstr.phroi ;saves phase roi information in case changes need to be undone
       ENDIF

    ;If roi is still being dragged then plot the window
      IF imstr.drag EQ 1 THEN BEGIN
       imstr.roi = [[imstr.roi(0),imstr.roi(1)],[x,y]]
;       ;Makes sure the max and min are in the right order
;        xmin=MIN(imstr.roi(0,*)) & ymin=MIN(imstr.roi(1,*))
;        xmax=MAX(imstr.roi(0,*)) & ymax=MAX(imstr.roi(1,*))
;        imstr.roi=[[xmin,ymin],[xmax,ymax]]
       ;Dynamically shows ROI during re-sizing; also dynamically updates ROI table
       WIDGET_CONTROL,imstr.roitype,GET_VALUE=roiindex
       set_roi,imstr,imstr.roi,roiindex        
      ENDIF
      
    ;Only respond on left button release
       IF event.release EQ 1 AND imstr.drag EQ 1 THEN BEGIN

       imstr.roi = [[imstr.roi(0),imstr.roi(1)],[x,y]] ;[[x0,y0],[x1,y1]]
       imstr.drag = 0

       ;Makes sure the max and min are in the right order
        xmin=MIN(imstr.roi(0,*)) & ymin=MIN(imstr.roi(1,*))
        xmax=MAX(imstr.roi(0,*)) & ymax=MAX(imstr.roi(1,*))
       imstr.roi=[[xmin,ymin],[xmax,ymax]]

       WIDGET_CONTROL,imstr.roitype,GET_VALUE=roiindex

          ;Checks to see if phase ROI really does want to be changed, since this will
          ;erase all stored phase information.
          IF imstr.whichimage NE 'IMAGE' AND imstr.roilist(roiindex) EQ 'Set Inversion ROI' THEN BEGIN
            junk=DIALOG_MESSAGE("Changing Inversion ROI can be performed in image view only")
            set_roi,imstr,imstr.oldphroi,roiindex
            GOTO,break_fw
          ENDIF

          IF N_ELEMENTS(imstr.phase) GT 1 AND imstr.roilist(roiindex) EQ 'Set Inversion ROI' THEN BEGIN
            answer=DIALOG_MESSAGE("Changing Inversion ROI will delete all existing inversion data. Do you wish to continue?",/QUESTION)
            IF answer EQ 'Yes' THEN BEGIN
              WIDGET_CONTROL,/HOURGLASS
              imstr1  =imagetostructure(imstr,'PHASE',[-1])
              imstr2  =imagetostructure(imstr1,'CUMULPHASE',[-1])
              imstr3  =imagetostructure(imstr2,'NORMALAMPL',[-1])
              newimstr=imagetostructure(imstr3,'AMPL',[-1])
              WIDGET_CONTROL,event.handler, SET_UVALUE=newimstr, /NO_COPY
              WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
            ENDIF ELSE BEGIN
              set_roi,imstr,imstr.oldphroi,roiindex
              GOTO,break_fw
            ENDELSE
          ENDIF


       set_roi,imstr,imstr.roi,roiindex
       
       ENDIF ;execute only on button release

    break_fw:
    WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
END
'InverseAbel':BEGIN
    WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
    IF imstr.phroi[0,0] EQ -1 THEN BEGIN
       junk=DIALOG_MESSAGE("Define Calculation Region")
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
           RETURN
    ENDIF
    IF imstr.loadedtarget EQ 0 THEN BEGIN
       junk=DIALOG_MESSAGE("Must load target data first")
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
           RETURN
    ENDIF
    
    IF imstr.xyroi[0,0] GE imstr.xyroi[0,1] OR $ ;xmin GE xmax
       imstr.xyroi[1,0] GE imstr.xyroi[1,1] OR $ ;ymin GE ymax
       imstr.phroi[0,0] GE imstr.phroi[0,1] OR $
       imstr.phroi[1,0] GE imstr.phroi[1,1] OR $
       imstr.phroi[0,0] LT 0 OR imstr.phroi[0,1] GE imstr.xsize*imstr.factor OR $
       imstr.phroi[1,0] LT 0 OR imstr.phroi[1,1] GE imstr.ysize*imstr.factor $
       THEN BEGIN
       junk=DIALOG_MESSAGE("Make sure xmin < xmax and ymin < ymax for ROI's")
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
       RETURN
    ENDIF

    WIDGET_CONTROL,/HOURGLASS

    WIDGET_CONTROL,imstr.mag_id,GET_VALUE=micpx  & micpx=FLOAT(micpx(0)) ;mics/pixel
    WIDGET_CONTROL,imstr.mag0_id,GET_VALUE=d0    & d0=FLOAT(d0(0)) ;pixels
    WIDGET_CONTROL,imstr.dx_id,GET_VALUE=dx      & dx=FLOAT(dx(0)) ;pixels
    WIDGET_CONTROL,imstr.dy_id,GET_VALUE=dy      & dy=FLOAT(dy(0)) ;pixels

    WIDGET_CONTROL,imstr.bpoly_id,GET_VALUE=Npoly      & Npoly=FLOAT(Npoly(0))  ;Polynomial order for background
    WIDGET_CONTROL,imstr.skew_id,GET_VALUE=skew  & skew=FLOAT(skew(0))  ;skew for density
    WIDGET_CONTROL,imstr.mix_id,GET_VALUE=Mix    & Mix=FLOAT(Mix(0));(mics) mix thickness

    WIDGET_CONTROL,imstr.slit_id,GET_VALUE=slit       &  slit=FLOAT(slit(0)) ;mic
    WIDGET_CONTROL,imstr.sloty0_id,GET_VALUE=sl_y0    & sl_y0=FLOAT(sl_y0(0));(mics) hohl slot position (y0=0=centered on capsule)
    WIDGET_CONTROL,imstr.slotdy_id,GET_VALUE=sl_dy    & sl_dy=FLOAT(sl_dy(0));(mics) hohl slot width

    WIDGET_CONTROL,imstr.R0_id,GET_VALUE=R0      & R0=FLOAT(R0(0))  ;(mics) outer shell edge position
    WIDGET_CONTROL,imstr.rho_id,GET_VALUE=rhom   & eta=FLOAT(rhom(0))  ;rho_max of shell
    WIDGET_CONTROL,imstr.del_id,GET_VALUE=del    & del=FLOAT(del(0));(mics) shell stdev thickness

;Npoly=4
;skew=-7.
;Mix=10.            ;Thickness of mix layer (in mics). Should be mix=10. for Omega49780

    ;This is already loaded during "Target Description" but needs to be loaded again in case mix is different
    ;tg=loadtarget(XKAPPAFILE=imstr.xrayfile,TDAT=Tgdata,MIX=mix) 
        
;    mur0=mu*rho0/10000. ;convert to inverse mics.
    micpx=micpx*dy  ;ensure that spatial scale takes into account dy bin size
;*************

    ;Rebin roi region
    phxmin=imstr.phroi[0,0]
    phxmax=phxmin+FIX((imstr.phroi[0,1]-phxmin+1)/dx)*dx-1 ;ensure that roi dimensions are
    phymin=imstr.phroi[1,0] 
    phymax=phymin+FIX((imstr.phroi[1,1]-phymin+1)/dy)*dy-1 ;integer multiples of dx & dy
    set_roi,imstr,[[phxmin,phymin],[phxmax,phymax]],1; change phase roi to be an integer multiple of dx
;*************
;  Check that d0 is within the y-range of the Invert ROI. If not, reject.
    IF d0 LT phymin OR d0 GT phymax THEN BEGIN
       junk=DIALOG_MESSAGE("Make sure d=0 is within y-range of Invert ROI")
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
       RETURN
    ENDIF
;*************
    ghxmin=imstr.ghroi[0,0] & ghxmax=ghxmin+FIX((imstr.ghroi[0,1]-ghxmin+1)/dx)*dx-1
    ghymin=imstr.ghroi[1,0] & ghymax=ghymin+FIX((imstr.ghroi[1,1]-ghymin+1)/dy)*dy-1
    IF ghxmin GE 0 AND ghxmax GE 0 AND ghymin GE 0 AND ghymax GE 0 THEN BEGIN
      set_roi,imstr,[[ghxmin,ghymin],[ghxmax,ghymax]],3 ; change phase roi to be an integer multiple of dx
    ENDIF
    
    im=imstr.image[phxmin:phxmax,phymin:phymax]   ;roi subimage
    
    nx=phxmax-phxmin+1 & ny=phymax-phymin+1
    snx=nx/dx & sny=ny/dy
    sub=REBIN(im,snx,sny)   ;rebinned roi image

    gnx=ghxmax-ghxmin+1 & gny=ghymax-ghymin+1
    IF gnx GE 1 AND gny GE 1 THEN BEGIN
      pblock=REBIN(FINDGEN(gnx)+ghxmin,gnx/dx)
      yblock=(REBIN(FINDGEN(gny)+ghymin,gny/dy)-d0)*micpx/dy
      block={p:pblock,y:yblock}
    ENDIF ELSE block=0
    
    px=phxmin+INDGEN(phxmax-phxmin+1)
    py=phymin+INDGEN(phymax-phymin+1)
    tt=time(imstr,px,py,dist)
    tsub=REBIN(tt,snx) ;binned time steps
    psub=REBIN(px,snx) ;binned pixel positions
    dsub=REBIN(dist,sny)
    x0=INTERPOL(FINDGEN(sny),dsub,0.) ;Find the exact position of x=0 in the binned vector

    ;x0=((d0-phymin)/dy) ;(pixels) capsule center; MUST be an integer for abel func to work
;************************
;Setup array of shell parameters with initial guesses
;    a=[eta,R0,del,x0,slit,micpx,mur0,rho0,sl_y0,sl_dy]
;    eta=eta*del ;iterate over rho*R for better stability (5/30/11)
;    rhomaxb=1.0 ;initial density guess for ablated material   -THIS IS STILL BEING TESTED
;    IF imstr.avmufit EQ 1 THEN BEGIN  ;modify initial density guess if fitting mu*rho
;      eta=eta*6. & rhomaxb=rhomaxb*6. 
;    ENDIF
    a=[R0,eta,del,x0,slit,micpx,sl_y0,sl_dy];,rhomaxb]
    
;    tg =LoadTarget(mix=mix) ;Transfer this step into xstreak_fit.pro
;    IF imstr.avmufit EQ 1 THEN fitmurho=1 ELSE fitmurho=0 
    
    dat=xstreak_basis(imstr,sub,tsub,psub,a,BLOCK=block,$
        MIX=mix,Npoly=Npoly,SKEW=skew,FITMURHO=imstr.avmufit);
    ;phase->inten; cumulative phase->bkg; amplitude->projected thickness; normalized ampl-> density

    ;Expand 2D arrays back to original size for viewing
;    inten=dat.inten     & inten=REBIN(inten,nx,ny,/SAMPLE)   
;    bkg=dat.bkg         & bkg=REBIN(bkg,nx,ny,/SAMPLE)       
;    OD=dat.OD           & OD=REBIN(OD,nx,ny,/SAMPLE)
;    density=dat.density & density=REBIN(density,nx,ny,/SAMPLE)
;    mu=dat.mu           & mu=REBIN(mu,nx,ny,/SAMPLE)
;**************************
;    imstr1  =imagetostructure(imstr,'PHASE',[-1])
;    imstr2  =imagetostructure(imstr1,'CUMULPHASE',[-1])
;    imstr3  =imagetostructure(imstr2,'NORMALAMPL',[-1])
;    newimstr=imagetostructure(imstr3,'AMPL',[-1])
;    WIDGET_CONTROL,event.handler, SET_UVALUE=newimstr, /NO_COPY
;    WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
;
;    imstr0  =imagetostructure(imstr,'DAT',dat)
;    imstr1  =imagetostructure(imstr0,'PHASE',inten)
;    imstr2  =imagetostructure(imstr1,'CUMULPHASE',bkg)
;    imstr3  =imagetostructure(imstr2,'AMPL',OD)
;    newimstr=imagetostructure(imstr3,'NORMALAMPL',density) ;mu
    newimstr  =imagetostructure(imstr,'DAT',dat)

    break_InverseAbel:
    WIDGET_CONTROL,event.handler, SET_UVALUE=newimstr, /NO_COPY
END
'ForwardAbel':BEGIN
    WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
    IF imstr.phroi(0,0) EQ -1 THEN BEGIN
       junk=DIALOG_MESSAGE("Define Calculation Region")
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
           RETURN
    ENDIF
    IF imstr.loadedtarget EQ 0 THEN BEGIN
       junk=DIALOG_MESSAGE("Must load target data first")
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
           RETURN
    ENDIF
    
    IF imstr.xyroi(0,0) GE imstr.xyroi(0,1) OR $ ;xmin GE xmax
       imstr.xyroi(1,0) GE imstr.xyroi(1,1) OR $ ;ymin GE ymax
       imstr.phroi(0,0) GE imstr.phroi(0,1) OR $
       imstr.phroi(1,0) GE imstr.phroi(1,1) OR $
       imstr.phroi(0,0) LT 0 OR imstr.phroi(0,1) GE imstr.xsize*imstr.factor OR $
       imstr.phroi(1,0) LT 0 OR imstr.phroi(1,1) GE imstr.ysize*imstr.factor $
       THEN BEGIN
       junk=DIALOG_MESSAGE("Make sure xmin < xmax and ymin < ymax for ROI's")
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
           RETURN
    ENDIF

    WIDGET_CONTROL,/HOURGLASS

    WIDGET_CONTROL,imstr.mag_id,GET_VALUE=micpx  & micpx=FLOAT(micpx(0)) ;mics/pixel
    WIDGET_CONTROL,imstr.mag0_id,GET_VALUE=d0    & d0=FLOAT(d0(0)) ;pixels
    WIDGET_CONTROL,imstr.dx_id,GET_VALUE=dx      & dx=FLOAT(dx(0)) ;pixels
    WIDGET_CONTROL,imstr.dy_id,GET_VALUE=dy      & dy=FLOAT(dy(0)) ;pixels

    WIDGET_CONTROL,imstr.bpoly_id,GET_VALUE=Npoly      & Npoly=FLOAT(Npoly(0))  ;Polynomial order for background
    WIDGET_CONTROL,imstr.skew_id,GET_VALUE=skew   & skew=FLOAT(skew(0))  ;skew for density
    WIDGET_CONTROL,imstr.mix_id,GET_VALUE=Mix    & Mix=FLOAT(Mix(0));(mics) mix thickness

    WIDGET_CONTROL,imstr.slit_id,GET_VALUE=slit       &  slit=FLOAT(slit(0)) ;mic
    WIDGET_CONTROL,imstr.sloty0_id,GET_VALUE=sl_y0    & sl_y0=FLOAT(sl_y0(0));(mics) hohl slot position (y0=0=centered on capsule)
    WIDGET_CONTROL,imstr.slotdy_id,GET_VALUE=sl_dy    & sl_dy=FLOAT(sl_dy(0));(mics) hohl slot width

    WIDGET_CONTROL,imstr.R0_id,GET_VALUE=R0      & R0=FLOAT(R0(0))  ;(mics) outer shell edge position
    WIDGET_CONTROL,imstr.rho_id,GET_VALUE=rhom   & eta=FLOAT(rhom(0))  ;rho_max of shell
    WIDGET_CONTROL,imstr.del_id,GET_VALUE=del    & del=FLOAT(del(0));(mics) shell stdev thickness

;Npoly=4
;skew=-7.
;Mix=10.            ;Thickness of mix layer (in mics). Should be mix=10. for Omega49780

    ;This is already loaded during "Target Description" but needs to be loaded again in case mix is different
    ;tg=loadtarget(XKAPPAFILE=imstr.xrayfile,TDAT=Tgdata,MIX=mix) 
        
;    mur0=mu*rho0/10000. ;convert to inverse mics.
    micpx=micpx*dy  ;ensure that spatial scale takes into account dy bin size
;*************

    ;Rebin roi region
    phxmin=imstr.phroi(0,0)
    phxmax=phxmin+FIX((imstr.phroi(0,1)-phxmin+1)/dx)*dx-1 ;ensure that roi dimensions are
    phymin=imstr.phroi(1,0) 
    phymax=phymin+FIX((imstr.phroi(1,1)-phymin+1)/dy)*dy-1 ;integer multiples of dx & dy
    set_roi,imstr,[[phxmin,phymin],[phxmax,phymax]],1; change phase roi to be an integer multiple of dx
;*************
;  Check that d0 is within the y-range of the Invert ROI. If not, reject.
    IF d0 LT phymin OR d0 GT phymax THEN BEGIN
       junk=DIALOG_MESSAGE("Make sure d=0 is within y-range of Invert ROI")
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
       RETURN
    ENDIF
;*************
    ghxmin=imstr.ghroi(0,0) & ghxmax=ghxmin+FIX((imstr.ghroi(0,1)-ghxmin+1)/dx)*dx-1
    ghymin=imstr.ghroi(1,0) & ghymax=ghymin+FIX((imstr.ghroi(1,1)-ghymin+1)/dy)*dy-1
    IF ghxmin GE 0 AND ghxmax GE 0 AND ghymin GE 0 AND ghymax GE 0 THEN BEGIN
      set_roi,imstr,[[ghxmin,ghymin],[ghxmax,ghymax]],3 ; change phase roi to be an integer multiple of dx
    ENDIF
    
    im=imstr.image(phxmin:phxmax,phymin:phymax)   ;roi subimage
    
    nx=phxmax-phxmin+1 & ny=phymax-phymin+1
    snx=nx/dx & sny=ny/dy
    sub=REBIN(im,snx,sny)   ;rebinned roi image

    gnx=ghxmax-ghxmin+1 & gny=ghymax-ghymin+1
    IF gnx GE 1 AND gny GE 1 THEN BEGIN
      pblock=REBIN(FINDGEN(gnx)+ghxmin,gnx/dx)
      yblock=(REBIN(FINDGEN(gny)+ghymin,gny/dy)-d0)*micpx/dy
      block={p:pblock,y:yblock}
    ENDIF ELSE block=0
    
    px=phxmin+INDGEN(phxmax-phxmin+1)
    py=phymin+INDGEN(phymax-phymin+1)
    tt=time(imstr,px,py,dist)
    tsub=REBIN(tt,snx) ;binned time steps
    psub=REBIN(px,snx) ;binned pixel positions
    dsub=REBIN(dist,sny)
    x0=INTERPOL(FINDGEN(sny),dsub,0.) ;Find the exact position of x=0 in the binned vector

    ;x0=((d0-phymin)/dy) ;(pixels) capsule center; MUST be an integer for abel func to work
;************************
;Setup array of shell parameters with initial guesses
;    a=[eta,R0,del,x0,slit,micpx,mur0,rho0,sl_y0,sl_dy]
;    eta=eta*del ;iterate over rho*R for better stability (5/30/11)
    rhomaxb=1.0 ;initial density guess for ablated material   -THIS IS STILL BEING TESTED
;    IF imstr.avmufit EQ 1 THEN BEGIN  ;modify initial density guess if fitting mu*rho
;      eta=eta*6. & rhomaxb=rhomaxb*6. 
;    ENDIF
    a=[R0,eta,del,x0,slit,micpx,sl_y0,sl_dy,rhomaxb]
    
;    tg =LoadTarget(mix=mix) ;Transfer this step into xstreak_fit.pro
;    IF imstr.avmufit EQ 1 THEN fitmurho=1 ELSE fitmurho=0 
    
    dat=xstreak_fit(imstr,sub,tsub,psub,a,BLOCK=block,INIT=1,$
        MIX=mix,Npoly=Npoly,SKEW=skew,FITMURHO=imstr.avmufit);
    ;phase->inten; cumulative phase->bkg; amplitude->projected thickness; normalized ampl-> density
 
    ;Expand 2D arrays back to original size for viewing
    inten=dat.inten     & inten=REBIN(inten,nx,ny,/SAMPLE)   
    bkg=dat.bkg         & bkg=REBIN(bkg,nx,ny,/SAMPLE)       
    OD=dat.OD           & OD=REBIN(OD,nx,ny,/SAMPLE)
    density=dat.density & density=REBIN(density,nx,ny,/SAMPLE)
    mu=dat.mu           & mu=REBIN(mu,nx,ny,/SAMPLE)
;**************************
;    imstr1  =imagetostructure(imstr,'PHASE',[-1])
;    imstr2  =imagetostructure(imstr1,'CUMULPHASE',[-1])
;    imstr3  =imagetostructure(imstr2,'NORMALAMPL',[-1])
;    newimstr=imagetostructure(imstr3,'AMPL',[-1])
;    WIDGET_CONTROL,event.handler, SET_UVALUE=newimstr, /NO_COPY
;    WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
;
;    imstr0  =imagetostructure(imstr,'DAT',dat)
;    imstr1  =imagetostructure(imstr0,'PHASE',inten)
;    imstr2  =imagetostructure(imstr1,'CUMULPHASE',bkg)
;    imstr3  =imagetostructure(imstr2,'AMPL',OD)
;    newimstr=imagetostructure(imstr3,'NORMALAMPL',density) ;mu
    newimstr  =imagetostructure(imstr,'DAT',dat)

    break_ForwardAbel:
    WIDGET_CONTROL,event.handler, SET_UVALUE=newimstr, /NO_COPY
END
'PlotNow':BEGIN
    WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
    IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
       junk=DIALOG_MESSAGE("No image loaded")
       GOTO,break_plotnow
    ENDIF

;***************Check that max > min for all roi's
    IF imstr.xyroi(0,0) GE imstr.xyroi(0,1) OR $ ;xmin GE xmax
       imstr.xyroi(1,0) GE imstr.xyroi(1,1) OR $ ;ymin GE ymax
       imstr.phroi(0,0) GE imstr.phroi(0,1) OR $
       imstr.phroi(1,0) GE imstr.phroi(1,1) OR $
       imstr.inroi(0,0) GE imstr.inroi(0,1) OR $
       imstr.inroi(1,0) GE imstr.inroi(1,1) OR $
       imstr.ghroi(0,0) GE imstr.ghroi(0,1) OR $
       imstr.ghroi(1,0) GE imstr.ghroi(1,1) $       
       THEN BEGIN
       junk=DIALOG_MESSAGE('Make sure xmin < xmax and ymin < ymax for all ROI')
       GOTO,break_plotnow
    ENDIF

    IF imstr.xyroi(0,0) EQ -1 THEN BEGIN
         junk=DIALOG_MESSAGE("Define Plot ROI Region")
         GOTO,break_plotnow
       ENDIF

    ;The following represent the x,y limits for the ROI to
    ;be used as the averaging region.
    xmin=imstr.xyroi(0,0) & ymin=imstr.xyroi(1,0)
    xmax=imstr.xyroi(0,1) & ymax=imstr.xyroi(1,1)

    ;The following are the x, y limits for the phase ROI
    phxmin=imstr.phroi(0,0) & phxmax=imstr.phroi(0,1)
    phymin=imstr.phroi(1,0) & phymax=imstr.phroi(1,1)


       IF imstr.whichimage EQ 'IMAGE' THEN BEGIN
          IF xmin LT 0 OR xmax GE imstr.xsize*imstr.factor OR $
            ymin LT 0 OR ymax GE imstr.ysize*imstr.factor THEN BEGIN
            junk=DIALOG_MESSAGE('Make sure (xmin,xmax,ymin,ymax) for Plot ROI are within image')
            GOTO,break_plotnow
          ENDIF  
         im=imstr.image[xmin:xmax,ymin:ymax]
       ENDIF ELSE BEGIN
         IF xmin GE phxmin AND xmax LE phxmax AND $
          ymin GE phymin AND ymax LE phymax THEN BEGIN
          CASE imstr.whichimage OF
           'PHASE':im=imstr.phase(xmin-phxmin:xmax-phxmin,ymin-phymin:ymax-phymin)
           'CUMULPHASE':im=imstr.cumulphase(xmin-phxmin:xmax-phxmin,ymin-phymin:ymax-phymin)
           'AMPL':im=imstr.ampl(xmin-phxmin:xmax-phxmin,ymin-phymin:ymax-phymin)
           'NORMALAMPL':im=imstr.normalampl(xmin-phxmin:xmax-phxmin,ymin-phymin:ymax-phymin)
          ENDCASE
         ENDIF ELSE BEGIN
          junk=DIALOG_MESSAGE("Plot ROI must be contained within Inversion ROI")
          GOTO,break_plotnow
         ENDELSE
       ENDELSE


    CASE imstr.plotlist(imstr.plotindex) OF
    'Row Profile':BEGIN;Plot row

        ;WIDGET_CONTROL,imstr.vimageID,GET_DRAW_VIEW=pt
       ;vxmin=pt(0)*imstr.factor
       ;vymin=pt(1)*imstr.factor

       ;IF xmin EQ xmax THEN BEGIN
       ; sz=SIZE(imstr.image)
       ; xmin=0 & xmax=sz(1)-1
       ;ENDIF
       avrow=REBIN(FLOAT(im),xmax-xmin+1,1)

       sigrow=FLTARR(xmax-xmin+1)
       FOR ii=0,xmax-xmin DO BEGIN
         sigrow(ii)=STDDEV(FLOAT(im(ii,*)))
       ENDFOR

       IF imstr.ylog EQ 1 THEN BEGIN
         avrow=ALOG10(avrow)
         sigrow=(ALOG10(avrow+sigrow)-ALOG10(avrow-sigrow))/2.
       w=WHERE(~FINITE(avrow)) & IF w(0) NE -1 THEN avrow(w)=1E-1
       w=WHERE(~FINITE(sigrow)) & IF w(0) NE -1 THEN sigrow(w)=1E-1
       ENDIF

       OUTPLOT,time(imstr,xmin+INDGEN(xmax-xmin+1)),avrow,imstr,ERROR=sigrow,$
         HIST=0,XTITL="Time (ns)",YTITL="Intensity",TITL='';,$
;         TITL="X="+STRCOMPRESS(STRING(xmin),/REMOVE_ALL)+$
;          ":"+STRCOMPRESS(STRING(xmax),/REMOVE_ALL)+$
;          ", Y="+STRCOMPRESS(STRING(ymin),/REMOVE_ALL)+$
;          ":"+STRCOMPRESS(STRING(ymax),/REMOVE_ALL)

    END
    'Column Profile':BEGIN;Plot column

       ;WIDGET_CONTROL,imstr.vimageID,GET_DRAW_VIEW=pt

       ;vxmin=pt(0)*imstr.factor
       ;vymin=pt(1)*imstr.factor

       ;IF ymin EQ ymax THEN BEGIN
       ; sz=SIZE(imstr.image)
       ; ymin=0 & ymax=sz(2)-1
       ;ENDIF

       avcol=REFORM(REBIN(FLOAT(im),1,ymax-ymin+1))
       sigcol=FLTARR(ymax-ymin+1)
       FOR ii=0,ymax-ymin DO BEGIN
         sigcol(ii)=STDDEV(FLOAT(im(*,ii)))
       ENDFOR

       IF imstr.ylog EQ 1 THEN BEGIN
         avcol=ALOG10(avcol)
         sigcol=(ALOG10(avcol+sigcol)-ALOG10(avcol-sigcol))/2.
       w=WHERE(~FINITE(avcol)) & IF w(0) NE -1 THEN avcol(w)=1E-1
       w=WHERE(~FINITE(sigcol)) & IF w(0) NE -1 THEN sigcol(w)=1E-1
       ENDIF


       t=time(imstr,xmin+INDGEN(xmax-xmin+1),ymin+INDGEN(ymax-ymin+1),dist)

       OUTPLOT,dist,avcol,imstr,HIST=0,ERROR=sigcol,$
         XTITL="Distance (mics)",YTITL="Intensity",TITL='';,$
;         TITL="X="+STRCOMPRESS(STRING(xmin),/REMOVE_ALL)+$
;          ":"+STRCOMPRESS(STRING(xmax),/REMOVE_ALL)+$
;          ", Y="+STRCOMPRESS(STRING(ymin),/REMOVE_ALL)+$
;          ":"+STRCOMPRESS(STRING(ymax),/REMOVE_ALL)

    END
    'Sweep Rate':BEGIN

        MinPix=imstr.MinPix & MaxPix=imstr.MaxPix
       cff=FLOAT(imstr.tpoly) & scff=FLOAT(imstr.spoly)

        px=FINDGEN(xmax-xmin+1)+xmin

       sweeprate=POLY(px,scff) ;Create the ps/pix vector

       wlo=WHERE(px LT MinPix)   & whi=WHERE(px GT MaxPix)
       IF wlo(0) NE -1 THEN BEGIN
         sweeprate(wlo)=POLY(MinPix,scff)
       ENDIF
       IF whi(0) NE -1 THEN BEGIN
         sweeprate(whi)=POLY(MaxPix,scff)
       ENDIF

;     t=time(imstr,px) & sweeprate=DERIV(px,t)

       OUTPLOT,px,sweeprate,imstr,$
        XTITL="Pixel",YTITL="Sweep Rate (ns/pix)",$
         TITL="Sweep Rate"

       break_sweep:
    END
;    'Time derivative':BEGIN ;Plot time derivative, with time delay given by time delay
;
;       WIDGET_CONTROL,imstr.vimageID,GET_DRAW_VIEW=pt
;
;       WIDGET_CONTROL,imstr.delay_id,GET_VALUE=tau_ps
;       WIDGET_CONTROL,imstr.a_id,GET_VALUE=a
;
;       tau_pix=FLOAT(tau_ps(0))/FLOAT(a(0))
;       IF tau_pix LT 2 THEN tau_pix=2.
;       IF tau_pix GE 3 THEN BEGIN
;         FOR ii=0,ymax-ymin DO im(*,ii)=SMOOTH(im(*,ii),tau_pix)
;       END
;
;       FOR ii=0,ymax-ymin DO BEGIN
;         im(*,ii)=SHIFT(im(*,ii),tau_pix/2.)-SHIFT(im(*,ii),-tau_pix/2.)
;         im(0:tau_pix/2.,ii)=im(tau_pix/2.+1)
;         im(xmax-xmin-tau_pix/2.:xmax-xmin,ii)=im(xmax-xmin-tau_pix/2.-1,ii)
;       ENDFOR
;
;       ac=REBIN(FLOAT(im),xmax-xmin+1,1)
;
;       sigac=FLTARR(xmax-xmin+1)
;       FOR ii=0,xmax-xmin DO BEGIN
;         sigac(ii)=STDDEV(FLOAT(im(ii,*)))
;       ENDFOR
;
;       OUTPLOT,time(imstr,xmin+INDGEN(xmax-xmin+1)),ac,imstr,$
;         ERROR=sigac,hist=1,$
;         XTITL='Time (ps)',YTITL='F(t-!4s!3/2)-F(t+!4s!3/2)',$
;         TITL="X="+STRCOMPRESS(STRING(xmin),/REMOVE_ALL)+$
;          ":"+STRCOMPRESS(STRING(xmax),/REMOVE_ALL)+$
;          ", Y="+STRCOMPRESS(STRING(ymin),/REMOVE_ALL)+$
;          ":"+STRCOMPRESS(STRING(ymax),/REMOVE_ALL)
;
;    END
;    'Edge Find':BEGIN
;       WIDGET_CONTROL,imstr.vimageID,GET_DRAW_VIEW=pt
;
;       WIDGET_CONTROL,imstr.delay_id,GET_VALUE=tau_ps
;       WIDGET_CONTROL,imstr.a_id,GET_VALUE=a
;
;       tau_pix=FLOAT(tau_ps(0))/FLOAT(a(0))
;       IF tau_pix LT 2 THEN tau_pix=2.
;       IF tau_pix GE 3 THEN BEGIN
;         FOR ii=0,ymax-ymin DO im(*,ii)=SMOOTH(im(*,ii),tau_pix)
;       END
;
;       edge=FLTARR(ymax-ymin+1)
;       err=FLTARR(ymax-ymin+1)
;       FOR ii=0,ymax-ymin DO BEGIN
;         im(*,ii)=ABS(SHIFT(im(*,ii),tau_pix/2.)-SHIFT(im(*,ii),-tau_pix/2.)  )
;         im(0:tau_pix/2.,ii)=im(tau_pix/2.+1)
;         im(xmax-xmin-tau_pix/2.:xmax-xmin,ii)=im(xmax-xmin-tau_pix/2.-1,ii)
;         ;edge(ii)=TOTAL(im(*,ii)*(FINDGEN(xmax-xmin+1)+xmin))/TOTAL(im(*,ii))
;         edge(ii)=MEAN(WHERE(im(*,ii) EQ MAX(im(*,ii))))+xmin
;         ;err(ii)=
;       ENDFOR
;
;       tt=time(imstr,edge,ymin+INDGEN(ymax-ymin+1),dist)
;       OUTPLOT,dist,edge,imstr,HIST=0,$
;         ;PSY=3,$
;         YTITL='Pixel #',XTITL='Distance (mic)',$
;         TITL="X="+STRCOMPRESS(STRING(xmin),/REMOVE_ALL)+$
;          ":"+STRCOMPRESS(STRING(xmax),/REMOVE_ALL)+$
;          ", Y="+STRCOMPRESS(STRING(ymin),/REMOVE_ALL)+$
;          ":"+STRCOMPRESS(STRING(ymax),/REMOVE_ALL)
;    END
    'Average Position':BEGIN
       IF imstr.phroi(0,0) EQ -1 THEN BEGIN
         junk=DIALOG_MESSAGE("Define Inversion Region")
         GOTO,break_position
       ENDIF
       IF SIZE(imstr.dat,/TYPE) NE 8 THEN BEGIN
        junk=DIALOG_MESSAGE("No Inversion data exists")
        GOTO,break_position
       ENDIF
       
       WIDGET_CONTROL,/HOURGLASS

       phxmin=imstr.phroi(0,0) & phxmax=imstr.phroi(0,1)
       phymin=imstr.phroi(1,0) & phymax=imstr.phroi(1,1)

       IF xmin GE phxmin AND xmax LE phxmax THEN BEGIN
         dat=imstr.dat
         tmin=time(imstr,xmin) & tmax=time(imstr,xmax)
         w=WHERE(dat.t GE tmin AND dat.t LE tmax)
         IF w(0) NE -1 THEN $
          OUTPLOT,dat.t(w),dat.R0(w),imstr,psy=4,$
            XTITL="Time (ns)",YTITL="Average Position (!9m!3m)",TITL=''
       ENDIF ELSE junk=DIALOG_MESSAGE("Plot ROI x-values must be contained within Phase ROI x-values")

       break_position:
    END
    'Velocity':BEGIN
       IF imstr.phroi(0,0) EQ -1 THEN BEGIN
         junk=DIALOG_MESSAGE("Define Inversion Region")
         GOTO,break_velocity
       ENDIF
       IF SIZE(imstr.dat,/TYPE) NE 8 THEN BEGIN
        junk=DIALOG_MESSAGE("No Inversion data exists")
        GOTO,break_velocity
       ENDIF
       
       WIDGET_CONTROL,/HOURGLASS

       phxmin=imstr.phroi(0,0) & phxmax=imstr.phroi(0,1)
       phymin=imstr.phroi(1,0) & phymax=imstr.phroi(1,1)

       IF xmin GE phxmin AND xmax LE phxmax THEN BEGIN
         dat=imstr.dat
         tmin=time(imstr,xmin) & tmax=time(imstr,xmax)
         w=WHERE(dat.t GE tmin AND dat.t LE tmax)
         IF w(0) NE -1 THEN $
          OUTPLOT,dat.t(w),dat.v(w),imstr,psy=4,$
            XTITL="Time (ns)",YTITL="Velocity (!9m!3m/ns)",TITL=''
       ENDIF ELSE junk=DIALOG_MESSAGE("Plot ROI x-values must be contained within Phase ROI x-values")

       break_velocity:
    END
    'Max Density':BEGIN ;Plot phase versus x
       IF imstr.phroi(0,0) EQ -1 THEN BEGIN
         junk=DIALOG_MESSAGE("Define Inversion Region")
         GOTO,break_maxdensity
       ENDIF
       IF SIZE(imstr.dat,/TYPE) NE 8 THEN BEGIN
        junk=DIALOG_MESSAGE("No Inversion data exists")
        GOTO,break_maxdensity
       ENDIF
       
       WIDGET_CONTROL,/HOURGLASS

       phxmin=imstr.phroi(0,0) & phxmax=imstr.phroi(0,1)
       phymin=imstr.phroi(1,0) & phymax=imstr.phroi(1,1)

       IF xmin GE phxmin AND xmax LE phxmax THEN BEGIN
         dat=imstr.dat
         tmin=time(imstr,xmin) & tmax=time(imstr,xmax)
         w=WHERE(dat.t GE tmin AND dat.t LE tmax)
         IF w(0) NE -1 THEN $
          OUTPLOT,dat.t(w),dat.maxdensity(w),imstr,psy=4,$
            XTITL="Time (ns)",YTITL="!9r!3!dmax!n (g/cm!u3!n)",TITL=''
       ENDIF ELSE junk=DIALOG_MESSAGE("Plot ROI x-values must be contained within Phase ROI x-values")

       break_maxdensity:

    END
    'Thickness':BEGIN ;Plot thickness of shell
       IF imstr.phroi(0,0) EQ -1 THEN BEGIN
         junk=DIALOG_MESSAGE("Define Inversion Region")
         GOTO,break_thickness
       ENDIF
       IF SIZE(imstr.dat,/TYPE) NE 8 THEN BEGIN
        junk=DIALOG_MESSAGE("No Inversion data exists")
        GOTO,break_thickness
       ENDIF
       
       WIDGET_CONTROL,/HOURGLASS

       phxmin=imstr.phroi(0,0) & phxmax=imstr.phroi(0,1)
       phymin=imstr.phroi(1,0) & phymax=imstr.phroi(1,1)

       IF xmin GE phxmin AND xmax LE phxmax THEN BEGIN
         dat=imstr.dat
         tmin=time(imstr,xmin) & tmax=time(imstr,xmax)
         w=WHERE(dat.t GE tmin AND dat.t LE tmax)
         IF w(0) NE -1 THEN $
          OUTPLOT,dat.t(w),dat.del(w),imstr,psy=4,$
            XTITL="Time (ns)",YTITL="Thickness (!9m!3m)",TITL=''
       ENDIF ELSE junk=DIALOG_MESSAGE("Plot ROI x-values must be contained within Phase ROI x-values")

       break_thickness:
    END
    'Areal Density':BEGIN
       IF imstr.phroi(0,0) EQ -1 THEN BEGIN
         junk=DIALOG_MESSAGE("Define Inversion Region")
         GOTO,break_rhoR
       ENDIF
       IF SIZE(imstr.dat,/TYPE) NE 8 THEN BEGIN
        junk=DIALOG_MESSAGE("No Inversion data exists")
        GOTO,break_rhoR
       ENDIF
       
       WIDGET_CONTROL,/HOURGLASS

       phxmin=imstr.phroi(0,0) & phxmax=imstr.phroi(0,1)
       phymin=imstr.phroi(1,0) & phymax=imstr.phroi(1,1)

       IF xmin GE phxmin AND xmax LE phxmax THEN BEGIN
         dat=imstr.dat
         tmin=time(imstr,xmin) & tmax=time(imstr,xmax)
         w=WHERE(dat.t GE tmin AND dat.t LE tmax)
         IF w(0) NE -1 THEN $
          OUTPLOT,dat.t(w),dat.rhoR(w),imstr,psy=4,$
            XTITL="Time (ns)",YTITL="!9r!3R (mg/cm!u2!n)",TITL=''
       ENDIF ELSE junk=DIALOG_MESSAGE("Plot ROI x-values must be contained within Phase ROI x-values")

       break_rhoR:

    END
    'Mass':BEGIN
       IF imstr.phroi(0,0) EQ -1 THEN BEGIN
         junk=DIALOG_MESSAGE("Define Inversion Region")
         GOTO,break_mass
       ENDIF
       IF SIZE(imstr.dat,/TYPE) NE 8 THEN BEGIN
        junk=DIALOG_MESSAGE("No Inversion data exists")
        GOTO,break_mass
       ENDIF
       
       WIDGET_CONTROL,/HOURGLASS

       phxmin=imstr.phroi(0,0) & phxmax=imstr.phroi(0,1)
       phymin=imstr.phroi(1,0) & phymax=imstr.phroi(1,1)

       IF xmin GE phxmin AND xmax LE phxmax THEN BEGIN
         dat=imstr.dat
         tmin=time(imstr,xmin) & tmax=time(imstr,xmax)
         w=WHERE(dat.t GE tmin AND dat.t LE tmax)
         IF w(0) NE -1 THEN $
          OUTPLOT,dat.t(w),dat.mass(w),imstr,psy=4,$
            XTITL="Time (ns)",YTITL="Mass (!9m!3g)",TITL=''
       ENDIF ELSE junk=DIALOG_MESSAGE("Plot ROI x-values must be contained within Phase ROI x-values")

       break_mass:

    END
    'Mass/Initial Mass':BEGIN
       IF imstr.phroi[0,0] EQ -1 THEN BEGIN
         junk=DIALOG_MESSAGE("Define Inversion Region")
         GOTO,break_mass_initialmass
       ENDIF
       IF SIZE(imstr.dat,/TYPE) NE 8 THEN BEGIN
        junk=DIALOG_MESSAGE("No Inversion data exists")
        GOTO,break_mass_initialmass
       ENDIF
       
       WIDGET_CONTROL,/HOURGLASS

       phxmin=imstr.phroi[0,0] & phxmax=imstr.phroi[0,1]
       phymin=imstr.phroi[1,0] & phymax=imstr.phroi[1,1]

       IF xmin GE phxmin AND xmax LE phxmax THEN BEGIN
         dat=imstr.dat
         tmin=time(imstr,xmin) & tmax=time(imstr,xmax)
         w=WHERE(dat.t GE tmin AND dat.t LE tmax)
         IF w(0) NE -1 THEN $
          OUTPLOT,dat.t(w),dat.mass(w)/MAX(dat.target.cmass),imstr,psy=4,$
            XTITL="Time (ns)",YTITL="Mass/Initial Mass",TITL=''
       ENDIF ELSE junk=DIALOG_MESSAGE("Plot ROI x-values must be contained within Phase ROI x-values")

       break_mass_initialmass:

    END

    ELSE:

    ENDCASE

    break_plotnow:
    WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
END

'Plus':BEGIN
WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
    IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
       junk=DIALOG_MESSAGE("No image loaded")
       GOTO,break_plus
    ENDIF
IF imstr.factor GT 1 THEN BEGIN
    image_size=SIZE(imstr.image)
    imstr.factor=imstr.factor/2.

    imstr.xsize=image_size(1)/imstr.factor
    imstr.ysize=image_size(2)/imstr.factor

    WIDGET_CONTROL,imstr.vimageID,DRAW_XSIZE=image_size(1)/imstr.factor
    WIDGET_CONTROL,imstr.vimageID,DRAW_YSIZE=image_size(2)/imstr.factor
ENDIF

wshow_images,imstr,'IMAGE'

set_roi,imstr,imstr.xyroi,0
set_roi,imstr,imstr.phroi,1
set_roi,imstr,imstr.inroi,2
set_roi,imstr,imstr.ghroi,3

break_plus:
WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
END

'Minus':BEGIN
WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
    IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
       junk=DIALOG_MESSAGE("No image loaded")
       GOTO,break_minus
    ENDIF
IF imstr.factor LT 30 THEN BEGIN
    image_size=SIZE(imstr.image)
    imstr.factor=imstr.factor*2.

    imstr.xsize=image_size(1)/imstr.factor
    imstr.ysize=image_size(2)/imstr.factor

    WIDGET_CONTROL,imstr.vimageID,DRAW_XSIZE=image_size(1)/imstr.factor
    WIDGET_CONTROL,imstr.vimageID,DRAW_YSIZE=image_size(2)/imstr.factor
ENDIF

wshow_images,imstr,'IMAGE'

set_roi,imstr,imstr.xyroi,0
set_roi,imstr,imstr.phroi,1
set_roi,imstr,imstr.inroi,2
set_roi,imstr,imstr.ghroi,3

break_minus:
WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
END
'Dattype':BEGIN
    WIDGET_CONTROL,event.handler,GET_UVALUE=imstr, /NO_COPY
    ;imstr.datlist=['Data','Reference','Grid','Bkg only']
    imstr.whichdat=imstr.datlist[event.value]
    break_dattype:
    WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
END
'No bkg':BEGIN
  WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
    imstr.bkgstate=event.select
  WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
END
;'Refbut':BEGIN
;  WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
;    imstr.refstate=event.select
;  WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
;END
;"Bkg only":BEGIN
;  WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
;    imstr.bkgonly=event.select
;  WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
;END
;'Negphase':BEGIN
;  WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
;    oldstate=imstr.negstate
;    IF event.select NE oldstate THEN BEGIN
;       IF N_ELEMENTS(imstr.cumulphase) GT 1 THEN $
;         imstr.cumulphase= -imstr.cumulphase
;
;       imstr.negstate=event.select
;    ENDIF
;  WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
;END
'ClrRoi':BEGIN
  WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
    oldstate=imstr.roiclr
    IF event.select NE oldstate THEN BEGIN
       imstr.roiclr=event.select
    ENDIF
  WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
END
'Ylog':BEGIN
  WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
    oldstate=imstr.ylog
    IF event.select NE oldstate THEN BEGIN
       imstr.ylog=event.select
    ENDIF
  WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
END
'AvMuFit':BEGIN
  WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
    oldstate=imstr.avmufit
    IF event.select NE oldstate THEN BEGIN
       imstr.avmufit=event.select
    ENDIF
  WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
END
'Contour':BEGIN
  WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
    
    IF SIZE(imstr.dat,/TYPE) NE 8 THEN BEGIN
      junk=DIALOG_MESSAGE("No Inversion data exists")
      GOTO,break_contour
    ENDIF
    
    oldstate=imstr.Rcont
    IF event.select NE oldstate THEN BEGIN
       imstr.Rcont=event.select
    ENDIF 
    IF imstr.Rcont EQ 1 AND N_ELEMENTS(imstr.dat) GT 0 THEN BEGIN
      WIDGET_CONTROL,imstr.mag_id,GET_VALUE=micpx  & micpx=FLOAT(micpx(0)) ;mics/pixel
      WIDGET_CONTROL,imstr.mag0_id,GET_VALUE=d0    & d0=FLOAT(d0(0)) ;pixels
      x=[imstr.dat.pix] & y_hi=d0+[imstr.dat.R0]/micpx & y_lo=d0-[imstr.dat.R0]/micpx
      x=[x,[REVERSE(x)]] & y=[y_hi,REVERSE(y_lo)]
      imstr.oimpolycont -> SetProperty, DATA=TRANSPOSE([[x],[y]])/imstr.factor,COLOR=[255,255,255],LINESTYLE=0
      imstr.oimwindow -> Draw,imstr.oimview
    ENDIF ELSE BEGIN
     imstr.oimpolycont -> SetProperty, DATA=[[0,0],[0,0]]/imstr.factor,$
      COLOR=[255,255,255],LINESTYLE=0
     imstr.oimwindow -> Draw,imstr.oimview
    ENDELSE
  break_contour:
  WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
END
'Plottype':BEGIN
WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
;   IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
;     junk=DIALOG_MESSAGE("No image loaded")
;     WIDGET_CONTROL,event.id,SET_DROPLIST_SELECT=0
;     GOTO,break_plottype
;   ENDIF

CASE imstr.plotlist(imstr.plotindex) OF
;    'Time derivative':BEGIN;autocorrelate
;       WIDGET_CONTROL,imstr.a_id,GET_VALUE=a
;       IF FLOAT(a(0)) EQ 0 THEN BEGIN
;         junk=DIALOG_MESSAGE("Time scale required")
;         WIDGET_CONTROL,event.id,SET_DROPLIST_SELECT=0
;         imstr.plotindex=0
;         GOTO,break_plottype
;       ENDIF
;    END
    'Freq spectrum' :BEGIN;freq spectrum
       IF imstr.phroi(0) EQ -1 THEN BEGIN
         junk=DIALOG_MESSAGE("Define an ROI first")
         WIDGET_CONTROL,event.id,SET_DROPLIST_SELECT=0
         imstr.plotindex=0
         GOTO,break_plottype
       ENDIF
    END
    'Phase Image':BEGIN;phase image
       IF imstr.phroi(0) EQ -1 THEN BEGIN
         junk=DIALOG_MESSAGE("Define an ROI first")
         WIDGET_CONTROL,event.id,SET_DROPLIST_SELECT=0
         imstr.plotindex=0
         GOTO,break_plottype
       ENDIF
    END
    'Phase Unwrap':BEGIN;phase unwrap
       IF imstr.phroi(0) EQ -1 THEN BEGIN
         junk=DIALOG_MESSAGE("Define a Phase ROI first")
         WIDGET_CONTROL,event.id,SET_DROPLIST_SELECT=0
         imstr.plotindex=0
         GOTO,break_plottype
       ENDIF
    END
    'Reflectivity':BEGIN
       IF imstr.inroi(0) EQ -1 THEN BEGIN
         junk=DIALOG_MESSAGE("Define phase and inten. ROI first")
         WIDGET_CONTROL,event.id,SET_DROPLIST_SELECT=0
         imstr.plotindex=0
         GOTO,break_plottype
       ENDIF
    END
    ELSE:
ENDCASE

imstr.plotindex=event.index

break_plottype:
WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
END
'Overtype':BEGIN
    WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
    IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
       junk=DIALOG_MESSAGE("No image loaded")
       WIDGET_CONTROL,event.id,SET_DROPLIST_SELECT=0
       GOTO,break_overtype
    ENDIF

    imstr.overindex=event.index

    break_overtype:
    WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
END
'Limtype':BEGIN
    WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY

    imstr.limindex=event.index

    break_limtype:
    WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
END
'Freetype':BEGIN
    WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
;    WIDGET_CONTROL,imstr.freetype_id,GET_VALUE=free  ;Which parameters to vary [R0,rho,del] - 1=vary,0=fix
;    print,free
    break_freetype:
    WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
END

ELSE: junk=DIALOG_MESSAGE("Program Error: Event User Value Not Found")
ENDCASE



END
