PRO savetif,event,bin=bin

       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
     IF N_ELEMENTS(imstr.image) LE 1 THEN BEGIN
         junk=DIALOG_MESSAGE("No image loaded")
         GOTO,break_save_tif
       ENDIF

       file=DIALOG_PICKFILE(PATH=imstr.pathsave,FILE="*.tif",$
         FILTER="*.tif",/WRITE,GET_PATH=path,GROUP=event.handler)

       IF file EQ '' THEN GOTO,break_save_tif

       imstr.pathsave=path
     OPENU,unit,imstr.pathfile,/GET_LUN
         PRINTF,unit,imstr.pathopen
         PRINTF,unit,imstr.pathsave
     CLOSE,unit
     FREE_LUN,unit

       TVLCT,r,g,b,/GET
;print,imstr.factor,imstr.xsize,imstr.ysize


       iname=WHERE(TAG_NAMES(imstr) EQ imstr.whichimage)
       im=imstr.(iname(0))

    ;*********************************
    IF imstr.roiclr EQ 1 THEN BEGIN
        xmin=imstr.inroi(0,0) & ymin=imstr.inroi(1,0)
        xmax=imstr.inroi(0,1) & ymax=imstr.inroi(1,1)
        IF xmax GT 0 AND xmin GT 0 AND xmax GT xmin AND $
           ymax GT 0 AND ymin GT 0 AND ymax GT ymin THEN BEGIN
        cmin=MIN(im(xmin:xmax,ymin:ymax))
        cmax=MAX(im(xmin:xmax,ymin:ymax))
        ENDIF ELSE BEGIN
           PRINT,"ROI is incompatible for color typing"
          cmin=MIN(im) & cmax=MAX(im)
        ENDELSE
    ENDIF ELSE BEGIN
        cmin=MIN(im) & cmax=MAX(im)
    ENDELSE
    ;**********************************
       szim=SIZE(im)
       IF KEYWORD_SET(bin) THEN BEGIN
        im=FLOAT(CONGRID(im,szim(1)/imstr.factor,szim(2)/imstr.factor))
       ENDIF ELSE BEGIN
        im=FLOAT(CONGRID(im,szim(1),szim(2)))
       ENDELSE
       
;       print,size(im)
      IF imstr.logplot EQ 1 THEN norm=(ALOG(im-cmin+1))/ALOG(cmax-cmin+2) $
      ELSE norm=(im-cmin)/(cmax-cmin)

        lo=WHERE(norm LE 0)
        hi=WHERE(norm GE 1)
        IF lo(0) NE -1 THEN norm(lo)=0
        IF hi(0) NE -1 THEN norm(hi)=1

       im8=BYTE(norm*255.)
       ;im=UINT(norm*65535.)
       
;       WRITE_TIFF,file,REVERSE(im8,2),red=r,green=g,blue=b
;       WRITE_TIFF,file,im8,red=r,green=g,blue=b

       IF KEYWORD_SET(bin) THEN BEGIN
        WRITE_TIFF,file,im8,red=r,green=g,blue=b ;easily seen with standard viewers
       ENDIF ELSE BEGIN
        WRITE_TIFF,file,im,/FLOAT   ;preserves original values
       ENDELSE       

;       im16=(norm*(LONG(2)^16-1))
;       WRITE_TIFF,file,REVERSE(im16,2),/SHORT;,red=r,green=g,blue=b

       break_save_tif:
       WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
       
RETURN
END