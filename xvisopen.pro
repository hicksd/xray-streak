PRO xvisopen,event,bkg=bkg,flatfield=flatfield

       WIDGET_CONTROL, event.handler, GET_UVALUE=imstr, /NO_COPY
       WIDGET_CONTROL,/HOURGLASS
       filename=DIALOG_PICKFILE(PATH=imstr.pathopen,GET_PATH=path,GROUP=event.handler)

       IF filename EQ '' THEN BEGIN
         WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
         RETURN
       ENDIF

       imstr.pathopen=path

     OPENU,unit,imstr.pathfile,/GET_LUN
         PRINTF,unit,imstr.pathopen
         PRINTF,unit,imstr.pathsave
     CLOSE,unit
     FREE_LUN,unit

       dot=RSTRPOS(filename,'.')
       IF dot EQ -1 THEN BEGIN
         junk=DIALOG_MESSAGE("File must have an extension")
         WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
         RETURN
       ENDIF
       ext=STRMID(filename,dot,STRLEN(filename))

       CASE STRUPCASE(ext) OF
         '.HDF':BEGIN
                ;imstr.datlist=['Data','Reference','Grid','Bkg only']
                CASE imstr.whichdat OF
                'Data':BEGIN
                 ref=0 & grid=0 & bkgonly=0
                END
                'Reference':BEGIN
                 ref=1 & grid=0 & bkgonly=0
                END
                'Grid':BEGIN
                 ref=0 & grid=1 & bkgonly=0
                END
                'Bkg only':BEGIN
                 ref=0 & grid=0 & bkgonly=1
                END
                ENDCASE
                
                image=omegahdf(filename,REF=ref,GRID=grid,$
                              BKGONLY=bkgonly,NOBKG=imstr.bkgstate,FILEATTR=fileattr)
                ;Incorporate any time calibration information in the hdf
                sz=SIZE(fileattr,/TYPE) ;If this is a structure then it contains useful info
                IF sz EQ 8 THEN BEGIN
                 result=DIALOG_MESSAGE("Time calibration is available. Do you wish to include?",/QUESTION)
                 IF result EQ 'Yes' THEN BEGIN
                   WIDGET_CONTROL,imstr.fidpx_id,SET_VALUE=0.0
                   WIDGET_CONTROL,imstr.fidtt_id,SET_VALUE=fileattr.t0;*1000.
                   imstr.tpoly=[0.,fileattr.a];*1000.
                   imstr.spoly=[fileattr.a,0.];*1000.
                 ENDIF
                ENDIF
               
               END
         '.IPL':image=read_ipl(filename)
         '.NU':image=read_nu(filename)
         '.TIF':BEGIN
                  image=read_tiff(filename)
                  sz=SIZE(image,/DIMENSIONS)
                  IF sz(0) EQ 3 THEN image=REFORM(REBIN(FLOAT(image),1,sz(1),sz(2)))
                END
         '.TIFF':BEGIN
                  image=read_tiff(filename)
                  sz=SIZE(image,/DIMENSIONS)
                  IF sz(0) EQ 3 THEN image=REFORM(REBIN(FLOAT(image),1,sz(1),sz(2)))
                END
         '.IMG':image=read_img(filename)
         '.BMP':image=read_bmp(filename)
         '.H5':image=NIF_h5_read(filename,bin=1)
         '.RAW':image=read_raw(filename)
         '.SAV':RESTORE,filename,/VERBOSE
         '.SPE':image=read_spe(filename)
         '.PNG':BEGIN
                 image=READ_PNG(filename)
                 sz=SIZE(image,/DIMENSIONS)
                 IF sz(0) EQ 3 THEN image=REFORM(REBIN(FLOAT(image),1,sz(1),sz(2)))
                END
         '.JPG':BEGIN
                 READ_JPEG,filename,image
                 sz=SIZE(image,/DIMENSIONS)
                 IF sz(0) EQ 3 THEN image=REFORM(REBIN(FLOAT(image),1,sz(1),sz(2)))
                END
         ELSE:BEGIN
          junk=DIALOG_MESSAGE("Readable file types: *.hdf, *.ipl, *.tif, *.nu, *.img, *.bmp,,*.raw,*.h5,*.spe")
          WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
          RETURN
         END
       ENDCASE

       IF image(0) EQ -999 THEN BEGIN
         WIDGET_CONTROL,event.handler, SET_UVALUE=imstr, /NO_COPY
         RETURN
       ENDIF

       IF NOT(KEYWORD_SET(bkg)) AND NOT(KEYWORD_SET(flatfield)) THEN $
       imstr.file=STRMID(filename,RSTRPOS(filename,'\')+1,STRLEN(filename))

       image=FLOAT(image)
       image_size=SIZE(image)

    PRINT,filename,image_size(1),image_size(2)

;    IF image_size(1) GT 2048 OR image_size(2) GT 2048 THEN BEGIN
;        nxs=STRING(image_size(1),FORMAT='(G6.4)') & nys=STRING(image_size(2),FORMAT='(G6.4)')
;        resp=DIALOG_MESSAGE("Large image: nx="+nxs+" ny="+nys+" Do you wish to bin by 2?",/QUESTION)
;        IF resp EQ 'Yes' THEN BEGIN
;            image=CONGRID(image,image_size(1)/2,image_size(2)/2)
;            image_size=SIZE(image)
;        ENDIF
;    ENDIF
    ;*********************Background image subtraction
       IF KEYWORD_SET(bkg) THEN BEGIN
           bkg=image & bkg_size=image_size & orig=imstr.image & orig_size=SIZE(orig)
           IF bkg_size(1) EQ orig_size(1) AND bkg_size(2) EQ orig_size(2) THEN BEGIN
            image=orig-bkg
           ENDIF ELSE BEGIN
        junk=DIALOG_MESSAGE("Dimensions of bkg image do not match that of original image")
        image=orig & image_size=orig_size
           ENDELSE
       ENDIF
    ;*********************Flatfield
       IF KEYWORD_SET(flatfield) THEN BEGIN
           ff=image & ff_size=image_size & orig=imstr.image & orig_size=SIZE(orig)
           IF ff_size(1) EQ orig_size(1) AND ff_size(2) EQ orig_size(2) THEN BEGIN
         wz=WHERE(ff EQ 0)
         IF wz(0) NE -1 THEN ff(wz)=1
            image=orig/ff
           ENDIF ELSE BEGIN
        junk=DIALOG_MESSAGE("Dimensions of bkg image do not match that of original image")
        image=orig & image_size=orig_size
           ENDELSE

       ENDIF
    ;*********************
       newimstr=imagetostructure(imstr,'IMAGE',image)

       newimstr.xsize=image_size(1)/newimstr.factor
       newimstr.ysize=image_size(2)/newimstr.factor

       WIDGET_CONTROL,newimstr.vimageID,DRAW_XSIZE=image_size(1)/newimstr.factor
       WIDGET_CONTROL,newimstr.vimageID,DRAW_YSIZE=image_size(2)/newimstr.factor

       wshow_images,newimstr,'IMAGE'
       newimstr.whichimage='IMAGE'

       newimstr.lastrot=0
       newimstr.lastshear=0

       WIDGET_CONTROL,newimstr.fname_id,SET_VALUE='File: '+newimstr.file
       WIDGET_CONTROL,event.handler, SET_UVALUE=newimstr, /NO_COPY
       ;This is the end of the section placing the new image into the
       ;storage structure


RETURN
END