FUNCTION omegahdf,filename,bkg,dat,ref=ref,grid=grid,nobkg=nobkg,bkgonly=bkgonly,fileattr=fileattr
;A simple function to read in an OMEGA hdf file
;and output the raw image data. No file information is
;attached at present
;Allow any number of fileattributes to be returned in a structure

; Open the HDF file:
sd_id = HDF_SD_START(filename)

indexNewAsbo=hdf_sd_nametoindex(sd_id,'IMAGE_DATA')

IF indexNewAsbo NE -1 THEN BEGIN
    ; Access the data SD in the HDF file:
    sds_id_DAT=HDF_SD_SELECT(sd_id, indexNewAsbo)
    ;Get the sd data
    HDF_SD_GETDATA, sds_id_DAT, data
    ; End access to any SD ids:
    HDF_SD_ENDACCESS, sds_id_DAT
    ; Close the file:
    HDF_SD_END, sd_id
    RETURN,data
ENDIF

indexXRF=hdf_sd_nametoindex(sd_id,'pds_image')
IF indexXRF NE -1 THEN BEGIN
    ; Access the data SD in the HDF file:
    sds_id_DAT=HDF_SD_SELECT(sd_id, indexXRF)
    ;Get the sd data
    HDF_SD_GETDATA, sds_id_DAT, data
    ; End access to any SD ids:
    HDF_SD_ENDACCESS, sds_id_DAT
    ; Close the file:
    HDF_SD_END, sd_id
    RETURN,data
ENDIF

indexTVS=hdf_sd_nametoindex(sd_id,'pdv_image')
IF indexTVS NE -1 THEN BEGIN
    ; Access the data SD in the HDF file:
    sds_id_DAT=HDF_SD_SELECT(sd_id, indexTVS)
    ;Get the sd data
    HDF_SD_GETDATA, sds_id_DAT, data
    ; End access to any SD ids:
    HDF_SD_ENDACCESS, sds_id_DAT
    ; Close the file:
    HDF_SD_END, sd_id
    RETURN,data
ENDIF

indexTVS=hdf_sd_nametoindex(sd_id,'fg_image')
IF indexTVS NE -1 THEN BEGIN
    ; Access the data SD in the HDF file:
    sds_id_DAT=HDF_SD_SELECT(sd_id, indexTVS)
    ;Get the sd data
    HDF_SD_GETDATA, sds_id_DAT, data
    ; End access to any SD ids:
    HDF_SD_ENDACCESS, sds_id_DAT
    ; Close the file:
    HDF_SD_END, sd_id
    RETURN,data
ENDIF

indexCID=hdf_sd_nametoindex(sd_id,'cid_foreground')
IF indexCID NE -1 THEN BEGIN
    ; Access the data SD in the HDF file:
    sds_id_DAT=HDF_SD_SELECT(sd_id, indexCID)
    ;Get the sd data
    HDF_SD_GETDATA, sds_id_DAT, data
    ; End access to any SD ids:
    HDF_SD_ENDACCESS, sds_id_DAT
    ;Obtain index for background
    indexCID=hdf_sd_nametoindex(sd_id,'cid_background')
    ; Access the data SD in the HDF file:
    sds_id_DAT=HDF_SD_SELECT(sd_id, indexCID)
    ;Get the sd data
    HDF_SD_GETDATA, sds_id_DAT, bkg
    ; End access to any SD ids:
    HDF_SD_ENDACCESS, sds_id_DAT
    ; Close the file:
    HDF_SD_END, sd_id
    RETURN,data-bkg
ENDIF

; Find the number of SD datasets and attributes
HDF_SD_FILEINFO, sd_id, datasets, attributes

CASE datasets OF
1:BEGIN
    ; Return the index of the 'Streak_array' dataset:
    indexDAT = HDF_SD_NAMETOINDEX(sd_id, 'Streak_array')
    ; Access the data SD in the HDF file:
    sds_id_DAT=HDF_SD_SELECT(sd_id, indexDAT)
    ;Get the sd data
    HDF_SD_GETDATA, sds_id_DAT, data
    ; End access to any SD ids:
    HDF_SD_ENDACCESS, sds_id_DAT
END
2:BEGIN
    ; Return the index of the 'Reference' dataset:
    indexREF = HDF_SD_NAMETOINDEX(sd_id, 'Reference')
    ; Access the reference SD in the HDF file:
    sds_id_ref=HDF_SD_SELECT(sd_id, indexREF)
    ;Get the sd data
    HDF_SD_GETDATA, sds_id_ref, refdat
    ; End access to any SD ids:
    HDF_SD_ENDACCESS, sds_id_ref
    
    ; Return the index of the 'Streak_array' dataset:
    indexDAT = HDF_SD_NAMETOINDEX(sd_id, 'Streak_array')
    ; Access the data SD in the HDF file
    sds_id_dat=HDF_SD_SELECT(sd_id, indexDAT)
    ;Get the sd data
    HDF_SD_GETDATA, sds_id_dat, data
    ; End access to any SD ids:
    HDF_SD_ENDACCESS, sds_id_dat
END
3:BEGIN
    ; Return the index of the 'Grid' dataset:
    indexGrid = HDF_SD_NAMETOINDEX(sd_id, 'Grid')
    ; Access the reference SD in the HDF file:
    sds_id_grid=HDF_SD_SELECT(sd_id, indexGRID)
    ;Get the sd data
    HDF_SD_GETDATA, sds_id_grid, griddat
    ; End access to any SD ids:
    HDF_SD_ENDACCESS, sds_id_grid
    
    ; Return the index of the 'Reference' dataset:
    indexREF = HDF_SD_NAMETOINDEX(sd_id, 'Reference')
    ; Access the reference SD in the HDF file:
    sds_id_ref=HDF_SD_SELECT(sd_id, indexREF)
    ;Get the sd data
    HDF_SD_GETDATA, sds_id_ref, refdat
    ; End access to any SD ids:
    HDF_SD_ENDACCESS, sds_id_ref
    
    ; Return the index of the 'Streak_array' dataset:
    indexDAT = HDF_SD_NAMETOINDEX(sd_id, 'Streak_array')
    ; Access the data SD in the HDF file
    sds_id_dat=HDF_SD_SELECT(sd_id, indexDAT)
    ;Get the sd data
    HDF_SD_GETDATA, sds_id_dat, data
    ; End access to any SD ids:
    HDF_SD_ENDACCESS, sds_id_dat

END
ELSE:
ENDCASE

;*********7/28/09 Check to see if time calibration is available
indexDAT   = HDF_SD_NAMETOINDEX(sd_id, 'Streak_array')
sds_id_DAT = HDF_SD_SELECT(sd_id, indexDAT)
HDF_SD_GETINFO, sds_id_DAT,natts=natts
IF natts GT 20 THEN GOTO,skipatt ;this is a cheap way of discriminating between files with/out cal info
t0_index   = HDF_SD_ATTRFIND(sds_id_DAT, 't0_of_first_px_in_ns') 
a_index    = HDF_SD_ATTRFIND(sds_id_DAT,'deltat_per_px_in_ns')
IF t0_index NE -1 AND a_index NE -1 THEN BEGIN
  HDF_SD_ATTRINFO, sds_id_DAT, t0_index, DATA=t0 
  HDF_SD_ATTRINFO, sds_id_DAT, a_index, DATA=a 
  fileattr={t0:t0,a:a} ;
ENDIF ELSE fileattr=0.
skipatt:
HDF_SD_ENDACCESS, sds_id_DAT
;*******************

; Close the file:
HDF_SD_END, sd_id

sz=SIZE(data,/STRUCTURE)
;Create a zero background if one doesn't yet exist. 
;This allows for compatibility with arrays which do not have backgrounds.
IF sz.n_dimensions EQ 2 THEN bkg=FLTARR(sz.dimensions(0),sz.dimensions(1)) $
 ELSE bkg=data(*,*,1) 

IF KEYWORD_SET(bkgonly) THEN BEGIN
;    bkg=data(*,*,1)
    sub=bkg
    GOTO,jump
ENDIF

IF KEYWORD_SET(ref) AND N_ELEMENTS(refdat) NE 0 THEN BEGIN
;    bkg=data(*,*,1)
    IF KEYWORD_SET(nobkg) THEN sub=refdat ELSE sub=refdat-bkg
    neg=WHERE(sub LT 0 OR sub GE 2L^16-1000)
    IF neg(0) NE -1 THEN sub(neg)=0
    GOTO,jump
ENDIF

IF KEYWORD_SET(grid) AND N_ELEMENTS(griddat) NE 0 THEN BEGIN
;    bkg=data(*,*,1)
    IF KEYWORD_SET(nobkg) THEN sub=griddat ELSE sub=griddat-bkg
    neg=WHERE(sub LT 0 OR sub GE 2L^16-1000)
    IF neg(0) NE -1 THEN sub(neg)=0
    GOTO,jump
ENDIF

;IF N_ELEMENTS(refdat) NE 0 AND KEYWORD_SET(ref) EQ 1 THEN ref=1 ELSE ref=0
;CASE ref OF
;0: BEGIN
;bkg=data(*,*,1) & dat=data(*,*,0)

IF KEYWORD_SET(nobkg) THEN sub=data(*,*,0) ELSE sub=data(*,*,0)-bkg
neg=WHERE(sub LT 0 OR sub GE 2L^16-1000)
IF neg(0) NE -1 THEN sub(neg)=0

;END
;1: BEGIN
;    bkg=data(*,*,1)
;    IF KEYWORD_SET(nobkg) THEN sub=refdat ELSE sub=refdat-bkg
;   neg=WHERE(sub LT 0 OR sub GE 2L^16-1000)
; IF neg(0) NE -1 THEN sub(neg)=0
;END
;ENDCASE


jump:
RETURN, sub

END