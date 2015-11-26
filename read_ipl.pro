
FUNCTION read_ipl,file

nbytes=128
junk=BYTARR(nbytes) ;initial 128 byte stream of junk present in some files
version=BYTARR(4) ;version # to indicate start of real file
width=LONG(1) ;width of data array
height=LONG(1) ;height of data array

OPENR,unit,file,/GET_LUN
 READU,unit,version
 offset=0
 IF STRING(version) NE '3.5a' THEN BEGIN ;test for start of real data
 	offset=128							;by checking for version #
 	POINT_LUN,unit,offset
 	READU,unit,version
	IF STRING(version) NE '3.5a' THEN BEGIN
		PRINT,"Can't find start of file"
		CLOSE,unit
		FREE_LUN,unit
		RETURN,-1
	ENDIF
 ENDIF

 POINT_LUN,unit,offset+6 & READU,unit,width
 POINT_LUN,unit,offset+10 & READU,unit,height
 POINT_LUN,unit,offset+2120

 width=SWAP_ENDIAN(width)	;data comes from a Mac so
 height=SWAP_ENDIAN(height)	;byte order needs to be swapped

 data=INTARR(width,height)
 READU,unit,data
 data=SWAP_ENDIAN(data)
CLOSE,unit
FREE_LUN,unit

;PRINT,STRING(version)
;PRINT,width,height

RETURN,data
END

PRO viewimage,image,file,ct=ct

file=dialog_pickfile(PATH="c:\hicks\vulcan\",GET_PATH=path,/MUST_EXIST)
IF file EQ '' THEN file="c:\hicks\25-001222-f-06.ipl"

data=read_ipl(file)
image=data

sdata=SIZE(data,/STRUCTURE)
width=sdata.dimensions(0) & height=sdata.dimensions(1)

factor=2
data=REBIN(data,width/factor,height/factor)

;range=MAX(data)-MIN(data)
;data=(data-MIN(data))/FLOAT(range)*256.

DEVICE,DECOMPOSED=0 ;allows truecolor displays to show color

IF KEYWORD_SET(ct) THEN BEGIN
LOADCT,2
ENDIF

WINDOW,/FREE,XSIZE=width/factor,YSIZE=height/factor,$
	TITLE=file
TVSCL,data

RETURN
END


PRO TEST
;n=10000
;PLOT,[0,N],[0,1],/NODATA
;FOR i=0,n DO BEGIN
;OPLOT,[i,i],[0,1],COLOR=i
;ENDFOR

image0=DIST(100,100)
image=BYTARR(3,100,100)
image(0,*,*)=image0
image(1,*,*)=image0
image(2,*,*)=image0

TV,image,true=1


RETURN
END

PRO writeastiff,im,im8

TVLCT,r,g,b,/GET

im=FLOAT(im)
im8=BYTE((im-MIN(im))/(MAX(im)-MIN(im))*255.)

WRITE_TIFF,"c:\windows\desktop\junk.tif",$
	REVERSE(im8,2),red=r,green=g,blue=b

RETURN
END