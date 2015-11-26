FUNCTION read_nu, file

comment=BYTARR(100) & exptime=0US & theader=0UL
orig_col=0US & orig_row=0US
ncol=0US & nrow=0US
bin_col=0US & bin_row=0US
gain=0US & read_rate=BYTARR(20) & im_type=BYTARR(10)


OPENR,unit,file,/GET_LUN

 POINT_LUN,unit,0 	& READU,unit,comment & comment = SWAP_ENDIAN(comment)
 POINT_LUN,unit,100 & READU,unit,exptime & exptime = SWAP_ENDIAN(exptime)
 POINT_LUN,unit,112	& READU,unit,theader & theader = SWAP_ENDIAN(theader)
 POINT_LUN,unit,116 & READU,unit,orig_col & orig_col = SWAP_ENDIAN(orig_col)
 POINT_LUN,unit,118 & READU,unit,orig_row & orig_row = SWAP_ENDIAN(orig_row)
 POINT_LUN,unit,120 & READU,unit,nrow & nrow = SWAP_ENDIAN(nrow)
 POINT_LUN,unit,122 & READU,unit,ncol & ncol = SWAP_ENDIAN(ncol)
 POINT_LUN,unit,124 & READU,unit,bin_col & bin_col = SWAP_ENDIAN(bin_col)
 POINT_LUN,unit,126 & READU,unit,bin_row & bin_row = SWAP_ENDIAN(bin_row)
 POINT_LUN,unit,128 & READU,unit,gain & gain = SWAP_ENDIAN(gain)
 POINT_LUN,unit,130 & READU,unit,read_rate & read_rate = SWAP_ENDIAN(read_rate)
 POINT_LUN,unit,150 & READU,unit,im_type & im_type = SWAP_ENDIAN(im_type)

 data=INTARR(ncol,nrow)

 POINT_LUN,unit,160 & READU,unit,data & data=SWAP_ENDIAN(data)

CLOSE,unit
FREE_LUN,unit

data=REVERSE(data,2)

RETURN,data
END

PRO viewimage,image,file,path=path,ct=ct

IF NOT(KEYWORD_SET(path)) THEN path="c:\hicks\janus\data\"
file=dialog_pickfile(PATH=path,GET_PATH=path,/MUST_EXIST)
IF file EQ '' THEN RETURN

data=read_nu(file)
image=data

sdata=SIZE(data,/STRUCTURE)
width=sdata.dimensions(0) & height=sdata.dimensions(1)

factor=2.
data=CONGRID(data,width/factor,height/factor,/INTERP)

DEVICE,DECOMPOSED=0 ;allows truecolor displays to show color

IF KEYWORD_SET(ct) THEN BEGIN
LOADCT,2
ENDIF

WINDOW,/FREE,XSIZE=width/factor,YSIZE=height/factor,$
	TITLE=file
TVSCL,data;,/ORDER

RETURN
END


PRO writeasjpg

COMMON COLORS, R_orig, G_orig, B_orig, R_curr, G_curr, B_curr


viewimage,image,file
livect,image

image=REVERSE(BYTE((FLOAT(image)-MIN(image))/MAX(image)*256),2)

WRITE_IMAGE, "c:\windows\desktop\image.jpg",'JPEG',image,$
	R_curr, G_curr, B_curr

;plot,r_curr
;oplot,g_curr
;oplot,b_curr,linestyle=2

RETURN
END

PRO tvimage,data=data

TVSCL,data,/ORDER

RETURN
END

PRO livect,image

XLOADCT,UPDATECALLBACK='tvimage',UPDATECBDATA=image

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


