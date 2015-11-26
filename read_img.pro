FUNCTION read_img,filename

header=INTARR(7)
OPENR,unit,filename,/GET_LUN
READU,unit,header

offset=64+header(1)

POINT_LUN,unit,offset
image=INTARR(header(2),header(3))

READU,unit,image

CLOSE,unit
FREE_LUN,unit

RETURN,image
END