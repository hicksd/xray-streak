FUNCTION subtime,px,cff,scff,MinPix,MaxPix
;Calculates the time at a given pixel position taking into account that
;regions outside MinPix and MaxPix should have a constant sweep rate; otherwise
;high order polynomial fits create problems outside the region where they were fit.

	t=POLY(px,cff)

	wlo=WHERE(px LT MinPix)	& whi=WHERE(px GT MaxPix)
	IF wlo(0) NE -1 THEN BEGIN
		smin=POLY(MinPix,scff) & tmin=POLY(MinPix,cff)
		t(wlo)=tmin+smin*(px(wlo)-MinPix)
	ENDIF
	IF whi(0) NE -1 THEN BEGIN
		smax=POLY(MaxPix,scff) & tmax=POLY(MaxPix,cff)
		t(whi)=tmax+smax*(px(whi)-MaxPix)
	ENDIF

RETURN,t
END

FUNCTION time,imstr,px,py,d
;Converts pixel positions into time and distance

    px=FLOAT(px)

;    WIDGET_CONTROL,imstr.a_id,GET_VALUE=a
;    WIDGET_CONTROL,imstr.b_id,GET_VALUE=b
;    WIDGET_CONTROL,imstr.c_id,GET_VALUE=c
    WIDGET_CONTROL,imstr.mag_id,GET_VALUE=mag
    WIDGET_CONTROL,imstr.mag0_id,GET_VALUE=d0
    WIDGET_CONTROL,imstr.fidpx_id,GET_VALUE=fidpx
    WIDGET_CONTROL,imstr.fidtt_id,GET_VALUE=fidtt
    MinPix=imstr.MinPix & MaxPix=imstr.MaxPix
	cff=FLOAT(imstr.tpoly) & scff=FLOAT(imstr.spoly)

    mag=FLOAT(mag(0)) & d0=FLOAT(d0(0))
    fidpx=FLOAT(fidpx(0)) & fidtt=FLOAT(fidtt(0))

;    K = fidtt - (a*fidpx + b*fidpx^2 + c*fidpx^3)
;    t = a*px + b*px^2 + c*px^3 + K

	K = fidtt - subtime(fidpx,cff,scff,MinPix,MaxPix)
    t = subtime(px,cff,scff,MinPix,MaxPix) + K



    IF N_ELEMENTS(py) GT 0 THEN BEGIN
       py=FLOAT(py)
       d=mag*(py-d0)
    ENDIF
RETURN,t
END