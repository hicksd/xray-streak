FUNCTION read_raw,filename

Nx=UINT(0) & Ny=UINT(0)

OPENR,unit,filename,/GET_LUN
    READU,unit,Nx & READU,unit,Ny
    image=FLTARR(Nx,Ny)
    READU,unit,image
CLOSE,unit
FREE_LUN,unit

RETURN,image
END