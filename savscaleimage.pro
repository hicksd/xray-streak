PRO savscaleimage,file,im,t,d

TVLCT,r,g,b,/GET
;im=omegahdf("c:\hicks\omega\data\010911-13\asbo1_s24557.hdf")
;im=read_tiff("c:\windows\desktop\test.tif")

sz=SIZE(im)
nx=sz(1) & ny=sz(2)

;d=FINDGEN(nx)
;t=FINDGEN(ny)

;DEVICE,DECOMPOSED=0
!P.FONT=1
SET_PLOT,'ps'
;!X.MARGIN=[3.5,1.] & !Y.MARGIN=[3.,1.]    
thk=5. & charsz=0.8
DEVICE,ENCAPSULATED=1,XSIZE=8,YSIZE=6,/CM,FILENAME=file,/COLOR,BITS_PER_PIXEL=16,/HELVETICA;,COLORS=255
TVLCT,r,g,b

aspect=FLOAT(ny)/FLOAT(nx) 
IF aspect LT 1 THEN  pos=[0.2,0.2,0.95,aspect*(0.95-0.2)] $
          ELSE pos=[0.2,0.2,(0.95-0.2)/aspect+0.2,0.95]

PLOT,[t(0),t(N_ELEMENTS(t)-1)]/1.e0,[d(0),d(N_ELEMENTS(d)-1)],$
    /NODATA,TICKLEN=-0.02,/XSTYLE,/YSTYLE,CHARSIZE=charsz,$
   XTITLE="!3Time (ns)",YTITLE="Distance (!Mmm)",COLOR=0,XTHICK=thk,YTHICK=thk,CHARTHICK=thk,$
    ;POSITION=pos ;maintains original pixel size (for 2D images)
    ;POSITION=[0.22,0.2,0.75,0.9] ;image aspect ratio fits a set rectangular image. (for streaks)
;    POSITION=[0.22,0.22,0.85,0.85] ;image aspect ratio fits a set rectangular image. (for streaks);ConA-long paper
    POSITION=[0.19,0.19,0.9,0.9] ;image aspect ratio fits a set rectangular image. (for streaks);Hyades ConA-long paper
    ;POSITION=[0.22,0.12,0.95,FLOAT(ny)/FLOAT(nx)*(0.95-0.12)+0.12] ;MAINTAINS image aspect ratio?????

;Get size of plot window in device pixels.
PX = !X.WINDOW * !D.X_VSIZE
PY = !Y.WINDOW * !D.Y_VSIZE
;Desired size of image in pixels.
SX = PX[1] - PX[0] + 1
SY = PY[1] - PY[0] + 1

;Display the image with its lower-left corner at
;the origin of the plot window and with its size
;scaled to fit the plot window.
TVSCL,im,XSIZE=SX,YSIZE=SY, PX[0], PY[0];,/ORDER

DEVICE,/CLOSE
SET_PLOT,'win'
!P.FONT=-1
!X.MARGIN=[10,3] & !Y.MARGIN=[4,2] 
RETURN
END