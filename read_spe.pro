; NAME: read_spe.pro (IDL 5.1)
;
; PURPOSE: Read image data in Princeton's WinView .spe format
;
; CATEGORY: Image data i/o
;
; CALLING SEQUENCE: Result = READ_SPE( [filename] )
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   filename            A string scalar containing the filename.
;                       If undefined, a fileselector box pops up.
;
; KEYWORD PARAMETERS:
; a) Informational output parameters
;     DATETIME          Date and time at which the image was taken (string).
;     NUMFRAMES         Number of frames
;     XDIM,YDIM         Dimensions of the frame(s).
;     XDIMET,YDIMDET    Dimensions of the detector.
;     SOFTWARE_VER      Version of WinView used to save the image (FLOAT).
;     VERSION           Well, something (STRING)
; b) Control parameters
;     QUIET       Set this keyword to suppress screen output during
;                 processing.
;
; OUTPUTS:
;   Result: A two or three-dimensional array containing the frames in
;           the file. If the file contains more than one frame, the
;           third dimension represents different frames. The type of Result
;           corresponds to the data type in the file.
;             Unsigned 16-bit integers are converted to long integers.
;
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS:
;   WinView32 does not seem to store the exposure time correctly.
;
; PROCEDURE:
;   First, the header section is read in with many READU commands.
;   Then, more READUs are used in a loop that reads in the frame(s).
;
; MODIFICATION HISTORY:
;   09.11.2001 Included keyw EXP_SEC (T.W.)
;   Created 08 Feb 1999 by T. Weitkamp, ESRF
;
;-


function read_spe, filename $
                   , DATETIME=datetime $
                   , EXP_SEC=exp_sec $
                   , XDIM=xdim $
                   , YDIM=ydim $
                   , NUMFRAMES=numframes $
                   , SOFTWARE_VER=software_ver $
                   , VERSION=version $
                   , QUIET=quiet

IF N_ELEMENTS(filename) EQ 0 THEN $
  filename = DIALOG_PICKFILE(TITLE="Select WinView file to read")
IF filename EQ '' THEN RETURN, 0

quiet = KEYWORD_SET(quiet)

;; ---------------- Define header data ------------------

dioden=0                        ;int
avgexp=0                        ;int
exposure=0                      ;int
xdimdet=0                       ;unsigned int
mode=0                          ;int
exp_sec=0.0                     ;float
asyavg=0
asyseq=0
ydimdet=0                       ;unsigned
date=STRING(255b+BYTARR(10))
ehour=0
eminute=0
noscan=0
fastacc=0
seconds=0
dettype=0
xdim=0                          ;unsigned
stdiode=0
nanox=0.0
calibdio=fltarr(10)
fastfile=STRING(255b+BYTARR(16))
asynen=0
datatype=0
calibnan=fltarr(10)
backgrndapplied=0
astdiode=0
minblk=0                        ;unsigned
numminblk=0                     ;unsigned
calibpol=dblarr(4)
adcrate=0                       ;unsigned
adctype=0                       ;unsigned
adcresolution=0                 ;unsigned
adcbitadjust=0                  ;unsigned
gain=0                          ;unsigned
exprem=STRING(255b+BYTARR(5,80))
geometric=0                     ;unsigned
xlabel=STRING(255b+BYTARR(16))
cleans=0                        ;unsigned
numskppercln=0                  ;unsigned
califile=STRING(255b+BYTARR(16))
bkgdfile=STRING(255b+BYTARR(16))
srccmp=0
ydim=0                          ;unsigned
scramble=0
lexpos=0L
lnoscan=0L
lavgexp=0L
stripfil=STRING(255b+BYTARR(16))
version=STRING(255b+BYTARR(16))
type=0
flatfieldapplied=0
spare=intarr(8)
kin_trig_mode=0
empty=BYTARR(702)
clkspd_us=0.0
hwaccumflag=0
storesync=0
blemishapplied=0
cosmicapplied=0
cosmictype=0
cosmicthreshold=0.0
numframes=0L
maxintensity=0.0
minintensity=0.0
ylabel=STRING(255b+BYTARR(16))
shuttertype=0                   ;unsigned
shuttercomp=0.0
readoutmode=0                   ;unsigned
windowsize=0                    ;unsigned
clkspd=0                        ;unsigned
interface_type=0                ;unsigned
ioadd1=0L                       ;unsigned
ioadd2=0L                       ;unsigned
ioadd3=0L                       ;unsigned
intlevel=0                      ;unsigned
gpibadd=0                       ;unsigned
controladd=0                    ;unsigned
controllernum=0                 ;unsigned
swmade=0                        ;unsigned
numroi=0
;
roiinfoblk=BYTARR(1632-1512)
;
flatfield=STRING(255b+BYTARR(120))
background=STRING(255b+BYTARR(120))
blemish=STRING(255b+BYTARR(120))
software_ver=0.0
userinfo=BYTARR(1000)
winview_id=0L
;
xcalib=BYTARR(3489-3000)
ycalib=BYTARR(3978-3489)
;
istring=STRING(255b+BYTARR(40))
empty3=BYTARR(80)
lastvalue=0

;
; -------------------- Read header ------------------------
;

OPENR, u, filename, /GET_LUN, /SWAP_IF_BIG_ENDIAN

READU,u,dioden
READU,u,avgexp
READU,u,exposure
READU,u,xdimdet
READU,u,mode
READU,u,exp_sec
READU,u,asyavg
READU,u,asyseq
READU,u,ydimdet
READU,u,date
READU,u,ehour
READU,u,eminute
READU,u,noscan
READU,u,fastacc
READU,u,seconds
READU,u,dettype
READU,u,xdim
READU,u,stdiode
READU,u,nanox
READU,u,calibdio
READU,u,fastfile
READU,u,asynen
READU,u,datatype
READU,u,calibnan
READU,u,backgrndapplied
READU,u,astdiode
READU,u,minblk
READU,u,numminblk
READU,u,calibpol
READU,u,adcrate
READU,u,adctype
READU,u,adcresolution
READU,u,adcbitadjust
READU,u,gain
READU,u,exprem
READU,u,geometric
READU,u,xlabel
READU,u,cleans
READU,u,numskppercln
READU,u,califile
READU,u,bkgdfile
READU,u,srccmp
READU,u,ydim
READU,u,scramble
READU,u,lexpos
READU,u,lnoscan
READU,u,lavgexp
READU,u,stripfil
READU,u,version
READU,u,type
READU,u,flatfieldapplied
READU,u,spare
READU,u,kin_trig_mode
READU,u,empty
READU,u,clkspd_us
READU,u,hwaccumflag
READU,u,storesync
READU,u,blemishapplied
READU,u,cosmicapplied
READU,u,cosmictype
READU,u,cosmicthreshold
READU,u,numframes
READU,u,maxintensity
READU,u,minintensity
READU,u,ylabel
READU,u,shuttertype
READU,u,shuttercomp
READU,u,readoutmode
READU,u,windowsize
READU,u,clkspd
READU,u,interface_type
READU,u,ioadd1
READU,u,ioadd2
READU,u,ioadd3
READU,u,intlevel
READU,u,gpibadd
READU,u,controladd
READU,u,controllernum
READU,u,swmade
READU,u,numroi
;
READU,u,roiinfoblk
;
READU,u,flatfield
READU,u,background
READU,u,blemish
READU,u,software_ver
READU,u,userinfo
READU,u,winview_id
;
READU,u,xcalib
READU,u,ycalib
;
READU,u,istring
READU,u,empty3
READU,u,lastvalue
;
; ----------------- Read Frame(s) --------------------
;
CASE datatype OF
    0: d=FLTARR(xdim,ydim,numframes) ; float
    1: d=LONARR(xdim,ydim,numframes) ; 32-bit integer
    2: d=INTARR(xdim,ydim,numframes) ; signed 16-bit integer
    3: d=INTARR(xdim,ydim,numframes) ; unsigned 16-bit integer
    ELSE: MESSAGE, 'Unknown datatype:'+STRING(datatype)
ENDCASE
oneframe=d[*,*,0]

FOR i=0, numframes-1 DO BEGIN
  IF NOT quiet THEN $
    PRINT, FORM='($,"Frame",i," of",i,"'+STRING(13b)+'")', i+1, numframes
  READU,u,oneframe
  d[*,*,i]=oneframe
ENDFOR
;
; Convert unsigned 16-bit integers to long integers
IF datatype EQ 3 THEN BEGIN
  d=LONG(d)
  negative=WHERE(d LT 0, count)
  IF count NE 0 THEN d[negative] = d[negative] + 65536
ENDIF

datetime = date + " " + STRTRIM(ehour,2) + ":" + STRTRIM(eminute,2)

CLOSE,u
FREE_LUN,u

return, d
end