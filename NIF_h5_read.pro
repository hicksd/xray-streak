FUNCTION NIF_h5_read,file,bin=bin

a=h5_browser(file,/DIALOG_READ)

IF SIZE(a,/type) NE 8 THEN BEGIN
    dat=-999
    GOTO,jump
ENDIF

dat=a._data

sz=SIZE(dat)
;Rebin the data if necessary.
IF KEYWORD_SET(bin) THEN BEGIN
    IF bin NE 1 THEN dat=CONGRID(dat,sz(1)/bin,sz(2)/bin)
ENDIF

jump:
RETURN,dat
END





FUNCTION NIF_h5_read_old,file,bin=bin,data=data

IF N_ELEMENTS(file) EQ 0 THEN file='c:\hicks\nif\visar_031108_test1new.h5'
group='/streaker1'

IF NOT(FILE_TEST(file)) THEN RETURN,0

;Available data sets (request is case sensitive):
;Background_1,Background_2,FullReference_1,FullReference_2
;Leg1_1,Leg1_2,Leg2_1,Leg2_2

IF KEYWORD_SET(data) THEN BEGIN
    datname=STRLOWCASE(data)
    CASE datname OF
       'background_1': dataname='Background_1'
       'background_2': dataname='Background_2'
       'fullreference_1':dataname='FullReference_1'
       'fullreference_2':dataname='FullReference_2'
       'leg1_1':dataname='Leg1_1'
       'leg1_2':dataname='Leg1_2'
       'leg2_1':dataname='Leg2_1'
       'leg2_2':dataname='Leg2_2'
    ELSE:dataname='FullReference_1'
    ENDCASE
ENDIF ELSE dataname='FullReference_1'

;Print the file, group, and data names as a check
PRINT,"File Name: ",file
PRINT,"Group Name: ",group
PRINT,"Dataset: ",dataname

;Open the h5 file and extract the dataset
file_id=h5f_open(file)
 group_id=h5g_open(file_id,group)
  dat_id=h5d_open(group_id,dataname)

   dat=h5d_read(dat_id)

  h5d_close,dat_id
 h5g_close,group_id
h5f_close,file_id

sz=SIZE(dat)

;Rebin the data if necessary.
IF KEYWORD_SET(bin) THEN BEGIN
    IF bin EQ 1 THEN dat=REBIN(dat,sz(1)/4,sz(2)/4) $
    ELSE dat=CONGRID(dat,sz(1)/bin,sz(2)/bin)
ENDIF

RETURN,dat
END