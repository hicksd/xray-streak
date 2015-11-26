PRO get_camera_data,file,datnames,camnames,alldata

line='' & ncams=0 & ndats=0
OPENR,unit,file,/GET_LUN
READF,unit,line
READF,unit,line

;READF,unit,ncams & READF,unit,ndats ;ncams=# of cameras (settings); ndats=# of entries per camera

;cam_names=STRARR(ncams)
;dat_names=STRARR(ndats)
;alldata=STRARR(ncams,ndats)

;aline=STRARR(ncams+1)

;READF,unit,line
;READS,line,aline,FORMAT='('+STRTRIM(STRING(ncams+1))+'A18)'
;cam_names=aline(1:ncams)

;FOR i=0,ndats-1 DO BEGIN
;   READF,unit,line
;   READS,line,aline,FORMAT='('+STRTRIM(STRING(ncams+1))+'A18)'
;   dat_names(i)=aline(0)
;   alldata(*,i)=aline(1:ncams)
;ENDFOR

READF,unit,line
aline=STRSPLIT(line,COUNT=N,/EXTRACT) & ncams=N-1 ;counts the number of camera names
cam_names=aline(1:ncams)

WHILE NOT(EOF(unit)) DO BEGIN
    ON_IOERROR,skip
    READF,unit,line & READS,line,aline,FORMAT='('+STRTRIM(STRING(ncams+1))+'A18)'
    IF SIZE(all_data,/TYPE) NE 0 THEN BEGIN
        all_data=[[all_data],[aline(1:ncams)]]
        dat_names=[dat_names,aline(0)]
    ENDIF ELSE BEGIN
       all_data=aline(1:ncams)
       dat_names=aline(0)
    ENDELSE
    skip:
ENDWHILE

CLOSE,unit
FREE_LUN,unit

datnames=dat_names & camnames=cam_names & alldata=all_data
RETURN
END

;PRO get_camera_data_old,file,dat_names,cam_names,alldata
;
;line='' & ncols=0 & nrows=0
;
;OPENR,unit,file,/GET_LUN
;READF,unit,line
;READF,unit,ncams & READF,unit,ndats ;ncams=# of cameras (settings); ndats=# of entries per camera
;
;cam_names=STRARR(ncams)
;dat_names=STRARR(ndats)
;alldata=STRARR(ncams,ndats)
;
;READF,unit,line
;s=STR_SEP(line,'\')
;
;FOR i=0,ndats-1 DO dat_names(i)=s(i+1)
;
;FOR j=0,ncams-1 DO BEGIN
;    READF,unit,line
;    s=STR_SEP(line,'\')
;    cam_names(j)=s(0)
;    FOR i=0,ndats-1 DO alldata(j,i)=s(i+1)
;ENDFOR
;
;alldata=FLOAT(alldata)
;
;CLOSE,unit
;FREE_LUN,unit
;
;RETURN
;END

PRO savecameradata,select
    junk=DIALOG_MESSAGE("Do you want to overwrite existing SC data file?",/DEFAULT_NO,/QUESTION)

    IF junk EQ 'Yes' THEN BEGIN
        ncams=N_ELEMENTS(select.cam_names)
        ndats=N_ELEMENTS(select.dat_names)
        WIDGET_CONTROL,select.table,GET_VALUE=alldata

        OPENW,unit,select.file,/GET_LUN
             PRINTF,unit,"STREAK CAMERA SETTINGS, "+" Created on "+SYSTIME(0)
             ;PRINTF,unit,ncams
             ;PRINTF,unit,ndats
             PRINTF,unit,' '
             PRINTF,unit,['Camera',select.cam_names],FORMAT='('+STRTRIM(STRING(ncams+1))+'A18)'
             PRINTF,unit,[TRANSPOSE(select.dat_names),STRING(alldata)],FORMAT='('+STRTRIM(STRING(ncams+1))+'A18)'
        CLOSE,unit
        FREE_LUN,unit
        select.newdata=0 ;current data table is now the same as saved file
    ENDIF
RETURN
END

PRO load_new_name_event,event
COMMON sc_info, select

WIDGET_CONTROL, event.id, GET_UVALUE = eventval

IF N_ELEMENTS(eventval) EQ 0 THEN RETURN
CASE eventval OF
    "OK":BEGIN
       WIDGET_CONTROL,event.handler,GET_UVALUE=bname_text
       WIDGET_CONTROL,bname_text,GET_VALUE=newname
       select.newdataname=STRCOMPRESS(newname,/REMOVE_ALL)
       WIDGET_CONTROL,event.handler,/DESTROY
    END
    "CANCEL":BEGIN
       select.newdataname=''
       WIDGET_CONTROL,event.handler,/DESTROY
     END
ELSE:junk=DIALOG_MESSAGE("No match to user value")
ENDCASE

RETURN
END

PRO load_new_name,event
;Creates a modal text widget to accept a new name for the current camera settings

bname = WIDGET_BASE(GROUP_LEADER=event.handler,/MODAL,TAB_MODE=1,UVALUE=0,$
       TITLE="Enter New Camera Name",EVENT_PRO='new_name_event')

bname1 = WIDGET_BASE(bname,/COLUMN)
bname2 = WIDGET_BASE(bname1,/ROW)
bname3 = WIDGET_BASE(bname1,/ROW)

bname_text = WIDGET_TEXT(bname2,/EDITABLE)
bname_ok = WIDGET_BUTTON(bname3,VALUE="OK",UVALUE="OK")
bname_cancel = WIDGET_BUTTON(bname3,VALUE="CANCEL",UVALUE="CANCEL")

WIDGET_CONTROL,bname,/REALIZE
WIDGET_CONTROL,bname,SET_UVALUE=bname_text

XMANAGER, 'load_new_name', bname, event='load_new_name_event', $
       GROUP_LEADER = event.handler

END

PRO camera_event,event
COMMON sc_info, select

WIDGET_CONTROL, event.id, GET_UVALUE = eventval

IF N_ELEMENTS(eventval) EQ 0 THEN RETURN
CASE eventval OF


    "TABLE0":BEGIN

    END
    "RESTORE":BEGIN
       WIDGET_CONTROL,select.table0,SET_VALUE=select.cdata
    END
    "CLEAR":BEGIN
      zero=FLTARR(N_ELEMENTS(select.cdata))
      zero(2)=1. ;mic/pixel fixed to 1.
      zero(4)=2048 ;last pixel for sweep calib
      zero(5)=1 ;linear sweep fixed to 1 ps/pix
      zero=TRANSPOSE(zero)
      WIDGET_CONTROL,select.table0,SET_VALUE=zero

    END
    "OK":BEGIN
;       IF select.top NE 0 THEN $ ;OR (select.bottom-select.top LT 6) THEN $
;         junk=DIALOG_MESSAGE("Make selection from the list of camera names") $
;       ELSE WIDGET_CONTROL, event.handler, /DESTROY
       ;WIDGET_CONTROL,
       IF select.newdata EQ 1 THEN BEGIN
        junk=DIALOG_MESSAGE("Do you wish to save new SC data?",/QUESTION)
        IF junk EQ 'Yes' THEN savecameradata,select
       ENDIF

       WIDGET_CONTROL,select.table0,GET_VALUE=current_data
       select.cdata=FLOAT(current_data)
       select.loadcurrent=1
       WIDGET_CONTROL, event.handler, /DESTROY
    END
    "CANCEL":BEGIN
       IF select.newdata EQ 1 THEN BEGIN
        junk=DIALOG_MESSAGE("Do you wish to save new SC data?",/QUESTION)
        IF junk EQ 'Yes' THEN savecameradata,select
       ENDIF

       select.loadcurrent=0
       WIDGET_CONTROL,event.handler,/DESTROY
    END
    "ADD COLUMN":BEGIN
;    ;s_cells:[ left, top, right, bottom ]
;       s_cells = WIDGET_INFO(select.table, /TABLE_SELECT)
;         ;Insert a new column
;     WIDGET_CONTROL,select.table,INSERT_COLUMNS=1
;     num=N_ELEMENTS(select.cam_names)
;     ;Append the new camera name to the structure and to the table
;
;     select=imagetostructure(select,'CAM_NAMES',new_cam_names)
;     WIDGET_CONTROL,select.table,COLUMN_LABELS=select.cam_names
;     ;Append the new camera data to 'alldata'
;     WIDGET_CONTROL,select.table,GET_VALUE=alldata
;     alldata=[current_data,alldata]
;     WIDGET_CONTROL,select.table,SET_VALUE=alldata
    END
    "DELETE COLUMN":BEGIN
       junk=DIALOG_MESSAGE("Do you want to delete the entire column of data?",/DEFAULT_NO,/QUESTION)

       IF junk EQ 'Yes' THEN BEGIN
        s_cells = WIDGET_INFO(select.table, /TABLE_SELECT)
        num=N_ELEMENTS(select.cam_names)
        camnum=INDGEN(num) & w=WHERE(camnum LT s_cells(0) OR camnum GT s_cells(2))
        new_cam_names=select.cam_names(w)
        select=imagetostructure(select,'CAM_NAMES',new_cam_names)
        WIDGET_CONTROL,select.table,/DELETE_COLUMNS ;alldata is automatically adjusted here.
        select.newdata=1 ;current data table is different to saved file
       ENDIF
    END
    "SAVE CAMERA DATA":BEGIN

       savecameradata,select

    END
    "LOAD CURRENT":BEGIN
    ;s_cells:[ left, top, right, bottom ]
        s_cells = WIDGET_INFO(select.table, /TABLE_SELECT)
        IF s_cells(0) NE s_cells(2) THEN BEGIN
          junk=DIALOG_MESSAGE("Selected cell(s) must be from a single column")
        ENDIF ELSE BEGIN
         WIDGET_CONTROL,select.table0,GET_VALUE=current_data
         current_data=FLOAT(current_data)
         WIDGET_CONTROL,select.table,GET_VALUE=alldata
           current_data(0,s_cells(1):s_cells(3))=alldata(s_cells(0),s_cells(1):s_cells(3))
           WIDGET_CONTROL,select.table0,SET_VALUE=current_data
         ENDELSE
    END
    "CREATE NEW":BEGIN
       WIDGET_CONTROL,select.table0,GET_VALUE=current_data
       current_data=FLOAT(current_data)
       load_new_name,event ;user inputs assigned name for current data
       IF select.newdataname EQ '' THEN BEGIN
         junk=DIALOG_MESSAGE("No camera name specified; cannot create a new data array.")
         GOTO,stopcreate
       ENDIF

       WIDGET_CONTROL,select.table,GET_VALUE=alldata
       ;Insert a new column
       WIDGET_CONTROL,select.table,INSERT_COLUMNS=1
       ;Append the new camera name to the structure and to the table
       new_cam_names=[select.newdataname,select.cam_names]
       select=imagetostructure(select,'CAM_NAMES',new_cam_names)
       WIDGET_CONTROL,select.table,COLUMN_LABELS=select.cam_names
       ;Append the new camera data to 'alldata'
       alldata=[current_data,alldata]
       WIDGET_CONTROL,select.table,SET_VALUE=alldata
       select.newdata=1 ;Data array is different to saved data file

       stopcreate:
    END
    "TABLE":BEGIN
;       IF event.type EQ 4 THEN BEGIN
;         select.left=event.sel_left & select.top=event.sel_top
;         select.right=event.sel_right & select.bottom=event.sel_bottom
;       IF select.top EQ 0 THEN BEGIN
;        WIDGET_CONTROL,select.table,GET_VALUE=alldata
;        WIDGET_CONTROL,select.table0,SET_VALUE=alldata(select.left,*)
;       ENDIF
;       END
    END

ELSE:junk=DIALOG_MESSAGE("No match to user value")
ENDCASE

RETURN
END


FUNCTION camera,imstr,event
COMMON sc_info, select

file_exist=FILE_TEST(imstr.camerafile)
IF file_exist THEN BEGIN
  sc_file=imstr.camerafile 
ENDIF ELSE BEGIN
  junk=DIALOG_MESSAGE("File "+imstr.camerafile+" cannot be found!")
  sc_file=DIALOG_PICKFILE(PATH=imstr.pathopen,GET_PATH=path,$
            GROUP=event.handler,TITLE="Please find the file xstreak_sc_settings.txt")
  IF sc_file EQ '' THEN BEGIN
    junk=DIALOG_MESSAGE("No camera file was selected.")
    GOTO,skip_camera
  ENDIF
  imstr.camerafile=sc_file
ENDELSE 

get_camera_data,sc_file,$
          dat_names,cam_names,alldata

sc = WIDGET_BASE(GROUP_LEADER=event.handler,/MODAL,TAB_MODE=1,UVALUE=0,$
       TITLE="Camera Settings",EVENT_PRO='camera_event')

tsc0 = WIDGET_BASE(sc,/ROW)
tsc1 = WIDGET_BASE(tsc0,/COLUMN)
lsc = WIDGET_BASE(tsc0,/COLUMN)
tsc2 = WIDGET_BASE(tsc0,/COLUMN)

;Compile values for "Current" Table
WIDGET_CONTROL,imstr.rot_id,GET_VALUE=rads
WIDGET_CONTROL,imstr.shr_id,GET_VALUE=shearr
WIDGET_CONTROL,imstr.mag_id,GET_VALUE=mag
MinPix=imstr.MinPix
MaxPix=imstr.MaxPix
cam_tpoly=imstr.tpoly(1:N_ELEMENTS(imstr.tpoly)-1)

current_data=TRANSPOSE([rads,shearr,mag,MinPix,MaxPix,cam_tpoly])

table0_id = WIDGET_TABLE(tsc1,COLUMN_LABELS=['Current'],$ ;
         ROW_LABELS=STRCOMPRESS(dat_names),ALL_EVENTS=0,EDITABLE=1,$
         VALUE=current_data,$
         UVALUE='TABLE0')

table_id = WIDGET_TABLE(tsc2,COLUMN_LABELS=cam_names,$
         ROW_LABELS=STRCOMPRESS(dat_names),/ALL_EVENTS,X_SCROLL_SIZE=5,Y_SCROLL_SIZE=N_ELEMENTS(dat_names),$
         VALUE=alldata,UVALUE='TABLE')

;instruct = WIDGET_LABEL(tsc2,FONT="HELVETICA*BOLD*18",$
;       VALUE="Click on camera name (at top) to select settings")
;
;clabel = WIDGET_LABEL(tsc1,FONT="HELVETICA*BOLD*18",$
;       VALUE="Edit current values")

rsc = WIDGET_BASE(tsc1,/ROW)
restoreID = WIDGET_BUTTON(rsc,VALUE="RESTORE",UVALUE="RESTORE")
clearID = WIDGET_BUTTON(rsc,VALUE="CLEAR",UVALUE="CLEAR")

;loadc_label = WIDGET_LABEL(lsc, VALUE="Loads highlighted cells to current")
loadc_id = WIDGET_BUTTON(lsc,VALUE="<-- LOAD TO CURRENT",UVALUE="LOAD CURRENT")
createnew_id = WIDGET_BUTTON(lsc,VALUE="CREATE AS NEW -->",UVALUE="CREATE NEW")
okID = WIDGET_BUTTON(lsc, VALUE="OK",UVALUE="OK")
cancelID = WIDGET_BUTTON(lsc,VALUE="CANCEL",UVALUE="CANCEL")

bsc = WIDGET_BASE(tsc2,/ROW)
;add_ID = WIDGET_BUTTON(bsc,VALUE="ADD COLUMN",UVALUE="ADD COLUMN")
del_ID = WIDGET_BUTTON(bsc,VALUE="DELETE COLUMN",UVALUE="DELETE COLUMN")
save_ID = WIDGET_BUTTON(bsc,VALUE="SAVE CAMERA DATA",UVALUE="SAVE CAMERA DATA")

WIDGET_CONTROL,sc,/REALIZE

select={left:-1,top:-1,right:-1,bottom:-1,$
       loadcurrent:0,file:sc_file,$
       table0:table0_id,table:table_id,$
       cdata:current_data,cam_names:cam_names,dat_names:dat_names,$
       newdataname:'',newdata:0 }

XMANAGER, 'camera', sc, event='camera_event', $
       GROUP_LEADER = event.handler


IF select.loadcurrent EQ 0 THEN RETURN,-1 ELSE RETURN,REFORM(select.cdata)

skip_camera:
RETURN,-1
END