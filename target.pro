PRO CreateTargetStructure
;Creates the initial target structure that will contain all future target information.
;Once TgStr is created this procedure should never (in theory) need to be run again.
;Once created, additional targets are saved to and deleted from this structure.

;Array of structures containing info about each target.
MaxN_Region=100 ;Maximum number of layers
MaxN_Mats=50 ;Maximum number of elements per layer
mat0 = REPLICATE({sym:'',A:0.,mu:0.,f:FLTARR(MaxN_Region)},MaxN_Mats)
TgStr0={Name:'',N_Region:0.,N_Mat:0.,En:0.,Rin:0.,dr:0.,$
        d_Region:FLTARR(MaxN_Region),rho_Region:FLTARR(MaxN_Region),mat:mat0} 
TgStr=REPLICATE(TgStr0,100) ;Number of different types of target

;Load a sample first target
i=0
TgStr(i).Name="CHOGe0.9/8.35/784/134"
TgStr(i).N_Region=4 ;Number of layers
TgStr(i).N_Mat=4 ;number of different elements
TgStr(i).En=8.35 ;keV, photon energy to be used
TgStr(i).Rin=784.15 ;mic, inner radius of capsule
TgStr(i).dr=0.5 ;mic, step size for tracking profile through target
TgStr(i).d_Region=[23.1,29.4,12.8,68.3];um, thickness of each region
TgStr(i).rho_Region=[1.02,1.07,1.05,1.02];g/cc, density of each region
;************Element details for each region
;f can be first listed in STOICHIOMETRIC form for convenience; automatically normalized afterwards
mat =[{sym:'C' ,A:12.0107,mu:0.,f:REPLICATE(0.4235,MaxN_Region)},$
      {sym:'H' ,A: 1.0079,mu:0.,f:REPLICATE(0.5715,MaxN_Region)},$
      {sym:'O' ,A:15.9994,mu:0.,f:REPLICATE(0.0050,MaxN_Region)},$
      {sym:'Ge',A:  72.64,mu:0.,f:[0.,0.0062,0.0039,0.,FLTARR(MaxN_Region-4)]}]
;*************
FOR j=0,TgStr(i).N_Mat-1 DO TgStr(i).mat(j)=mat(j)


Last_Target=0 ;Index of which was the last target to have been used 
SAVE,TgStr,Last_Target,FILENAME="c:\Hicks\xstreak\Initial_Target_Structure.sav"

END



FUNCTION savetarget,select
;This inserts current new target information into the Target Structure (TgStr) 
;and saves it to disk. Returns 0 if not saved, 1 if saved successfully.

;MaxN_Region=100 & MaxN_Mats=50
;mat0 = REPLICATE({sym:'',A:0.,mu:0.,f:FLTARR(MaxN_Region)},MaxN_Mats)
;TgStr0={Name:'',N_Region:0.,N_Mat:0.,En:0.,Rin:0.,dr:0.,$
;        d_Region:FLTARR(MaxN_Region),rho_Region:FLTARR(MaxN_Region),mat:mat0} 

;WIDGET_CONTROL,select.targetlist_id,GET_VALUE=Target_Names
WIDGET_CONTROL,select.name_id,GET_VALUE=name & name=STRCOMPRESS(STRING(name),/REMOVE_ALL) 
WIDGET_CONTROL,select.energy_id,GET_VALUE=En 
WIDGET_CONTROL,select.Rmin_id,GET_VALUE=Rin  
WIDGET_CONTROL,select.step_id,GET_VALUE=dr   
WIDGET_CONTROL,select.table_id,GET_VALUE=tabledata 

IF name EQ '' THEN BEGIN
  junk=DIALOG_MESSAGE("Target must have a name")
  RETURN,0
ENDIF

TgStr=select.TgStr

;*******Find place to put this target into TgStr.

;Check if name already exists
exist=STRCMP(REPLICATE(name,N_ELEMENTS(TgStr.Name)),TgStr.Name,/FOLD_CASE) 
IF MAX(exist) NE 0 THEN BEGIN ;If name exists
;  answer=DIALOG_MESSAGE("This target name already exists. Overwrite?",/QUESTION)
;  IF answer EQ 'Yes' THEN BEGIN
    w=MAX(WHERE(exist GT 0)) ;take maximum position with the same name
    Target_Index=w 
;  ENDIF ELSE RETURN,0
ENDIF ELSE BEGIN ;If name doesn't already exist
  w=MAX(WHERE(TgStr(*).name NE '')) ;highest filled position
  IF w EQ -1 OR w GE N_ELEMENTS(TgStr)-1 THEN BEGIN
    junk=DIALOG_MESSAGE("Target Structure Full; Delete a target first.",/ERROR)
    RETURN,0
  ENDIF ELSE BEGIN
    w=w+1 ;First vacant spot in TgStr
    Target_Index=w
  ENDELSE
ENDELSE

;*************** Fill the identified target structure
TgStr0=TgStr(Target_Index) ;Define TgStr0 as the single target structure to be filled

TgStr0.Name=name
TgStr0.En=FLOAT(En)
TgStr0.Rin=FLOAT(Rin)
TgStr0.dr=FLOAT(dr)

column_names=select.column_names  & row_names=select.row_names
N_Mat=N_ELEMENTS(WHERE(column_names NE ''))-2 & N_Region=N_ELEMENTS(WHERE(row_names NE ''))-3

TgStr0.N_Mat=N_Mat & TgStr0.N_Region=N_Region

;i=Target_index
;name=TgStr(i).name ;Target name
;En=TgStr(i).En ;keV. x-ray photon energy
;Rin=TgStr(i).Rin ;inner radius of capsule
;dr=TgStr(i).dr ;step size for target profile
;d_Region=TgStr(i).d_Region;thickness of each region 
;rho_Region=TgStr(i).rho_Region;density of each region
;N_Region=TgStr(i).N_Region ;Number of regions
;N_Mat=TgStr(i).N_Mat ;Number of material elements
;mat=TgStr(i).mat
;*************
;column_names=STRCOMPRESS(['Thickness(um)','Density(g/cc)','Element'+STRING(INDGEN(N_Mat)+1)],/REMOVE_ALL)
;row_names=STRCOMPRESS(['Element','Atomic Mass','Opacity(cm2/g)','Region'+STRING(INDGEN(N_Region)+1)],/REMOVE_ALL)
;tabledata=STRARR(2+20,15)

TgStr0.mat(0:N_Mat-1).sym=tabledata(2:N_Mat-1+2,0);=STRCOMPRESS(STRING(mat(0:N_Mat-1).sym),/REMOVE_ALL)
TgStr0.mat(0:N_Mat-1).A  =tabledata(2:N_Mat-1+2,1);=STRCOMPRESS(STRING(mat(0:N_Mat-1).A),/REMOVE_ALL)
TgStr0.mat(0:N_Mat-1).mu =tabledata(2:N_Mat-1+2,2);=STRCOMPRESS(STRING(mat(0:N_Mat-1).mu),/REMOVE_ALL)
FOR i=0,N_Region-1 DO BEGIN
  TgStr0.d_Region(i)=FLOAT(tabledata(0,i+3));=STRCOMPRESS(STRING(d_Region(i)),/REMOVE_ALL)
  TgStr0.rho_Region(i)=FLOAT(tabledata(1,i+3));=STRCOMPRESS(STRING(rho_Region(i)),/REMOVE_ALL)
  TgStr0.mat(0:N_Mat-1).f(i)=tabledata(2:N_Mat-1+2,i+3);=STRCOMPRESS(STRING(mat(0:N_Mat-1).f(i)),/REMOVE_ALL)
ENDFOR

TgStr(Target_Index)=TgStr0 ;Fill the appropriate target element with this new target
select=imagetostructure(select,'TGSTR',TgStr) ;Put new structure back into memory

w=WHERE(TgStr(*).name NE '')
IF w(0) NE -1 THEN Target_Names=TgStr(w).name ELSE Target_Names='' ;list of target names
WIDGET_CONTROL,select.targetlist_id,SET_VALUE=Target_Names ; updates list of target names

select.Target_index=Target_Index
select.highlighted_target=Target_Index
WIDGET_CONTROL,select.targetlist_id,SET_LIST_SELECT=Target_Index

Last_Target=Target_Index
SAVE,TgStr,Last_Target,FILENAME=select.file

select.unsaved_data=0

skip_savetarget:
RETURN,1
END



FUNCTION Load_New_Target
COMMON tg_info, select
;This loads the new target indexed by Target_Index out of the existing structure TgStr

TgStr=select.TgStr
Target_Index=select.Target_Index

w=WHERE(TgStr(*).name NE '')
IF w(0) NE -1 THEN Target_Names=TgStr(w).name ELSE Target_Names='' ;list of target names
N_targets=N_ELEMENTS(Target_Names) ;number of targets
IF Target_Index LT 0 OR Target_Index GE N_Targets THEN Target_Index=0

;************Extract info for current target
i=Target_index
name=TgStr(i).name ;Target name
En=TgStr(i).En ;keV. x-ray photon energy
Rin=TgStr(i).Rin ;inner radius of capsule
dr=TgStr(i).dr ;step size for target profile
d_Region=TgStr(i).d_Region;thickness of each region 
rho_Region=TgStr(i).rho_Region;density of each region
N_Region=TgStr(i).N_Region ;Number of regions
N_Mat=TgStr(i).N_Mat ;Number of material elements
mat=TgStr(i).mat
;*************

column_names=STRCOMPRESS(['Thickness(um)','Density(g/cc)','Element'+STRING(INDGEN(N_Mat)+1)],/REMOVE_ALL)
row_names=STRCOMPRESS(['Element','Atomic Mass','Opacity(cm2/g)','Region'+STRING(INDGEN(N_Region)+1)],/REMOVE_ALL)
tabledata=STRARR(2+20,15)

tabledata(2:N_Mat-1+2,0)=STRCOMPRESS(STRING(mat(0:N_Mat-1).sym),/REMOVE_ALL)
tabledata(2:N_Mat-1+2,1)=STRCOMPRESS(STRING(mat(0:N_Mat-1).A),/REMOVE_ALL)
tabledata(2:N_Mat-1+2,2)=STRCOMPRESS(STRING(mat(0:N_Mat-1).mu),/REMOVE_ALL)
FOR i=0,N_Region-1 DO BEGIN
  tabledata(0,i+3)=STRCOMPRESS(STRING(d_Region(i)),/REMOVE_ALL)
  tabledata(1,i+3)=STRCOMPRESS(STRING(rho_Region(i)),/REMOVE_ALL)
  tabledata(2:N_Mat-1+2,i+3)=STRCOMPRESS(STRING(mat(0:N_Mat-1).f(i)),/REMOVE_ALL)
ENDFOR
;*******************

WIDGET_CONTROL,select.targetlist_id,SET_VALUE=Target_Names
WIDGET_CONTROL,select.name_id,SET_VALUE=name
WIDGET_CONTROL,select.energy_id,SET_VALUE=En
WIDGET_CONTROL,select.Rmin_id,SET_VALUE=Rin
WIDGET_CONTROL,select.step_id,SET_VALUE=dr
WIDGET_CONTROL,select.table_id,COLUMN_LABELS=column_names
WIDGET_CONTROL,select.table_id,ROW_LABELS=row_names
WIDGET_CONTROL,select.table_id,SET_VALUE=tabledata

d=d_region(0:N_Region-1) & rho=rho_region(0:N_Region-1)
tThick=TOTAL(d) & Rmax=Rin+tThick & trhoR=TOTAL(1.E-4*d*rho*1.E3)
RR=Rin+[0.,TOTAL(d,/CUMULATIVE)] & mass=0.
FOR i=0,N_Region-1 DO BEGIN
;  print,(4./3.)*!pi*rho(i)*(RR(i+1)^3-RR(i)^3)*1.E-6,rho(i),RR(i+1)-RR(i) 
  mass=mass+(4./3.)*!pi*rho(i)*(RR(i+1)^3-RR(i)^3)*1.E-6 
ENDFOR

WIDGET_CONTROL,select.Rmax_id,SET_VALUE="Outer Radius = "+STRING(RMax,FORMAT='(F6.1)')+" um"
WIDGET_CONTROL,select.TotalThick_id,SET_VALUE="Total Thickness = "+STRING(tThick,FORMAT='(F6.1)')+" um"
WIDGET_CONTROL,select.Mass_id,SET_VALUE="Total Mass = "+STRING(Mass,FORMAT='(F6.1)')+" ug"
WIDGET_CONTROL,select.rhoR_id,SET_VALUE="Areal Density = "+STRING(trhoR,FORMAT='(F6.1)')+" mg/cm2" 


select=imagetostructure(select,'COLUMN_NAMES',column_names)
select=imagetostructure(select,'ROW_NAMES',row_names)

select.highlighted_target=Target_Index
WIDGET_CONTROL,select.targetlist_id,SET_LIST_SELECT=Target_Index

RETURN,1
END 

PRO target_event,event
COMMON tg_info, select

WIDGET_CONTROL, event.id, GET_UVALUE = eventval

IF N_ELEMENTS(eventval) EQ 0 THEN RETURN
CASE eventval OF


;    "TABLE0":BEGIN
;
;    END
;    "RESTORE":BEGIN
;       WIDGET_CONTROL,select.table0,SET_VALUE=select.cdata
;    END
;    "CLEAR":BEGIN
;      zero=FLTARR(N_ELEMENTS(select.cdata))
;      zero(2)=1. ;mic/pixel fixed to 1.
;      zero(4)=2048 ;last pixel for sweep calib
;      zero(5)=1 ;linear sweep fixed to 1 ps/pix
;      zero=TRANSPOSE(zero)
;      WIDGET_CONTROL,select.table0,SET_VALUE=zero
;
;    END
    "OK":BEGIN
       IF select.unsaved_data  EQ 1 THEN BEGIN
        junk=DIALOG_MESSAGE("Target data has changed. Do you wish to save?",/QUESTION)
        IF junk EQ 'Yes' THEN a=savetarget(select)
       ENDIF ELSE BEGIN
        TgStr=select.TgStr & Last_Target=select.target_index
        SAVE,TgStr,Last_Target,FILENAME=select.file
       ENDELSE
       select.loadcurrent=1 ;ensures that currently loaded target gets sent to xstreak
       WIDGET_CONTROL, event.handler, /DESTROY
    END
    "CANCEL":BEGIN
       IF select.unsaved_data  EQ 1 THEN BEGIN
        junk=DIALOG_MESSAGE("Target data has changed. Do you wish to save?",/QUESTION)
        IF junk EQ 'Yes' THEN a=savetarget(select)
       ENDIF

       select.loadcurrent=0 ;no new target gets sent to xstreak
       WIDGET_CONTROL,event.handler,/DESTROY
    END
    "NAME":BEGIN
      select.unsaved_data=1
    END
    "ENERGY":BEGIN
      select.unsaved_data=1
    END
    "RMIN":BEGIN
      select.unsaved_data=1
    END
    "STEP":BEGIN
      select.unsaved_data=1
    END
    "INSERT ROW":BEGIN
     Target_Index=select.Target_Index
     select.TgStr(Target_Index).N_Region=select.TgStr(Target_Index).N_Region+1
     a=Load_New_Target()
     select.unsaved_data=1 ;current data table is different to saved file
    END
    "INSERT COLUMN":BEGIN
     Target_Index=select.Target_Index
     select.TgStr(Target_Index).N_Mat=select.TgStr(Target_Index).N_Mat+1
     a=Load_New_Target()
     select.unsaved_data=1 ;current data table is different to saved file
    END
    "DELETE ROW":BEGIN
       WIDGET_CONTROL,select.table_id,GET_VALUE=table ;find size of table
       sz=SIZE(table) & ncols=sz(1) & nrows=sz(2)
       s_cells = WIDGET_INFO(select.table_id, /TABLE_SELECT) ;find highlighted region
       
       ;Only proceed with column delection if an entire single column is highlighted
       IF s_cells(1) EQ s_cells(3) AND s_cells(0) EQ 0 AND s_cells(2) EQ ncols-1 $
        AND s_cells(1) NE 0 AND s_cells(1) NE 1 AND s_cells(1) NE 2 THEN BEGIN
        answer=DIALOG_MESSAGE("Do you want to delete the entire row of data?",/DEFAULT_NO,/QUESTION)
        IF answer EQ 'Yes' THEN BEGIN
          num=N_ELEMENTS(select.row_names)
          rownum=INDGEN(num) & w=WHERE(rownum LT s_cells(1) OR rownum GT s_cells(3))
          new_row_names=select.row_names(w)
          select=imagetostructure(select,'ROW_NAMES',new_row_names)
          WIDGET_CONTROL,select.table_id,/DELETE_ROWS ;tabledata is automatically adjusted here.
          select.unsaved_data=1 ;current data table is different to saved file
        ENDIF
       ENDIF ELSE BEGIN
        junk=DIALOG_MESSAGE("Click on Region# row to delete that region")
       ENDELSE
    END
    "DELETE COLUMN":BEGIN
       WIDGET_CONTROL,select.table_id,GET_VALUE=table 
       sz=SIZE(table) & ncols=sz(1) & nrows=sz(2)
       s_cells = WIDGET_INFO(select.table_id, /TABLE_SELECT)
       
       ;Only proceed with column delection if an entire single column is highlighted
       IF s_cells(0) EQ s_cells(2) AND s_cells(1) EQ 0 AND s_cells(3) EQ nrows-1 $
       AND s_cells(0) NE 0 AND s_cells(0) NE 1 THEN BEGIN
        answer=DIALOG_MESSAGE("Do you want to delete the entire column of data?",/DEFAULT_NO,/QUESTION)
        IF answer EQ 'Yes' THEN BEGIN
          num=N_ELEMENTS(select.column_names)
          colnum=INDGEN(num) & w=WHERE(colnum LT s_cells(0) OR colnum GT s_cells(2))
          new_column_names=select.column_names(w)
          select=imagetostructure(select,'COLUMN_NAMES',new_column_names)
          WIDGET_CONTROL,select.table_id,/DELETE_COLUMNS ;tabledata is automatically adjusted here.
          select.unsaved_data=1 ;current data table is different to saved file
        ENDIF
       ENDIF ELSE BEGIN
        junk=DIALOG_MESSAGE("Click on top of Element# column to delete that element")
       ENDELSE
    END
    "SAVE TARGET":BEGIN
       a=savetarget(select)
       a=Load_New_Target() ;Run this to update Target List
    END
;    "LOAD TARGET":BEGIN
;     IF select.unsaved_data THEN BEGIN
;      junk=DIALOG_MESSAGE("Target data has changed. Save as new?",/QUESTION)
;      IF junk EQ 'Yes' THEN BEGIN
;        a=savetarget(select)
;        IF NOT(a) THEN GOTO,skip_load ;not successfully saved so skip load
;      ENDIF
;     ENDIF
;     select.Target_Index=select.highlighted_target
;     a=Load_New_Target()
;     skip_load:
;    END
    "MOVE TARGET UP":BEGIN
      TgStr=select.TgStr & Target_Index=select.Target_Index
      IF Target_Index EQ 0 THEN GOTO,skip_targetup
      TgStr_Temp=TgStr(Target_Index-1)
      TgStr(Target_Index-1)=TgStr(Target_Index)
      TgStr(Target_Index)=TgStr_Temp
      select.TgStr=TgStr & select.Target_Index=Target_Index-1 & select.highlighted_target=Target_Index-1
      a=Load_New_Target()
      TgStr=select.TgStr & Last_Target=select.target_index
      SAVE,TgStr,Last_Target,FILENAME=select.file
      
      skip_targetup:
    END
    "MOVE TARGET DOWN":BEGIN
      TgStr=select.TgStr & Target_Index=select.Target_Index 
      Target_Names=TgStr.name & N_targets=N_ELEMENTS(WHERE(Target_Names NE '')) ;number of targets with names
      IF Target_Index EQ N_Targets-1 THEN GOTO,skip_targetdown
      TgStr_Temp=TgStr(Target_Index+1)
      TgStr(Target_Index+1)=TgStr(Target_Index)
      TgStr(Target_Index)=TgStr_Temp
      select.TgStr=TgStr & select.Target_Index=Target_Index+1 & select.highlighted_target=Target_Index+1
      a=Load_New_Target()
      TgStr=select.TgStr & Last_Target=select.target_index
      SAVE,TgStr,Last_Target,FILENAME=select.file
      
      skip_targetdown:
    END
    "DELETE TARGET":BEGIN
      junk=DIALOG_MESSAGE("Permanently delete currently selected target?",/QUESTION)
      IF junk EQ 'Yes' THEN BEGIN
        TgStr=select.TgStr
        Target_Index=select.Target_Index
        Target_Names=TgStr.name
        N_targets=N_ELEMENTS(Target_Names) ;number of target positions
        NamesNum=INDGEN(N_Targets)
        w=WHERE(NamesNum LT Target_Index OR NamesNum GT Target_Index);find all other targets
        IF w(0) NE -1 THEN New_TgStr=TgStr(w) ELSE New_TgStr=TgStr
        WIDGET_CONTROL,select.targetlist_id,SET_VALUE=New_TgStr.name
        select=imagetostructure(select,'TGSTR',New_TgStr)
;        w=WHERE(TgStr(*).name NE '')
;        IF w(0) NE -1 THEN Target_Names=TgStr(w).name ELSE Target_Names='' ;list of target names
        ;IF Target_Index LT 0 OR Target_Index GE N_Targets THEN Target_Index=0        

        
        ;Reset index to zero
        select.target_index=0
        select.highlighted_target=0
        
        ;WIDGET_CONTROL,select.targetlist_id,SET_LIST_SELECT=0
        
        a=Load_New_Target()
        TgStr=select.TgStr & Last_Target=select.target_index
        SAVE,TgStr,Last_Target,FILENAME=select.file

;        IF target_index GT N_ELEMENTS(N_targets)-2 THEN $
;          select.target_index=N_ELEMENTS(N_targets)-2
;        IF target_index LT 0 THEN select.target_index=0
        
        ;a=Load_New_Target()
;        select.unsaved_data=0

      ENDIF
    END
    "CALC OPACITY":BEGIN
     
     IF select.unsaved_data EQ 1 THEN BEGIN
      junk=DIALOG_MESSAGE("Must save target changes before opacities can be calculated.")
      GOTO,skip_opacity
     ENDIF
     
     xrayfile=select.xrayfile
     IF NOT(FILE_TEST(xrayfile)) THEN BEGIN
      junk=DIALOG_MESSAGE("xKappa.sav not found! Cannot calculate opacities.")
      GOTO,skip_opacity
     ENDIF 
     
     RESTORE,xrayfile ;structure XRAY now contains Henke data
     TgStr0=select.tgstr(select.Target_Index)

     FOR i=0,TgStr0.N_Mat-1 DO BEGIN
      sym=STRCOMPRESS(TgStr0.mat(i).sym,/REMOVE_ALL)
      En=TgStr0.En
      exist=STRCMP(REPLICATE(sym,N_ELEMENTS(xray.symb)),xray.symb,/FOLD_CASE) ;finds position of that
      w=MIN(WHERE(exist GT 0))                                                ;element in XRAY
      IF w(0) NE -1 THEN TgStr0.mat(i).mu=GetMu(XRAY,sym,En) ELSE $
        junk=DIALOG_MESSAGE('Warning! Henke data only available for '+STRJOIN(xray.symb,' '))
     ENDFOR
     select.tgstr(select.Target_Index)=TgStr0
     a=Load_New_Target()
     TgStr=select.TgStr & Last_Target=select.target_index
     SAVE,TgStr,Last_Target,FILENAME=select.file
     
     skip_opacity:
    END
;    "CREATE NEW":BEGIN
;       WIDGET_CONTROL,select.table0,GET_VALUE=current_data
;       current_data=FLOAT(current_data)
;       load_new_name,event ;user inputs assigned name for current data
;       IF select.newdataname EQ '' THEN BEGIN
;         junk=DIALOG_MESSAGE("No target name specified; cannot create a new data array.")
;         GOTO,stopcreate
;       ENDIF
;
;       WIDGET_CONTROL,select.table,GET_VALUE=tabledata
;       ;Insert a new column
;       WIDGET_CONTROL,select.table,INSERT_COLUMNS=1
;       ;Append the new target name to the structure and to the table
;       new_column_names=[select.newdataname,select.column_names]
;       select=imagetostructure(select,'column_names',new_column_names)
;       WIDGET_CONTROL,select.table,COLUMN_LABELS=select.column_names
;       ;Append the new target data to 'tabledata'
;       tabledata=[current_data,tabledata]
;       WIDGET_CONTROL,select.table,SET_VALUE=tabledata
;       select.newdata=1 ;Data array is different to saved data file
;
;       stopcreate:
;    END
    "TARGETLIST":BEGIN
    ;IF event.clicks EQ 1 THEN select.highlighted_target=event.index
    
    IF select.unsaved_data THEN BEGIN
      junk=DIALOG_MESSAGE("Target data has changed. Save as new?",/QUESTION)
      IF junk EQ 'Yes' THEN BEGIN
       a=savetarget(select) 
       IF NOT(a) THEN BEGIN
        WIDGET_CONTROL,select.targetlist_id,SET_LIST_SELECT=select.target_index
        GOTO,skip_targetlist ;not successfully saved so skip load
       ENDIF
      ENDIF ELSE select.unsaved_data=0
    ENDIF
    
    select.highlighted_target=event.index
    select.Target_Index=event.index
    a=Load_New_Target()
        
    skip_targetlist:
    END
    "TABLE":BEGIN
    IF event.type EQ 0 THEN BEGIN ;Change the state to "unsaved" if a single character is typed in
      select.unsaved_data=1
;      WIDGET_CONTROL,select.table_id,GET_VALUE=tabledata 
    ENDIF
    END

ELSE:junk=DIALOG_MESSAGE("No match to user value")
ENDCASE

RETURN
END


FUNCTION target,imstr,event
COMMON tg_info, select

fonttype='';'HELVETICA*BOLD*14'

tg_file=FILE_SEARCH(imstr.targetfile)
IF tg_file(0) EQ '' THEN BEGIN
  junk=DIALOG_MESSAGE("File "+imstr.targetfile+" cannot be found!")
  tg_file(0)=DIALOG_PICKFILE(PATH=imstr.pathopen,GET_PATH=path,$
            GROUP=event.handler,TITLE="Please find xstreak_target.sav")
  IF tg_file(0) EQ '' THEN BEGIN
    junk=DIALOG_MESSAGE("No target file was selected.")
    GOTO,skip_target
  ENDIF
  imstr.targetfile=tg_file(0)
ENDIF ;Make sure to get target name correctly

;get_target_data,tg_file,$
;          row_names,column_names,tabledata

sc = WIDGET_BASE(GROUP_LEADER=event.handler,/MODAL,TAB_MODE=1,UVALUE=0,$
       TITLE="Target Settings",EVENT_PRO='target_event')

tsc0 = WIDGET_BASE(sc,/ROW)
tsc1 = WIDGET_BASE(tsc0,/COLUMN)
lsc = WIDGET_BASE(tsc0,/COLUMN)
tsc2 = WIDGET_BASE(tsc0,/COLUMN)
tsc2a0= WIDGET_BASE(tsc2,/COLUMN)
tsc2a = WIDGET_BASE(tsc2,/ROW)
tsc2a1= WIDGET_BASE(tsc2,/ROW)
tsc2b = WIDGET_BASE(tsc2,/ROW)
tsc2c = WIDGET_BASE(tsc2,/ROW)
tsc2c1 =WIDGET_BASE(tsc2c,/COLUMN)
tsc2c2 =WIDGET_BASE(tsc2c,/COLUMN)
tsc2c3 =WIDGET_BASE(tsc2c,/COLUMN)
tsc2c4 =WIDGET_BASE(tsc2c,/COLUMN)
tsc2c5 =WIDGET_BASE(tsc2c,/COLUMN)

;****OLD STUFF MARKED FOR DELETION
;Compile values for "Current" Table
;WIDGET_CONTROL,imstr.rot_id,GET_VALUE=rads
;WIDGET_CONTROL,imstr.shr_id,GET_VALUE=shearr
;WIDGET_CONTROL,imstr.mag_id,GET_VALUE=mag
;MinPix=imstr.MinPix
;MaxPix=imstr.MaxPix
;cam_tpoly=imstr.tpoly(1:N_ELEMENTS(imstr.tpoly)-1)

;current_data=TRANSPOSE([rads,shearr,mag,MinPix,MaxPix,cam_tpoly])
;table0_id = WIDGET_TABLE(tsc1,COLUMN_LABELS=['Available Targets'],$ ;
;         /NO_ROW_HEADERS,ALL_EVENTS=0,EDITABLE=1,Y_SCROLL_SIZE=10,$
;         VALUE=STRARR(1,100),$
;         UVALUE='TABLE0')

tglabel_id = WIDGET_LABEL(tsc1,VALUE="Available targets",FONT=fonttype)
targetlist_id=WIDGET_LIST(tsc1,VALUE='',FONT=fonttype,XSIZE=30,YSIZE=10,$
                          UNITS=1,SCR_YSIZE=2.8,SCR_XSIZE=1.5,$
                          UVALUE='TARGETLIST')
;loadtarget_id = WIDGET_BUTTON(tsc1,VALUE="LOAD TARGET",UVALUE="LOAD TARGET")
MoveUptarget_id   = WIDGET_BUTTON(tsc1,VALUE="MOVE TARGET UP",UVALUE="MOVE TARGET UP")
MoveDowntarget_id = WIDGET_BUTTON(tsc1,VALUE="MOVE TARGET DOWN",UVALUE="MOVE TARGET DOWN")
deletetarget_id   = WIDGET_BUTTON(tsc1,VALUE="DELETE TARGET",UVALUE="DELETE TARGET")
;**************************

;targetlist_id=WIDGET_DROPLIST(tsc2a0,VALUE='',FONT=fonttype,UVALUE='Targetlist',UNITS=1,XSIZE=2.0)

name_id  =CW_FIELD(tsc2a1, TITLE='Target Name',XSIZE=30,FONT=fonttype,/COLUMN,/ALL_EVENTS,UVALUE='NAME',VALUE='')
energy_id=CW_FIELD(tsc2a1, TITLE='Energy (keV)',XSIZE=11,FONT=fonttype,/COLUMN,/ALL_EVENTS,UVALUE='ENERGY',VALUE='')
Rmin_id  =CW_FIELD(tsc2a1, TITLE='Inner Radius (um)',XSIZE=11,FONT=fonttype,/COLUMN,/ALL_EVENTS,UVALUE='RMIN',VALUE='')
step_id  =CW_FIELD(tsc2a1, TITLE='Step Size (um)',XSIZE=11,FONT=fonttype,/COLUMN,/ALL_EVENTS,UVALUE='STEP',VALUE='')
WIDGET_CONTROL,step_id,SENSITIVE=0

table_id = WIDGET_TABLE(tsc2b,COLUMN_LABELS='',$
         ROW_LABELS='',ALL_EVENTS=0,EDITABLE=1,X_SCROLL_SIZE=7,Y_SCROLL_SIZE=10,$
         VALUE=FLTARR(22,15),UVALUE='TABLE')

;***********************
;rsc = WIDGET_BASE(tsc1,/ROW)
;restoreID = WIDGET_BUTTON(rsc,VALUE="RESTORE",UVALUE="RESTORE")
;clearID = WIDGET_BUTTON(rsc,VALUE="CLEAR",UVALUE="CLEAR")
;
;
;;loadc_label = WIDGET_LABEL(lsc, VALUE="Loads highlighted cells to current")
;loadc_id = WIDGET_BUTTON(lsc,VALUE="<-- LOAD TO CURRENT",UVALUE="LOAD CURRENT")
;createnew_id = WIDGET_BUTTON(lsc,VALUE="CREATE AS NEW -->",UVALUE="CREATE NEW")

;bsc = WIDGET_BASE(lsc,/COLUMN)
;add_ID = WIDGET_BUTTON(bsc,VALUE="ADD COLUMN",UVALUE="ADD COLUMN")

okID = WIDGET_BUTTON(tsc2c1, VALUE="OK",UVALUE="OK")
cancelID = WIDGET_BUTTON(tsc2c1,VALUE="CANCEL",UVALUE="CANCEL")

insR_ID = WIDGET_BUTTON(tsc2c2,VALUE="INSERT REGION",UVALUE="INSERT ROW")
insC_ID = WIDGET_BUTTON(tsc2c2,VALUE="INSERT ELEMENT",UVALUE="INSERT COLUMN")

delR_ID = WIDGET_BUTTON(tsc2c3,VALUE="DELETE REGION",UVALUE="DELETE ROW")
delC_ID = WIDGET_BUTTON(tsc2c3,VALUE="DELETE ELEMENT",UVALUE="DELETE COLUMN")

save_ID = WIDGET_BUTTON(tsc2c4,VALUE="SAVE TARGET",UVALUE="SAVE TARGET")
CalcOpacity_ID = WIDGET_BUTTON(tsc2c4,FONT=fonttype,VALUE="Calculate Opacities",UVALUE="CALC OPACITY")

Rmax_id = WIDGET_LABEL(tsc2c5,/ALIGN_LEFT,/DYNAMIC_RESIZE,VALUE="Outer Radius (um)=")
TotalThick_id = WIDGET_LABEL(tsc2c5,/ALIGN_LEFT,/DYNAMIC_RESIZE,VALUE="Total Thickness (um)=")
Mass_id = WIDGET_LABEL(tsc2c5,/ALIGN_LEFT,/DYNAMIC_RESIZE,VALUE="Total Mass (um)=")
rhoR_id = WIDGET_LABEL(tsc2c5,/ALIGN_LEFT,/DYNAMIC_RESIZE,VALUE="Areal Density (mg/cm2)=")


;*************Get full Target Structure array
RESTORE,tg_file(0) ;Restores TgStr and Last_Target

select={$;left:-1,top:-1,right:-1,bottom:-1,$
       loadcurrent:0,file:tg_file,xrayfile:imstr.xrayfile,$
       highlighted_target:Last_Target,$
;       table0:table0_id,$
;       cdata:current_data,$
       unsaved_data:0,$ 
       column_names:'',row_names:'',$
;       newdataname:'',newdata:0,$
       targetlist_id:targetlist_id,name_id:name_id,$
       energy_id:energy_id,Rmin_id:Rmin_id,step_id:step_id,$
       Rmax_id:Rmax_id,TotalThick_id:TotalThick_id,Mass_id:Mass_id,rhoR_id:rhoR_id,$
       table_id:table_id,TgStr:TgStr,Target_Index:Last_Target        }

a=Load_New_Target() ;This loads the last target that was used
;***************

WIDGET_CONTROL,sc,/REALIZE

XMANAGER, 'target', sc, event='target_event', $
       GROUP_LEADER = event.handler

;help,/struct,select.TgStr(select.Target_Index)
;print,select.Target_Index


IF select.loadcurrent EQ 0 THEN RETURN,-1 ELSE RETURN,select.TgStr(select.Target_Index)

skip_target:

RETURN,-1
END

;PRO get_target_data,file,row_names,column_names,tabledata
;
;line='' & ncams=0 & ndats=0
;OPENR,unit,file,/GET_LUN
;READF,unit,line
;READF,unit,line
;
;;READF,unit,ncams & READF,unit,ndats ;ncams=# of targets (settings); ndats=# of entries per target
;
;;column_names=STRARR(ncams)
;;row_names=STRARR(ndats)
;;tabledata=STRARR(ncams,ndats)
;
;;aline=STRARR(ncams+1)
;
;;READF,unit,line
;;READS,line,aline,FORMAT='('+STRTRIM(STRING(ncams+1))+'A18)'
;;column_names=aline(1:ncams)
;
;;FOR i=0,ndats-1 DO BEGIN
;;   READF,unit,line
;;   READS,line,aline,FORMAT='('+STRTRIM(STRING(ncams+1))+'A18)'
;;   row_names(i)=aline(0)
;;   tabledata(*,i)=aline(1:ncams)
;;ENDFOR
;
;READF,unit,line
;aline=STRSPLIT(line,COUNT=N,/EXTRACT) & ncams=N-1 ;counts the number of target names
;column_names=aline(1:ncams)
;
;WHILE NOT(EOF(unit)) DO BEGIN
;    ON_IOERROR,skip
;    READF,unit,line & READS,line,aline,FORMAT='('+STRTRIM(STRING(ncams+1))+'A18)'
;    IF SIZE(all_data,/TYPE) NE 0 THEN BEGIN
;        all_data=[[all_data],[aline(1:ncams)]]
;        row_names=[row_names,aline(0)]
;    ENDIF ELSE BEGIN
;       all_data=aline(1:ncams)
;       row_names=aline(0)
;    ENDELSE
;    skip:
;ENDWHILE
;
;CLOSE,unit
;FREE_LUN,unit
;
;row_names=row_names & column_names=column_names & tabledata=all_data
;RETURN
;END

;PRO get_target_data_old,file,row_names,column_names,tabledata
;
;line='' & ncols=0 & nrows=0
;
;OPENR,unit,file,/GET_LUN
;READF,unit,line
;READF,unit,ncams & READF,unit,ndats ;ncams=# of targets (settings); ndats=# of entries per target
;
;column_names=STRARR(ncams)
;row_names=STRARR(ndats)
;tabledata=STRARR(ncams,ndats)
;
;READF,unit,line
;s=STR_SEP(line,'\')
;
;FOR i=0,ndats-1 DO row_names(i)=s(i+1)
;
;FOR j=0,ncams-1 DO BEGIN
;    READF,unit,line
;    s=STR_SEP(line,'\')
;    column_names(j)=s(0)
;    FOR i=0,ndats-1 DO tabledata(j,i)=s(i+1)
;ENDFOR
;
;tabledata=FLOAT(tabledata)
;
;CLOSE,unit
;FREE_LUN,unit
;
;RETURN
;END

;PRO load_new_name_event,event
;COMMON tg_info, select
;
;WIDGET_CONTROL, event.id, GET_UVALUE = eventval
;
;IF N_ELEMENTS(eventval) EQ 0 THEN RETURN
;CASE eventval OF
;    "OK":BEGIN
;       WIDGET_CONTROL,event.handler,GET_UVALUE=bname_text
;       WIDGET_CONTROL,bname_text,GET_VALUE=newname
;       select.newdataname=STRCOMPRESS(newname,/REMOVE_ALL)
;       WIDGET_CONTROL,event.handler,/DESTROY
;    END
;    "CANCEL":BEGIN
;       select.newdataname=''
;       WIDGET_CONTROL,event.handler,/DESTROY
;     END
;ELSE:junk=DIALOG_MESSAGE("No match to user value")
;ENDCASE
;
;RETURN
;END
;
;PRO load_new_name,event
;;Creates a modal text widget to accept a new name for the current target settings
;
;bname = WIDGET_BASE(GROUP_LEADER=event.handler,/MODAL,TAB_MODE=1,UVALUE=0,$
;       TITLE="Enter New target Name",EVENT_PRO='new_name_event')
;
;bname1 = WIDGET_BASE(bname,/COLUMN)
;bname2 = WIDGET_BASE(bname1,/ROW)
;bname3 = WIDGET_BASE(bname1,/ROW)
;
;bname_text = WIDGET_TEXT(bname2,/EDITABLE)
;bname_ok = WIDGET_BUTTON(bname3,VALUE="OK",UVALUE="OK")
;bname_cancel = WIDGET_BUTTON(bname3,VALUE="CANCEL",UVALUE="CANCEL")
;
;WIDGET_CONTROL,bname,/REALIZE
;WIDGET_CONTROL,bname,SET_UVALUE=bname_text
;
;XMANAGER, 'load_new_name', bname, event='load_new_name_event', $
;       GROUP_LEADER = event.handler
;
;END