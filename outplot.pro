PRO OUTPLOT,xxx,yyy,imstr,xtitl=xtitl,ytitl=ytitl,titl=titl,xstyl=xstyl,$
         error=error,hist=hist,psy=psy,newplot=newplot,overplot=overplot,color=color,$
         ystyl=ystyl,xrang=xrang,yrang=yrang,charsiz=charsiz

    szx=SIZE(xxx) & szy=SIZE(yyy)
    IF szx(0) EQ 1 THEN xx=xxx ELSE xx=xxx(*,0)
    IF szy(0) EQ 1 THEN yy=yyy ELSE yy=yyy(*,0)
    
    fin=FINITE(yy) & w=WHERE(fin EQ 0) & IF w[0] NE -1 THEN yy[w]=0.0
     
    IF imstr.limlist(imstr.limindex) EQ 'Auto Limits' AND $
       imstr.overlist(imstr.overindex) EQ 'New Plot' THEN BEGIN

         xr = [MIN(xx),MAX(xx)]
         IF xr(0) EQ xr(1) THEN BEGIN
         	xr(1)=xr(0)+1
		    xr(0)=xr(0)-1
		 ENDIF
         IF NOT(KEYWORD_SET(error)) THEN BEGIN
         	yr = [MIN(yy)*0.95,MAX(yy)*1.05]
			IF yr(0) EQ yr(1) THEN BEGIN
				yr(1)=yr(0)+1
				yr(0)=yr(0)-1
			ENDIF
         ENDIF ELSE yr = [MIN(yy-error),MAX(yy+error)]

         imstr.xyvals(4,0)=xr(0) & imstr.xyvals(4,1)=xr(1)
         imstr.xyvals(4,2)=yr(0) & imstr.xyvals(4,3)=yr(1)
         WIDGET_CONTROL,imstr.xytable,SET_VALUE=imstr.xyvals

    ENDIF

    IF imstr.limlist(imstr.limindex) EQ 'Manual Limits' OR $
       imstr.overlist(imstr.overindex) EQ 'Overplot' THEN BEGIN

         xr = [FLOAT(imstr.xyvals(4,0)),FLOAT(imstr.xyvals(4,1))]
         yr = [FLOAT(imstr.xyvals(4,2)),FLOAT(imstr.xyvals(4,3))]

         xwithin=WHERE(xx GE xr(0) AND xx LE xr(1))

;         IF xwithin(0) EQ -1 THEN BEGIN
;          junk=DIALOG_MESSAGE('Plot ROI lies outside plot limits')
;          GOTO,break_outplot
;         ENDIF ELSE BEGIN
;          xx=xx(xwithin) & yy=yy(xwithin)
;         ENDELSE
    ENDIF

    xrange=xr(1)-xr(0)
    yrange=yr(1)-yr(0)

    newornot=imstr.overlist(imstr.overindex)          ;This can be used to plot multiple plots, 
    IF KEYWORD_SET(overplot) THEN newornot='Overplot' ;even if "New Plot" is set
    IF KEYWORD_SET(newplot) THEN newornot='New Plot' ;Sets a new plot, regardless 
    IF KEYWORD_SET(color) THEN clr=color ELSE clr=0   ;fixes color of plot
    
    CASE newornot OF
       'New Plot':BEGIN
         FOR pp=1,imstr.nplots-1 DO BEGIN
          imstr.oplmodel -> Remove,imstr.oplplot(pp)
          OBJ_DESTROY,imstr.oplplot(pp)
         ENDFOR

         IF KEYWORD_SET(psy) THEN BEGIN
          oplsymbol=OBJ_NEW('IDLgrSymbol')
          oplsymbol -> SetProperty,DATA=psy,$
            SIZE=[FLOAT(xrange)/150.,FLOAT(yrange)/150.]
          oplot0=imstr.oplplot(0) & oplot0->SetProperty, SYMBOL=oplsymbol,LINESTYLE=6 & imstr.oplplot(0)=oplot0
         ENDIF ELSE BEGIN
          oplot0=imstr.oplplot(0) & oplot0->SetProperty, SYMBOL='',LINESTYLE=0 & imstr.oplplot(0)=oplot0
         ENDELSE
         
         oplot0=imstr.oplplot(0) & oplot0-> SetProperty, DATAX=xx, DATAY=yy,$
          MIN_VALUE=yr(0),MAX_VALUE=yr(1),HISTOGRA=hist,COLOR=clr & imstr.oplplot(0)=oplot0

         imstr.nplots = 1

         IF KEYWORD_SET(error) THEN BEGIN
          imstr.oplplot(1) = OBJ_NEW('IDLgrPlot')
          imstr.oplmodel -> Add, imstr.oplplot(1)
          imstr.oplplot(2) = OBJ_NEW('IDLgrPlot')
          imstr.oplmodel -> Add, imstr.oplplot(2)

          oplot1=imstr.oplplot(1) & oplot1-> SetProperty, DATAX=xx, DATAY=yy+error,$
              MIN_VALUE=yr(0),MAX_VALUE=yr(1),LINESTYLE=1,COLOR=clr & imstr.oplplot(1)=oplot1
          oplot2=imstr.oplplot(2) & oplot2->SetProperty, DATAX=xx, DATAY=yy-error,$
              MIN_VALUE=yr(0),MAX_VALUE=yr(1),LINESTYLE=1,COLOR=clr & imstr.oplplot(2)=oplot2

          imstr.nplots=3
         ENDIF

         xtl = 0.02 * (yrange)
         ytl = 0.02 * (xrange)

         imstr.oplxaxis0 -> GetProperty, TITLE=oxtitle
         imstr.oplyaxis0 -> GetProperty, TITLE=oytitle
         imstr.oplxaxis1 -> GetProperty, TITLE=otitle

         oxtitle -> SetProperty, STRINGS=xtitl
         oytitle -> SetProperty, STRINGS=ytitl
         otitle -> SetProperty, STRINGS=titl

         imstr.oplxaxis0 -> SetProperty, RANGE=xr, TICKLEN=xtl, LOCATION=[xr(0),yr(0),0.]
         imstr.oplxaxis1 -> SetProperty, RANGE=xr, TICKLEN=xtl, LOCATION=[xr(0),yr(1),0.]
         imstr.oplyaxis0 -> SetProperty, RANGE=yr, TICKLEN=ytl, LOCATION=[xr(0),yr(0),0.]
         imstr.oplyaxis1 -> SetProperty, RANGE=yr, TICKLEN=ytl, LOCATION=[xr(1),yr(0),0.]

         imstr.oplmodel -> Reset
         imstr.oplmodel -> Scale, 1.0/xrange,1.0/yrange,1.0
         imstr.oplmodel -> Translate, -xr(0)/xrange, -yr(0)/yrange, 0.

         imstr.oplxaxis0 -> GetProperty, TICKTEXT=oXTickText
         oxtitle -> SetProperty, RECOMPUTE_DIMENSIONS=2
         oXTickText -> SetProperty, RECOMPUTE_DIMENSIONS=2

         imstr.oplYaxis0 -> GetProperty, TICKTEXT=oYTickText
         oytitle -> SetProperty, RECOMPUTE_DIMENSIONS=2
         oYTickText -> SetProperty, RECOMPUTE_DIMENSIONS=2

         imstr.oplxaxis1 -> GetProperty, TICKTEXT=ox1TickText
         otitle -> SetProperty, RECOMPUTE_DIMENSIONS=2
         ox1TickText -> SetProperty, RECOMPUTE_DIMENSIONS=2,STRINGS=''

         imstr.oplwindow -> Draw, imstr.oplview
       END
       'Overplot':BEGIN
         IF imstr.nplots GE N_ELEMENTS(imstr.oplplot) THEN BEGIN
          junk=DIALOG_MESSAGE('Cannot have more than '$
           +STRCOMPRESS(STRING(imstr.nplots),/REMOVE_ALL)+' plots')
          GOTO, break_overplot
         END
         IF imstr.nplots EQ 0 THEN BEGIN
          junk=DIALOG_MESSAGE('Create a new plot first before overplotting')
          WIDGET_CONTROL,imstr.overtype,SET_DROPLIST_SELECT=0
          imstr.overindex=0
          GOTO, break_overplot
         END

         imstr.oplplot(imstr.nplots) = OBJ_NEW('IDLgrPlot')
         imstr.oplmodel -> Add, imstr.oplplot(imstr.nplots)

;         IF KEYWORD_SET(psy) THEN BEGIN
;          oplsymbol=OBJ_NEW('IDLgrSymbol',psy)
;          imstr.oplplot(imstr.nplots) -> SetProperty, SYMBOL=oplsymbol, LINESTYLE=6
;         ENDIF
         IF KEYWORD_SET(psy) THEN BEGIN
          oplsymbol=OBJ_NEW('IDLgrSymbol')
          imstr.oplplot[0] -> GetProperty,XRANGE=xrang,YRANGE=yrang
          xrange=xrang[1]-xrang[0] & yrange=yrang[1]-yrang[0]
          oplsymbol -> SetProperty,DATA=psy,SIZE=[FLOAT(xrange)/100.,FLOAT(yrange)/100.],FILL_COLOR=clr,/FILLED
          oplotn=imstr.oplplot(imstr.nplots) & oplotn->SetProperty, SYMBOL=oplsymbol,LINESTYLE=6
          imstr.oplplot(imstr.nplots)=oplotn
         ENDIF ELSE BEGIN
          oplotn=imstr.oplplot(imstr.nplots) & oplotn->SetProperty, SYMBOL='',LINESTYLE=0
          imstr.oplplot(imstr.nplots)=oplotn
         ENDELSE
         
         oplotn=imstr.oplplot(imstr.nplots) & oplotn-> SetProperty, DATAX=xx, DATAY=yy,$
          MIN_VALUE=yr(0),MAX_VALUE=yr(1),HISTOGRA=hist,COLOR=clr & imstr.oplplot(imstr.nplots)=oplotn
           
         ;imstr.oplplot(0) -> GetProperty, XRANGE=xr0,YRANGE=yr0
         ;imstr.oplplot(imstr.nplots) -> SetProperty, DATAX=xx,DATAY=yy, $
         ;     MIN_VALUE=yr(0),MAX_VALUE=yr(1),HISTOGRA=hist, LINESTYLE=0,COLOR=clr

         imstr.nplots=imstr.nplots+1

         IF KEYWORD_SET(error) THEN BEGIN
          imstr.oplplot(imstr.nplots) = OBJ_NEW('IDLgrPlot')
          imstr.oplmodel -> Add, imstr.oplplot(imstr.nplots)
          oplotn=imstr.oplplot(imstr.nplots) & oplotn-> SetProperty, DATAX=xx, DATAY=yy+error,$
              MIN_VALUE=yr(0),MAX_VALUE=yr(1),LINESTYLE=1,COLOR=clr & imstr.oplplot(imstr.nplots)=oplotn
          imstr.nplots=imstr.nplots+1

          imstr.oplplot(imstr.nplots) = OBJ_NEW('IDLgrPlot')
          imstr.oplmodel -> Add, imstr.oplplot(imstr.nplots)
          oplotn=imstr.oplplot(imstr.nplots) & oplotn-> SetProperty, DATAX=xx, DATAY=yy-error,$
              MIN_VALUE=yr(0),MAX_VALUE=yr(1),LINESTYLE=1,COLOR=clr & imstr.oplplot(imstr.nplots)=oplotn
          imstr.nplots=imstr.nplots+1
         ENDIF



         imstr.oplwindow -> Draw, imstr.oplview

         break_overplot:
       END
    ENDCASE

    ;PLOT,xx,yy(*,0),XTITLE=xtitl,YTITLE=ytitl,TITLE=titl,XSTYLE=xstyl,$
    ;  YSTYLE=ystyl,CHARSIZE=charsiz,PSYM=psy,XRANGE=xrang,YRANGE=YRANG
    ; OPLOT,xx,yy(*,0),PSYM=psy




break_outplot:
RETURN
END

