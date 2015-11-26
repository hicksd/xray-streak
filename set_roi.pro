PRO set_roi,imstr,roi,roiindex
;This takes a 2x2 roi array, loads it into the specific
;roi structure given by roiindex, and draws that roi.


		imroi=[ [roi(0,0),roi(1,0)],$
		 		[roi(0,0),roi(1,1)],$
		 		[roi(0,1),roi(1,1)],$
		 		[roi(0,1),roi(1,0)],$
		 		[roi(0,0),roi(1,0)]]/imstr.factor

	;Set ROI for plot/phase/intensity and draw appropriate boundary on image
	  CASE roiindex OF
	;	CASE imstr.roilist(roiindex) OF
	;	'Set Plot ROI':BEGIN
	  0:BEGIN
			imstr.xyroi=roi
			imstr.oimpolyxyroi -> SetProperty, DATA=imroi,$
					COLOR=[255,255,255],LINESTYLE=0
		END
;		'Set Phase ROI':BEGIN
    1:BEGIN
			imstr.phroi=roi
			imstr.oimpolyphroi -> SetProperty, DATA=imroi,$
					COLOR=[255,255,255],LINESTYLE=2
		END
;		'Set Intensity ROI':BEGIN
    2:BEGIN
			imstr.inroi=roi
			imstr.oimpolyinroi -> SetProperty, DATA=imroi, $
					COLOR=[255,255,255], LINESTYLE=3
		END
;    'Set Ghost ROI':BEGIN
    3:BEGIN
      imstr.ghroi=roi
      imstr.oimpolyghroi -> SetProperty, DATA=imroi, $
          COLOR=[255,255,255], LINESTYLE=4
    END		
		ENDCASE

		imstr.xyvals=[$	
		  [imstr.xyroi(0,0),imstr.phroi(0,0),imstr.inroi(0,0),imstr.ghroi(0,0),imstr.xyvals(4,0)],$
			[imstr.xyroi(0,1),imstr.phroi(0,1),imstr.inroi(0,1),imstr.ghroi(0,1),imstr.xyvals(4,1)],$
			[imstr.xyroi(1,0),imstr.phroi(1,0),imstr.inroi(1,0),imstr.ghroi(1,0),imstr.xyvals(4,2)],$
			[imstr.xyroi(1,1),imstr.phroi(1,1),imstr.inroi(1,1),imstr.ghroi(1,1),imstr.xyvals(4,3)] ]

		WIDGET_CONTROL,imstr.xytable,SET_VALUE=imstr.xyvals

		imstr.oimwindow -> Draw,imstr.oimview
END