Function DISC_deWarp, Data_File, Warp_file
;; This function de-warps DISC data given the raw streak camera image and the appropriate Warp correction
;; Inputs are the full file locations of the Data file and Warp Correction file (from NSTEC)
;; This function subtracts the pre shot background from the raw image then interpolates the data from its spatial and temporal locations to the de-warped spatial and temporal positions.
;; The output of this function is the 4200 X 4200 de-warped streak image
;; This function also writes a tiff file of the de-warped image in the current directory 

;Syntax:
;Data = Disc_dewarp(Dialog_pickfile(title = 'Pick the data file', filter = '*.h5'), Dialog_pickfile(title = 'Pick the Warp file', filter = '*.h5'))

;I downloaded the warp correction file and broke it apart in IDL.  I saved it as an idl save file.
;If you use IDL, you can restore this file, read in the raw image, and apply the warp correction directly with the following IDL commands (paste into the command line)
;
;Read in the raw image into IDL
;                X=h5_browser(pickfile()
;If you load the image with the default name, it is called data._data
;
;Save the file attached to this email to a location on your computer
;
;Enter the IDL commands:
;restore,filename=pickfile() 
;this is where you select the image you just saved to disk (DISC_3_warp_correction_variables.dat)
;
;triangulate,lin_xpos,lin_ypos,tr,b
;warp_correct=interpolate(data._data,trigrid(lin_xpos,lin_ypos,nonlin_xpos,tr,[1,1],[0,0,4199,4199]),trigrid(lin_xpos,lin_ypos,nonlin_ypos,tr,[1,1],[0,0,4199,4199]))
;
;the warp-corrected image is then called warp_correct

sig='DATA/SHOT/DATA'
bg='DATA/PRESHOT_01/DATA'
file_id=h5f_open(warp_file)
dataset_id1=h5d_open(file_id,'DATA/LINEAR_XPOS/DATA')
linxpos=h5d_read(dataset_id1)
h5d_close,dataset_id1

dataset_id2=h5d_open(file_id,'DATA/LINEAR_YPOS/DATA')
linypos=h5d_read(dataset_id2)
h5d_close,dataset_id2

dataset_id3=h5d_open(file_id,'DATA/NONLINEAR_XPOS/DATA')
nonlinxpos=h5d_read(dataset_id3)
h5d_close,dataset_id3

dataset_id4=h5d_open(file_id,'DATA/NONLINEAR_YPOS/DATA')
nonlinypos=h5d_read(dataset_id4)
h5d_close,dataset_id4
 
file_id=h5f_open(data_file)
dataset_id5=h5d_open(file_id,sig)
DISC_shot_raw=float(h5d_read(dataset_id5))
h5d_close,dataset_id5

file_id=h5f_open(data_file)
dataset_id6=h5d_open(file_id,bg)
DISC_shot_bg=float(h5d_read(dataset_id6))
h5d_close,dataset_id6
DISC_shot=DISC_shot_raw - DISC_shot_bg
triangulate,linxpos,linypos,tr,b
DISC_shot_warp_correct=interpolate(DISC_shot,trigrid(linxpos,linypos,nonlinxpos,tr,[1,1],[0,0,4199,4199]),trigrid(linxpos,linypos,nonlinypos,tr,[1,1],[0,0,4199,4199]))

WRITE_TIFF,'N120227_dewarp.tif',DISC_shot_warp_correct,/float
;WRITE_TIFF,'N120227BgSub.tif',DISC_shot,/float

Return, DISC_shot_warp_correct
END