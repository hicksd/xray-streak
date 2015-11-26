FUNCTION RotateFlip,imstr,action


CASE action OF
    'Right 90':image=TRANSPOSE(REVERSE(imstr.image))
    'Left 90':image=REVERSE(TRANSPOSE(imstr.image))
    'Flip Horizontal':image=REVERSE(imstr.image)
    'Flip Vertical':image=TRANSPOSE(REVERSE(TRANSPOSE(imstr.image)))
END

newimstr=imagetostructure(imstr,'IMAGE',image)

image_size=SIZE(image)
newimstr.xsize=image_size(1)/newimstr.factor
newimstr.ysize=image_size(2)/newimstr.factor

WIDGET_CONTROL,newimstr.vimageID,DRAW_XSIZE=image_size(1)/newimstr.factor
WIDGET_CONTROL,newimstr.vimageID,DRAW_YSIZE=image_size(2)/newimstr.factor

wshow_images,newimstr,'IMAGE'
newimstr.whichimage='IMAGE'

newimstr.lastrot=0 & newimstr.lastshear=0


RETURN,newimstr
END