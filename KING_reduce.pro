@headfits.pro
@fxposit.pro
@mrd_hread.pro
@sxpar.pro
@gettok.pro
@mrdfits.pro
@fxpar.pro
@valid_num.pro
@mrd_skip.pro
@writefits.pro
@check_fits.pro
@fxaddpar.pro
@sxdelpar.pro
@uniq.pro
@resistant_mean.pro
@mean.pro
@array_indices.pro
@mkhdr.pro
@sxaddpar.pro
@get_date.pro
@daycnv.pro
@cgimage.pro
@image_dimensions.pro
@cgdefcharsize.pro
@setdefaultvalue.pro
@cgdefaultcolor.pro
@cggetcolorstate.pro
@cgerase.pro
@cgsetcolorstate.pro
@cgcolor.pro
@strsplit.pro
@cgsnapshot.pro
@cgcolor24.pro
@clipscl.pro
@cgresizeimage.pro
@congrid.pro
@julday.pro
@sixlin.pro
@proceeding_text.pro
@fixpix_mod.pro
@dist_circle.pro
@cgscalevector.pro
@fpufix.pro
@cghistoplot.pro
@convert_to_type.pro
@cgcheckforsymbols.pro
@correl_optimize.pro
@correl_images.pro
@corrmat_analyze.pro
@mpfit2dpeak.pro
@mpfit.pro
@mpfit2dfun.pro
@shift_sub.pro
@scale_image_am.pro
@cgplot.pro
@cgbitget.pro
@colorsareidentical.pro
@la_cosmic_MOD.pro
@sky.pro
@mmm.pro
@reverse.pro
@asinh.pro
@poly.pro
@fftshift.pro

@KING_Bias.pro
@KING_BPmap_SkyFlat.pro
@KING_Dark.pro
; @KING_BPmap.pro
; @KING_SkyFlat.pro
@KING_reduce_Science.pro
@KING_cosmic_ray.pro
@KING_stack_and_add.pro



pro KING_reduce, display=display

path = '/home/amueller/work/KING/20180222/AFGL5173_180.0_Ha/Raw/'
;path = '/home/amueller/work/KING/20180213/Gaia_180.0_R/Raw/'
;path = '/home/amueller/work/KING/20171101/M57_180.0_R/Raw/'
;path = '/home/amueller/work/KING/20171019/M57_300.0_V/Raw/'
;path = '/home/amueller/work/KING/20171122/Gaia_180.0_R/Raw/'
;path = '/home/amueller/work/KING/20171122/Gaia_200.0_R/Raw/'

nodark = 1	;1: don't use darks
stackmethod = 'gauss'
;stackmethod = 'correl'

;=========================================================================

pathr = path+'../Reduced/'
spawn, 'mkdir -p '+path+'../Reduced/'

;===============================================================================================

;Bias
KING_Bias, path

;===============================================================================================

;DARK
if (nodark ne 1) then KING_Dark, path

;===============================================================================================

;BPmap and SkyFlat
KING_BPmap_SkyFlat, path, display, nodark

;===============================================================================================

;reduce Science
;nodark = 0
KING_reduce_Science, path, nodark, sout

;===============================================================================================
;sout = 'Gaia_180.000s_R_reduced.fits'
;stack and add
KING_stack_and_add, pathr, sout, stackmethod

;===============================================================================================

;sout = 'Gaia_180.000s_R_reduced_shift.fits'
;cosmic ray correction
KING_cosmic_ray, pathr, sout

;===============================================================================================

print, ''
print, 'Reduction finished!'
print, ''

end
