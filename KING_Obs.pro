;this script has to be inside the observers directory, e.g. '/home/obs70/data/amueller/'
;offsets: kappa70 -> kappa_gui
;focus V: +62
;focus R: +47
;for focus choose ~6mag star, go 10 steps within a range of 50, start at +35, 1sec
;object moves 0.15px/min on the detector -> max. DIT 6min


@caldat.pro
@mrdfits.pro
@fxposit.pro
@mrd_hread.pro
@fxpar.pro
@gettok.pro
@valid_num.pro
@mrd_skip.pro
@sxaddpar.pro
@writefits.pro
@check_fits.pro
@sxdelpar.pro
@sxpar.pro
@headfits.pro
@julday.pro
@sky.pro
@mmm.pro
@asinh.pro
@cgimage.pro
@image_dimensions.pro
@cgdefcharsize.pro
@setdefaultvalue.pro
@cgdefaultcolor.pro
@cgerase.pro
@cgcolor.pro
@cgsnapshot.pro
@cgcolor24.pro
@cgresizeimage.pro
@cgplot.pro
@cgcheckforsymbols.pro
@colorsareidentical.pro
@array_indices.pro
@mpfit2dfun.pro
@mpfit.pro
@poly_fit.pro
@mpfit2dpeak.pro
; @cggetcolorstate.pro
; @cgsetcolorstate.pro
@strsplit.pro
@congrid.pro
; @cgbitget.pro
@convert_to_type.pro
; @linspace.pro
@poly.pro
@mean.pro
@mkhdr.pro
@get_date.pro
@daycnv.pro
@fxaddpar.pro
@fxparpos.pro


function get_mid_jd, hdr

  exptime = double(sxpar(hdr, 'EXPTIME'))
  dateobs = sxpar(hdr, 'DATE-OBS')
  t1 = double(strmid(dateobs, 0, 4))
  t2 = double(strmid(dateobs, 5, 2))
  t3 = double(strmid(dateobs, 8, 2))
  t4 = double(strmid(dateobs, 11, 2))
  t5 = double(strmid(dateobs, 14, 2))
  t6 = double(strmid(dateobs, 17, 2))
  mid_jd = julday(t2,t3,t1,t4,t5,t6)+(exptime/2.)/86400.d0

  return, mid_jd

end

function linspace, start, stop, num, step=step, bincentres=bincentres
  compile_opt idl2

  if n_elements(num) eq 0 then $
     if keyword_set(bincentres) then $
        num = 50 $
     else $
        num = 51

  if (size(start, /type) eq 5) or (size(stop, /type) eq 5) then begin
      integers = dindgen(num)   ; Double precision</span>

      num = double(num)
  endif else begin
      integers = findgen(num)   ; Floating point </span>
      num = float(num)
  endelse

  if keyword_set(bincentres) then begin
      step = (stop - start) / num
      start = start + step / 2
      stop = stop - step / 2
  endif else $
     step = (stop - start) / (num - 1.)

  linspace_array = start + integers * step

  return, linspace_array

end

pro hak, NoPrompt=noprompt

  while Get_KBRD(0) do junk=1
  if Keyword_Set(noprompt) then begin
    junk = Get_KBRD(1)
  endif else begin
    Print, 'Hit any key to continue...'
    junk = Get_KBRD(1)
    if StrUpCase(junk) eq 'Q' then return
    print, 'Continuing...
  endelse

  while Get_KBRD(0) do junk=1

end

function scale_image_am, im

  w = where(finite(im) and (im ne 0.), goodcount)
  sky, im[w], skymode, skysig, /silent
  mn = skymode - (2.*skysig) > min(im)
  mx = asinh((max(im)-mn)/skysig)
  scaled_im = bytscl(asinh((im-mn)/skysig),min=0,max=mx,top=246)+8

  return, scaled_im

end

pro time_name, name=name, stamp=stamp, night=night, fileflag=fileflag

  if (fileflag eq 0) then caldat, double(systime(/julian, /utc))-0.5, mo, dy, yr
  if (fileflag eq 1) then caldat, double(systime(/julian, /utc)), mo, dy, yr
  yr = strcompress(yr,/rem)
  if (mo lt 10.) then mo = '0'+strcompress(mo,/rem) else mo = strcompress(mo,/rem)
  if (dy lt 10.) then dy = '0'+strcompress(dy,/rem) else dy = strcompress(dy,/rem)

  caldat, double(systime(/julian, /utc)), t1, t2, t3, hr, mn, sec
  if (hr lt 10.) then hr = '0'+strcompress(hr,/rem) else hr = strcompress(hr,/rem)
  if (mn lt 10.) then mn = '0'+strcompress(mn,/rem) else mn = strcompress(mn,/rem)
  sec = round(sec)
  if (sec lt 10.) then sec = '0'+strcompress(sec,/rem) else sec = strcompress(sec,/rem)

  name = 'KING.'+yr+'-'+mo+'-'+dy+'T'+hr+'_'+mn+'_'+sec+'.fits'
  stamp = yr+'-'+mo+'-'+dy+'T'+hr+'_'+mn+'_'+sec
  night = yr+mo+dy

end

pro write_night_log, night, path, fn, catg, object, selfilt

  dum = file_search(path+night+'.log', count=n)
  if (n ne 1) then begin
  
    openw, lun, path+night+'.log', width=2000, /get_lun, /append
    
      printf, lun, '                         File        CATG         OBJECT            MID_JD     DIT FILTER   T_CCD'
      printf, lun, '-------------------------------------------------------------------------------------------------'
    
    close, lun
    free_lun, lun

  endif

  hdr = headfits(path+fn, exten=0, /silent)
  mid_jd = get_mid_jd(hdr)
  exptime = sxpar(hdr, 'EXPTIME')
  ccdtemp = sxpar(hdr, 'CCDTEMP')

  openw, lun, path+night+'.log', width=2000, /get_lun, /append
  
    printf, lun, fn, catg, object, mid_jd, exptime, selfilt, ccdtemp, format='(a29, a12, a15, f18.8, f8.1, a7, f8.2)'
  
  close, lun
  free_lun, lun

end

pro show_image, path, file, object, flag, binfact, quest

  im = mrdfits(path+file,0,hdr,/fscale,/silent)

  if (flag eq 0) then begin
  
    window, 0, xs=1000, ys=1000, title=file+ '/ '+object
    cgimage, scale_image_am(im)

    maxv = sxpar(hdr ,'DATAMAXI')
    minv = sxpar(hdr ,'DATAMINI')
    meanv = sxpar(hdr ,'DATAMEAN')
    print, ''
    print, 'min value: '+strcompress(round(minv),/rem)+' / max value: '+strcompress(round(maxv),/rem)+' / mean value: '+strcompress(round(meanv),/rem)

    if (quest eq '3') then begin
    
      window, 2, xs=500, ys=500
      cghistoplot, im, xr=[0,65540], xtitle='Counts', ytitle='Number of pixels'
      oplot, 60.d3*[1,1], !y.crange, linestyle=2
    
    endif
    
  endif
  
  if (flag eq 1) then begin
  
    window, 0, xs=1000, ys=1000, title=file+ '/ '+object
    cgimage, scale_image_am(im), /axis

    read, 'Enter DEC [deg]: ', dec
    
    print, 'Click on target you want to center!'
    cursor, x, y, /data
    
    npix = 2048./binfact
    xc = 2048./binfact/2. & yc = xc
    
    res = binfact*15./5600.*206.265d0	;arcsec/px, 15um pixelsize, f=5.6m
    ;RA offset
    dx = xc-x
    dxsec = (dx*res/15.)*cos(dec*!dtor)
    if (dxsec lt 0.) then da_min = ceil(dxsec/60.) else da_min = uint(dxsec/60.)
    da_sec = round(dxsec mod 60.)
    
    ;DEC offset
    dy = y-yc
    dysec = dy*res
    if (dysec lt 0.) then dd_min = ceil(dysec/60.) else dd_min = uint(dysec/60.)
    dd_sec = round(dysec mod 60.)
    
    print, ''
    print, 'Offset in R.A.: '+strcompress(da_min,/rem)+' min '+strcompress(da_sec,/rem)+' sec'
    print, 'Offset in DEC : '+strcompress(dd_min,/rem)+' min '+strcompress(dd_sec,/rem)+' sec'

  endif
  
end

function moffat_slope, x, y, p

  widx  = abs(p[1]) > 1e-20 & widy  = abs(p[2]) > 1e-20 
  xp    = x-p[3]            & yp    = y-p[4]
  theta = p[5]
  c  = cos(theta) & s  = sin(theta)
  u = ( (xp * (c/widx) - yp * (s/widx))^2 + $
                (xp * (s/widy) + yp * (c/widy))^2 )

  moffat = p[0] / (u + 1)^p[6]

  slope = p[7]*x+p[8]*y

  fit = moffat+slope

  return, fit

end


pro KING_Obs

;lat: 8deg 43.42min
;long: +49deg 23.72min
;alt: 615m

;=========================================================================

;define main directory for the night

;base = '/home/obs70/data/amueller/'
base = '/home/amueller/Downloads/amueller/'
spawn, 'mkdir -p '+base


time_name, night=night
path = base+night+'/'
spawn, 'mkdir -p '+path

;TEST TEST TEST
; path = '/home/amueller/Downloads/test/amueller/20171019/'

; cd, path

prime = '"'
observer = 'Andre Mueller'
binfact = '2'

filter = ['U', 'B', 'V', 'R', 'I', 'Ha', 'SIIc', 'SII', 'OIII', 'OII']
print, ''
for i=0,n_elements(filter)-1 do print, strcompress(i+1,/rem)+'  '+filter[i]
read, 'Select filter: ', fquest
selfilt = filter[fquest-1]

;=========================================================================

;selection menu

menujump:
quest = ''
print, ''
print, '0: Reset electronics'
print, '1: Bias'
print, '2: Dark'
print, '3: Flat Field'
print, '4: Acquisition'
print, '5: Focus'
print, '6: Science'
print, 't: Test exposure'
print, 'c: Cal Checker'
print, '9: Exit'
print, 'B: Set Binning (default: 2)'
print, 'f: Reanalyze Focus'
print, 'w: Change Filter'
read, 'What to do?: ', quest

;=========================================================================

;Electronic reset

if (quest eq '0') then begin

  spawn, 'ccdread -R -I -C 1'

  goto, menujump

endif

;=========================================================================

;Bias

if (quest eq '1') then begin

  catg = 'CALIB'
  object = 'BIAS'
  selfiltB = 'none'

  read, 'Number of Bias frames: ', nquest

  fn_all = strarr(nquest)
  if (binfact eq 1) then im_all = dblarr(2100,2048,nquest)
  if (binfact eq 2) then im_all = dblarr(1050,1024,nquest)

  for i=0,nquest-1 do begin

    print, ''
    print, '******************************************************'
    print, 'Exposure '+strcompress(i+1,/rem)+' / '+strcompress(round(nquest),/rem)
  
    time_name, name=fn, fileflag=1
    fn_all[i] = fn
    spawn, 'ccdread -C 0 -b '+binfact+' -g 5 -T -f '+fn+' -E 0 -o '+prime+object+prime+' -u '+prime+observer+prime+' -O -C 1'

    im = mrdfits(fn,0,hdr,/fscale,/silent)
    im_all[*,*,i] = im
    sxaddpar, hdr, 'CATG', catg
    sxaddpar, hdr, 'FILTER', selfiltB
    mid_jd = get_mid_jd(hdr)
    sxaddpar, hdr, 'MID_JD', mid_jd, format='(f16.8)'
    writefits, fn, im, hdr
    spawn, 'mv '+fn+' '+path

    write_night_log, night, path, fn, catg, object, selfiltB
    
    show_image, path, fn, object, 0, binfact, quest

  endfor

  ;create temporarly master bias
  if (nquest gt 1) then bias = median(im_all, dim=3, /even) else bias = im_all
  writefits, path+'tmpBias.fits', bias
  
  goto, menujump

endif

;=========================================================================

;Dark

if (quest eq '2') then begin

  catg = 'CALIB'
  object = 'DARK'
  selfiltD = 'none'

  read, 'Number of Dark frames: ', nquest
  expquest = ''
  read, 'Exposure time [sec]: ', expquest

  for i=0,nquest-1 do begin

    print, ''
    print, '******************************************************'
    print, 'Exposure '+strcompress(i+1,/rem)+' / '+strcompress(round(nquest),/rem)
  
    time_name, name=fn, fileflag=1
    spawn, 'ccdread -C 0 -b '+binfact+' -g 5 -T -f '+fn+' -E '+expquest+' -o '+prime+object+prime+' -u '+prime+observer+prime+' -O -C 1'

    im = mrdfits(fn,0,hdr,/fscale,/silent)
    sxaddpar, hdr, 'CATG', catg
    sxaddpar, hdr, 'FILTER', selfiltD
    mid_jd = get_mid_jd(hdr)
    sxaddpar, hdr, 'MID_JD', mid_jd, format='(f16.8)'
    writefits, fn, im, hdr
    spawn, 'mv '+fn+' '+path
    
    write_night_log, night, path, fn, catg, object, selfiltD
    
    show_image, path, fn, object, 0, binfact, quest

  endfor

  goto, menujump

endif

;=========================================================================

;Flat Field

if (quest eq '3') then begin

  print, ''
  print, 'Current filter: ', selfilt
  print, 'Start telescope guiding!'
  print, 'Set Focus, press enter when done!'
  hak
  
  catg = 'CALIB'
  object = 'SKYFLAT'

  read, 'Number of flat field frames: ', nquest
  expquest = ''
  read, 'Exposure time [sec]: ', expquest

  for i=0,nquest-1 do begin

    print, ''
    print, '******************************************************'
    print, 'Exposure '+strcompress(i+1,/rem)+' / '+strcompress(round(nquest),/rem)

    movequest = ''
    if (i gt 0 and i lt nquest) then begin

      print, ''
      read, 'Move telescope in RA/DEC (f=finished): ', movequest
      if (movequest ne 'f') then stop
    
    endif
  
    time_name, name=fn, fileflag=1
    spawn, 'ccdread -C 0 -b '+binfact+' -g 5 -T -f '+fn+' -e '+expquest+' -o '+prime+object+prime+' -u '+prime+observer+prime+' -O -C 1'

    im = mrdfits(fn,0,hdr,/fscale,/silent)
    sxaddpar, hdr, 'CATG', catg
    sxaddpar, hdr, 'FILTER', selfilt
    mid_jd = get_mid_jd(hdr)
    sxaddpar, hdr, 'MID_JD', mid_jd, format='(f16.8)'
    writefits, fn, im, hdr
    spawn, 'mv '+fn+' '+path
    
    write_night_log, night, path, fn, catg, object, selfilt
    
    show_image, path, fn, object, 0, binfact, quest
    
  endfor

  goto, menujump

endif

;=========================================================================

;Acquisition

if (quest eq '4') then begin

  print, ''
  print, 'Current filter: ', selfilt
  
  catg = 'ACQUISITION'
  object = catg

  ;read, 'Number of flat field frames: ', nquest
  nquest = 1
  expquest = ''
  read, 'Exposure time [sec]: ', expquest

  for i=0,nquest-1 do begin
  
    print, ''
    print, '******************************************************'
    print, 'Exposure '+strcompress(i+1,/rem)+' / '+strcompress(round(nquest),/rem)

    time_name, name=fn, fileflag=1
    spawn, 'ccdread -C 0 -b '+binfact+' -g 5 -T -f '+fn+' -e '+expquest+' -o '+prime+object+prime+' -u '+prime+observer+prime+' -O -C 1'

    im = mrdfits(fn,0,hdr,/fscale,/silent)
    sxaddpar, hdr, 'CATG', catg
    sxaddpar, hdr, 'FILTER', selfilt
    mid_jd = get_mid_jd(hdr)
    sxaddpar, hdr, 'MID_JD', mid_jd, format='(f16.8)'
    writefits, fn, im, hdr
    spawn, 'mv '+fn+' '+path
    
    write_night_log, night, path, fn, catg, object, selfilt
    
    show_image, path, fn, object, 1, binfact, quest

  endfor

  goto, menujump

endif

;=========================================================================

;Focus sequence

if (quest eq '5') then begin

  print, ''
  print, 'Current filter: ', selfilt

  CATG = 'CALIB'
  object = 'FOCUS'

  expquest = ''
  read, 'Exposure time [sec] (1.5s for 6mag): ', expquest
  read, 'Focus steps (10, start at +30): ', nsteps

  fpos = dblarr(nsteps)	;focus position as displayed on terminal
  fnall = strarr(nsteps)

  for i=0,nsteps-1 do begin
  
    print, ''
    print, '******************************************************'
    print, 'Exposure '+strcompress(i+1,/rem)+' / '+strcompress(round(nsteps),/rem)

    read, 'New focus value: ', fval
    fpos[i] = fval
    time_name, name=fn, fileflag=1
    fnall[i] = path+fn
    spawn, 'ccdread -C 0 -b '+binfact+' -g 5 -T -f '+fn+' -e '+expquest+' -o '+prime+object+prime+' -u '+prime+observer+prime+' -O -C 1'

    im = mrdfits(fn,0,hdr,/fscale,/silent)
    sxaddpar, hdr, 'CATG', catg
    sxaddpar, hdr, 'FILTER', selfilt
    mid_jd = get_mid_jd(hdr)
    sxaddpar, hdr, 'MID_JD', mid_jd, format='(f16.8)'
    writefits, fn, im, hdr
    spawn, 'mv '+fn+' '+path
    
    write_night_log, night, path, fn, catg, object, selfilt
    
    show_image, path, fn, object, 0, binfact, quest

  endfor

  time_name, stamp=tstamp, fileflag=1
  fnout = path+'FocusSequence_'+selfilt+'_'+tstamp+'.sav'
  save, fnall, fpos, nsteps, filename=fnout

  ;-----------------------------------------------------------------------

  ;Analysis
  
  bias = mrdfits(path+'tmpBias.fits',0,/silent)
  
  restore, fnout, /v
  if (binfact eq 1) then im = dblarr(2100,2048,nsteps)
  if (binfact eq 2) then im = dblarr(1050,1024,nsteps)
  for i=0,nsteps-1 do im[*,*,i] = mrdfits(fnall[i],0,hdr,/fscale,/silent)-bias
  naxis1 = sxpar(hdr ,'NAXIS1')
  naxis2 = sxpar(hdr ,'NAXIS2')

  radius = 25
  cutim = dblarr(2.*radius+1, 2.*radius+1, nsteps)

  medim = median(im,dim=3,/even)
  window, 1, xs=naxis1, ys=naxis2
  cgimage, scale_image_am(medim), /axis

  print, ''
  print, 'Select Star'
  ;manually select star for focus evaluation
  cursor, x, y, 3, /data
  x = ceil(x) & y = ceil(y)

  tmp_cutim = medim[x-radius:x+radius, y-radius:y+radius]
  mx = max(tmp_cutim, location)
  idxmax = array_indices(tmp_cutim, location)
  xc = x+(idxmax[0]-radius)
  yc = y+(idxmax[1]-radius)
  
  spos = dblarr(2,nsteps)	;fiitted stellar position
  samp = dblarr(nsteps)	;fitted amplitude
  sampe = dblarr(nsteps)
  sfwhm = dblarr(2,nsteps)	;fitted FWHM
  sfwhme = dblarr(2,nsteps)
  
  for i=0,nsteps-1 do begin

    ;final cutted region
    cutim[*,*,i] = im[xc-radius:xc+radius, yc-radius:yc+radius,i]

    ;fit stellar image
;     xa = (dindgen(2.*radius+1)+1.d0)#(dblarr(2.*radius+1.)+1.)
;     ya = (dblarr(2.*radius+1)+1.)#(dindgen(2.*radius+1.)+1.d0)
;     ;		amplitude		fwhm	fwhm	xc	yc	rot	power	slope x,y
;     estimates = [mx-median(cutim[*,*,i]), 10., 10., radius, radius, 0.1, 1., 1., 1.]
;     A = mpfit2dfun('moffat_slope', xa, ya, reform(cutim[*,*,i]), dumerr, estimates, dof=dof, chisq=chisq, perror=perror, yfit=yfit, /quiet);, weights=weights, parinfo=pi
;     spos[*,i] = A[3:4]
;     samp[i] = A[0]
;     sfwhm[*,i] = A[1:2]

    yfit = mpfit2dpeak(reform(cutim[*,*,i]), A, /GAUSSIAN, perror=perror, /tilt, /quiet);/circular, /quiet
    
    spos[*,i] = A[4:5]
    samp[i] = A[1]
    sampe[i] = perror[1]
    sfwhm[*,i] = A[2:3]
    sfwhme[*,i] = perror[2:3]
  endfor

  ra = poly_fit(fpos, samp, 2, yfit=yfit, measure_errors=sampe)
  rx = poly_fit(fpos, reform(sfwhm[0,*]), 2, yfit=yfit, measure_errors=reform(sfwhme[0,*]))  
  ry = poly_fit(fpos, reform(sfwhm[1,*]), 2, yfit=yfit, measure_errors=reform(sfwhme[1,*]))

  ra = reform(ra) & rx = reform(rx) & ry = reform(ry)
  min0 = -1.*ra[1]/(2*ra[2])
  min1 = -1.*rx[1]/(2*rx[2])
  min2 = -1.*ry[1]/(2*ry[2])
  
  focus = strcompress(round(mean([min1,min2])),/rem)
  
  xa = linspace(fpos[0], fpos[nsteps-1],100)
  
  
  window, 1, xs=700, ys=700, title='Focus'
  !p.multi=[0,1,2]
  plot, fpos, samp, psym=2, /yn, charsize=2.5, ytitle='Amplitude'
    oplot, xa, poly(xa, ra)
  plot, fpos, sfwhm[0,*], /nodata, yr=[0.5,3.5], yst=1, charsize=2.5, ytitle='FWHM', xtitle='Focus position', title='Focus = '+strcompress(focus,/rem)
    oplot, fpos, sfwhm[0,*], psym=2, color=cgcolor('red')
    oplot, xa, poly(xa, rx), color=cgcolor('red')
    oplot, fpos, sfwhm[1,*], psym=2, color=cgcolor('green')
    oplot, xa, poly(xa, ry), color=cgcolor('green')
;     legend, ['Focus: '+strcompress(focus,/rem)], /top, /left, margon=0, box=0, charsize=2
  !p.multi=[0,1,0]

  print, 'Min. values: ', min0, min1, min2
  ;print, 'Average: ', median([min0,min1,min2],/even)
  print, '*******************'
  print, 'Focus position: ', focus
  print, '*******************'
  
  goto, menujump

endif

;=========================================================================

;Science

if (quest eq '6') then begin

  print, ''
  print, 'Current filter: ', selfilt

  catg = 'SCIENCE'

  id = ''
  read, 'Object ID: ', id
  object = id
  read, 'Number of frames: ', nquest
  expquest = ''
  read, 'Exposure time [sec]: ', expquest

  for i=0,nquest-1 do begin
  
    print, ''
    print, '******************************************************'
    print, 'Exposure '+strcompress(i+1,/rem)+' / '+strcompress(round(nquest),/rem)

    time_name, name=fn, fileflag=1
    spawn, 'ccdread -C 0 -b '+binfact+' -g 5 -T -f '+fn+' -e '+expquest+' -o '+prime+object+prime+' -u '+prime+observer+prime+' -O -C 1'

    im = mrdfits(fn,0,hdr,/fscale,/silent)
    sxaddpar, hdr, 'CATG', catg
    sxaddpar, hdr, 'FILTER', selfilt
    mid_jd = get_mid_jd(hdr)
    sxaddpar, hdr, 'MID_JD', mid_jd, format='(f16.8)'
    writefits, fn, im, hdr
    spawn, 'mv '+fn+' '+path
    
    write_night_log, night, path, fn, catg, object, selfilt
    
    show_image, path, fn, object, 0, binfact, quest


  endfor

  goto, menujump

endif

;=========================================================================

;Test exposure

if (quest eq 't') then begin

  print, ''
  print, 'Current filter: ', selfilt

  catg = 'TEST'
  object = 'TEST'

  expquest = ''
  read, 'Exposure time [sec]: ', expquest
  
  time_name, name=fn, fileflag=1
  spawn, 'ccdread -C 0 -b '+binfact+' -g 5 -T -f '+fn+' -e '+expquest+' -o '+prime+object+prime+' -u '+prime+observer+prime+' -O -C 1'

  im = mrdfits(fn,0,hdr,/fscale,/silent)
  sxaddpar, hdr, 'CATG', catg
  sxaddpar, hdr, 'FILTER', selfilt
  mid_jd = get_mid_jd(hdr)
  sxaddpar, hdr, 'MID_JD', mid_jd, format='(f16.8)'
  writefits, fn, im, hdr
  spawn, 'mv '+fn+' '+path
    
  write_night_log, night, path, fn, catg, object, selfilt
    
  show_image, path, fn, object, 0, binfact, quest

  goto, menujump

endif

;=========================================================================

;Check for calibration files

if (quest eq 'c') then begin

  print, ''

  readcol, path+night+'.log', file, catg, object, jd, dit, filter, format='a,a,a,d,a,a', /silent

  st = dit+'_'+filter

  ;check BIAS
  idx = where(catg eq 'BIAS')
  if (idx[0] ne -1) then print, '*** BIAS frams are missing ***'


  ;check science
  idx = where(catg eq 'SCIENCE')
  if (idx[0] ne -1) then begin

    sdit = dit[idx]
    sfilter = filter[idx]

    usdit = sdit[uniq(sdit, sort(sdit))]
    usfilter = sfilter[uniq(sfilter, sort(sfilter))]

    ;darks for science
    for i=0,n_elements(usdit)-1 do begin
      idx = where(dit eq usdit[i] and catg eq 'CALIB' and object eq 'DARK')
      if (idx[0] eq -1) then print, '*** DARK frames with DIT '+usdit[i]+' sec are missing ***'
    endfor

    ;flats for science
    for i=0,n_elements(usfilter)-1 do begin
      idx = where(filter eq usfilter[i] and catg eq 'CALIB' and object eq 'SKYFLAT')
      if (idx[0] eq -1) then print, '*** Flat field frames with filter '+usfilter[i]+' are missing ***'
    endfor

  endif

  ;check darks for flat field
  idx = where(catg eq 'CALIB' and object eq 'SKYFLAT')
  if (idx[0] ne -1) then begin

    fdit = dit[idx]
    ffilter = filter[idx]

    ufdit = fdit[uniq(fdit, sort(fdit))]
    uffilter = ffilter[uniq(ffilter, sort(ffilter))]

    ;darks for flat fields
    for i=0,n_elements(ufdit)-1 do begin
      idx = where(dit eq ufdit[i] and catg eq 'CALIB' and object eq 'DARK')
      if (idx[0] eq -1) then print, '*** DARK frames with DIT '+ufdit[i]+' sec are missing ***'
    endfor

  endif

  goto, menujump
   
endif

;=========================================================================

;binning factor

if (quest eq 'B') then begin

  binfact = ''
  read, 'Binning factor (default: 2): ', binfact

  goto, menujump
  
endif

;=========================================================================

;Reanalysis of focus sequence

if (quest eq 'f') then begin

  bias = mrdfits(path+'tmpBias.fits',0,/silent)

  file = file_search(path+'FocusSequence_*.sav')
  for i=0,n_elements(file)-1 do print, strcompress(i+1, /rem)+' '+file[i]
  read, 'Select focus sequence: ', self
  file = file[self-1]

  restore, file, /v

  if (binfact eq 1) then im = dblarr(2100,2048,nsteps)
  if (binfact eq 2) then im = dblarr(1050,1024,nsteps)
  for i=0,nsteps-1 do im[*,*,i] = mrdfits(fnall[i],0,hdr,/fscale,/silent)-bias
  naxis1 = sxpar(hdr ,'NAXIS1')
  naxis2 = sxpar(hdr ,'NAXIS2')

  radius = 25
  cutim = dblarr(2.*radius+1, 2.*radius+1, nsteps)

  medim = median(im,dim=3,/even)
  window, 1, xs=naxis1, ys=naxis2
  cgimage, scale_image_am(medim), /axis

  print, ''
  print, 'Select Star'
  ;manually select star for focus evaluation
  cursor, x, y, 3, /data
  x = ceil(x) & y = ceil(y)

  tmp_cutim = medim[x-radius:x+radius, y-radius:y+radius]
  mx = max(tmp_cutim, location)
  idxmax = array_indices(tmp_cutim, location)
  xc = x+(idxmax[0]-radius)
  yc = y+(idxmax[1]-radius)
  
  spos = dblarr(2,nsteps)	;fiitted stellar position
  samp = dblarr(nsteps)	;fitted amplitude
  sampe = dblarr(nsteps)
  sfwhm = dblarr(2,nsteps)	;fitted FWHM
  sfwhme = dblarr(2,nsteps)
  
  for i=0,nsteps-1 do begin

    ;final cutted region
    cutim[*,*,i] = im[xc-radius:xc+radius, yc-radius:yc+radius,i]

    yfit = mpfit2dpeak(reform(cutim[*,*,i]), A, /GAUSSIAN, perror=perror, /tilt, /quiet);/circular, /quiet

    spos[*,i] = A[4:5]
    samp[i] = A[1]
    sampe[i] = perror[1]
    sfwhm[*,i] = A[2:3]
    sfwhme[*,i] = perror[2:3]

  endfor

  ra = poly_fit(fpos, samp, 2, yfit=yfit, measure_errors=sampe)
  rx = poly_fit(fpos, reform(sfwhm[0,*]), 2, yfit=yfit, measure_errors=reform(sfwhme[0,*]))
  ry = poly_fit(fpos, reform(sfwhm[1,*]), 2, yfit=yfit, measure_errors=reform(sfwhme[1,*]))

  ra = reform(ra) & rx = reform(rx) & ry = reform(ry)
  min0 = -1.*ra[1]/(2*ra[2])
  min1 = -1.*rx[1]/(2*rx[2])
  min2 = -1.*ry[1]/(2*ry[2])
  
  focus = strcompress(round(mean([min1,min2])),/rem)
  
  xa = linspace(fpos[0], fpos[nsteps-1],100)
  window, 1, xs=700, ys=700, title='Focus'
  !p.multi=[0,1,2]
  plot, fpos, samp, psym=2, /yn, charsize=2.5, ytitle='Amplitude'
    oplot, xa, poly(xa, ra)
  plot, fpos, sfwhm[0,*], /nodata, yr=[0.5,3.5], yst=1, charsize=2.5, ytitle='FWHM', xtitle='Focus position', title='Focus = '+strcompress(focus,/rem)
    oplot, fpos, sfwhm[0,*], psym=2, color=cgcolor('red')
    oplot, xa, poly(xa, rx), color=cgcolor('red')
    oplot, fpos, sfwhm[1,*], psym=2, color=cgcolor('green')
    oplot, xa, poly(xa, ry), color=cgcolor('green')
  !p.multi=[0,1,0]

  print, 'Min. values: ', min0, min1, min2
  ;print, 'Average: ', median([min0,min1,min2],/even)
  print, '*******************'
  print, 'Focus position: ', focus
  print, '*******************'
  
  goto, menujump

endif

;=========================================================================

;change filter

if (quest eq 'w') then begin

  print, ''
  for i=0,n_elements(filter)-1 do print, strcompress(i+1,/rem)+'  '+filter[i]
  read, 'Select filter: ', fquest
  selfilt = filter[fquest-1]

  goto, menujump

endif

;=========================================================================

if (quest eq '9') then return


stop
end
