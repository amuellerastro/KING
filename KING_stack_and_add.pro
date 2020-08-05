pro KING_stack_and_add, path, sout, method

print, ''
print, 'IMAGE REGISTRATION'
print, ''

im = mrdfits(path+sout,0,/silent)
fnout = strmid(sout, 0, strpos(sout, '.', /reverse_search))

nframes = n_elements(im[0,0,*])
dim = n_elements(im[*,0,0])

;=========================================================================

radius = 7.

window, 0, xs=1000, ys=1000
device, cursor_standard=2
guess_off = dblarr(2,nframes)
xc = dblarr(nframes)
yc = dblarr(nframes)
fitpos = dblarr(2,nframes)

print, 'Click on same star!'

for i=0,nframes-1 do begin
  
  print, strcompress(i+1,/rem)+' / '+strcompress(nframes)
  
  cgimage, scale_image_am(im[*,*,i]), /axis
  cursor, tx, ty, 3, /data
  x = tx
  y = ty

  x = ceil(x) & y = ceil(y)

  tmp_cutim = im[x-radius:x+radius, y-radius:y+radius, i]

  mx = max(tmp_cutim, location)
  idxmax = array_indices(tmp_cutim, location)

  xc[i] = x+(idxmax[0]-radius)
  yc[i] = y+(idxmax[1]-radius)

  cutim = im[xc[i]-radius:xc[i]+radius,yc[i]-radius:yc[i]+radius,i]
  
  if (i eq 0) then begin
  
    guess_off[0,i] = 0.
    guess_off[1,i] = 0.
    
  endif else begin
  
    guess_off[0,i] = xc[0]-xc[i]
    guess_off[1,i] = yc[0]-yc[i]
  
  endelse
  
  ;-----------------------------------------------------------------------
  
  ;fit Gaussian
  
  ;xa = (dindgen(2*radius+1)+1.d0)#(dblarr(2.*radius+1)+1)
  ;ya = (dblarr(2.*radius+1)+1)#(dindgen(2*radius+1)+1.d0)
  xa = dindgen(n_elements(cutim[*,0]))
  ya = dindgen(n_elements(cutim[0,*]))
  
  dumerr = xa
  dumerr[*,*] = 1.

  ;		amplitude		fwhm	fwhm	xc	yc	rot	power
  estimates = [median(cutim), mx-median(cutim), 1.5, 1.5, radius, radius, 0.0]
  
  pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},7)
  pi[1].limited(0) = 1
  pi[1].limits(0) = 1.0d0
  
  yfit = mpfit2dpeak(cutim, A, xa, ya, estimates=estimates, dof=dof, chisq=chisq, perror=perror, sigma=sigma, /Gaussian, /tilt, /quiet, parinfo=pi)
  
  fitpos[0,i] = A[4]-radius+xc[i]
  fitpos[1,i] = A[5]-radius+yc[i]

  ;-----------------------------------------------------------------------  
  
endfor
;restore, 'tmp.sav', /v

;=========================================================================

if (nframes gt 1) then begin

  ;=======================================================================
  
  ;Method 1: correlate images
  
  if (method eq 'correl') then begin

    stmp = im
    for i=0,nframes-1 do begin

      stmp[*,*,i] = shift_sub(im[*,*,i], guess_off[0,i], guess_off[1,i])

    endfor

    im = stmp
  
    tmp_im = im;[400:800,400:800,*]
    master = tmp_im[*,*,0]

    off1 = dblarr(2,nframes)
    off2 = off1

    xyshift = 7	;odd number
    shift_im = im

    for i=0,nframes-1 do begin

  ;     if (i eq 0) then correl_optimize, master, reform(tmp_im[*,*,i]), xoff, yoff, Init_factor=1 $
  ;     else 
      correl_optimize, master, reform(tmp_im[*,*,i]), xoff, yoff;, xoff_init = guess_off[0,i], yoff_init = guess_off[1,i];, Init_factor=1
      off1[*,i] = [xoff,yoff]
      
      corrmatrix = correl_images(master, reform(tmp_im[*,*,i]), xoffset_b=xoff, yoffset_b=yoff, xshift=xyshift, yshift=xyshift, /monitor)

      xa = dindgen(n_elements(corrmatrix[*,0]))
      ya = dindgen(n_elements(corrmatrix[0,*]))
      yfit = mpfit2dpeak(corrmatrix, A, xa, ya, /gauss, dof=dof, chisq=chisq, perror=perror, sigma=sigma, /quiet)

      off2[*,i] = A[4:5]

      shift_im[*,*,i] = shift_sub(im[*,*,i], off1[0,i]-(xyshift-off2[0,i]), off1[1,i]-(xyshift-off2[1,i]))

    endfor

  endif
    
  ;=======================================================================
  
  ;Method 2: Fit Gaussian

  if (method eq 'gauss') then begin
  
    shift_im = im
    dim = n_elements(im[*,0,0])
    
    for i=0,nframes-1 do begin

      ;shift_im[*,*,i] = shift_sub(im[*,*,i], xc[0]-fitpos[0,i], yc[0]-fitpos[1,i])
      
      wframe = 200.	;width of frame
      tmp = dblarr(dim+2.*wframe, dim+2.*wframe)
      tmp[wframe:wframe+dim-1, wframe:wframe+dim-1] = im[*,*,i]
      stmp = fftshift(tmp, xc[0]-fitpos[0,i], yc[0]-fitpos[1,i])
      shift_im[*,*,i] = stmp[(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1,(dim+2.*wframe)/2.-dim/2:(dim+2.*wframe)/2.+dim/2.-1]

    endfor

  endif
  
  ;=======================================================================
  
  fn = path+fnout+'_shift.fits'
  snout = fn
  writefits, fn, shift_im

  fn = path+fnout+'_shift_total.fits'
  writefits, fn, total(shift_im,3)
  
  fn = path+fnout+'_shift_mean.fits'
  writefits, fn, mean(shift_im,dim=3)

  fn = path+fnout+'_shift_product.fits'
  writefits, fn, product(shift_im,3)/max(product(shift_im,3))

  fn = path+fnout+'_shift_median.fits'
  writefits, fn, median(shift_im,dim=3,/even)
  
;   resm_im = dblarr(n_elements(shift_im[*,0,0]), n_elements(shift_im[0,*,0]))
;   for xx=0,n_elements(shift_im[*,0,0])-1 do begin
;     for yy=0,n_elements(shift_im[0,*,0])-1 do begin
; 
;       resistant_mean, shift_im[xx,yy,*], 1.5, val, sigma, num_rej, goodvec=good, badvec=bad
;       resm_im[xx,yy] = val
; 
;     endfor
;   endfor
;   writefits, path+fnout+'_shift_resistant_mean.fits', resm_im
  

endif

end
