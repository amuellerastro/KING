pro KING_track_SSobject

path = '/home/amueller/work/KING/20171122/Gaia_180.0_R/Reduced/'
file = path+'Gaia_180.000s_R_reduced_reduced_shift.fits'
fn_out1 = path+'manual_shift_stack_add.fits'
fn_out2 = path+'manual_shift_stack_mean.fits'

im = mrdfits(file,0,/silent)
sim = im

nframes = n_elements(im[0,0,*])


window, 0, xs=1200, ys=1200
device, cursor_standard=2

print, 'Click on same object!'

radius = 5.
for i=0,nframes-1 do begin

  print, 'Image ', i+1
  ;cgimage, scale_image_am(im[*,*,i]), /axis
  cgimage, im[*,*,i], /axis, stretch=1, minvalue=1200, maxvalue=1400
  cursor, tx, ty, 3, /data
  x = tx
  y = ty

  x = ceil(x) & y = ceil(y)

  tmp_cutim = im[x-radius:x+radius, y-radius:y+radius]

  mx = max(tmp_cutim, location)
  idxmax = array_indices(tmp_cutim, location)

  xc = x+(idxmax[0]-radius)
  yc = y+(idxmax[1]-radius)

  if (i eq 0) then begin
  
    xref = xc
    yref = yc
  
  endif
  
  sim[*,*,i] = shift(im[*,*,i], xref-xc, xref-yc)  

endfor

writefits, fn_out1, total(sim,3)
writefits, fn_out2, mean(sim,dim=3)


stop
end
