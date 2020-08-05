pro KING_cosmic_ray, path, sout

print, ''
print, 'COSMIC RAY DETECTION'
print, ''

im = mrdfits(path+sout,0,/silent)
fnout = strmid(sout, 0, strpos(sout, '.', /reverse_search))

nframes = n_elements(im[0,0,*])

im_corr = im
for i=0,nframes-1 do begin

  print, 'Frame '+strcompress(i+1,/rem)+' / '+strcompress(nframes,/rem)
  sky, im[200:600,200:600,i], skymode, skysig, /silent
  tmp = im[*,*,i]-skymode
  la_cosmic_MOD, tmp, outim, '.', skyval=skymode, gain=5., readn=2., sigclip=2.

  im_corr[*,*,i] = outim
  
endfor


writefits, path+fnout+'_cosmic.fits', im_corr

if (nframes gt 1) then begin

  fn = path+fnout+'_cosmic_total.fits'
  writefits, fn, total(im_corr,3)
    
  fn = path+fnout+'_cosmic_mean.fits'
  writefits, fn, mean(im_corr,dim=3)

  fn = path+fnout+'_cosmic_product.fits'
  writefits, fn, product(im_corr,3)/max(product(im_corr,3))

  fn = path+fnout+'_cosmic_median.fits'
  writefits, fn, median(im_corr,dim=3,/even)
    
  resm_im = dblarr(n_elements(im_corr[*,0,0]), n_elements(im_corr[0,*,0]))
  for xx=0,n_elements(im_corr[*,0,0])-1 do begin
  
    for yy=0,n_elements(im_corr[0,*,0])-1 do begin

      resistant_mean, im_corr[xx,yy,*], 1.5, val, sigma, num_rej, goodvec=good, badvec=bad
      resm_im[xx,yy] = val

    endfor
    
  endfor
  writefits, path+fnout+'_cosmic_resistant_mean.fits', resm_im

endif

  
stop
end
