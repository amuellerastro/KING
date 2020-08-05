pro KING_Dark, path

print, ''
print, 'DARK'
print, ''

fileo = file_search(path+'KING.*fits', count=nfiles)

catg = strarr(nfiles)
object = strarr(nfiles)
dit = strarr(nfiles)
filter = strarr(nfiles)

st = strarr(nfiles)

for i=0,nfiles-1 do begin

  hdr = headfits(fileo[i], exten=0, /silent)
  object[i] = strcompress(sxpar(hdr,'OBJECT'),/rem)
  dit[i] = strcompress(sxpar(hdr,'EXPTIME'),/rem)
  catg[i] = strcompress(sxpar(hdr,'CATG'),/rem)
  filter[i] = strcompress(sxpar(hdr,'FILTER'),/rem)
  st[i] = dit[i]+'_'+filter[i]

endfor


idx = where(catg eq 'CALIB' and object eq 'DARK')
fileo = fileo[idx]
catg = catg[idx]
object = object[idx]
dit = dit[idx]
filter = filter[idx]
st = st[idx]

idx = uniq(st, sort(st))
ust = st[idx]
udit = dit[idx]
nu = n_elements(ust)

bias = mrdfits(path+'../Reduced/Bias.fits', 0, /silent)

for xx=0,n_elements(udit)-1 do begin

  idx = where(dit eq udit[xx])
  
  nframes = n_elements(idx)

  hdrraw = headfits(fileo[idx[0]], exten=0, /silent)
  naxis1 = float(sxpar(hdrraw, 'NAXIS1'))
  naxis2 = float(sxpar(hdrraw, 'NAXIS2'))

  im = dblarr(naxis1, naxis2, nframes)
  for i=0,nframes-1 do im[*,*,i] = mrdfits(fileo[idx[i]], 0, hdr, /silent, /fscale)

  if (naxis1 gt naxis2) then im = im[0:naxis2-1,*,*]
  if (naxis1 lt naxis2) then im = im[*,0:naxis1-1,*]

  for i=0,nframes-1 do im[*,*,i] = im[*,*,i]-bias

  ;if (n_elements(im[0,0,*]) gt 1) then medim = median(im, dim=3, /even) else medim = im
  if (n_elements(im[0,0,*]) gt 1) then medim = mean(im, dim=3) else medim = im
  
  writefits, path+'../Reduced/MasterDark_'+udit[xx]+'s.fits', medim, hdrraw

endfor


end
