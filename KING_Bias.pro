pro KING_Bias, path

print, ''
print, 'BIAS'
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

idx = where(catg eq 'CALIB' and object eq 'BIAS')
file = fileo[idx]
catg = catg[idx]
object = object[idx]
dit = dit[idx]
filter = filter[idx]
st = st[idx]

nframes = n_elements(file)

hdrraw = headfits(file[0], exten=0, /silent)
naxis1 = float(sxpar(hdrraw, 'NAXIS1'))
naxis2 = float(sxpar(hdrraw, 'NAXIS2'))

im = dblarr(naxis1, naxis2, nframes)
for i=0,nframes-1 do im[*,*,i] = mrdfits(file[i], 0, hdr, /silent, /fscale)

if (naxis1 gt naxis2) then im = im[0:naxis2-1,*,*]
if (naxis1 lt naxis2) then im = im[*,0:naxis1-1,*]

medim = median(im, dim=3, /even)

writefits, path+'../Reduced/Bias.fits', medim, hdrraw


end
