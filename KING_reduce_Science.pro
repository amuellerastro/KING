pro KING_reduce_Science, path, nodark, sout

print, ''
print, 'SCIENCE'
print, ''

fileo = file_search(path+'KING.*fits', count=nfiles)

catg = strarr(nfiles)
object = strarr(nfiles)
dit = strarr(nfiles)
filter = strarr(nfiles)
jd = dblarr(nfiles)

st = strarr(nfiles)

for i=0,nfiles-1 do begin

  hdr = headfits(fileo[i], exten=0, /silent)
  object[i] = strcompress(sxpar(hdr,'OBJECT'),/rem)
  dit[i] = strcompress(sxpar(hdr,'EXPTIME'),/rem)
  catg[i] = strcompress(sxpar(hdr,'CATG'),/rem)
  filter[i] = strcompress(sxpar(hdr,'FILTER'),/rem)
  jd[i] = double(sxpar(hdr,'MID_JD'))
  st[i] = dit[i]+'_'+filter[i]

endfor


idx = where(catg eq 'SCIENCE')
fileo = fileo[idx]
catg = catg[idx]
object = object[idx]
dit = dit[idx]
filter = filter[idx]
st = st[idx]
jd = jd[idx]

idx = uniq(st, sort(st))
ust = st[idx]
udit = dit[idx]
ufilter = filter[idx]
nu = n_elements(ust)

bias = mrdfits(path+'../Reduced/Bias.fits', 0, /silent)
if (nodark ne 1) then dark = mrdfits(path+'../Reduced/MasterDark_'+udit+'s.fits', 0, /silent)
bp = mrdfits(path+'../Reduced/BadPixelMap.fits', 0, /silent)
ff = mrdfits(path+'../Reduced/MasterFlat.fits', 0, /silent)

nframes = n_elements(fileo)

hdrraw = headfits(fileo[idx[0]], exten=0, /silent)
naxis1 = float(sxpar(hdrraw, 'NAXIS1'))
naxis2 = float(sxpar(hdrraw, 'NAXIS2'))
im = dblarr(naxis1, naxis2, nframes)

for i=0,nframes-1 do im[*,*,i] = mrdfits(fileo[i], 0, hdr, /silent, /fscale)

if (naxis1 gt naxis2) then im = im[0:naxis2-1,*,*]
if (naxis1 lt naxis2) then im = im[*,0:naxis1-1,*]
naxis1 = n_elements(im[*,0,0])
naxis2 = n_elements(im[0,*,0])

if (naxis1 eq 1024) then begin

 im = im[512-400:512+399,512-400:512+399,*]
 ff = ff[512-400:512+399,512-400:512+399,*]
 bias = bias[512-400:512+399,512-400:512+399]
 bp = bp[512-400:512+399,512-400:512+399]
 if (nodark ne 1) then dark = dark[512-400:512+399,512-400:512+399]
  
endif else begin

  print, 'Adopt for different binning!'
  stop

endelse


corim = im
corim2 = corim
for i=0,nframes-1 do begin

  if (nodark ne 1) then corim[*,*,i] = (im[*,*,i]-dark-bias)/ff else corim[*,*,i] = (im[*,*,i]-bias)/ff
  tmp = corim[*,*,i]
  fixpix_mod, tmp, bp, outim, npix=24, /weight, /silent
  corim[*,*,i] = outim

  proceeding_text,loop=nframes, i=i, prompt='> Frame   '+string(i+1,form='(I4)')
  
endfor

sout = object[0]+'_'+udit+'s_'+ufilter+'_reduced.fits'
writefits, path+'../Reduced/'+sout, corim, hdrraw
save, jd, filename=path+'../Reduced/JD.sav'

end
