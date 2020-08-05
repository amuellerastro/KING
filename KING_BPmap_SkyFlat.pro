pro KING_BPmap_SkyFlat, path, display, nodark

fileo = file_search(path+'KING.*fits', count=nfiles)

catgo = strarr(nfiles)
objecto = strarr(nfiles)
dito = strarr(nfiles)
filtero = strarr(nfiles)


sto = strarr(nfiles)

for i=0,nfiles-1 do begin

  hdr = headfits(fileo[i], exten=0, /silent)
  objecto[i] = strcompress(sxpar(hdr,'OBJECT'),/rem)
  dito[i] = strcompress(sxpar(hdr,'EXPTIME'),/rem)
  catgo[i] = strcompress(sxpar(hdr,'CATG'),/rem)
  filtero[i] = strcompress(sxpar(hdr,'FILTER'),/rem)
  sto[i] = dito[i]+'_'+filtero[i]

endfor

;=========================================================================

print, ''
print, 'BAD PIXEL MAP'
print, ''

;1st iteration
;Skyflat without BP map

idx = where(catgo eq 'CALIB' and objecto eq 'SKYFLAT')
fileo = fileo[idx]
catg = catgo[idx]
object = objecto[idx]
dit = dito[idx]
filter = filtero[idx]
st = sto[idx]

idx = uniq(st, sort(st))
ust = st[idx]
udit = dit[idx]
nu = n_elements(ust)

bias = mrdfits(path+'../Reduced/Bias.fits', 0, /silent)
if (nodark ne 1) then dark = mrdfits(path+'../Reduced/MasterDark_'+udit+'s.fits', 0, /silent)

nframes = n_elements(fileo)

hdrraw = headfits(fileo[idx[0]], exten=0, /silent)
naxis1 = float(sxpar(hdrraw, 'NAXIS1'))
naxis2 = float(sxpar(hdrraw, 'NAXIS2'))
ff = dblarr(naxis1, naxis2, nframes)

jd = dblarr(nframes)
for i=0,nframes-1 do begin

  ff[*,*,i] = mrdfits(fileo[i], 0, hdr, /silent, /fscale)
  if (i eq 0) then hdrraw = hdr
  dateobs = sxpar(hdr, 'DATE-OBS')
  t1 = double(strmid(dateobs, 0, 4))
  t2 = double(strmid(dateobs, 5, 2))
  t3 = double(strmid(dateobs, 8, 2))
  t4 = double(strmid(dateobs, 11, 2))
  t5 = double(strmid(dateobs, 14, 2))
  t6 = double(strmid(dateobs, 17, 2))
  exptime = sxpar(hdr,'EXPTIME')
  jd[i] = julday(t2,t3,t1,t4,t5,t6)+(exptime/2.)/86400.d0
  
endfor  
  
if (naxis1 gt naxis2) then ff = ff[0:naxis2-1,*,*]
if (naxis1 lt naxis2) then ff = ff[*,0:naxis1-1,*]

if (nodark ne 1) then for i=0,nframes-1 do ff[*,*,i] = ff[*,*,i]-bias-dark else for i=0,nframes-1 do ff[*,*,i] = ff[*,*,i]-bias


nx = n_elements(ff[*,0,0])
ny = nx

slope = dblarr(nx,ny)
for xx=0,nx-1 do begin

  for yy=0,ny-1 do begin

    sixlin, jd, ff[xx,yy,*], a, siga, b, sigb
    slope[xx,yy] = b[0]

  endfor
  
  proceeding_text,loop=(nx), i=xx, prompt='> x direction   '+string(xx+1,form='(I4)')
  
endfor

slope = slope/median(slope)

;=========================================================================

if (nodark ne 1) then begin

  ;tmp = max([double(udit), double(dito)],idxmax)
  ;if (idxmax eq 0) then maxdit = udit else maxdit = dito
  maxdit = max(double(dito),idxmax)
  dark = mrdfits(path+'../Reduced/MasterDark_'+dito[idxmax]+'s.fits', 0, /silent)
  dark = dark-bias
  
endif

bp1 = intarr(nx,ny)
bp1[*,*] = 1	;ratio FF BP
bp2 = bp1	;slope flux BP
bp3 = bp1	;linear fit FF BP
if (nodark ne 1) then bp4 = bp1	;dark BP

ff2 = ff[*,*,0]
ff1 = ff[*,*,nframes-1]


resistant_mean, slope, 2., t1, sigma, num_rej, goodvec=good, badvec=bad
idxgood = array_indices(bp2, good)
for i=0L,n_elements(idxgood[0,*])-1 do bp2[idxgood[0,i], idxgood[1,i]] = 0

if (nodark ne 1) then begin

  resistant_mean, dark, 3., t1, sigma, num_rej, goodvec=good, badvec=bad
  idxgood = array_indices(bp4, good)
  for i=0L,n_elements(idxgood[0,*])-1 do bp4[idxgood[0,i], idxgood[1,i]] = 0

endif

if (nodark ne 1) then bp = bp2+bp4 else bp = bp2

idx0 = where(bp gt 0.)
idx = array_indices(bp, idx0)
for i=0L,n_elements(idx[0,*])-1 do bp[idx[0,i], idx[1,i]] = 1

writefits, path+'../Reduced/'+'BadPixelMap.fits', bp

if keyword_set(display) then begin

  window, 2, xs=nx, ys=ny
  cgimage, bp, stretch=2

endif

print, ''
print, 'Percantage of bad Pixel: ', 100.*double(n_elements(where(bp eq 1.)))/double(n_elements(bp))
print, ''

;=========================================================================

print, ''
print, 'SKY FLAT'
print, ''

;re-compute flat field with new bad pixel map

slope = dblarr(nx,ny)
for xx=0,nx-1 do begin

  for yy=0,ny-1 do begin

    if (bp[xx,yy] ne 1.) then begin

      sixlin, jd, ff[xx,yy,*], a, siga, b, sigb
      slope[xx,yy] = b[0]
  
    endif else begin
    
      slope[xx,yy] = 0.
    
    endelse
  
  
  endfor
  
  proceeding_text,loop=(nx), i=xx, prompt='> x direction   '+string(xx+1,form='(I4)')
  
endfor

;median slope
idx = where(bp eq 0.d0)
msl = median(slope[idx], /even)

ff_final = slope/msl

; tmp = ff_final
; fixpix_mod, tmp, bp, outim, npix=24, /weight, /silent
; ff_final = outim

writefits, path+'../Reduced/'+'MasterFlat.fits', ff_final, hdrraw

if keyword_set(display) then begin

  window, 0
  !p.multi=[0,1,0]

  plot, slope, psym=3, /yn, xst=1
  !p.multi=[0,1,0]

  window, 2, xs=nx, ys=ny
  cgimage, ff_final, stretch=2

  window, 1, xs=500, ys=500
  cghistoplot, ff_final, xr=[0.5,1.2], bin=0.003, yr=[0,5.d4]

endif

end
