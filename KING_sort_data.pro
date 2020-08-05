@readcol.pro
@remchar.pro
@gettok.pro
@strsplit.pro
@strnumber.pro
@uniq.pro
@headfits.pro
@fxposit.pro
@mrd_hread.pro
@sxpar.pro

pro KING_sort_data

base = '/home/amueller/work/KING/'
;night = '20180222'
night = ''
read, 'Night YYYMMDD: ', night

nightpath = base+night+'/'
readcol, nightpath+night+'.log', file, catg, object, jd, dit, filter, ccdtemp, format='a,a,a,d,a,a,d', /silent

st = dit+'_'+filter

idx = where(catg eq 'SCIENCE')

sfile = file[idx]
sst = st[idx]
osst = object[idx]+'_'+sst
sdit = dit[idx]
sfilter = filter[idx]
uosst = osst[uniq(osst, sort(osst))]
usst = sst[uniq(sst, sort(sst))]
usdit = sdit[uniq(sdit, sort(sdit))]
usfilter = sfilter[uniq(sfilter, sort(sfilter))]

path = strarr(n_elements(uosst))
for xx=0,n_elements(uosst)-1 do begin

  path0 = nightpath+uosst[xx]+'/'
  spawn, 'mkdir -p '+path0
  path[xx] = path0+'Raw/'
  spawn, 'mkdir -p '+path[xx]

  idx = where(osst eq uosst[xx])
  for i=0,n_elements(idx)-1 do spawn, 'cp '+nightpath+sfile[idx[i]]+' '+path[xx]

  idx = where(catg eq 'CALIB' and object eq 'BIAS')
  for i=0,n_elements(idx)-1 do spawn, 'cp '+nightpath+file[idx[i]]+' '+path[xx]

  pos1 = strpos(uosst[xx], '_')
  pos2 = strpos(uosst[xx], '_', /reverse_search)
  tdit = strmid(uosst[xx], pos1+1, pos2-pos1-1)
  tfilt = strmid(uosst[xx], pos2+1, strlen(uosst[xx])-pos2-1)

  ;darks for science
  idx = where(tdit eq dit and catg eq 'CALIB' and object eq 'DARK')
  if (idx[0] ne -1) then for i=0,n_elements(idx)-1 do spawn, 'cp '+nightpath+file[idx[i]]+' '+path[xx]

  ;flats for science
  idxff = where(tfilt eq filter and catg eq 'CALIB' and object eq 'SKYFLAT')
  for i=0,n_elements(idxff)-1 do spawn, 'cp '+nightpath+file[idxff[i]]+' '+path[xx]

  ;darks for flats
  pos1 = strpos(st[idxff], '_')
  tfdit = strarr(n_elements(idxff))
  for i=0,n_elements(idxff)-1 do tfdit[i] = strmid(st[idxff[i]], 0, pos1[i])
  for i=0,n_elements(tfdit)-1 do begin

    idx = where(dit eq tfdit[i] and catg eq 'CALIB' and object eq 'DARK')
    if (idx[0] ne -1) then for j=0,n_elements(idx)-1 do spawn, 'cp '+nightpath+file[idx[j]]+' '+path[xx]

  endfor

endfor

;create overview
;doing this in separate loop to not overwrite variables

for xx=0,n_elements(uosst)-1 do begin

  file = file_search(path[xx]+'KING*fits', count=n)

  out = file
  pos1 = strpos(file, '/KING.', /reverse_search)
  for i=0,n_elements(file)-1 do out[i] = strmid(file[i], pos1[i]+1, strlen(file[i])-pos1[i]-1)


  catg = strarr(n_elements(file))
  object = strarr(n_elements(file))
  dit = dblarr(n_elements(file))
  filter = strarr(n_elements(file))

  for i=0,n-1 do begin

    hdr = headfits(file[i], exten=0, /silent)
    catg[i] = strcompress(sxpar(hdr, 'CATG'),/rem)
    object[i] = strcompress(sxpar(hdr, 'OBJECT'),/rem)
    dit[i] = sxpar(hdr, 'EXPTIME')
    filter[i] = strcompress(sxpar(hdr, 'FILTER'),/rem)

  endfor

  openw, lun, path[xx]+'files.txt', width=1400, /get_lun

    printf, lun, '                         File        CATG         OBJECT     DIT FILTER'
    printf, lun, '-----------------------------------------------------------------------'

    for i=0,n_elements(file)-1 do printf, lun, out[i], catg[i], object[i], dit[i], filter[i], format='(a29, a12, a15, f8.1, a7)'

  close, lun
  free_lun, lun

endfor


stop
end
