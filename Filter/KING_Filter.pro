pro KING_Filter

tr = dblarr(180,6)

fn = 'Johnson_U'
readcol, 'UG1.txt', l1, t1, format='d,d'
readcol, 'BG40.txt', l2, t2, format='d,d'

wl = l1
tr[*,0] = t1*t2


fn = 'Johnson_B'
readcol, 'GG395.txt', l1, t1, format='d,d'
readcol, 'BG25.txt', l2, t2, format='d,d'
readcol, 'BG39.txt', l3, t3, format='d,d'
tr[*,1] = t1*t2*t3

fn = 'Johnson_V'
readcol, 'GG495.txt', l1, t1, format='d,d'
readcol, 'BG39.txt', l2, t2, format='d,d'
tr[*,2] = t1*t2

fn = 'Johnson_R'
readcol, 'OG570.txt', l1, t1, format='d,d'
readcol, 'KG3.txt', l2, t2, format='d,d'
tr[*,3] = t1*t2

fn = 'Johnson_I'
readcol, 'RG9.txt', l1, t1, format='d,d'
tr[*,4] = t1

fn = 'SDSS_z'
readcol, 'RG830.txt', l1, t1, format='d,d'
tr[*,5] = t1

window, 0, xs=1200, ys=600
!p.thick=3
plot, wl, tr[*,0], /nodata, xtitle='Wavelenth [nm]', ytitle='Transmission', xr=[300,1100], xst=1, yr=[0,1], yst=1, charsize=2
oplot, wl, tr[*,0], color=cgcolor('purple')
oplot, wl, tr[*,1], color=cgcolor('blue')
oplot, wl, tr[*,2], color=cgcolor('green')
oplot, wl, tr[*,3], color=cgcolor('orange')
oplot, wl, tr[*,4], color=cgcolor('red')
oplot, wl, tr[*,5], color=cgcolor('violet red')
xyouts, 350,0.9,'U', charsize=3, alignment=0.5, color=cgcolor('purple')
xyouts, 425,0.9,'B', charsize=3, alignment=0.5, color=cgcolor('blue')
xyouts, 475,0.9,'V', charsize=3, alignment=0.5, color=cgcolor('green')
xyouts, 600,0.9,'R', charsize=3, alignment=0.5, color=cgcolor('orange')
xyouts, 800,0.9,'I', charsize=3, alignment=0.5, color=cgcolor('red')
xyouts, 1000,0.9,'SDSS z', charsize=3, alignment=0.5, color=cgcolor('violet red')

oplot, 372.*[1,1], !y.crange
xyouts, 372,0.8,'OII', charsize=3, alignment=0
oplot, 501.*[1,1], !y.crange
xyouts, 501,0.8,'OIII', charsize=3, alignment=0
oplot, 656.*[1,1], !y.crange
xyouts, 656,0.8,'Halpha', charsize=3, alignment=0
oplot, 673.*[1,1], !y.crange
xyouts, 673,0.7,'SII', charsize=3, alignment=0
oplot, 700.*[1,1], !y.crange
xyouts, 700,0.6,'SII cont.', charsize=3, alignment=0

stop
end
