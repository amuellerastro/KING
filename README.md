
Weather:
http://calsky.com amueller / Ipot2k01
https://www.meteoblue.com/en/weather/forecast/meteogramfive/k%c3%b6nigstuhl_germany_2885754
http://meteox.com/forecastloop.aspx?type=4&continent=europa
http://www.myweather2.com/forecastcloud/player.aspx?fc=3
https://en.sat24.com/en/de/visual
https://www.lsw.uni-heidelberg.de/wetter/?lang=de
http://www.wetter.com/deutschland/heidelberg/koenigstuhl/DE0004329012.html

KING
https://svn.mpia.de/trac/gulli/king/
https://svn.mpia.de/trac/gulli/king/wiki/WikiStart/FP30Doc
https://svn.mpia.de/trac/gulli/king/blog

Image processing
http://billsnyderastrophotography.com/?page_id=1767
https://free-astro.org/index.php?title=Image_Processing:Main#Siril

Comets:
http://news.astronomie.info/sky201710/kometen.html
https://www.calsky.com/cs.cgi/Comets/1
https://theskylive.com/comets

Gaia or other spacecrafts:
https://ssd.jpl.nasa.gov/horizons.cgi#top
http://projects.familie-steinel.de/stellarium-comet-jpl/

https://astronomy.tools/calculators/ccd

https://www.eso.org/sci/observing/tools/standards/Landolt.html

Pointing model:
kappa70 -> kappa_gui
offsets written in /opt/kappa70/data/models.dat (last 2 parameters are offsets in sec/arcsec)

To Do:
Webcam check
GAL 192.16-3.82
M15 with planetary nebulae
R Mon (Hubble variable nebula)
Gaia, GRB, SN, transiting planets

GPS:
N 49° 23' 43.20"
E 8° 43' 25.20"
altitude: 615m

calsky:
49°23'44.01" N, 8°43'25.15" E
Google Earth
49°23'43.48"N, 8°43'25.35"E

Focus:
-----
V: +62
R: +47
Ha: +41

Observation and reduction
KING_Obs
KING_sort_data
KING_reduce
KING_OrbitCorrect    provide coordinates of moving object
KING_track_SSobject    manual selection of moving object

root pwd: xFresco5

CCD failure
-----------
"Reading out...ERROR: unable to allocate Memory for Buffer1 (dataport=0) size 0x0
ERROR: No buffer allocated: all datasize were zero !"
-Clean restart of PC or as root: echo 1 > /proc/sys/vm/drop_caches
If images are black then restart of electronics:
1. ccdstop
2. switch off Instrumentierung at the console
3. wait 10 seconds
4. switch on Instrumentierung
5. ccdstart


Verbesserungen:
==============

webcam oder zeitaufloesende Kamera
neuer Flat Field Schirm
Lichtschutz um Schreibtisch und 180deg gedreht

# KING
# KING
# KING
