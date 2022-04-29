
;; NAME:
;   lonlat2utm
;
; PURPOSE:
;
;      Convert longitude and latitude from WGS 84 to UTM (WGS84) from
;       http://www.uwgb.edu/dutchs/UsefulData/UTMFormulas.HTM
;
; AUTHOR:
;      Céline Danilo celine.danilo@gmail.com
;
;
; INPUTS:
;
;       long : longitude array in degrees decimal
;       lat : latitude array in degrees decimal
; OUTPUTS :
;       utm : a structure including East and North component and the
;       number of the UTM zone


function lonlat2utm, long, lat ; lon lat en degre

a = double(6378137) ; rayon à l'équateur expression en m
b = double(6356752.3142) ; rayon au pole en m 
k0 = 0.9996

e = sqrt(1.-b^2/a^2) ; excentricité de la Terre
e_p2 = (e*a/b)^2
eccSquared = e^2.

;//Make sure the longitude is between -180.00 .. 179.9
p = where(abs(Long) gt 180 )
if p[0] ne -1 then LongTemp = (Long+180.)-round((Long+180.)/360.)*360.-180. else LongTemp = Long; // -180.00 .. 179.9;
p= -1

LatRad = Lat*!dtor               ;
LongRad = LongTemp*!dtor         ;

ZoneNumber = ceil((temporary(LongTemp) + 180.)/6.); +1 ; cic normalement égale à 22

LongOrigin = (ZoneNumber-1)*6 - 180 + 3    ; //+3 puts origin in middle of zone
LongOriginRad = LongOrigin *!dtor        ;

;//compute the UTM Zone from the latitude and longitude

eccPrimeSquared = (eccSquared)/(1.-eccSquared) ;

N = a/sqrt(1-eccSquared*sin(LatRad)*sin(LatRad))    ;
T = tan(LatRad)*tan(LatRad)                         ;
C = eccPrimeSquared*cos(LatRad)*cos(LatRad)         ;
AA = cos(LatRad)*(LongRad-LongOriginRad)             ;

LonOrigin = 0
LonOriginRad =0
Longrad = 0

M = a*((1.-eccSquared/4.-3.*eccSquared*eccSquared/64.- 5.*eccSquared*eccSquared*eccSquared/256.)*LatRad $
       -(3.*eccSquared/8.+ 3.*eccSquared*eccSquared/32.+45.*eccSquared*eccSquared*eccSquared/1024.)*sin(2.*LatRad)$
       +(15.*eccSquared*eccSquared/256. + 45.*eccSquared*eccSquared*eccSquared/1024.)*sin(4.*LatRad) $
       -(35.*eccSquared*eccSquared*eccSquared/3072.)*sin(6.*LatRad));

UTMNorthing = double(k0*(M+N*tan(LatRad)*(AA^2/2.+(5.-T+9.*C+4.*C*C)*AA^4/24. $
                                          + (61.-58.*T+T*T+600.*C-330.*eccPrimeSquared)*AA^6/720.))) ;
M = 0

UTMEasting = double(k0*N*(AA+(1.-T+C)*AA^3/6. $
                          + (5.-18.*T+T*T+72*C-58.*eccPrimeSquared)*AA^5/120.)$
                          + 500000.0) ;


AA = 0
C = 0
N = 0
T = 0

if(Lat[0] lt 0) then UTMNorthing += 10000000.0       ; //10000000 meter offset for southern hemisphere

ec_zone = stddev(ZoneNumber)
if ec_zone  eq 0 then ZoneNumber = mean(ZoneNumber)

utm = CREATE_STRUCT('E', UTMEasting, 'N', UTMNorthing,'zone',ZoneNumber) 

return, utm

end
