function improv_prof_via_deriv, prof, periode

  nazim = n_elements(prof)
  
  nb_periode = nazim/periode
  deriv_unif = fltarr(nazim)
  deriv_prof = fltarr(nazim)
  
  p = where(prof ne 0)
  if p[0] ne -1 and n_elements(p) gt 3 then deriv_prof[p] = deriv(prof[p])
  if p[0] eq -1 or periode eq 0 or n_elements(p) le 3  then return, prof
  p =-1
  
  cc_reg = 0.85;0.85
  cc_fai = 0.75;0.75
  
  affi = 0
  
  p = where(deriv_prof ne 0)
  if p[0] ne -1 then min_deriv_init = min(deriv_prof[p])
  if p[0] ne -1 then max_deriv_init = max(deriv_prof[p])
  p = -1
  
  ip = 0
  ind0 = findgen(periode)
  ind1 = findgen(periode)+periode
  indf = findgen(periode)+nazim-periode
  indf_1 = findgen(periode)+nazim-periode*2
  
  ccm = correlate(prof[ind0],prof[ ind1])
  if ccm ge cc_reg then begin
    p = where(deriv_prof[ind0] ne 0. and deriv_unif[ ind0] eq 0 )
    if p[0] ne -1 then deriv_unif[ ind0[p]] = deriv_prof[ ind0[p]]
    p = -1
    p = where(deriv_prof[ind1] ne 0. and deriv_unif[ ind1] eq 0 )
    if p[0] ne -1 then deriv_unif[ ind1[p]] = deriv_prof[ ind1[p]]
    p = -1
  endif
  
  ccp = correlate(prof[ indf],prof[ indf_1])
  if ccp gt cc_reg  then begin
    p = where( deriv_unif[ indf] eq 0 and deriv_prof[ indf] ne 0 )
    if p[0] ne -1 then deriv_unif[ indf[p]] = deriv_prof[ indf[p]]
    p = where( deriv_unif[ indf_1] eq 0 and deriv_prof[ indf_1] ne 0 )
    if p[0] ne -1 then deriv_unif[ indf_1[p]] = deriv_prof[ indf_1[p]]
  endif
  
  for ip = 1L, nb_periode-2 do begin
    ccm = correlate( reform(prof[ ip*periode:(ip+1)*periode-1]), reform(prof[ (ip+1)*periode:(ip+2)*periode-1]) )
    ccp = correlate( reform(prof[ (ip-1)*periode:ip*periode-1]), reform(prof[ ip*periode:(ip+1)*periode-1]) )
    cce = correlate( reform(prof[ (ip-1)*periode:ip*periode-1]), reform(prof[ (ip+1)*periode:(ip+2)*periode-1]) )
    
    if affi then print, ip, ccm, ccp, cce
    
    ind0 = findgen(periode)+periode*ip
    ind1 = findgen(periode)+periode*(ip+1)
    ind_1 = findgen(periode)+periode*(ip-1)
    
    if ccm ge cc_reg then begin
      p = where(deriv_prof[ind0] ne 0. and deriv_unif[ ind0] eq 0 )
      if p[0] ne -1 then deriv_unif[ind0[p]] = deriv_prof[ind0[p]]
      if affi then print, 'ccm > cc_reg :'
      if p[0] ne -1 and affi then print, ind0[p]
      p=-1
      p = where(deriv_prof[ind1] ne 0. and deriv_unif[ind1] eq 0.)
      if p[0] ne -1 then deriv_unif[ind1[p]] = deriv_prof[ind1[p]]
      if affi then print, 'ccm > cc_reg :'
      if p[0] ne -1 and affi then print, ind1[p]
      p=-1
    endif
    
    if ccp ge cc_reg then begin
      p = where(deriv_prof[ind0] ne 0. and deriv_unif[ind0] eq 0)
      if p[0] ne -1 then deriv_unif[ind0[p]] = deriv_prof[ind0[p]]
      if affi then print, 'ccp > cc_reg : '
      if p[0] ne -1 and affi then print, ind0[p]
      p=-1
      p = where(deriv_prof[ind_1] ne 0. and deriv_unif[ind_1] eq 0.)
      if p[0] ne -1 then deriv_unif[ind_1[p]] = deriv_prof_1[ind_1[p]]
     if affi then print, 'ccp > cc_reg'
      if p[0] ne -1 and affi then print, ind_1[p]
      p=-1
    endif
    
  endfor
  
  
  decal_ind = round(periode/2)
  for ip = 1L, nb_periode-3 do begin
    ccm = correlate(reform(prof[ ip*periode+decal_ind:(ip+1)*periode-1+decal_ind]), reform(prof[ (ip+1)*periode+decal_ind:(ip+2)*periode-1+decal_ind]))
    ccp = correlate(reform(prof[ ip*periode+decal_ind:(ip+1)*periode-1+decal_ind]), reform(prof[ (ip-1)*periode+decal_ind:ip*periode-1+decal_ind]))
    cce = correlate(reform(prof[ (ip-1)*periode+decal_ind:ip*periode-1+decal_ind]), reform(prof[ (ip+1)*periode+decal_ind:(ip+2)*periode-1+decal_ind]))
    
    if affi then print, ' decal'
    if affi then print, ip, ccm, ccp, cce
  ;  stop
    
    ind0 = findgen(periode)+periode*ip+decal_ind
    ind1 = findgen(periode)+periode*(ip+1)+decal_ind
    ind_1 = findgen(periode)+periode*(ip-1)+decal_ind
    
    if ccm ge cc_reg  then begin
      p = where(deriv_prof[ind0] ne 0. and deriv_unif[ind0] eq 0. )
      if p[0] ne -1 then deriv_unif[ind0[p]] = deriv_prof[ind0[p]]
      if affi then print, 'ccm decal > cc_reg'
      if p[0] ne -1 and affi then print, ind0[p]
      p=-1
      p = where(deriv_prof[ind1] ne 0. and deriv_unif[ind1] eq 0.)
      if p[0] ne -1 then deriv_unif[ind1[p]] = deriv_prof[ind1[p]]
      if affi then print,  'ccm decal > cc_reg'
      if p[0] ne -1 and affi then print, ind1[p]
      p=-1
    endif
    
    if ccp ge cc_reg then begin
      p = where(deriv_prof[ind0] ne 0. and deriv_unif[ind0] eq 0)
      if p[0] ne -1 then deriv_unif[ind0[p]] = deriv_prof[ind0[p]]
      if affi then print, 'ccp decal > cc_reg '
      if p[0] ne -1 and affi then print, ind0[p]
      p=-1
      p = where(deriv_prof[ind_1] ne 0. and deriv_unif[ind_1] eq 0.)
      if p[0] ne -1 then deriv_unif[ind_1[p]] = deriv_prof[ind_1[p]]
      if affi then print, 'ccp decal > cc_reg '
      if p[0] ne -1 and affi then print, ind_1[p]
      p=-1      
    endif ;else print, 'C''est quoi ce border ?!'
    
  endfor
  
  for ip = 1L, nb_periode-2 do begin
  
    ind0 = findgen(periode)+periode*ip
    ind1 = findgen(periode)+periode*(ip+1)
    ind_1 = findgen(periode)+periode*(ip-1)
    
    ccm = correlate( prof[ ind0], prof[ ind1] )
    ccp = correlate( prof[ ind_1], prof[ ind0] )
    cce = correlate( prof[ ind_1], prof[ ind1] )
    
     if affi then print, ip, ccm, ccp, cce
    
    if ccm ge cc_fai and ccp ge cc_fai then begin
      p = where(deriv_unif[ind0] eq 0 and deriv_prof[ ind0] ne 0 and deriv_prof[ ind1] ne 0 and deriv_prof[ ind_1] ne 0)
      if p[0] ne -1 then deriv_unif[ind0[p]] = ( deriv_prof[ ind0[p]] + deriv_prof[ ind1[p]] + deriv_prof[ ind_1[p]] )/3
      if affi then print, 'ccm et ccp > cc_fai'
      if p[0] ne -1 and affi then print,  ind0[p]
      p= -1
    endif
    
    if cce ge cc_reg then begin
      p =where(deriv_unif[ ind0] eq 0 and deriv_prof[ ind_1] ne 0 and deriv_prof[ ind1] ne 0)
      if p[0] ne -1 then deriv_unif[ ind0[p]] = (deriv_prof[ ind_1[p]]+deriv_prof[ ind1[p]])/2.
      if affi then print, 'cce > cc_reg '
      if p[0] ne -1 and affi then print, ind0[p]
      p=-1
    endif
    
    if ccm ge cc_fai then begin
      p = where( deriv_unif[ ind0] eq 0 and  deriv_prof[ ind0] ne 0)
      if p[0] ne -1 then deriv_unif[ ind0[p]] = deriv_prof[ ind0[p]]
      if affi then print,'ccm > cc_fai'
      if p[0] ne -1 and affi then print, [ ind0[p]]
      p=-1
      p = where( deriv_unif[ ind1] eq 0 and  deriv_prof[ ind1] ne 0)
      if p[0] ne -1 then deriv_unif[ ind1[p]] = deriv_prof[ ind1[p]]
      if affi then print, 'ccm > cc_fai'
      if p[0] ne -1 and affi then print, [ ind1[p]]
      p= -1
      
    endif
    
    if ccp ge cc_fai then begin
      p = where(deriv_unif[ ind0] eq 0 and deriv_prof[ ind0] ne 0 )
      if p[0] ne -1 then deriv_unif[ ind0[p]] = deriv_prof[ ind0[p]]
      if affi then print, 'ccp > cc_fai'
      if p[0] ne -1 and affi then print, ind0[p]
      p =-1
      p = where(deriv_unif[ ind_1] eq 0 and deriv_prof[ ind_1] ne 0 )
      if p[0] ne -1 then deriv_unif[ ind_1[p]] = deriv_prof[ ind_1[p]]
      if affi then print, 'ccp > cc_fai'
      if p[0] ne -1 and affi then print, ind_1[p]
      p=-1
    endif
    
    if cce gt cc_fai then begin
      p = where(deriv_unif[ ind0] eq 0 and deriv_prof[ ind_1] ne  deriv_prof[ind1] ne 0 )
      if p[0] ne -1 then deriv_unif[ ind0[p]] = (deriv_prof[ ind_1[p]]+deriv_prof[ind1[p]])/2
      if affi then print,  'cce > cc_fai'
      if p[0] ne -1 and affi then print,ind0[p]
      p=-1
    endif
  endfor
  
  pp = where( deriv_unif eq 0 and deriv_prof ne 0 )
  if pp[0] ne -1 then begin
    ;stop
    disc = pp-shift(pp, -1)
    ;print, n_elements(pp)
    for ii = 0L, n_elements(pp)-1 do begin
      if pp[ii] ge periode+1 and pp[ii] le nazim-1-periode-1 then begin
        if deriv_unif[pp[ii]+periode] ne 0 and deriv_unif[pp[ii]-periode] ne 0 then deriv_unif[pp[ii]] = (deriv_unif[pp[ii]-periode]+deriv_unif[pp[ii]+periode])/2; (deriv_unif[pp[ii]-1]+deriv_unif[pp[ii]+1])/2
      endif
      if pp[ii] lt periode+1 then  begin
        if deriv_unif[pp[ii]+periode] ne 0 then  deriv_unif[pp[ii]] = deriv_unif[pp[ii]+periode]
      endif
      if pp[ii] gt nazim-1-periode-1 then begin
        if deriv_unif[pp[ii]-periode] ne 0 then  deriv_unif[pp[ii]] = deriv_unif[pp[ii]-periode]
      endif
    endfor
    
  endif
  
  deriv_unif = fill_gap_one_point_inavector(deriv_unif)
  
  pn = where(deriv_unif ne 0)
  if pn[0] ne -1 then min_deriv_unif = min(deriv_unif[pn]) else min_deriv_unif = 0
  if pn[0] ne -1 then max_deriv_unif = max(deriv_unif[pn]) else max_deriv_unif = 0
  
  ;deriv_unif = fill_gap_one_point_inavector(deriv_unif)
  ;for ii = periode, nazim-1-periode do begin ;
  ;	if abs(deriv_unif[ii]-deriv_unif[ii-periode]) gt min_deriv and abs(deriv_unif[ii]-deriv_unif[ii+periode]) gt min_deriv then deriv_unif[ii] = (deriv_unif[[ii]-periode]+deriv_unif[ii+periode])/2
  ; endfor
  
  pz = where(deriv_unif eq 0)
  pn = where(deriv_unif ne 0)
  deriv_unif_smooth = fltarr(nazim)
  deriv_smooth = deriv_unif
  if pn[0] ne -1 and n_elements(pn) gt nazim/3 then deriv_unif_smooth[pn] = smooth(deriv_unif[pn],3)
  if pz[0] ne -1 then deriv_unif_smooth[pz] = 0
  
  ;; continuite de l'integrale
  disc = pn-shift(pn, -1)
  qq = where(abs(disc) gt 1)
  
  ;; attention la periode n'est pas vraiment constante le long de l'azimut.... on ne peut pas utiliser une derivee ideale à partir de la periode
  ;moy_deriv = fltarr(periode)
  ;for ip=0L, periode-1 do begin
  ;  pp = where(deriv_unif[  ip:(nazim-1):periode] ne 0 and abs(prof[ ip:(nazim-1):periode]) gt 0.000001 )
  ;  if pp[0] ne -1 then moy_deriv[ip] = mean((deriv_unif[  ip:(nazim-1):periode])[pp]) else stop
  ;endfor
  ;	moy_deriv_long = fltarr(nazim)
  ;	for ip = 0L, nb_periode-1 do moy_deriv_long[ip*periode:(ip+1)*periode-1]=moy_deriv[*]
  
  integral = fltarr(nazim)
  if pn[0] ne -1 and qq[0] ne -1 then begin
    integral[pn[0]] = deriv_unif[pn[0]]+prof[pn[0]]
    ii=1
    while ii le n_elements(pn)-2 do begin
      while pn[ii]-pn[ii-1] eq 1 and ii le n_elements(pn)-2 do begin
        integral[pn[ii]] = integral[pn[ii]-1] +deriv_unif[ pn[ii]]
        ii = ii+1
      endwhile
      if pn[ii]-pn[ii-1] gt 1 then begin
        integral[pn[ii]] = deriv_unif[pn[ii]]+prof[pn[ii]]
      ;		print, 'Cas de discontinuite : ', pn[ii], pn[ii-1]
      endif
      ii = ii+1
    endwhile
  endif
  
  if affi then print, 'max_deriv_unif, min_deriv_unif : ', max_deriv_unif, min_deriv_unif, min_deriv_init, max_deriv_init
  
  p = where(integral ne 0 and prof ne 0 and deriv_prof lt max_deriv_unif and deriv_prof gt min_deriv_unif)
  if p[0] ne -1  then const_integration = mean(integral[p] - prof[p]) else const_integration = 0
  ;print, const_integration
  p=-1
  
  test = fltarr(nazim)
  p = where(prof ne 0 and abs(deriv_unif-deriv_prof) eq 0 ); or deriv_prof lt max_deriv_unif and deriv_prof gt min_deriv_unif ))
  if p[0] ne -1 then test[p] = prof[p]
  p= -1
  
  p = where((prof eq 0 or abs(deriv_unif-deriv_prof) gt 0.) and integral ne 0 );deriv_unif ne deriv_prof)
  if p[0] ne -1 then begin
    for ii = 0, n_elements(p)-1 do begin
      test[p[ii]] = integral[p[ii]] -const_integration
    ;if p[ii] ge 0 and p[ii] lt nazim-2 then begin
    ;if deriv_unif_smooth[p[ii]] ne 0 then test[p[ii]+1] = deriv_unif_smooth[p[ii]]+test[p[ii]]
    ;endif
      
    ;if p[ii] le nazim and p[ii] gt 0 then begin
      
    ;if deriv_unif_smooth[p[ii]] ne 0 then test[p[ii]-1] = test[p[ii]]-deriv_unif_smooth[p[ii]]
    ;endif
    endfor
  endif
  
  prof_new = prof*0
  p = where( prof ne 0 and deriv_prof eq deriv_unif); abs(integral-prof) lt 0.5)
  if p[0] ne -1 then prof_new[p] = prof[p]
  p = where( integral ne 0 and ( prof eq 0 or abs(deriv_prof-deriv_unif) gt 0.0001) ););abs(integral-prof) ge 0.5) )
  if p[0] ne -1 then prof_new[p] = integral[p]-const_integration
  p = -1
  
;  window, 6
;  plot, deriv_prof, col=-1, back=0, ytitle='Deriv profil', psym=-1
;  ;oplot, smooth(deriv(prof[ *]),3), col=200
;  oplot, deriv_unif, col=100, psym=-1
;  ;oplot, deriv(prof_new), col= 200
;  ;oplot, moy_deriv_long, col=100
;  
;  window,5
;  plot, prof, col=-1, back=0, ytitle='profil initial'
;  ; oplot, smooth(prof[ *],3), col=50
;  oplot, integral-const_integration, col=100
;  oplot, prof_new, col=200, thick=2
;  ;oplot, test, col=240
;  
;  stop
;  
  
  return, prof_new
  
end
