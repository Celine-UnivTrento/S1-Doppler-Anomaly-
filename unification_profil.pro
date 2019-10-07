function unification_profil, profils, std_profils, num_profils, deriv_num_profils, num_min 

  ss = size(profils)
  nx  = ss[2]
  nprofil = ss[1]
  
  profil_unifie = fltarr(nx)
  max_prfls = max(profils)
  min_prfld = min(profils)

  max_impose = 20.

  if max_prfls gt max_impose then print, ' ''From unification_profil.pro '' Carefull variation larger than expected. '

  p = where(abs(profils) gt 0.0001) 
  if p[0] ne -1 and n_elements(p) gt 3 then begin 
     std_prfls = stddev(profils[p])
     mean_prfls= mean(profils[p])
     median_prfls = median(profils[p])
  endif else begin
     std_prfls = stddev(profils)
     mean_prfls= mean(profils)
     median_prfls = median(profils)
  endelse

  max_profil_tolerat = min([median_prfls+3*std_prfls, max_impose]) 
  max_std_tolerat = std_prfls*2.
  
  print, ' ''From unification_profil.pro '' Extremum obtained for unification (profil and std_profil) ' , max_profil_tolerat ,max_std_tolerat
  
  for kim = 0L, nx-1 do begin  
     p = where( abs(profils[*,kim]) lt max_profil_tolerat and $
                abs(profils[*,kim]) gt 0.0001 and $                                
                num_profils[*,kim] ge num_min and $
                abs(deriv_num_profils[*,kim]) lt 20 )
     
     if n_elements(p) eq 1 and p[0] ne -1 then profil_unifie[kim] = profils[p,kim] 
     
     if n_elements(p) gt 1 then begin
        mean_dop = median(profils[p,kim])
        std_dop = stddev(profils[p,kim])
        pp = where(abs(profils[p,kim]-mean_dop) lt std_dop*3.)
        if pp[0] ne -1 and n_elements(pp) gt 1 then begin
           std = stddev(profils[p[pp],kim])
           if std le max_std_tolerat then profil_unifie[kim] = mean(profils[p[pp],kim])     
           
       endif
        if  pp[0] ne -1 and n_elements(pp) eq 1 then profil_unifie[kim] = profils[p[pp],kim]
     endif
  endfor
    
  return, profil_unifie
  
end
