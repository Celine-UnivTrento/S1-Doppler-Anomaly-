function rad_profil_std, array, flag_array, flag_value, array_def_val, std_array, max_std=max_std

  ss =size(array)
   
  nrad = ss[1]
  nazim = ss[2]
  
  rad_prof = fltarr(nazim)
  rad_prof_std = fltarr(nazim)
  rad_prof_num = fltarr(nazim)
    
  if n_elements(max_std) eq 0 then max_std = 50

  for kim = 0L, nazim-1 do begin
     
     pf = where(flag_array[*, kim] eq flag_value and abs(std_array[*,kim]) lt max_std and array[*, kim] ne array_def_val)
    
     if pf[0] ne -1 and n_elements(pf) gt 1 then begin
   
       std_rvl_f = stddev(array[ pf,kim])	      
       mean_rvl_f =  mean(array[ pf,kim])
       median_rvl_f = median(array[ pf,kim])
       
       ;; on cherche a enlever les extremum de retrodiffusion et de vitesse radiale
       
       p_calcul_dopt =  where(flag_array[ *,kim] eq flag_value and $
                              abs(array[*,kim]-mean_rvl_f) lt std_rvl_f*2. and $
                              abs(std_array[*,kim]) lt max_std  and array[*, kim] ne array_def_val)  
       
       if n_elements(p_calcul_dopt) gt 3 then begin
          rad_prof[ kim] = median(array[ p_calcul_dopt, kim ])
          rad_prof_num[ kim] = n_elements(p_calcul_dopt)
          rad_prof_std[ kim] = stddev(array[p_calcul_dopt, kim ])
          
       endif
       
       if abs(rad_prof[kim]) lt 0.0000001 then rad_prof[ kim] = 0.
       
       p_calcul_dopt =-1
       
    endif
     pf =-1
     
  endfor
                                
  return, {rad_prof: rad_prof, rad_prof_std:rad_prof_std, rad_prof_num:rad_prof_num} 
end
