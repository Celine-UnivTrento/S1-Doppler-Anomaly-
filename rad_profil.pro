function rad_profil, array, flag_array, flag_value, array_def_val, ancilla_array, min_nrcs = min_nrcs, max_nrcs = max_nrcs

  ss =size(array)
  ;array_def_val = -999.00
  
  nrad = ss[1]
  nazim = ss[2]
  
  rad_prof = fltarr(nazim)
  rad_prof_std = fltarr(nazim)
  rad_prof_num = fltarr(nazim)
  
  mean_nrcs= fltarr(nazim)
  median_nrcs= fltarr( nazim)
  std_nrcs = fltarr( nazim)
  
  if n_elements(min_nrcs) eq 0 then min_nrcs = - 5
  if n_elements(max_nrcs) eq 0 then max_nrcs = - 0.5

  for kim = 0L, nazim-1 do begin
  
    pf = where(flag_array[*, kim] eq flag_value and alog(ancilla_array[*,kim]) gt min_nrcs and alog(ancilla_array[*,kim]) lt max_nrcs); and array[*, kim] ne array_def_val)
    
    if pf[0] ne -1 and n_elements(pf) gt 1 then begin
      mean_ncrs_f = mean(ancilla_array[ pf,kim])
      median_ncrs_f = median(ancilla_array[pf,kim])
      std_ncrs_f = stddev(ancilla_array[ pf,kim])
      std_rvl_f = stddev(array[ pf,kim])	      
      mean_rvl_f =  mean(array[ pf,kim])
      median_rvl_f = median(array[ pf,kim])

      ; on cherche a enlever les extremum de retrodiffusion et de vitesse radiale
      p = where( flag_array[*, kim] eq flag_value and abs(ancilla_array[*,kim]-median_ncrs_f) lt std_ncrs_f*2. and $
		 alog(ancilla_array[*,kim]) gt min_nrcs and alog(ancilla_array[*,kim]) lt max_nrcs );and abs(array[ *, kim]-median_rvl_f) lt std_rvl_f*1.)
      ;p = where( flag_array[*, kim] eq flag_value and (ancilla_array[*,kim]-median_ncrs_f) lt std_ncrs_f)
      if p[0] ne -1 and n_elements(p) gt 1 then begin
        mean_nrcs[kim]   = mean(ancilla_array[ p,kim])
        median_nrcs[kim] = median(ancilla_array[p,kim])
        std_nrcs[ kim]   = stddev(ancilla_array[ p,kim])
        mean_radvelec    = median(array[ p , kim ])
        std_radvelec     = stddev(array[ p , kim ])
        
        ;p_calcul_dopt =  where(flag_array[ *,kim] eq flag_value and $
        ;  abs(ancilla_array[*,kim]-median_nrcs[kim]) lt std_nrcs[kim]*2 and $
        ;  array[ *, kim] ne 0. )
        p_calcul_dopt =  where(flag_array[ *,kim] eq flag_value and $
          abs(ancilla_array[*,kim]-median_nrcs[kim]) lt std_nrcs[kim]*1. and $
          abs(array[*,kim]-mean_radvelec) lt std_radvelec*1. )  
                              
        if (p_calcul_dopt[0] ne -1 and  n_elements(p_calcul_dopt) gt 3) then begin
          rad_prof[ kim] = median(array[ p_calcul_dopt, kim ])
          rad_prof_num[ kim] = n_elements(p_calcul_dopt)
          rad_prof_std[ kim] = stddev(array[p_calcul_dopt, kim ])
          
;          if abs(mean_radvelec - rad_prof[kim]) gt 0.2 then begin
;            print, mean_radvelec , rad_prof[kim], n_elements(p_calcul_dopt), n_elements(p)
;            stop
;          endif
          
        endif
        
        if abs(rad_prof[kim]) lt 0.0000001 then rad_prof[ kim] = 0.
        
        p_calcul_dopt =-1
        
      endif
      p = -1
    
    endif
    pf =-1
    
  endfor
  

  return, {rad_prof : rad_prof, rad_prof_std:rad_prof_std, rad_prof_num:rad_prof_num, median_nrcs : median_nrcs}
  
end
