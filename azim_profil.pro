function azim_profil, array, flag_array, flag_value, def_val, std_array, lim_std

  ss =size(array)
  nrad = ss[1]
  nazim = ss[2]  
  ss = 0
  
  azim_prof = fltarr(nrad)
  azim_prof_std = fltarr(nrad)
  azim_prof_num = fltarr(nrad)
   
  for jim = 0L, nrad-1 do begin
  
    p = where(flag_array[jim,*] eq flag_value and array[jim,*] ne def_val and std_array[jim,*] lt lim_std)
   
    if p[0] ne -1 then begin
      dop_median_jim = median(array[jim, p])
      dop_std_jim = stddev(array[jim,p])
      pp = where(flag_array[jim,*] eq flag_value and $
                 array[jim,*] ne def_val and $
                 std_array[jim,*] lt lim_std and $
                 abs(array[jim,*]-dop_median_jim) lt 2*dop_std_jim )
   
      if pp[0] ne -1 then begin
        azim_prof[jim] = mean(array[jim, pp])
        azim_prof_std[jim] = stddev(array[jim,pp])
        azim_prof_num[jim] = n_elements(pp)       
      endif
      pp = -1
    endif
    p = -1
    
  endfor
 
  return, {azim_prof: azim_prof, azim_prof_std:azim_prof_std, azim_prof_num:azim_prof_num}
  
  
end
