function lissage_profil, vector_brut, vector_brut_num, num_min, vector_flag, affi=affi, nb_pts_smooth_min= nb_pts_smooth_min, fill_gap=fill_gap

  if n_elements(affi) eq 0 then affi =0
  if n_elements(fill_gap) eq 0 then fill_gap =0
  
  ss = size(vector_brut)
  nx = ss[1]
  
  p = where(vector_flag eq 1 and vector_brut_num ge num_min and abs(vector_brut) gt 0.000001 )

  ; verif rapide
  mmo = 0.
  mme = 1.
  if p[0] ne -1 then begin 
    mme = mean(vector_brut[p])
    mmo = median(vector_brut[p])
    mms = stddev(vector_brut[p])
  endif 

  p = where(vector_flag eq 1 and vector_brut_num ge num_min and abs(vector_brut) gt 0.000001 and abs(vector_brut-mmo) lt 20./abs(mme-mmo))
  ; parametre de lissage fixé suivant le nombre de points present utilisable pour augmenter un maximum le lissage final
  ;if  n_elements(nb_pts_smooth_min) eq 0 then  nb_pmin_smooth= nx/3. else  nb_pmin_smooth=nb_pts_smooth_min
  ;nb_pts_smooth = round(n_elements(p)/3)

  if n_elements(nb_pts_smooth_min) eq 0 then  nb_pmin_smooth= 20. else  nb_pmin_smooth=nb_pts_smooth_min

  nb_pts_smooth = min([round(n_elements(p)/3), nb_pmin_smooth])
  vector_smooth = fltarr(nx)
  
  if affi then begin
    window, 3
    plot, vector_brut, col=0, back=-1, psym=2
  endif

  if affi then print,' nb_pts_smooth : ', nb_pts_smooth
  
  if p[0] ne -1 and  nb_pts_smooth gt nb_pmin_smooth/3 then begin
    if affi then oplot, p, vector_brut[p], col=50, psym=-1
      
    if p[n_elements(p)-1]-p[0]+1 eq n_elements(p) then begin

       if affi then print, 'points consecutifs' 
       ;; premier lissage le plus large possible
       vector_smooth[p[0]:p[n_elements(p)-1]-1] = smooth( vector_brut[p[0]:p[n_elements(p)-1]-1], nb_pts_smooth )
       
       if affi then oplot, vector_smooth, col=-50
            
       pp = where(vector_smooth ne 0)
       
       ;; reglage au niveau des extremites des points selectionnes
       ;; la limite correspond à la taille du precedent lissage /2
       
       ind_smooth_lim = ceil(n_elements(pp)/3)+1	       
       
       if ind_smooth_lim/2 gt 1 then begin        
          polyd = poly_fit( p[0]+findgen(ind_smooth_lim/2)+floor(ind_smooth_lim/8), $
                            vector_smooth[p[0]+floor(ind_smooth_lim/8):( p[0]+ind_smooth_lim/2+floor(ind_smooth_lim/8)-1)], 1 )
          vector_smooth[p[0]:(p[0]+ind_smooth_lim/2-1)] = polyd[1]*(findgen(ind_smooth_lim/2)+p[0])+polyd[0]
       endif else vector_smooth[ p[0]:(p[0]+ind_smooth_lim/2-1)] =0
      
       if affi then oplot, vector_smooth, col=-70
      
       pf = p[n_elements(p)-1]
       if pf ge ind_smooth_lim/2 and ind_smooth_lim/2 gt 1 then begin
          polyf = poly_fit( reverse(pf-findgen(ind_smooth_lim/2)-floor(ind_smooth_lim/8)), $
                            vector_smooth[pf-(ind_smooth_lim/2-1)-floor(ind_smooth_lim/8): pf-floor(ind_smooth_lim/8)], 1 )
          vector_smooth[pf-(ind_smooth_lim/2):pf ] = polyf[1]*(reverse(pf-findgen(ind_smooth_lim/2+1)))+polyf[0]                     
      endif
       
       if affi then oplot, vector_smooth, col=-80
       
       ;; lisasge prenant en compte les extrapolations sur les cotes
       ppp = where(vector_smooth ne 0)
       if n_elements(ppp) ge 16. then vector_smooth[ ppp] = smooth(vector_smooth[ ppp], round(n_elements(ppp)/8)-1)
       if affi then oplot, vector_smooth, col=-90
       ppp= -1                              
       pp = -1
       
       pi = p[0]
       pf= p[n_elements(p)-1]-1
       if  pf-pi gt nb_pmin_smooth then  vector_smooth[pi:pf] = smooth(vector_smooth[pi:pf], nb_pmin_smooth)
       
       smooth_indi = nb_pmin_smooth-2
       while  smooth_indi gt 4 do begin        
          vector_smooth[pi:pf] = smooth(vector_smooth[pi:pf],smooth_indi)
          smooth_indi = smooth_indi -2
       endwhile
       
       if affi then oplot, vector_smooth, col=-80      
       if affi then stop  
      
    endif else begin  ;;; les points selectionnes ne sont pas tous consecutifs
       
       interm = fltarr(nx)
       interm2 = fltarr(nx)
       interm[p] = vector_brut[p]   
       
       ;; remplissage si eventuellement un pas manquant (moyennes des 2 autres valeurs
       interm = fill_gap_one_point_inavector(interm)
       
      ;; afin d'eviter les petits bouts de données tous seuls dans les coins qui introduisent des erreurs lors du lissage suivant
       pp = where(abs(interm) gt 0.000001)
       q = abs(shift(pp,-1)-pp)
       p_disc = where(q gt 1)
       
       n_disc = n_elements(p_disc)
       
       val_q = where(q lt 4 and q gt 1)
       ;print, 'Val q : ', n_elements(val_q)
       if fill_gap eq 0 and n_disc lt 5 and n_elements(val_q) lt 5 then begin
          ;print, 'Use the fill_gap_several_points_inavector in lissage pour de petits lissages'
          interm = fill_gap_several_points_inavector(interm)
          interm2 = fill_gap_several_points_inavector(interm2)
          pp = where(abs(interm) gt 0.000001)
          q = abs(shift(pp,-1)-pp)
          p_disc = where(q gt 1)
          n_disc = n_elements(p_disc)
       endif
       
       if fill_gap eq 1 then begin
          ;print, 'Use the fill_gap_several_points_inavector in lissage pour d''eventuels longs lissages - fill_gap'         
          interm = fill_gap_several_points_inavector(interm)
          
          pa1 = p[0]
          pa2 = floor(nb_pmin_smooth/2)+1+p[0] 
          pb1 = p[n_elements(p)-1]-1 -floor(nb_pmin_smooth/2)
          pb2 = p[n_elements(p)-1]-1 
          
          interm2[pa1:pb2] = smooth( interm[pa1:pb2], nb_pmin_smooth )       
          interm2[pa1:pa2] = interm[pa1: pa2]
          interm2[pb1:pb2] = interm[pb1: pb2]
          dinterm2 = fltarr(nx)
          dinterm2[pa1:pb2] = deriv(interm2[pa1: pb2])         
          mean_ds = mean(dinterm2[pa1: pb2])
          std_ds = stddev(dinterm2[pa1: pb2])
          ;;print, 'meazn et STd  : ',mean_ds, std_ds
          ps = where(abs(dinterm2-mean_ds) gt std_ds or interm eq 0)
          if ps[0] ne -1 then interm2[ps] = 0.
          ps = -1
          
          diff = abs(interm2-interm)
          mean_diff = mean(diff[p])
          std_diff = stddev(diff[p])
          
         ps = where(abs(diff-mean_diff) gt std_diff and interm ne 0)
         if ps[0] ne -1 then interm2[ps] = interm[ps] 
         ps = -1
                  
         if affi then oplot, interm2, col=-80
         if affi then stop
         interm = fill_gap_several_points_inavector(interm2)
        
         pp = where(abs(interm) gt 0.000001)
         q = abs(shift(pp,-1)-pp)
         p_disc = where(q gt 1)
         n_disc = n_elements(p_disc)
      endif else print, 'No use of fill_gap n_disc :', n_disc
       
      pi = pp[0]
      for i_disc = 0L, n_disc-1 do begin
         
        pf = pp[p_disc[i_disc]]        
        n_tibou = pf-pi+1
        if affi then print, 'n tibou : ', n_tibou
        if affi then print, 'Bornes consideres : ', pi, pf
        if pi gt pf then stop
        
        if n_tibou gt nb_pmin_smooth then begin           
           ind_smooth_lim = min([round(n_tibou/6)*2+1, nb_pmin_smooth])
           
           vector_smooth[pi:pf] = smooth(interm[pi:pf], ind_smooth_lim)
           
           if affi then oplot, vector_smooth, col=-50
           
           if ind_smooth_lim/2 gt 1 then begin 
              polyd = poly_fit( pi+findgen(ind_smooth_lim/2)+floor(ind_smooth_lim/8), $
                                vector_smooth[pi+floor(ind_smooth_lim/8):( pi+ind_smooth_lim/2+floor(ind_smooth_lim/8)-1)], 1 )
              vector_smooth[ pi:(pi+ind_smooth_lim/2-1)] = polyd[1]*(findgen(ind_smooth_lim/2)+pi)+polyd[0]
           endif          
           
           if affi then oplot, vector_smooth, col=-60
           
           if ind_smooth_lim/2 gt 1 then begin
              polyf = poly_fit( reverse(pf-findgen(floor(ind_smooth_lim/2))-floor(ind_smooth_lim/8)), $
                                vector_smooth[pf-floor((ind_smooth_lim/2)-1)-floor(ind_smooth_lim/8): pf-floor(ind_smooth_lim/8)], 1 )           
              vector_smooth[pf-floor(ind_smooth_lim/2):pf ] = polyf[1]*(reverse(pf-findgen(floor(ind_smooth_lim/2)+1)))+polyf[0]
           endif
           
           if affi then oplot, vector_smooth, col=-70
           
           vector_smooth[pi:pf] = smooth(vector_smooth[pi:pf], nb_pmin_smooth)
           
           if affi then print, 'Nbre de pts de lissage deuxieme option :', ind_smooth_lim, nb_pmin_smooth
           
           smooth_indi = ind_smooth_lim-2
           while smooth_indi gt 4 do begin        
              vector_smooth[pi:pf] = smooth(vector_smooth[pi:pf],smooth_indi)
              smooth_indi = smooth_indi - 2        
           endwhile
           
           if affi then oplot, vector_smooth, col=-80         
           
        endif
        
        pi = pp[p_disc[i_disc]]+q[p_disc[i_disc]]
      
      endfor
      
    endelse
    
  endif
        
  if affi then oplot, vector_smooth, col=-100, psym=-1
  if affi then stop
  
  p = where(vector_smooth ne 0 and vector_brut ne 0)
  if p[0] ne -1 then begin
     diff_std = stddev(vector_smooth[p]-vector_brut[p])
     diff_mean = mean(vector_smooth[p]-vector_brut[p])
     pp = where(abs( (vector_smooth-vector_brut)-diff_mean) gt diff_std)
     if pp[0] ne -1 then vector_smooth[pp] = 0.
     pp = -1
  endif
  
  vector_smooth = fill_gap_one_point_inavector(vector_smooth)
  if fill_gap eq 1 then vector_smooth = fill_gap_several_points_inavector(vector_smooth)
  
  dsmooth = fltarr(nx)
  p = where(vector_smooth ne 0)
  if p[0] ne -1 then begin
     dsmooth[p] = deriv(vector_smooth[p])         
     mean_ds = mean(dsmooth[p])
     std_ds = stddev(dsmooth[p])
     ;print, 'meazn et STd  : ',mean_ds, std_ds
     ps = where(abs(dsmooth-mean_ds) gt 3*std_ds or vector_smooth eq 0)
     if ps[0] ne -1 then vector_smooth[ps] = 0.
     ps = -1          
  endif

  if affi then oplot, vector_smooth, col=240, thick=2, psym=-1
  if affi then stop

  return, vector_smooth
end
