
;; procedure read_dop_nc
;;    read the OCN-RVL file from SENTINEL-1A and B 
;;    unbias the Doppler anomaly
;;    INPUT : - but the data directory has to be specified at the
;;               first line. data are supposed to be in the OCN*.SAFE format 
;;    OUTPUT : fichiers tif, fichiers png
;;    Author :C. Danilo, juin 2018
;;
;;    1 - Plot the Wind information l.68
;;    2 - First step correction : calculation of azimuthal profile l.125
;;    3 - Second step correction correction : calculation of
;;        reference radial profile (subswath, ground, ocean, etc...) l.306
;;    4 - Last part of the correction - how to deal with the
;;        difference reference profile l.453
;;    5 - Last part : transformation from Doppler anomaly to radial
;;        velocity -  printing l.854


pro read_dop_nc
  
  dir = 'name_of_your_directory';; directory where are SAFE files
  date = '' ;; initialisation of date to an empty string
  
  dirocn = file_search(dir+'*OCN*'+date+'*.SAFE')
  
  nb_ocn = n_elements(dirocn)
  print, 'Nb of OCN files :', nb_ocn, dir
  
  wrange=[0,20]
  nrange=[-10,-1]
  rrange = [-10,10]
  dcrange= [-60.,60.]
  rvrange= [-3.,3.]
  
  flag_terre = 100
  flag_mer = 0
  def_val_rvl = -999.00
  num_dop_min = 10
  max_dc_std = 4.

  dop_max_tolerated = 80.
  err_ano = 15.
  
  xsize = 800
  pos_cb = [0.15, 0.90, 0.95, 0.95]
  
  suf ='_0'

  for i_ocn =0L,  nb_ocn-1 do begin
    info_safe = read_safe_s1(dirocn[i_ocn])
    status = -1
    
    if info_safe.check ne 0 then begin

      title=info_safe.product+' '+strmid(info_safe.starttime, 0,16)+' '+info_safe.swath
      
      time_str = strmid(info_safe.starttime, 0,4)+strmid(info_safe.starttime, 5,2)+$
                 strmid(info_safe.starttime, 8,2)+strmid(info_safe.starttime, 10,3) 

      status = -1

      read_netcdf, info_safe.file[0], data, attributes, status

      if status eq 0 then begin 
       
         srvl = size(data.rvldcobs)   

         ;; 1 - Plot the Wind information
         if total(data.rvldcgeo) ne 0. and srvl[2] gt 1 and srvl[3] gt 1 and srvl[0] gt 2 then begin
            minlon = min( data.owilon, /nan)
            minlat = min( data.owilat, /nan)
            maxlon = max( data.owilon, /nan)
            maxlat = max( data.owilat, /nan)
            
            pta = lonlat2utm(minlon, minlat)
            ptb = lonlat2utm(maxlon, maxlat)
            if pta.zone eq ptb.zone then s_ratio = (ptb.n-pta.n)/(ptb.e-pta.e)
            if pta.zone eq ptb.zone-1 then s_ratio = (ptb.n-pta.n)/(ptb.e+500000.-pta.e)
            if pta.zone eq ptb.zone-2 then s_ratio = (ptb.n-pta.n)/(ptb.e+500000.*2-pta.e)
            if pta.zone eq ptb.zone-3 then s_ratio = (ptb.n-pta.n)/(ptb.e+500000.*3-pta.e) 
            
            l_ratio = (maxlat-minlat)/(maxlon-minlon)
            ysize = xsize*l_ratio
            
            sca = 0.05
            pas = 20 
            A = FINDGEN(17) * (!PI*2/16.)
            USERSYM, COS(A)/2, SIN(A)/2, /FILL

            p = where(data.owilon ne 9999 and data.owiwindspeed gt 0. and data.owiwindspeed lt wrange[1] )
            if p[0] ne -1 then begin
               window, 0, xsize = xsize, ysize = ysize
               scatter, data.owilon[p], data.owilat[p], data.owiwindspeed[p], /cbar, zrange=wrange, back=-1, cltable =40, $
                        xtitle='Long. [deg.]', ytitle='Lat. [deg.]', psym=8, title=title, cbtitle='Wind speed [m/s]', symsize=1, thick= 1, $
                        xrange=[minlon, maxlon], yrange=[minlat, maxlat]
               arrow, data.owilon[p[0]:p[n_elements(p)-1]:pas], data.owilat[p[0]:p[n_elements(p)-1]:pas], $
                      data.owilon[p[0]:p[n_elements(p)-1]:pas]+cos(-data.owiwinddirection[p[0]:p[n_elements(p)-1]:pas]*!dtor-!pi/2)*sca, $
                      data.owilat[p[0]:p[n_elements(p)-1]:pas]+sin(-data.owiwinddirection[p[0]:p[n_elements(p)-1]:pas]*!dtor-!pi/2)*sca, $
                      /data, hsize = -0.3, col = 0, thick = 1                   
            endif
 
            p = where(data.owilon ne 9999 and data.owiwindspeed gt 0. )
            if p[0] ne -1 then begin 
               window, 1, xsize = xsize, ysize = ysize
               scatter, data.owilon[p], data.owilat[p], data.owiecmwfwindspeed[p], /cbar, zrange=wrange, back=-1, cltable =40, $
                        xtitle='Long. [deg.]', ytitle='Lat. [deg.]', psym=8, title=title, cbtitle='Wind speed ECMWF [m/s]', symsize=1, thick= 1,$
                        xrange=[minlon, maxlon], yrange=[minlat, maxlat]

               arrow, data.owilon[p[0]:p[n_elements(p)-1]:pas], data.owilat[p[0]:p[n_elements(p)-1]:pas], $
                      data.owilon[p[0]:p[n_elements(p)-1]:pas]+cos(-data.owiecmwfwinddirection[p[0]:p[n_elements(p)-1]:pas]*!dtor-!pi/2)*sca, $
                      data.owilat[p[0]:p[n_elements(p)-1]:pas]+sin(-data.owiecmwfwinddirection[p[0]:p[n_elements(p)-1]:pas]*!dtor-!pi/2)*sca, $
                /data, hsize = -0.3, col = 0, thick = 2
               
               pp = where(data.rvllandcoverage eq flag_terre)   
               if pp[0] ne -1 then oplot, data.rvllon[pp], data.rvllat[pp], col=50, psym=8, thick=2
               
               ;filename=dirocn[i_ocn]+'/OWI_'+time_str+'ecmwfwindspeed.png'
               ;saveimage, filename, PNG='PNG'
            endif    
    
      nrad = srvl[2]
      nazim = srvl[3]
      nbwin = srvl[1]

      ;; 2 - First correction : calculation of azimuthal profile
      print, ' "From read_dop_nc " -- Radial correction -- '
      ;; Variable initialisation
      ;; 'terre' means ground
      ;; 'mer' means sea
      doppler_terre = fltarr(nbwin, nazim)
      num_doppler_terre = fltarr(nbwin, nazim)
      std_doppler_terre= fltarr(nbwin, nazim)
     
      doppler_mer = fltarr(nbwin, nazim)
      num_doppler_mer = fltarr(nbwin, nazim)
      std_doppler_mer= fltarr(nbwin, nazim)
      
      mean_ncrs_mer = fltarr(nbwin, nazim)
      median_ncrs_mer= fltarr(nbwin, nazim)
      std_ncrs_mer = fltarr(nbwin, nazim)
      
      rvl_azim_terre = fltarr(nbwin, nrad)
      rvl_azim_terre_std = fltarr(nbwin, nrad)
      rvl_azim_terre_num = fltarr(nbwin, nrad)
      
      rvl_azim_mer = fltarr(nbwin, nrad)
      rvl_azim_mer_std = fltarr(nbwin, nrad)
      rvl_azim_mer_num = fltarr(nbwin, nrad)
      
      nrcs_azim_mer = fltarr(nbwin, nrad)
      nrcs_azim_mer_std = fltarr(nbwin, nrad)
      nrcs_azim_mer_num = fltarr(nbwin, nrad)
      
      rvl_rad = fltarr(nbwin, nazim)
      rvl_rad_std = fltarr(nbwin, nazim)
            
      for iim = 0L, nbwin-1 do begin
      
         dc_to_treat = data.rvldcobs-data.rvldcmiss-data.rvldcgeo

         ;; Doppler azimuthal correction on ground
         res = azim_profil( reform(dc_to_treat[iim, *, *]), reform(data.rvllandcoverage[iim, *, *]), flag_terre, def_val_rvl, reform(data.rvldcobsstd[iim, *, *]), max_dc_std)
         rvl_azim_terre[iim, *] = res.azim_prof
         rvl_azim_terre_num[iim, *] = res.azim_prof_num
         rvl_azim_terre_std[iim, *] = res.azim_prof_std
         res = 0
        
         ;; Idem on sea
         res = azim_profil( reform(dc_to_treat[iim, *, *]), reform(data.rvllandcoverage[iim, *, *]), flag_mer, def_val_rvl, reform(data.rvldcobsstd[iim, *, *]), max_dc_std)
         rvl_azim_mer[iim, *] = res.azim_prof
         rvl_azim_mer_num[iim, *] = res.azim_prof_num
         rvl_azim_mer_std[iim, *] = res.azim_prof_std
         res = 0

         ;; NRCS azimuthal correction on sea
         res = azim_profil( reform(data.rvlnrcs[iim, *, *]), reform(data.rvllandcoverage[iim, *, *]), flag_mer, def_val_rvl, reform(data.rvldcobsstd[iim, *, *]), max_dc_std)
         nrcs_azim_mer[iim, *] = res.azim_prof
         nrcs_azim_mer_num[iim, *] = res.azim_prof_num
         nrcs_azim_mer_std[iim, *] = res.azim_prof_std
         res = 0
        
         p = where(abs(rvl_azim_terre[iim,*]) lt 0.00001)
         if p[0] ne -1 then rvl_azim_terre[iim,p] = 0.
         p = where(abs(rvl_azim_mer[iim,*]) lt 0.00001)
         if p[0] ne -1 then rvl_azim_mer[iim,p] = 0.

         ;; calculation of some stats
         for kim = 0L, nazim-1 do begin
            p = where(data.rvldcobs[iim, *, kim] ne def_val_rvl and data.rvldcmiss[iim, *, kim] ne def_val_rvl and $
                      data.rvlnrcs[iim, *, kim] ne def_val_rvl and abs(data.rvldcobsstd[iim, *, kim]) lt max_dc_std )
            if p[0] ne -1 then begin
               rvl_rad[iim, kim] = mean(dc_to_treat[iim, p, kim])
               rvl_rad_std[iim, kim] = stddev(dc_to_treat[iim, p, kim])
            endif
         endfor
        
      endfor
      
      p = where(rvl_azim_terre ne 0)
      if p[0] ne -1 then begin 
	mean_rvl_azim_terre = mean(rvl_azim_terre[p])
        std_rvl_azim_terre = stddev(rvl_azim_terre[p])
      endif else begin
        mean_rvl_azim_terre = 0.
        std_rvl_azim_terre = 0
      endelse

      ;; looking for biased border points 
      rad_bon=fltarr(nbwin, nrad)+1
      for iim = 0L, nbwin-1 do $
         for jim=0L, nrad-1 do begin
            nbc = 0
            nb0 = 0
            for kim =0L, nazim-1 do begin
              if abs(dc_to_treat[iim, jim, kim]-rvl_rad[iim, kim]) lt 3*rvl_rad_std[iim, kim] then nb0 = nb0+1
              if abs(dc_to_treat[iim, jim, kim]-rvl_rad[iim, kim]) gt 3*rvl_rad_std[iim, kim] or data.rvldcobs[iim, jim, kim] eq 0 then nbc = nbc+1
            endfor
            if nbc gt 10 then rad_bon[iim, jim] = 0
         endfor
           
      ;; looking for a smoother solution in the radial direction correction
      rvl_azim_terre_smooth = fltarr(nbwin, nrad)
      rvl_azim_mer_smooth = fltarr(nbwin, nrad)
      nrcs_azim_mer_smooth = fltarr(nbwin, nrad)
      liss = 10
      
      for iim = 0L, nbwin-1 do begin
      
        rvl_azim_terre_smooth[iim, *] = lissage_profil(reform(rvl_azim_terre[iim,*]),reform(rvl_azim_terre_num[iim,*]), num_dop_min ,reform(rad_bon[iim,*]),/fill_gap, nb_pts_smooth_min= liss)
        rvl_azim_mer_smooth[iim, *]   = lissage_profil(reform(rvl_azim_mer[iim,*]),reform(rvl_azim_mer_num[iim,*]), num_dop_min,reform(rad_bon[iim,*]), /fill_gap, nb_pts_smooth_min= liss)
        nrcs_azim_mer_smooth[iim, *]  = lissage_profil(reform(nrcs_azim_mer[iim,*]),reform(nrcs_azim_mer_num[iim,*]), num_dop_min ,reform(rad_bon[iim,*]), /fill_gap, nb_pts_smooth_min= liss)
        
        ;; looking for a profil on the entire width of the sub-swath
        largeur_totale = n_elements(where(rad_bon[iim,*] eq 1))
        largeur_terre = n_elements(where(rvl_azim_terre_smooth[iim, *] ne 0))
        largeur_mer = n_elements(where(rvl_azim_mer_smooth[iim, *] ne 0))

        print, ' "From read_dop_nc " -- swath width in pixels (total, on seas , on ground ) ', largeur_totale, largeur_terre, largeur_mer
        
        if largeur_totale le largeur_terre then print, 'Radial correction optimized for the sub-swath'
        
        if largeur_totale gt largeur_terre then begin
          print, 'Radial correction non optimized for the sub-swath'
          
          diff_rvl_terre_mer =  fltarr(nrad)
          
          p = where(rvl_azim_terre_smooth[iim, *] ne 0 and rvl_azim_mer_smooth[iim, *] ne 0 and rad_bon[iim,*] eq 1 and nrcs_azim_mer_smooth[iim, *] ne 0.)
          
          if p[0] ne -1 and n_elements(p) gt 2 then begin
            diff_rvl_terre_mer[p] = rvl_azim_terre_smooth[iim, p]-rvl_azim_mer_smooth[iim, p]
            
            p_mer = where(rvl_azim_mer_smooth[iim, *] ne 0 and rad_bon[iim,*] eq 1 and nrcs_azim_mer_smooth[iim, *] ne 0.)
            
            if stddev(diff_rvl_terre_mer[p]) lt 5. and stddev(alog(nrcs_azim_mer_smooth[iim, p])) lt 1. and $
              stddev(alog(nrcs_azim_mer_smooth[iim, p_mer])) lt 1.  then begin
              
              p_mer_plus = where(rvl_azim_terre_smooth[iim, *] eq 0 and rvl_azim_mer_smooth[iim, *] ne 0 and rad_bon[iim,*] eq 1 and nrcs_azim_mer_smooth[iim, *] ne 0.)
              if p_mer_plus[0] ne -1 then begin

                 loadct, 39
                 device, decomposed = 0
                 window, 2
                 plot, rvl_azim_terre_smooth[iim, *], col=-1, back=0, xtitle='Azim', ytitle='RVL azim terre et terre elargi avec la composante maritime'
                
                 print, 'Use of maritime component !!'
                 
                 rvl_azim_terre_smooth[iim, p_mer_plus] = rvl_azim_mer_smooth[iim, p_mer_plus]+mean(diff_rvl_terre_mer[p])
                 oplot, rvl_azim_terre_smooth[iim, *], col=-1 
                 rvl_azim_terre_smooth[iim,*] = lissage_profil(reform(rvl_azim_terre_smooth[iim,*]),reform(rvl_azim_mer_num[iim,*]+rvl_azim_terre_num[iim,*]),num_dop_min ,reform(rad_bon[iim,*]), /fill_gap)                             

                oplot, rvl_azim_terre_smooth[iim, *], psym=-2, col=50
	              
              endif else print, 'No commune points between  azim_terre nul and azim_mer no nul.'
              p_mer_plus = -1
              p_mer = -1
            endif else begin
              print, 'Strong variations between land and sea (STD diff, STD NRCS (composante commune), STD NRCS (composante sea))', $
                stddev(diff_rvl_terre_mer[p]), stddev(alog(nrcs_azim_mer_smooth[iim, p])),stddev(alog(nrcs_azim_mer_smooth[iim, p_mer]))
                
              window, 2
              plot, rvl_azim_terre_smooth[iim, *], psym=-1, col=-1, back=0
              oplot, rvl_azim_mer_smooth[iim, *], col=50
              
            endelse            
            p = -1
          endif else print, 'No commune values between land and sea profiles'
        
        endif
        
      endfor
      
      loadct, 39
      device, decomposed = 0
      
      window, 3
      plot,  findgen(nrad*3), rvl_azim_terre[*,*], yrange=[min([mean_rvl_azim_terre-3*std_rvl_azim_terre,min(rvl_azim_terre)]),max([mean_rvl_azim_terre+3*std_rvl_azim_terre,max(rvl_azim_terre)]) ], $
                col=0, back =-1, title=title, ytitle='Mean Doppler on the ground [Hz]', /nodata, /xstyle
     
      for iim = 0L, nbwin-1 do oplot, findgen(nrad)+iim*nrad,rvl_azim_terre[iim,*], col=50+50*iim
      for iim = 0L, nbwin-1 do oplot, findgen(nrad)+iim*nrad, rvl_azim_terre_smooth[iim,*], col=50+50*iim, thick=2
      
      filename=dirocn[i_ocn]+'/Correc_radial_'+time_str+suf+'.png'
      saveimage, filename, PNG='PNG'       

      ;; 3 - Correction in the azimuthal correctio - calculation of
      ;;     radial profiles
      
      radvelec = fltarr(nbwin, nrad, nazim)
      
      ;; if enough points on the ground : profiles on ground is used as a reference
      p = where(rvl_azim_terre_smooth ne 0)
      ref_terre = 0

      if n_elements(p) gt 3. then begin
      ;; if n_elements(p) gt nrad*nbwin/3 then begin
        ref_terre = 1
        for iim = 0L, nbwin-1 do $
           for kim = 0L, nazim-1 do begin
           p = where(rvl_azim_terre_smooth[iim,*] ne 0 and $
                     rad_bon[iim, *] eq 1 and $
                     data.rvldcobs[iim, *, kim] ne def_val_rvl and $
                     data.rvldcmiss[iim, *, kim] ne def_val_rvl and $
                     abs(data.rvldcobsstd[iim, *, kim]) lt max_dc_std)
           if p[0] ne -1 then radvelec[iim, p, kim] = dc_to_treat[iim, p, kim]-rvl_azim_terre_smooth[iim,p]
           p=-1
        endfor
        
        ;; Check 
        rvlec_azim_terre = fltarr(nbwin, nrad)
        rvlec_azim_terre_smo = fltarr(nbwin, nrad)  
        rvlec_azim_terre_petit_lissage = fltarr(nbwin*nrad)
        rvlec_azim_terre_gd_lissage = fltarr(nbwin*nrad)   
        
        for iim = 0L, nbwin-1 do begin  
           res = azim_profil(reform(radvelec[iim, *, *]),reform(data.rvllandcoverage[iim, *, *]), flag_terre, def_val_rvl, reform(data.rvldcobsstd[iim, *, *]), max_dc_std)
           rvlec_azim_terre[iim, *] = res.azim_prof
           res = 0
           rvlec_azim_terre_smo[iim, *] = lissage_profil(reform(rvlec_azim_terre[iim,*]),reform(rvl_azim_terre_num[iim,*]), num_dop_min ,reform(rad_bon[iim,*]), /fill_gap, nb_pts_smooth_min= liss)
           rvlec_azim_terre_petit_lissage[iim*nrad:(iim+1)*nrad-1] = rvlec_azim_terre_smo[iim, *]            
        endfor    
        
        rvlec_azim_terre_gd_lissage = lissage_profil(rvlec_azim_terre_petit_lissage, rvlec_azim_terre_petit_lissage*0+num_dop_min, num_dop_min , rvlec_azim_terre_petit_lissage*0+1, /fill_gap, nb_pts_smooth_min= liss)       
          
      endif else begin  
      
         print, ' No ground reference - Use the maritime component '        
         ref_terre = 0
         for iim = 0L, nbwin-1 do $
            for kim = 0L, nazim-1 do begin
            p = where(rvl_azim_mer_smooth[iim,*] ne 0 and rad_bon[iim, *] eq 1)
            if p[0] ne -1 then radvelec[iim, p, kim] = dc_to_treat[iim, p, kim]-rvl_azim_mer_smooth[iim,p]                
            p=-1
            
         endfor
         title = title+' -- no ground reference --'
      endelse
      p = -1      
          
      ;; Into the azimuthal direction
      for iim=0L, nbwin-1 do begin

        res = rad_profil_std(reform(radvelec[iim, *,*]), reform(data.rvllandcoverage[iim, *, *]), flag_terre, 0., $
                             reform(data.rvldcobsstd[iim, *,*]), max_std =  max_dc_std/2. )

        doppler_terre[iim, *]= res.rad_prof
        num_doppler_terre[iim, *] = res.rad_prof_num
        std_doppler_terre[iim, *] = res.rad_prof_std
        res =0
        
        res = rad_profil_std(reform(radvelec[iim, *,*]), reform(data.rvllandcoverage[iim, *, *]), flag_mer, 0., $
                             reform(data.rvldcobsstd[iim, *,*]), max_std =  max_dc_std/2. )

        doppler_mer[iim, *]= res.rad_prof
        num_doppler_mer[iim, *] = res.rad_prof_num
        std_doppler_mer[iim, *] = res.rad_prof_std
        res = 0
      endfor
      
      deriv_num_dop_terre = fltarr(nbwin, nazim)
      deriv_num_dop_mer = fltarr(nbwin, nazim)
      for  iiw = 0L, nbwin-1 do begin
        deriv_num_dop_terre[iiw,*] = deriv(num_doppler_terre[iiw,*])
        deriv_num_dop_mer[iiw,*] = deriv(num_doppler_mer[iiw,*])
      endfor    
     
      doppler_terre_treated = fltarr(nbwin, nazim)
      doppler_mer_treated = fltarr(nbwin, nazim)
      doppler_terre_treated_verif = fltarr(nbwin, nazim)
      doppler_mer_treated_verif= fltarr(nbwin, nazim)
      trend_dop_mer = fltarr(nbwin, nazim)
      
      for iiw = 0L, nbwin-1 do begin
      
        ;;Terre
        p = where(num_doppler_terre[iiw,*] ge num_dop_min and doppler_terre[iiw,*] ne 0 )
        if p[0] ne -1 then doppler_terre_treated[iiw,p] = doppler_terre[iiw,p]
        p= -1
        
        doppler_terre_treated[iiw,*] = fill_gap_one_point_inavector(doppler_terre_treated[iiw,*]) 
        
        ;; Verification
        p = where(doppler_terre_treated[iiw,*] ne 0)
        if p[0] ne -1 then for kim = 0L, n_elements(p)-1 do begin
          pt = where(reform(data.rvllandcoverage[iiw, *, p[kim]]) eq flag_terre)
          if pt[0] ne -1 then doppler_terre_treated_verif[iiw, p[kim]] = mean(radvelec[iiw, pt,p[kim]]-doppler_terre_treated[iiw,p[kim]])
          pt =-1
        endfor
        p = -1
        
        p =where( abs(doppler_terre_treated_verif[iiw, *]) gt err_ano)
        if p[0] ne -1 then doppler_terre_treated[iiw,p] = fltarr(n_elements(p))
        p=-1
        
        doppler_terre_treated[iiw,*] = fill_gap_one_point_inavector(reform(doppler_terre_treated[iiw,*]))
        
        ;;Mer
        p = where(num_doppler_mer[iiw,*] ge num_dop_min  and doppler_mer[iiw,*] ne 0) 
        if p[0] ne -1 then doppler_mer_treated[iiw,p] = doppler_mer[iiw,p]
        p= -1
        
        doppler_mer_treated[iiw,*] = fill_gap_one_point_inavector((reform(doppler_mer_treated[iiw,*])))
        
        vect_flag = fltarr(nazim)
        p = where(doppler_mer_treated[iiw,*] ne 0 and  abs(doppler_mer_treated[iiw,*]) lt dop_max_tolerated)
        if p[0] ne -1 then vect_flag[p] =1
        p=-1     
        
	; ne pas jeter bebe avec l''eau du bain
        trend_dop_mer[iiw,*] = lissage_profil(reform(doppler_mer_treated[iiw,*]), reform(num_doppler_mer[iiw,*]), num_dop_min, vect_flag, nb_pts_smooth_min= 45, /fill_gap)
        
        p_t = where(trend_dop_mer[iiw,*] ne 0 and abs(doppler_mer_treated[iiw,* ]) gt 0.00001 )
        p_z = where(trend_dop_mer[iiw,*] eq 0 or  abs(doppler_mer_treated[iiw,* ]) lt 0.00001 )
        if p_t[0] ne -1 then doppler_mer_treated[iiw,p_t ] = doppler_mer_treated[iiw,p_t ]-trend_dop_mer[iiw,p_t]
        if p_z[0] ne -1 then doppler_mer_treated[iiw,p_z ] = fltarr(n_elements(p_z))
        p_t = -1
        p_z = -1          
             
        doppler_mer_treated[iiw,*] = fill_gap_one_point_inavector(reform(doppler_mer_treated[iiw,*]))  
      endfor
     
      loadct, 39
      device, decomposed = 0
      
      window, 4, xsize = xsize
      plot,  doppler_terre_treated[0,*], psym=-5, color=0, /nodata, back=-1,yrange=rrange*2., title=title, xrange=[0, nazim-1], /ystyle
      for iim=0L, nbwin-1 do  oplot, findgen(nazim), doppler_terre_treated[iim,*], color=50+50*iim, psym=-1;
      for iim=0L, nbwin-1 do  oplot, findgen(nazim), doppler_mer_treated[iim,*], color=50+50*iim, psym=-4;     
      for iim=0L, nbwin-1 do  oplot, findgen(nazim), trend_dop_mer[iim,*], color=50+50*iim, linestyle=4;     
      pointille, [0,nazim-1], rrange*2, col=0

      ;; 4 - Last part of the correction - how to deal with the
      ;;     difference reference profile
      rvl_corrected_2 = fltarr(nbwin, nrad, nazim)
      rvl_corrected_1 = fltarr(nbwin, nrad, nazim)
      doppler_win_test = fltarr(nbwin, nazim)
      doppler_win_test_terre = fltarr(nbwin, nazim)
      doppler_win_test_mer = fltarr(nbwin, nazim)
      
      for iim = 0L, nbwin-1 do begin
        ;; rvl_corrected_1
        p = where(abs(doppler_mer_treated[iim, *]) gt 0.00001 )
        if p[0] ne -1 then for jim=0L, nrad-1 do begin 
          if rad_bon[iim, jim] eq 0  then rvl_corrected_1[iim,jim,p] = fltarr(n_elements(p)) else rvl_corrected_1[iim,jim,p] = radvelec[iim, jim, p]-doppler_mer_treated[iim, p]
         
          for kim = 0L, nazim-1 do begin
            if (radvelec[iim, jim, kim] eq 0 or data.rvlnrcs[iim,jim,kim] eq def_val_rvl or data.rvldcobs[iim, jim, kim] eq def_val_rvl) then rvl_corrected_1[iim,jim,kim] = 0
          endfor
        endfor
        p=-1
        
        p = where(abs(doppler_terre_treated[iim, *]) gt 0.00001 )
        if p[0] ne -1 then for jim=0L, nrad-1 do begin
           if rad_bon[iim, jim] eq 0  then rvl_corrected_1[iim,jim,p] = fltarr(n_elements(p)) $
           else rvl_corrected_1[iim,jim,p] = radvelec[iim, jim, p]-doppler_terre_treated[iim, p]
          for kim = 0L, nazim-1 do begin
            if (radvelec[iim, jim, kim] eq 0 or data.rvlnrcs[iim,jim,kim] eq def_val_rvl or data.rvldcobs[iim, jim, kim] eq def_val_rvl) then rvl_corrected_1[iim,jim,kim] = 0
          endfor
        endfor
        p=-1

        ;; rvl_corrected_2 lissage dans la dir.azimutale
        for jim=0L, nrad-1 do begin
           p = where(rvl_corrected_1[iim, jim, *] ne 0 and radvelec[iim, jim, *] ne 0 and data.rvlnrcs[iim, jim, *] ne def_val_rvl and $        
                     abs(data.rvldcobsstd[iim,  jim, *]) lt max_dc_std and $
                     data.rvllandcoverage[iim,  jim, *] eq flag_mer)
           if n_elements(p) gt 10 then rvl_corrected_2[iim,jim,p] = smooth( rvl_corrected_1[iim,jim, p], 10 )
           p=-1

           p = where(rvl_corrected_1[iim, jim, *] ne 0 and radvelec[iim, jim, *] ne 0 and data.rvlnrcs[iim, jim, *] ne def_val_rvl and $ 
                     abs(data.rvldcobsstd[iim,  jim, *]) lt max_dc_std and $
                     data.rvllandcoverage[iim,  jim, *] eq flag_terre)
           if n_elements(p) gt 10 then rvl_corrected_2[iim,jim,p] = smooth( rvl_corrected_1[iim,jim, p], 10 )
           p=-1
        endfor
   
        ;; doppler_win_test_mer
        for kim=0L, nazim-1 do begin
          p  = where(radvelec[iim,*, kim] ne 0 and rvl_corrected_2[iim, *, kim] ne 0 and $
                     data.rvlnrcs[iim,*, kim] ne def_val_rvl and $                   
                     abs(data.rvldcobsstd[iim,  *, kim]) lt max_dc_std and $
                     data.rvllandcoverage[iim, *, kim] eq flag_mer)
          
          if p[0] ne -1 then doppler_win_test_mer[iim, kim] = median(radvelec[iim,p, kim]-rvl_corrected_2[iim,p, kim])
          p=-1       
                
          ;; doppler_win_test_terre     
          p  = where( radvelec[iim,*, kim] ne 0 and rvl_corrected_2[iim, *, kim] ne 0 and  $ 
                      data.rvlnrcs[iim,*, kim] ne def_val_rvl and $                     
                      abs(data.rvldcobsstd[iim,  *, kim]) lt max_dc_std and $
                      data.rvllandcoverage[iim,  *, kim] eq flag_terre )            
          if p[0] ne -1 then doppler_win_test_terre[iim, kim] = median(radvelec[iim,p, kim]-rvl_corrected_2[iim,p, kim])
          p=-1
          
        endfor 

        doppler_win_test_terre[iim,*] = fill_gap_one_point_inavector(reform(doppler_win_test_terre[iim,*]))     
        doppler_win_test_mer[iim,*] = fill_gap_one_point_inavector(reform(doppler_win_test_mer[iim,*]))     

        vect_flagm = fltarr(nazim)
        p = where(doppler_win_test_mer[iim,*] ne 0 and  abs(doppler_win_test_mer[iim,*]) lt dop_max_tolerated)
        if p[0] ne -1 then vect_flagm[p] =1
        p=-1
                  
        trend_mer_0 = lissage_profil(reform(doppler_win_test_mer[iim,*]), reform(num_doppler_mer[iim,*]), num_dop_min, vect_flagm, nb_pts_smooth_min= 45, /fill_gap)
                   
        p_t = where(trend_mer_0 ne 0 and abs(doppler_win_test_mer[iim,* ]) gt 0.00001 )
        p_z = where(trend_mer_0 eq 0 or  abs(doppler_win_test_mer[iim,* ]) lt 0.00001 )

        if p_t[0] ne -1 then doppler_win_test_mer[iim,p_t ] = doppler_win_test_mer[iim,p_t ]-trend_mer_0[p_t]
        if p_z[0] ne -1 then doppler_win_test_mer[iim,p_z ] = fltarr(n_elements(p_z)) 
               
        p_t = -1
        p_z = -1
        
        trend_mer_0 = 0
        doppler_win_test_mer[iim,*] = fill_gap_one_point_inavector(reform(doppler_win_test_mer[iim,*])) 
      endfor
    
      doppler_unif_treated = fltarr(nazim)      
      doppler_unif_terre = fltarr(nazim)
      doppler_unif_mer = fltarr(nazim)
      doppler_unif_tmw = fltarr(nbwin, nazim)
      
      doppler_terre_et_mer = [doppler_win_test_terre, doppler_win_test_mer]
      num_doppler_terre_et_mer = [num_doppler_terre, num_doppler_mer]
      std_doppler_terre_et_mer = [std_doppler_terre, std_doppler_mer]
      deriv_num_terre_et_mer = [deriv_num_dop_terre, deriv_num_dop_mer]
      
      ;; Looking for a signal period
      periode_possible = findgen(10)+16
      periode_terre_et_mer = fltarr(2*nbwin)
      max_correl = fltarr(2*nbwin)
      for iiw = 0L, 2*nbwin-1 do begin
        p = where((doppler_terre_et_mer[iiw,*]) ne 0 )
        if p[0] ne -1 and n_elements(p) gt max(periode_possible) then $
           correlate_val = a_correlate(reform(doppler_terre_et_mer[iiw,*]), periode_possible) else correlate_val = 0
        p = -1
        max_correl[iiw] = max(correlate_val)
        pper  =where( correlate_val eq max_correl[iiw])
        periode_terre_et_mer[iiw] = periode_possible[pper]
        pper = -1
      endfor      
      
      pper = where(max_correl gt 0.40 )
      if pper[0] ne -1 then periode_signal = round(total(max_correl[pper]*periode_terre_et_mer[pper])/total(max_correl[pper])) $
      else periode_signal = 0
      if stddev(periode_terre_et_mer) lt 1. and periode_signal eq 0 then periode_signal = round(mean(periode_terre_et_mer))
     
      doppler_unif_treated = unification_profil(doppler_terre_et_mer, std_doppler_terre_et_mer, num_doppler_terre_et_mer, deriv_num_terre_et_mer, num_dop_min )
      doppler_unif_terre = unification_profil(doppler_win_test_terre, std_doppler_terre, num_doppler_terre, deriv_num_dop_terre, num_dop_min )
      doppler_unif_mer = unification_profil(doppler_win_test_mer, std_doppler_mer, num_doppler_mer, deriv_num_dop_mer, num_dop_min )
      
      doppler_unif_treated = fill_gap_one_point_inavector(doppler_unif_treated)
      if periode_signal eq 0 then begin
        p = where((doppler_unif_treated) ne 0 )
        if p[0] ne -1 and n_elements(p) gt max(periode_possible) then correlate_val = a_correlate((doppler_unif_treated), periode_possible) else correlate_val = 0
        p = -1        
        max_cor = max(correlate_val,pper)
        if max_cor gt 0.1 then periode_signal = periode_possible[pper]
        pper = -1
      
      endif
            
      if periode_signal ne 21 then  print, '"From read_dop_nc " -- Signal Period : ', periode_signal

      doppler_unif_tw=fltarr(nbwin, nazim)
      doppler_unif_mw=fltarr(nbwin, nazim)
      doppler_unif_w = fltarr(nbwin, nazim)  
      
      if periode_signal gt 0 then begin
              
         for iiw = 0L, nbwin-1 do begin         
            doppler_unif_mw[iiw,*] = improv_prof_via_deriv( reform(doppler_win_test_mer[iiw,*]), periode_signal)
            doppler_unif_mw[iiw,*] = fill_gap_one_point_inavector(doppler_unif_mw[iiw,*])
            doppler_unif_tw[iiw,*] = improv_prof_via_deriv( reform(doppler_win_test_terre[iiw,*]), periode_signal)                            
            doppler_unif_tw[iiw,*] = fill_gap_one_point_inavector(doppler_unif_tw[iiw,*])   
         endfor
         
      endif else begin
         print, ' "From read_dop_nc " -- no evident signal period'
         doppler_unif_mw = doppler_win_test_mer
         doppler_unif_tw = doppler_win_test_terre        
      endelse                 
     
     for iim=0L, nbwin-1 do begin
        p = where(doppler_unif_tw[iim,*] ne 0)
        if n_elements(p) gt 40 then begin
          mmm = stddev(doppler_unif_tw[iim,p])
          pp = where( abs(doppler_unif_tw[iim,*]) lt 3*mmm)
          if n_elements(pp) gt 2 then res_trend_terre = poly_fit(pp,(doppler_unif_tw[iim,pp]), 1) else res_trend_terre = [0.,0.]
          if pp[0] ne -1 then oplot, pp, doppler_unif_tw[iim,pp], psym=2, thick=3 , col=50*iim+50  
          
        endif else res_trend_terre=[0.,0.]
        p = -1

        p = where(doppler_unif_mw[iim,*] ne 0)
        if p[0] ne -1 then doppler_unif_mw[iim,p] = doppler_unif_mw[iim,p]+res_trend_terre[1]*p+res_trend_terre[0]                         
        p = -1

        ;; raccord Mer et Terre 
        p_z = where(doppler_unif_mw[iim,*] eq 0)
        p = where(doppler_unif_mw[iim,*] ne 0 and doppler_unif_tw[iim,*] ne 0)
        if n_elements(p) gt 10 then begin
           mean_diff_terre_mer = mean(doppler_unif_mw[iim,p]-doppler_unif_tw[iim,p])
           doppler_unif_mw[iim,*] = doppler_unif_mw[iim,*]-mean_diff_terre_mer
           print, '"From read_dop_nc " -- RACCORD - Difference terre et mer [Hz] ', mean_diff_terre_mer
        endif
        p=-1
        if p_z[0] ne -1 then doppler_unif_mw[iim,p_z] = 0 
     
     endfor
  
     doppler_tm_mw = [doppler_unif_tw, doppler_unif_mw]
     doppler_unif_w = unification_profil(doppler_tm_mw, std_doppler_terre_et_mer*0+1, num_doppler_terre_et_mer*0+num_dop_min+1, deriv_num_terre_et_mer*0+4, num_dop_min )     
     doppler_unif = fill_gap_one_point_inavector(doppler_unif_w)
     
     doppler_unif =  improv_prof_via_deriv( reform(doppler_unif), periode_signal)      
     rvl_corrected0 = fltarr(nbwin, nrad, nazim)
     rvl_corrected1 = fltarr(nbwin, nrad, nazim)
     rvl_corrected2 = fltarr(nbwin, nrad, nazim)            
     rvl_corrected3 = fltarr(nbwin, nrad, nazim)
     
     for iim=0L, nbwin-1 do begin      
        ;; rvl_corrected0
        p = where(abs(doppler_unif_mw[iim, *]) gt 0.00001 )
        if p[0] ne -1 then for jim=0L, nrad-1 do begin
        
           if rad_bon[iim, jim] eq 0  then rvl_corrected0[iim,jim,p] = fltarr(n_elements(p)) $
           else rvl_corrected0[iim,jim,p] = radvelec[iim, jim, p]-doppler_unif_mw[iim, p]
        
           for kim = 0L, nazim-1 do begin
              if (radvelec[iim, jim, kim] eq 0 or data.rvlnrcs[iim,jim,kim] eq def_val_rvl or data.rvldcobs[iim, jim, kim] eq def_val_rvl) then rvl_corrected0[iim,jim,kim] = 0
           endfor
        endfor
        p=-1
        
        p = where(abs(doppler_unif_tw[iim, *]) gt 0.00001 )
        if p[0] ne -1 then for jim=0L, nrad-1 do begin
           if rad_bon[iim, jim] eq 0  then rvl_corrected0[iim,jim,p] = fltarr(n_elements(p)) else rvl_corrected0[iim,jim,p] = radvelec[iim, jim, p]-doppler_unif_tw[iim, p]
           for kim = 0L, nazim-1 do begin
              if (radvelec[iim, jim, kim] eq 0 or data.rvlnrcs[iim,jim,kim] eq def_val_rvl or data.rvldcobs[iim, jim, kim] eq def_val_rvl) then rvl_corrected0[iim,jim,kim] = 0
           endfor
        endfor
        p=-1
        
        ;; rvl_corrected1
        p = where(doppler_unif[*] ne 0 )
        if p[0] ne -1 then for jim=0L, nrad-1 do begin
           if rad_bon[iim, jim] eq 0 then rvl_corrected1[iim,jim,p] = fltarr(n_elements(p)) else rvl_corrected1[iim,jim,p] = radvelec[iim, jim, p]-doppler_unif[p]
                    
           for kim = 0L, nazim-1 do begin
              if (radvelec[iim, jim, kim] eq 0 or data.rvlnrcs[iim,jim,kim] eq def_val_rvl or data.rvldcobs[iim, jim, kim] eq def_val_rvl) then rvl_corrected1[iim,jim,kim] = 0          
           endfor
        endfor
        p = -1
     endfor
      
     ;; rvl_corrected2
     rvl_corrected2 = rvl_corrected0
     for iim=0L, nbwin-1 do begin
        for jim=0L, nrad-1 do begin
           p = where( abs(rvl_corrected2[iim,jim,*]) lt 0.0001 and abs(rvl_corrected1[iim,jim,*]) gt 0.0001)
           if p[0] ne -1  then  rvl_corrected2[iim,jim,p] = rvl_corrected1[iim,jim,p]
           p = -1
           for kim=0L, nazim-1 do begin
              if (data.rvlnrcs[iim,jim,kim] eq def_val_rvl or $
                  radvelec[iim, jim, kim] eq 0 or $
                  data.rvldcobs[iim, jim, kim] eq def_val_rvl $
                  or rad_bon[iim, jim] eq 0 and rad_bon [iim, jim] eq 0) $
              then rvl_corrected2[iim,jim,kim] = 0
           endfor
        endfor        
     endfor
           
     rvl_corrected3 = rvl_corrected2
     doppler_win_verifie = fltarr(nbwin, nazim)
     doppler_win_ver_terre = fltarr(nbwin, nazim)
     doppler_win_ver_mer = fltarr(nbwin, nazim)
     trend_dop_mer_fin = fltarr(nbwin, nazim)    
     doppler_win_ver_tm = fltarr(nbwin, nazim)
     
     for iim=0L, nbwin-1 do begin
        
        for jim=0L, nrad-1 do rvl_corrected3[iim,jim, *] = fill_gap_several_points_inavector( rvl_corrected3[iim,jim,*] )
        
        ; verification terre
        ; si valeur non nul sur Terre alors corrections faites
        for kim=0L, nazim-1 do begin
           rvl_corrected3[iim,*,kim] = fill_gap_several_points_inavector( rvl_corrected3[iim, *,kim])
          
           mean_terre_verif = 0
           std_terre_verif = 0
           num_terre_verif = 0
           
           p = where( data.rvllandcoverage[iim, *, kim] eq flag_terre and $
                      data.rvldcobs[iim,*,kim] ne def_val_rvl and $
                      abs(data.rvldcobsstd[iim,  *, kim]) lt max_dc_std  and $
                      rad_bon[iim,*] eq 1 )        

           if n_elements(p) gt 3 then begin
              median_terre_verif = median(rvl_corrected3[iim,p, kim])
              std_terre_verif = stddev(rvl_corrected3[iim,p, kim])
              
              pp = where( data.rvllandcoverage[iim, *, kim] eq flag_terre and $
                          abs(rvl_corrected3[iim,*, kim]-median_terre_verif) lt std_terre_verif*2. and $
                          data.rvldcobs[iim,*,kim] ne def_val_rvl and $
                          abs(data.rvldcobsstd[iim,  *, kim]) lt max_dc_std  and $
                          rad_bon[iim,*] eq 1) 
                 
              if pp[0] ne -1 then begin
                 median_terre_verif2 = median(rvl_corrected3[iim, pp, kim])                                  
                 if n_elements(pp) gt 2. then rvl_corrected3[iim,*, kim] = rvl_corrected3[iim,*, kim]-median_terre_verif2             
                
              endif             
              pp= -1           
           endif
           p=-1
        endfor
        
        ;; lissage pour une meilleure definition des profils azimutaux de scalloping
        for jim=0L, nrad-1 do begin
           rvl_corrected3[iim,jim,*] = fill_gap_several_points_inavector( rvl_corrected3[iim,jim, *] )
           p = where( radvelec[iim, jim, *] ne 0 and data.rvlnrcs[iim, jim, *] ne def_val_rvl and $
                      abs(data.rvldcobsstd[iim, jim, *]) lt max_dc_std  and $
                      data.rvllandcoverage[iim, jim, *] eq flag_mer  ) ; rvl_corrected3[iim, jim, *] ne 0 and
           
           if p[0] ne -1 and rad_bon[iim,jim]  eq 1 then begin 
              nb_el = min([periode_signal*8, 42]) ;10
              
              if n_elements(p) gt nb_el then rvl_corrected3[iim,jim,p] = smooth( rvl_corrected3[iim,jim, p], nb_el )
          
              nb_el = 10
              if n_elements(p) gt nb_el then rvl_corrected3[iim,jim,p] = smooth( rvl_corrected3[iim,jim, p], nb_el )

              nb_el = 3
              if n_elements(p) gt nb_el then rvl_corrected3[iim,jim,p] = smooth( rvl_corrected3[iim,jim, p], nb_el )
           endif                  
        endfor
              
        ;; nouveaux profils azimutaux a partir de radvelec
        for kim=0L, nazim-1 do begin
           p  = where(rvl_corrected3[iim, *, kim] ne 0 and $
                      radvelec[iim,*, kim] ne 0 and $
                      data.rvlnrcs[iim,*, kim] ne def_val_rvl and $
                      abs(data.rvldcobsstd[iim,  *, kim]) lt max_dc_std )
            
           if p[0] ne -1 then doppler_win_verifie[iim, kim] = median(radvelec[iim,p, kim]-rvl_corrected3[iim,p, kim])
           p=-1
           
           p1  = where( radvelec[iim,*, kim] ne 0 and data.rvlnrcs[iim,*, kim] ne def_val_rvl and $
                        abs(rvl_corrected3[iim, *, kim]) gt 0.0001 and $
                        abs(data.rvldcobsstd[iim,  *, kim]) lt max_dc_std  and $ 
                        data.rvllandcoverage[iim, *, kim] eq flag_terre)
           
           if n_elements(p1) gt 2 then doppler_win_ver_terre[iim, kim] = median(radvelec[iim,p1, kim]-rvl_corrected3[iim,p1, kim])
          
           p2  = where( radvelec[iim,*, kim] ne 0 and data.rvlnrcs[iim,*, kim] ne def_val_rvl and $
                        abs(rvl_corrected3[iim, *, kim]) gt 0.0001 and $
                        abs(data.rvldcobsstd[iim,  *, kim]) lt max_dc_std  and $
                        data.rvllandcoverage[iim, *, kim] eq flag_mer)
            
           if n_elements(p2) gt 2 then doppler_win_ver_mer[iim, kim] = median(radvelec[iim,p2, kim]-rvl_corrected3[iim,p2, kim])          

           p1 = -1
           p2 = -1
        endfor
        
        for kim=0L, nazim-1 do begin
           p =where( abs(doppler_win_ver_terre[iim, *]) gt err_ano)
           if p[0] ne -1 then doppler_win_ver_terre[iim,p] = fltarr(n_elements(p))
           p=-1

           p =where( abs(doppler_win_ver_mer[iim, *]) gt err_ano)
           if p[0] ne -1 then doppler_win_ver_mer[iim,p] = fltarr(n_elements(p))
           p=-1
           
           p1  = where(radvelec[iim,*, kim] ne 0 and data.rvlnrcs[iim,*, kim] ne def_val_rvl and $ 
                       abs(data.rvldcobsstd[iim,  *, kim]) lt max_dc_std  and $ 
                       data.rvllandcoverage[iim, *, kim] eq flag_terre  and abs(rvl_corrected3[iim, *, kim]) gt 0.0001 )
           
           p2  = where(radvelec[iim,*, kim] ne 0 and data.rvlnrcs[iim,*, kim] ne def_val_rvl and $ 
                       abs(data.rvldcobsstd[iim,  *, kim]) lt max_dc_std  and $
                       data.rvllandcoverage[iim, *, kim] eq flag_mer and abs(rvl_corrected3[iim, *, kim] ) gt 0.0001 ) 
           
           if n_elements(p2) gt n_elements(p1) then doppler_win_ver_tm[iim, kim] = doppler_win_ver_mer[iim, kim]
           if n_elements(p2) le n_elements(p1) then doppler_win_ver_tm[iim, kim] = doppler_win_ver_terre[iim, kim]
           
           p1 = -1
           p2 = -1
        endfor
        
     endfor
    
     doppler_unif_ver = unification_profil(doppler_win_ver_tm, std_doppler_terre_et_mer*0+1, num_doppler_terre_et_mer*0+num_dop_min+1, deriv_num_terre_et_mer*0+4, num_dop_min )   

     for iim=0L, nbwin-1 do begin
        for kim=0L, nazim-1 do begin
           if abs(doppler_win_ver_tm[iim, kim]-doppler_unif_ver[kim]) gt err_ano then doppler_win_ver_tm[iim, kim] = 0.
        endfor
     endfor
   
     
     loadct, 39
     device, decomposed = 0
      
     window, 8, xsize = xsize
     plot, doppler_win_verifie[0,*], /nodata, col=0, back=-1, yrange=rrange*1.3 ,ytitle='Scalloping [Hz]', xtitle='Nb of pixels in the azimuth direction', title=title, /xstyle, /ystyle
     for iim=0L, nbwin-1 do oplot, doppler_win_ver_tm[iim,*], thick=2, col=50*iim+50
     
     pointille, [0,nazim-1], rrange*1.3, col=0   
     filename=dirocn[i_ocn]+'/Correc_azim_'+time_str+suf+'.png'
     saveimage, filename, PNG='PNG'
     
     rvl_corrected4= fltarr(nbwin, nrad, nazim)     
     for iim=0L, nbwin-1 do begin
        
        for jim=0L, nrad-1 do begin
           p = where(doppler_win_ver_tm[iim,*] ne 0 and radvelec[iim,jim,*] ne 0 and data.rvlnrcs[iim,jim,*] ne def_val_rvl )
           if p[0] ne -1 then rvl_corrected4[iim,jim,p] = radvelec[iim,jim,p]-doppler_win_ver_tm[iim,p] 
           p= -1
        endfor
        
        p = where(rad_bon[iim,*] eq 0)
        if p[0] ne -1 then rvl_corrected4[iim,p,*] = 0     
        p=-1
     endfor
      
     print, '"From read_dop_nc " -- radar wavelength = 5.405 GHz'
     p = where(data.rvllon ne def_val_rvl)
     print, ' "From read_dop_nc " -- incidence angle [deg.] (min et max) : ', min(data.rvlincidenceangle[p]), max(data.rvlincidenceangle[p])

     ;; 5 - Last part : transformation from Doppler anomaly to radial
     ;;     velocity -  printing 
     
     vrad2 = fltarr(nbwin, nrad, nazim)
     vrad2[p] = -2.99792/(5.405*10.)*rvl_corrected4[p]/sin(data.rvlincidenceangle[p]*!dtor)/2.
                                      
     p = where(data.rvllon ne def_val_rvl and abs(data.rvldcobsstd) lt max_dc_std )  
     if p[0] ne -1 then begin 
        print, '"From read_dop_nc " -- retrieved radial velocities [m/s] (min et max) : ', min(vrad2[p]), max(vrad2[p]), mean(vrad2[p]),median(vrad2[p])
        print, '"From read_dop_nc " -- retrieved Doppler Anomalie [Hz] (min, max, mean and median) : ', min(rvl_corrected0[p]), max(rvl_corrected0[p]), mean(rvl_corrected0[p]), median(rvl_corrected0[p])
     endif
     p=-1
     
     A = FINDGEN(17) * (!PI*2/16.)    
     USERSYM, COS(A)/2, SIN(A)/2, /FILL
  
     tif = 1
 
     p = where(data.rvllat ne def_val_rvl  and rvl_corrected4 ne 0. and data.rvldcobs ne def_val_rvl )
     if p[0] ne -1 then begin
        minlon = min(data.rvllon[p], /nan)
        minlat = min( data.rvllat[p], /nan)
        maxlon = max( data.rvllon[p], /nan)
        maxlat = max( data.rvllat[p], /nan)
        
        pta = lonlat2utm(minlon, minlat)
        ptb = lonlat2utm(maxlon, maxlat)
        ptc = lonlat2utm(maxlon, minlat)
        ptd = lonlat2utm(minlon, maxlat)
        ys_test1 = ptb.n-ptc.n
        ys_test2 = ptd.n-pta.n
	if pta.zone eq ptc.zone then xs_test1 = ptc.e-pta.e else xs_test1 = ptc.e+500000.-pta.e
	if pta.zone eq ptc.zone then xs_test2 = ptb.e-ptd.e else xs_test2 = ptb.e+500000.-ptd.e 

        ptl = lonlat2utm(data.rvllon[0, 0, nazim-1],data.rvllat[0, 0, nazim-1] )
        ptm = lonlat2utm(data.rvllon[0, nrad-1, nazim-1],data.rvllat[0, nrad-1, nazim-1] )
        ptn = lonlat2utm(data.rvllon[0, nrad-1, 0],data.rvllat[0, nrad-1, 0] )
        pto = lonlat2utm(data.rvllon[0, 0, 0],data.rvllat[0, 0, 0] )
   
        if pta.zone eq ptb.zone then s_ratio = (ptb.n-pta.n)/(ptb.e-pta.e)
        if pta.zone eq ptb.zone-1 then s_ratio = (ptb.n-pta.n)/(ptb.e+500000.-pta.e)
        l_ratio = (maxlat-minlat)/(maxlon-minlon)
        ysize = xsize*l_ratio           
        
        limits = [minlon, minlat, maxlon, maxlat]
        GS = [0.001, 0.001]
        gtcitationgeokey = 'WGS 84 lonlat EPSG 4326'
        TRIANGULATE, data.rvllon[p], data.rvllat[p], tr, b
        geotiff_out = {modeltiepointtag:[0, 0, 0, minlon, maxlat, 0 ], modelpixelscaletag:[GS,0], GTCITATIONGEOKEY:gtcitationgeokey}        
         
        window, 9, xsize = xsize, ysize = ysize
        scatter, data.rvllon[p], data.rvllat[p], float(vrad2[p]), /cbar, zrange=rvrange, back=-1, cltable=41, $
          xtitle='Long. [deg.]', ytitle='Lat [deg.]', psym=8, title=title, cbtitle='RVL corrected final [m/s]', symsize=1, thick= 1, $
          xrange=[minlon, maxlon], yrange=[minlat, maxlat]

        pp = where(data.rvllandcoverage eq flag_terre)   
        if pp[0] ne -1 then oplot, data.rvllon[pp], data.rvllat[pp], col=50, psym=8, thick=2
        
        filename=dirocn[i_ocn]+'/RVL_corrected_fin_'+time_str+suf+'.png'
        saveimage, filename, PNG='PNG'
        
        if tif then begin           
           print, 'Tiff - RVL corrected : ', n_elements(pp)        
           rvl_cor_0_int = trigrid( data.rvllon[p], data.rvllat[p], float(vrad2[p]), tr, GS, limits, missing =  0.) 
           write_tiff, dirocn[i_ocn]+'/RVL_corrected_fin'+time_str+'.tif', reverse(rvl_cor_0_int ,2), geotiff=geotiff_out, /float           
        endif        
     endif
     p=-1
     
     p = where(data.rvllat ne def_val_rvl and $
               data.rvldcmiss ne def_val_rvl and $
               data.rvldcobs ne def_val_rvl and $
               data.rvlnrcs ne def_val_rvl and $
               alog(data.rvlnrcs) gt -8 ) 
     
     if p[0] ne -1 then begin
        
        window, 10, xsize = xsize, ysize = ysize
        scatter, data.rvllon[p], data.rvllat[p], alog(data.rvlnrcs[p]), /cbar, zrange=nrange, back=-1, cltable = -1, $
                 xtitle='Long. [deg.]', ytitle='Lat. [deg.]', psym=8, title=title, cbtitle='Normalized Radar Cross Section [dB]', $
                 symsize=1, thick= 1, $
                 xrange=[minlon, maxlon], yrange=[minlat, maxlat]
        
        filename=dirocn[i_ocn]+'/Nrcs_coloc_'+time_str+'.png'
        saveimage, filename, PNG='PNG'
        
        if tif then begin   
           
           limits = [minlon, minlat, maxlon, maxlat]
           GS = [0.001, 0.001]
           gtcitationgeokey = 'WGS 84 lonlat EPSG 4326'
           
           TRIANGULATE, data.rvllon[p], data.rvllat[p], tr, b
           geotiff_out = {modeltiepointtag:[0, 0, 0, minlon, maxlat, 0 ], modelpixelscaletag:[GS,0], GTCITATIONGEOKEY:gtcitationgeokey} 
                                
           miss_val = mean(float(alog(data.rvlnrcs[p])))
           rvl_cor_2_int = trigrid( data.rvllon[p], data.rvllat[p], float(alog(data.rvlnrcs[p])), tr, GS, limits, missing = Nan) ;else stop
           print, 'Tiff - NRCS : ', n_elements(pp)
           write_tiff, dirocn[i_ocn]+'/RVL_nrcs_'+time_str+'.tif', reverse(rvl_cor_2_int ,2), geotiff=geotiff_out, /float
           pp = -1
        endif
        
        pp = where(data.rvllat ne def_val_rvl and $
                   data.rvldcmiss ne def_val_rvl and $
                   data.rvldcobs ne def_val_rvl and $
                   data.rvllandcoverage ne flag_terre)
        vrad_init = float(-2.99792/(5.405*10.)*(dc_to_treat)/sin(data.rvlincidenceangle*!dtor)/2.)
        if pp[0] ne -1 then mean_terre_init = mean(vrad_init[pp])
        if pp[0] ne -1 then std_terre_init = stddev(vrad_init[pp])
        pp= -1        
        
        window, 15, xsize = xsize, ysize = ysize
        scatter, data.rvllon[p], data.rvllat[p], vrad_init[p], /cbar, $
                 zrange=[mean_terre_init-(rvrange[1]-rvrange[0])/2, mean_terre_init+(rvrange[1]-rvrange[0])/2 ], $
                 back=-1, cltable =41, $
                 xtitle='Long. [deg.]', ytitle='Lat. [deg.]', psym=8, title=title, cbtitle='RVL intial [m/s]', symsize=1, thick= 1, $
                 xrange=[minlon, maxlon], yrange=[minlat, maxlat]
        filename=dirocn[i_ocn]+'/RVL_init_'+time_str+'.png'
        saveimage, filename, PNG='PNG'
        
        if tif then begin
           pp = where(data.rvllat ne def_val_rvl)
           
           if pp[0] ne -1 then rvl_cor_3_int = trigrid( data.rvllon[pp], data.rvllat[pp], float(vrad_init[pp]), tr, GS, limits, missing =  0.) else stop
           print, 'Tiff - INIT : ', n_elements(pp)
           write_tiff, dirocn[i_ocn]+'/RVL_init_'+time_str+'.tif', reverse(rvl_cor_3_int ,2), geotiff=geotiff_out, /float
           pp = -1
        endif
                
     endif
      p=-1
      
  endif else print, 'Pb with file SAFE'   
      endif else print, 'STATUS netcdf file - no good', status
  
   endif else print,'No safe file in ', dirocn[i_ocn]
        
 endfor

end

