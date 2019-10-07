
;; read_safe_s1
;;    read the S1 format
;;    INPUT : SAFE file
;;    OUTPUT : following information on the file :
;;              file: filef, nb_file:nb_file, lut_sigma:a.lut_sigma, geo_loc_grid:a.geo_loc_grid, $
;;              mission:a.mission, product:a.product, polar:a.polar, $
;;              mode:a.mode, swath:a.swath, starttime:start_time, stoptime:stop_time, pass:pass, check:1
;;    Auteurs :C. Danilo, juin 2018


function read_safe_s1, dir

  file_safe = file_search(dir, '*.safe')
  
  quiet = 1
  
  IF n_elements(file_safe) eq 0 THEN stop $
  ELSE print,  'Read Safe file : '+file_safe
  
  str = ''
  
  if not(quiet) then print, 'Lecture Fichier SAFE '
  
  openr, lun, file_safe[0], /get_lun
  readf, lun, str

  while eof(lun) NE 1 do BEGIN
  
    if (strpos(str, '<informationPackageMap>'))[0] ne -1 then BEGIN
      if not(quiet) then print, '--> Lecture informationPackageMap'
      READF, lun, str
      pos_txt_deb = strpos(str, 'textInfo')
      pos_txt_fin = strpos(str, 'dmdID')
      info_txt = strmid(str, pos_txt_deb+10, pos_txt_fin-pos_txt_deb-12)
      pos_dmdID_fin = strpos(str, 'pdiID')
      
      if not(quiet) then print, '     Text Info : ', info_txt
   
    endif
    
    if (strpos(str, 'metadataObject ID="platform"'))[0] ne -1 then BEGIN
      if not(quiet) then print, '--> Lecture Plateforme'
      WHILE  (strpos(str, '</metadataObject>'))[0] eq -1 do  BEGIN
        READF, lun, str
        IF (strpos(str, '<safe:familyName'))[0] ne -1 then $
          capteur_famille = extrait_entre_inegalite(str)
          
        if (strpos(str, '<safe:number'))[0] ne -1 then $
          capteur_number = extrait_entre_inegalite(str)
          
        if (strpos(str, '<safe:extension>'))[0] ne -1 then $
          BEGIN
          
          platform_extension = ''
          platform_extension_name = ''
          WHILE (strpos(str, '</safe:extension>'))[0] eq -1 do BEGIN
            IF (stregex(str, '>[^*]+<'))[0] NE -1 THEN begin
              platform_extension = [ platform_extension, $
                extrait_entre_inegalite(str) ]
              platform_extension_name = [ platform_extension_name, $
                extrait_champ_avant_inegalite(str)]
            endif
            READF, lun, str
          ENDWHILE
        endif
      ENDWHILE
      if not(quiet) then begin
        print,  '    Famille du capteur : ', capteur_famille+' '+capteur_number
        print,  '    Information en plus : ', platform_extension_name
        print,  '     - - champ associés :', platform_extension
      endif
      
    endif
    
    if (strpos(str, 'metadataObject ID="measurementOrbitReference"'))[0] ne -1 then BEGIN
       if not(quiet) then  print, '--> Lecture measurement orbit'
      WHILE  (strpos(str, '</metadataObject>'))[0] eq -1 do BEGIN
      
        if (strpos(str, '<safe:extension>'))[0] ne -1 then BEGIN
        
          orbit_extension = ''
          orbit_extension_name = ''
          WHILE  (strpos(str, '</safe:extension>'))[0] eq -1 do BEGIN
          
            IF (stregex(str, '>[^*]+<',  /fold_case))[0] NE -1 THEN BEGIN
              orbit_extension = [ orbit_extension, $
                extrait_entre_inegalite(str) ]
                
              orbit_extension_name = [ orbit_extension_name, $
                extrait_champ_avant_inegalite(str)]
            ENDIF
            READF, lun, str
          ENDWHILE
          
        ENDIF
        READF, lun, str
      ENDWHILE
      if not(quiet) then begin
        print,'    Information en plus : ', orbit_extension_name
        print,'     - - champ associés :', orbit_extension
      endif
    endif
    
    if (strpos(str, '<metadataObject ID="generalProductInformation"'))[0] ne -1 then BEGIN
       if not(quiet) then  print, '--> Lecture generalProductInformation'
      geneprodinfo_extension = ''
      geneprodinfo_extension_name = ''
      
      WHILE  (strpos(str, '</metadataObject>'))[0] eq -1 do BEGIN
      
        IF (stregex(str, '>[^*]+<',  /fold_case))[0] NE -1 THEN BEGIN
          geneprodinfo_extension = [ geneprodinfo_extension, $
            extrait_entre_inegalite(str) ]
            
          geneprodinfo_extension_name = [ geneprodinfo_extension_name, $
            extrait_champ_avant_inegalite(str)]
        ENDIF
        READF, lun, str
      ENDWHILE
      if not(quiet) then begin
        print, '    Information en plus : ', geneprodinfo_extension_name
        print, '     - - champ associés :', geneprodinfo_extension
      endif
    endif
    
    if (strpos(str, '<metadataObject ID="acquisitionPeriod'))[0] ne -1 then BEGIN
      if not(quiet) then print, '--> Lecture acquisitionPeriode'
      
      WHILE  (strpos(str, '</metadataObject>'))[0] eq -1 do begin
        IF (strpos(str, '<safe:startTime>'))[0] ne -1 then $
          start_time = extrait_entre_inegalite(str)
        IF (strpos(str, '<safe:stopTime>'))[0] ne -1 then $
          stop_time = extrait_entre_inegalite(str)
        READF, lun, str
      ENDWHILE
      if not(quiet) then begin
        print, '    Start time : ',  start_time
        print, '    Stop time : ',  stop_time
      endif
    endif
    
    if (strpos(str, '<metadataObject ID="measurementFrameSet'))[0] ne -1 then BEGIN
      if not(quiet) then print, '--> Lecture measurementFrameset'
      
      WHILE (strpos(str, '</metadataObject>'))[0] eq -1 do begin
        if  (strpos(str, '<gml:coordinates>'))[0] ne -1 then $
          corner_coord = extrait_entre_inegalite(str)
        READF, lun, str
        
      endwhile
      if not(quiet) then print, '    Position geographiques des 4 coins de l''image : ', corner_coord
      
    ENDIF
    
    if (strpos(str, '<dataObjectSection>'))[0] NE -1 THEN BEGIN
      if not(quiet) then print,  '--> Lecture, dataObjectSection'
      
      file_loc = ''
      WHILE  (strpos(str, '</dataObjectSection>'))[0] eq -1 do BEGIN
      
        IF  (strpos(str, '<fileLocation'))[0] NE -1 THEN $
          file_loc = [ file_loc, $
          extrait_file_loc(str) ]
          
        READF, lun, str
      ENDWHILE
      
      if not(quiet) then print, '    Fichiers dependants a l''acquisition : ', file_loc
      
    endif
    
    READF, lun, str
      
  ENDWHILE
  
  free_lun, lun
  close, lun, /all
  
  file = file_loc[where(strmatch(file_loc, './measurement/*') EQ 1)]
  
  print, 'Fichiers de donnees ', file
  nb_file = n_elements(file)
  
  filef = strarr(nb_file)
  
  FOR i = 0, nb_file-1 do begin
    file_tiff = file_search(dir,file[i])
    res = file_test(file_tiff)
    
    IF res EQ 1 THEN BEGIN
    
      restiff  = query_tiff(file_tiff, info)
      
      if restiff eq 1 then begin
      
        pos_deb = strpos(file_tiff, 'measurement')
        pos_fin = strpos(file_tiff, '.tiff')
        file_name_tmp = strmid(file_tiff, pos_deb+12, pos_fin-pos_deb-12)
        
        file_geoloc = file_search(dir, 'annotation/*'+file_name_tmp+'*.xml')
        file_calib = file_search(dir, $
          'annotation/calibration/calibration-*'+file_name_tmp+'*.xml')
          
        IF file_calib NE '' THEN lut_sigma0 = read_calibration(file_calib) else stop
        
        IF file_geoloc NE '' and lut_sigma0.line[0] ne -1 THEN begin
          geo_loc_grid = read_geoloc_grid(file_geoloc)
          
          pos_pass = where(strmatch(orbit_extension_name,'pass'))
          if pos_pass[0] ne -1 and n_elements(pos_pass) eq 1 then pass=orbit_extension[pos_pass[0]] else pass=''
          
          if i eq 0 then a = {lut_sigma:lut_sigma0, geo_loc_grid:geo_loc_grid, $
            mission:lut_sigma0.mission, product:lut_sigma0.product, polar:lut_sigma0.polar, $
            mode:lut_sigma0.mode, swath:lut_sigma0.swath} else $
            b= {lut_sigma:lut_sigma0, geo_loc_grid:geo_loc_grid, $
            mission:lut_sigma0.mission, product:lut_sigma0.product, polar:lut_sigma0.polar, $
            mode:lut_sigma0.mode, swath:lut_sigma0.swath}
                           
            
        endif
        
      endif else begin
        print, 'Fichier lu mais pas fichier tiff'
        return,{file: file_tiff, check:1, product:geneprodinfo_extension[8], $
          mode:platform_extension[1], swath:orbit_extension[1], starttime:start_time, stoptime:stop_time}
      ENDelse
      
    endif else begin
      print, 'Erreur lecture de fichier'
      print, file_tiff
      return,{check: 0}
    endelse
    filef[i] = file_tiff
  endfor
 
  
  if nb_file eq 1 then $
     return, {file: filef, nb_file:nb_file, lut_sigma:a.lut_sigma, geo_loc_grid:a.geo_loc_grid, $
              mission:a.mission, product:a.product, polar:a.polar, $
              mode:a.mode, swath:a.swath, starttime:start_time, stoptime:stop_time, pass:pass, check:1}
    
  if nb_file eq 2 then begin
 
     return, {file: filef, nb_file:nb_file, lut_sigma:a.lut_sigma, geo_loc_grid:a.geo_loc_grid, $
              mission:a.mission, product:a.product, polar:a.polar, $
              mode:a.mode, swath:a.swath, starttime:start_time, stoptime:stop_time, pass:pass, $
              check:1,lut_sigma2:b.lut_sigma, geo_loc_grid2:b.geo_loc_grid, $
              mission2:b.mission, product2:b.product, polar2:b.polar, $
              mode2:b.mode, swath2:b.swath}
              
  endif

end
