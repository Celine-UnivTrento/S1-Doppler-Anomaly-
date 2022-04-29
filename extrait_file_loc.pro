FUNCTION extrait_file_loc,  chaine
  pos_href = strpos(chaine, 'href')
  pos_ine_fin = strpos(chaine, '>')
  return, strmid(chaine, pos_href+6, pos_ine_fin-pos_href-8)
end
