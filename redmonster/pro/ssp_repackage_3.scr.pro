;
; ssp_repackage_3.scr.pro
;
; repackage the theoretical SSPs from Charlie  Conroy
; into a FITS file.
;
; 3rd iteration, 3 sept 2014

; I THINK you have to say ".com read_spec" before you run this script!

;infile = 'SSP_Padova_RRLIB_Kroupa_Z0.0190.out.spec'
;infiles = file_search('/data/SSP/SSP_*.spec')
infiles = file_search('/data/FSPS/SSP_*.spec')
nf = n_elements(infiles)

for j = 0, nf - 1 do begin & $
   print, j & $
   junk = read_spec1(infiles[j]) & $
   specimage = junk.spec & $
   lambda = junk[0].lambda & $
   im_junk = {lambda: lambda, spec: specimage} & $
   sub_junk = struct_selecttags(junk, except_tags=['LAMBDA', 'SPEC']) & $
   spos = strpos(infiles[j], '.', /reverse_search) & $
   froot = strmid(infiles[j], 0, spos) & $
   ofile = froot + '.fits' & $
   mwrfits, im_junk, ofile, /create & $
   mwrfits, sub_junk, ofile & $
endfor

