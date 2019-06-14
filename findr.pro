function findr, nys, sector = sector
;;; rind the DFB radius R. based on two-point observation values.
;;; nys: positive means evening side, negative means evening side
;;; sector: must contain "dusk" or "dawn". default is dusk
;; r values. change here.
r_dusk_evening = 2.1
r_dusk_morning = 1.3
r_dawn_evening = 0.5
r_dawn_morning = 1.1
;; no need to change below
if strmatch(sector, 'dawn*') then begin
	r_evening = r_dawn_evening
	r_morning = r_dawn_morning
endif else begin
	r_evening = r_dusk_evening
	r_morning = r_dusk_morning
endelse

r_values = nys
r_values[*] = 0.
for i = 0, n_elements(nys)-1 do begin
	if nys[i] gt 0. then r_values[i] = r_evening else r_values[i] = r_morning	
endfor

return, r_values
end
