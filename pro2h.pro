function pro2h, value, h, h_out, quantity = quantity, method = method, cu_type = cu_type, b_type = b_type
;; project current to h = any location (can be other things than h)
;; quantity: specify type of quantity
;; cu_type: the method used to calculate the current, 'l' or 'db'
;; type: specify method of projection
if ~keyword_set(quantity) then quantity = 'i'
if ~keyword_set(cu_type) then cu_type = 'l'
if ~keyword_set(method) then method = 'linear_med_abs'

;;; make sure the vectors have the same lengths
if (n_elements(value) ne n_elements(h)) or (n_elements(value) ne n_elements(h_out)) or (n_elements(h_out) ne n_elements(h)) then begin
	print, 'PRO2H: Numbers not agree!'
	stop
endif

;;; projection settings
if strcmp(b_type, 'bin') then begin
	if strcmp(cu_type, 'l') then begin
		if strcmp(quantity, 'i') then begin
			hcntrs_more = [-0.875000, -0.625000, -0.375000, -0.125000, 0., 0.125000, 0.375000, 0.625000, 0.875000]
			med_more = [9.15901e+006, 8.30704e+006, 7.71190e+006, 6.05853e+006, 0., -5.28624e+006, -7.33338e+006, -8.32665e+006, -9.16148e+006]
			hcntrs_abs = [0.125000, 0.375000, 0.625000, 0.875000]
			med_abs = [7.80077e+006, 8.39430e+006, 8.84371e+006, 9.16148e+006]
		endif
		if strcmp(quantity, 'j') then begin
			hcntrs_abs = [0.125000, 0.375000, 0.625000, 0.875000]
			med_abs = [11.7778, 8.45031, 7.93155, 7.38024]
		endif
	endif
	
	if strcmp(cu_type, 'db') then begin
		if strcmp(quantity, 'i') then begin
			;;; for only good dbxbin
			;hcntrs_more = [-0.875000, -0.625000, -0.375000, -0.125000, 0., 0.125000, 0.375000, 0.625000, 0.875000]
			;med_more = [9.78891e+006, 8.52720e+006, 8.15212e+006, 6.40621e+006, 0., -5.50418e+006, -7.88000e+006, -8.70876e+006, -9.89348e+006]
			;hcntrs_abs = [0.1, 0.3, 0.5, 0.7, 0.9]
			;med_abs = [8.01381e+006, 7.88000e+006, 8.70876e+006, 9.54937e+006, 1.00439e+007]
			;;; for good dbxbuse and good boutxbin
			hcntrs_abs = [0.1, 0.3, 0.5, 0.7, 0.9]
			med_abs = [5.51559e+006, 6.24579e+006, 7.84819e+006, 9.24917e+006, 9.91777e+006]
		endif
		if strcmp(quantity, 'j') then begin
			hcntrs_abs = [0.1, 0.3, 0.5, 0.7, 0.9]
			med_abs = [12.3273, 8.75507, 7.03429, 8.80635, 8.35022]
		endif
	endif
endif

if strcmp(b_type, 'bave') then begin
	if strcmp(cu_type, 'l') then begin
		if strcmp(quantity, 'i') then begin
			;;; for only good dbxbin
			hcntrs_abs = [0.125000, 0.375000, 0.625000, 0.875000]
			med_abs = [9.10232e+006, 9.92999e+006, 9.88949e+006, 9.06804e+006]
		endif
	endif
	
	if strcmp(cu_type, 'db') then begin
		if strcmp(quantity, 'i') then begin
			;;; for only good dbxbin
			;hcntrs_abs = [0.125000, 0.375000, 0.625000, 0.875000]
			;med_abs = [9.29492e+006, 1.00527e+007, 9.94263e+006, 1.00590e+007]
			;;; for good dbxbuse and good boutxbin
			hcntrs_abs = [0.125000, 0.375000, 0.625000, 0.875000]
			med_abs = [8.46460e+006, 9.59267e+006, 9.80297e+006, 1.00590e+007]
		endif
	endif
endif
 
;;; begin projecting
if strcmp(method, 'linear') then begin
	results_out = value*h_out/h
endif

if strmatch(method, '*_med*') then begin
	if strcmp(method, 'linear_med_more') then begin
		value_nominal = interpol(med_more, hcntrs_more, h)
		value_nominal_pro = interpol(med_more, hcntrs_more, h_out)
	endif
	if strcmp(method, 'linear_med_abs') then begin
		value_nominal = interpol(med_abs, hcntrs_abs, abs(h))
		value_nominal_pro = interpol(med_abs, hcntrs_abs, abs(h_out))
	endif
	results_out = value*value_nominal_pro/value_nominal
endif

return, results_out
end
