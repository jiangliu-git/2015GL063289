pro main_model_plot
thm_init
@folders

;;;;;;;; Current values
;if MLT<MLT_seperate then I_out = 0.15MA, I_in = 0.05MA
;if MLT>MLT_seperate then I_out = 0.05MA, I_in = 0.15MA

;;;;;;; Values
r_dfb = 1. ;; in RE
space_dfb = 0.5 ;; in RE
MLT_seperate = -1 ;; the seperator of dawn/dusk sector of the magnetotail, in MLT

;;;;;;;; Locations
;; There are four lines of DFBs (Lines 1-4).

;; Xn is the X location for Line n.
X = dblarr(4) ;; in RE
X[0] = -8
for i = 1,3 do begin
	X[i] = X[0]-i*(r_dfb+space_dfb)
endfor

;; Y locations for the for lines n=1-4
Y_seperator = X*tan(MLT_seperate*15./180*!pi) ;; Y of the seperator of dawn/dusk sector
;; for each line, DFB centers:
Y_centers_dusk_1 = r_dfb+0.5*space_dfb
Y_centers_dawn_1 = -(r_dfb+0.5*space_dfb)
Y_centers_dusk_2 = 3*(r_dfb+0.5*space_dfb)
Y_centers_dawn_2 = -3*(r_dfb+0.5*space_dfb)
store_data, 'Y_centers_line1', data = Y_seperator[0]+[Y_centers_dusk_1, Y_centers_dawn_1]
store_data, 'Y_centers_line2', data = Y_seperator[1]+[Y_centers_dusk_2, Y_centers_dusk_1, Y_centers_dawn_1, Y_centers_dawn_2]
store_data, 'Y_centers_line3', data = Y_seperator[2]+[Y_centers_dusk_2, Y_centers_dusk_1, Y_centers_dawn_1, Y_centers_dawn_2]
store_data, 'Y_centers_line4', data = Y_seperator[3]+[Y_centers_dusk_1, Y_centers_dawn_1]
Y_lines = ['Y_centers_line1','Y_centers_line2','Y_centers_line3','Y_centers_line4']

;============= Plot ==================================
;;;;;;;;;;;; plot current locations
;;; size of currents, in times of r_dfb
cur_small = 0.09
cur_large = cur_small*sqrt(3)
xsize = 4
ysize = 3
;; locations
xrange = [-14.99, -5.01]
yrange = [9.99,-3.99]
Y_sepplot = xrange*tan(MLT_seperate*15./180*!pi) ;; Y of the seperator of dawn/dusk sector
;; start plot
popen, pic_folder+'/model_fac'
print_options,xsize=xsize, ysize=ysize;; use this for single plot
plot, Y_sepplot, xrange, yrange = xrange, xrange = yrange, ytitle = 'X [R!dE!n]', xtitle = 'Y [R!dE!n]', title = '', xstyle = 1, ystyle = 1, /iso, line = 1
for i = 0, 3 do begin
	get_data, Y_lines[i], data = y_centers
	y_in = y_centers-0.5*r_dfb
	y_out = y_centers+0.5*r_dfb
	x_this = replicate(X[i], n_elements(y_centers))
	for j = 0, n_elements(x_this)-1 do begin
		if y_centers[j] lt Y_seperator[i] then begin
			size_in = cur_large
			size_out = cur_small
		endif
		if y_centers[j] gt Y_seperator[i] then begin
			size_in = cur_small
			size_out = cur_large
		endif
		draw_circle, y_in[j], x_this[j], size_in*r_dfb, /fill, color_fill = 2
		draw_circle, y_out[j], x_this[j], size_out*r_dfb, /fill, color_fill = 6
	endfor
endfor
xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), !y.crange[0]+0.85*(!y.crange[1]-!y.crange[0]), '(a)', charsize = 1.2, /data
xyouts, 1.4, -6, '23MLT', charsize = 0.7, /data
xyouts, !x.crange[0]+0.5*(!x.crange[1]-!x.crange[0]), !y.crange[0]+1.03*(!y.crange[1]-!y.crange[0]), 'Equatorial locations of model wedgelets', align = 0.5
pclose

;;;;;;;;;;; plot results
;; values
ytitle = delta_letter+'B [nT]'
xsize = 6
ysize = 5
;; position setting for the location plot
n_panels = 2
left_margin = 0.15
right_margin = 0.1
top_margin = 0.1
bot_margin = 0.15
vspace = 0.01
positions = panel_positions([1, n_panels], lr_margins = [left_margin, right_margin], bt_margins = [bot_margin, top_margin], space = [0., vspace], height = height)

;; load data
restore, 'dmag_ground.sav'
store_data, 'ground', data={mlt:mlt_ground, dmag:dmag}
restore, 'dmag_goes.sav'
store_data, 'geo', data={mlt:mlt_goes, dmag:dmag}
locations = ['ground', 'geo']
abc = ['(b)', '(c)']
symsize = 0.6

popen, pic_folder+'/model_db'
print_options,xsize=xsize, ysize=ysize;; use this for single plot
for i = 0, n_elements(locations)-1 do begin
	if i eq n_elements(locations)-1 then begin
		xtitle = 'MLT [hr]'
		xticknames = ''
	endif else begin
		xtitle = ''
		xticknames = replicate(' ', 50)
	endelse
	if i eq 0 then begin
		title = 'Change of B caused by model wedgelets'
	endif else begin
		title = ''
	endelse
	get_data, locations[i], data = data
	mlt = data.mlt
	mlt[0] = abs(mlt[0])
	dmag = data.dmag
	plot, mlt, dmag[*,0], /noerase, title = title, xtitle = xtitle, xrange = [-12., 12.], xstyle = 1, yrange = [-59.9, 59.9], ystyle = 1, ytitle = ytitle, position = positions[*, n_elements(locations)-1-i], xtickname = xticknames, /nodata
	oplot, mlt, dmag[*,0], thick = l_thick, line = 5
	oplot, mlt, dmag[*,1], thick = l_thick
	;oplot, mlt, dmag[*,0], psym = 4, symsize = symsize
	;oplot, mlt, dmag[*,1], psym = 5, symsize = symsize
	oplot, !x.crange, [0,0], line = 1
	xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), !y.crange[0]+0.8*(!y.crange[1]-!y.crange[0]), abc[i], charsize = 1.2
	if i eq 0 then begin
		xyouts, -5, 30, 'D', /data
		xyouts, -0.2, 40, 'H', /data
		xyouts, !x.crange[0]+0.95*(!x.crange[1]-!x.crange[0]), !y.crange[0]+0.8*(!y.crange[1]-!y.crange[0]), 'Ground!c45!Uo!n Latitude', charsize = 0.8, alignment = 1
	endif
	if i eq 1 then begin
		xyouts, 0, -25, 'D', /data
		xyouts, -2, 40, 'H', /data
		xyouts, !x.crange[0]+0.95*(!x.crange[1]-!x.crange[0]), !y.crange[0]+0.8*(!y.crange[1]-!y.crange[0]), 'GEO!c10!Uo!n Latitude', charsize = 0.8, alignment = 1
	endif
endfor
pclose
    
stop
end
