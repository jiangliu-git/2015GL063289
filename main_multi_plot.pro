pro main_multi_plot
;; plot things from main_multi
thm_init

computer = 'I:'
;computer = '/home/jliu'

@folders

;;; choose method for normal direction
;method_suffix = 'minvar'
method_suffix = 'binxbout'

;;;; whether to translate to the propagate direction
v_suffix = '' ;;; X direction
;v_suffix = '_vperpX'
;v_suffix = '_vperp2X'
;v_suffix = '_vxyX'
;v_suffix = '_vdhtX'

;;;; which velocity to take when computing the expansion/contraction
;vc_use = 'vperp' ;; given by df_thick, experienced interpolation
vc_use = 'vperp2' ;; exactly the vperp over the few points of the DF.
;vc_use = 'vdht'

;;;; decide whether use events with large dy only
;large_dy = 'yes'
large_dy = 'no'

;;;; how to judge expand or contract
;judge_reshape_para = 'dvy'
;judge_reshape_para = 'angle'
judge_reshape_para = 'both'

;;; time range to get t_out (minus this seconds)
t_out_suf = '' ;; default: 15s
;t_out_suf = '_0'

list_suffix = '_earthward_df_'+method_suffix+t_out_suf

x = datain_simple(save_folder+'/x'+list_suffix+'.dat', dim = 2, type = 'double')
y = datain_simple(save_folder+'/y'+list_suffix+'.dat', dim = 2, type = 'double')
nfront = datain_simple(save_folder+'/nfront'+list_suffix+'.dat', dim = 2, type = 'double')
dxyz = datain_simple(save_folder+'/dxyz'+list_suffix+'.dat',  dim = 3, type = 'double')
dnangle = transpose(datain_simple(save_folder+'/dnangle'+list_suffix+'.dat', dim = 1, type = 'double'))*180./!pi
dfront = transpose(datain_simple(save_folder+'/dfront'+list_suffix+v_suffix+'.dat',  dim = 1, type = 'double'))
dnfront = transpose(datain_simple(save_folder+'/dnfront'+list_suffix+v_suffix+'.dat',  dim = 1, type = 'double'))
r_convex = transpose(datain_simple(save_folder+'/r_convex'+list_suffix+v_suffix+'.dat', dim = 1, type = 'double'))
i_convex = transpose(datain_simple(save_folder+'/i_convex'+list_suffix+v_suffix+'.dat', dim = 1, type = 'long'))
i_concave = transpose(datain_simple(save_folder+'/i_concave'+list_suffix+v_suffix+'.dat',dim = 1, type = 'long'))
i_expand = transpose(datain_simple(save_folder+'/i_expand'+list_suffix+v_suffix+'_'+vc_use+'_judge_'+judge_reshape_para+'.dat', dim = 1, type = 'long')) ;; of convex
i_contra = transpose(datain_simple(save_folder+'/i_contra'+list_suffix+v_suffix+'_'+vc_use+'_judge_'+judge_reshape_para+'.dat', dim = 1, type = 'long')) ;; of convex

rho = sqrt(x^2+y^2)
mlt = asin(-y/rho)/(15./180*!pi)
x_convex	= x(*, i_convex)
x_concave	= x(*, i_concave)
y_convex	= y(*, i_convex)
mlt_convex = mlt(*, i_convex)
nfront_convex	= nfront(*, i_convex)
dfront_convex	= dfront(i_convex)
dfront_concave	= dfront(i_concave)
dnangle_convex	= dnangle(i_convex)
dnangle_concave	= dnangle(i_concave)
dnfront_convex	= dnfront(i_convex)
x_expand	= x_convex(*, i_expand)
x_contra	= x_convex(*, i_contra)
dfront_expand	= dfront_convex(i_expand)
dfront_contra	= dfront_convex(i_contra)
dnangle_expand	= dnangle_convex(i_expand)
dnangle_contra	= dnangle_convex(i_contra)
dnfront_expand	= dnfront_convex(i_expand)

help, where(abs(dnangle_expand) gt 26)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; general quantities ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;dxyz_abs = abs(dxyz)
;rstat, dxyz_abs(0,*), x_med, x_low, x_high
;rstat, dxyz_abs(1,*), y_med, y_low, y_high
;rstat, dxyz_abs(2,*), z_med, z_low, z_high
;print, 'All dual observations:'
;print, 'dX: '+'lowq='+strcompress(string(x_low),/remove)+' med='+strcompress(string(x_med),/remove)+' highq='+strcompress(string(x_high),/remove)
;print, 'dY: '+'lowq='+strcompress(string(y_low),/remove)+' med='+strcompress(string(y_med),/remove)+' highq='+strcompress(string(y_high),/remove)
;print, 'dZ: '+'lowq='+strcompress(string(z_low),/remove)+' med='+strcompress(string(z_med),/remove)+' highq='+strcompress(string(z_high),/remove)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; results with requirements ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;dy_req = 0.5
;i_smalldy = where(abs(dfront) gt dy_req, n_smalldy)
;i_smalldy_convex = where(abs(dfront_convex) gt dy_req, n_smalldy_convex)
;print, n_smalldy
;print, n_smalldy_convex
;print, double(n_smalldy_convex)/double(n_smalldy)
;
;i_smalldy_expand = where(abs(dfront_expand) gt dy_req, n_smalldy_expand)
;i_smalldy_contra = where(abs(dfront_contra) gt dy_req, n_smalldy_contra)
;print, n_smalldy_expand
;print, n_smalldy_contra
;print, double(n_smalldy_expand)/double(n_smalldy_expand+n_smalldy_contra)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; x dependence of radius and convex/concave ratio ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;; set the range of events
;;;xbin_bounds = [-15., -12., -11, -10, -9, -8, -6] ;; highest resolution, used for ratio_convex, 15s.
;;xbin_bounds = [-15., -11., -8.5, -6] ;; for radius, 15s 
;
;
;store_data, 'convex_ratio', data = {xbin_bounds:[-15., -12., -11, -10, -9, -8, -6], x_abc:0.8, y_abc:0.1, yrange:[0.,1.]}
;store_data, 'radius', data = {xbin_bounds:[-15., -11., -8.5, -6], x_abc:0.8, y_abc:0.8, yrange:[0., 2.9]}
;store_data, 'width', data = {xbin_bounds:[-15., -11., -8.5, -6], x_abc:0.8, y_abc:0.8, yrange:[0., 5.8]}
;;things = ['convex_ratio', 'radius']
;things = ['convex_ratio', 'width']
;
;;;; plot parameters
;xrange_plot = [-6., -15.]
;abc = ['(b)', '(c)']
;left_margin = 0.15
;right_margin = 0.01
;top_margin = 0.05
;bot_margin = 0.05
;vspace = 0.013
;rate_bar = 0.1
;n_panels = n_elements(things)
;
;positions = panel_positions([1, n_panels], lr_margins = [left_margin, right_margin], bt_margins = [top_margin, bot_margin], space = [0., vspace], height = height)
;
;popen, pic_folder+'/convex_radius'
;print_options,xsize=2.6,ysize=3.5
;for i_plot = 0, n_panels-1 do begin
;	if i_plot eq n_panels-1 then begin
;	    xticknames = ''
;		xtitle = 'X [R!dE!n]'
;	endif else begin
;	    xticknames = replicate(' ', 59)
;		xtitle = ''
;	endelse
;	get_data, things(i_plot), data = thing
;	xbin_bounds = thing.xbin_bounds
;	;; arrays for results
;	n_concave_arr = intarr(n_elements(xbin_bounds)-1)
;	n_concave_arr(*) = !values.f_nan
;	n_convex_arr = n_concave_arr
;	n_expand_arr = n_convex_arr
;	n_contract_arr = n_convex_arr
;	r_mean_arr = fltarr(n_elements(xbin_bounds)-1)+!values.f_nan
;	r_std_arr = r_mean_arr
;	r_stat = fltarr(n_elements(xbin_bounds)-1, 3)+!values.f_nan ;; lower quartile, median, higher quartile
;	dfront_stat = r_stat
;	dfront_concave_stat = r_stat
;	dfront_convex_stat = r_stat
;	dfront_expand_stat = r_stat
;	dfront_contract_stat = r_stat
;	for i_dis = 0, n_elements(xbin_bounds)-2 do begin
;		xrange = xbin_bounds(i_dis:i_dis+1)
;	
;		i_this = where((x(0,*) gt xrange(0)) and (x(0,*) lt xrange(1)) and (x(1,*) gt xrange(0)) and (x(1,*) lt xrange(1)), j_this)
;		if j_this gt 0 then begin
;			dfront_this = dfront[i_this]
;			dnangle_this = dnangle[i_this]
;			rstat, abs(dfront_this), med_dfront, hing1, hing2
;			dfront_stat[i_dis, *] = [[hing1], [med_dfront], [hing2]]
;		endif
;	
;		i_this_convex = where((x_convex(0,*) gt xrange(0)) and (x_convex(0,*) lt xrange(1)) and (x_convex(1,*) gt xrange(0)) and (x_convex(1,*) lt xrange(1)), j_this_convex)
;		n_convex_arr(i_dis) = j_this_convex
;		if j_this_convex gt 0 then begin
;			r_this = r_convex(i_this_convex)
;			if strcmp(large_dy, 'yes') then begin
;				dfront_convex_this = dfront_convex(i_this_convex)
;				i_small_dy = where(abs(dfront_convex_this) gt 0.5, n_small_dy)
;				if n_small_dy gt 0 then r_bin = r_this(i_small_dy) else r_bin = !values.f_nan
;			endif else r_bin = r_this
;			rstat, r_bin, med_r, hing1, hing2
;;			if i_plot eq 1 then stop
;			print, '-- Median radius:'
;			print, med_r
;			r_stat(i_dis, *) = [[hing1], [med_r], [hing2]]
;		endif
;	
;		i_this_concave = where((x_concave(0,*) gt xrange(0)) and (x_concave(0,*) lt xrange(1)) and (x_concave(1,*) gt xrange(0)) and (x_concave(1,*) lt xrange(1)), j_this_concave)
;		n_concave_arr(i_dis) = j_this_concave
;	
;		i_this_expand = where((x_expand(0,*) gt xrange(0)) and (x_expand(0,*) lt xrange(1)) and (x_expand(1,*) gt xrange(0)) and (x_expand(1,*) lt xrange(1)), j_this_expand)
;		n_expand_arr(i_dis) = j_this_expand
;	
;		i_this_contra = where((x_contra(0,*) gt xrange(0)) and (x_contra(0,*) lt xrange(1)) and (x_contra(1,*) gt xrange(0)) and (x_contra(1,*) lt xrange(1)), j_this_contra)
;		n_contract_arr(i_dis) = j_this_contra
;	endfor
;	;;;;;; plot the variations
;	xbincntrs = 0.5*(xbin_bounds(1:*)+xbin_bounds(0:-2))
;	if strcmp(things(i_plot), 'convex_ratio') then begin
;		;;; convex ratio
;		usersym, 1.*[-1,0,1,0,-1], 1.*[0,1,0,-1,0], thick = l_thick;, /fill ;; diamond
;		ratios = double(n_convex_arr)/double(n_convex_arr+n_concave_arr)*100
;		plot, xbincntrs, ratios, xrange = xrange_plot, xtitle = xtitle, ytitle = 'Convex/All Percentage', position = positions(*,-(i_plot+1)), /noerase, xtickname = xticknames, xstyle = 1, thick = l_thick-0.4
;		oplot, xbincntrs, ratios, psym = 8
;		xyouts, xbincntrs, ratios-0.08*(!y.crange(1)-!y.crange(0)), strcompress(string(n_convex_arr), /remove)+'/'+strcompress(string(n_convex_arr+n_concave_arr), /remove), align = 0.5, charsize = 0.5
;		;makepng, pic_folder+'/ratio_convex'
;	endif
;	;;; expand ratio
;	;plot, xbincntrs, double(n_expand_arr)/double(n_expand_arr+n_contract_arr), xrange = [-6., -16.], xtitle = 'X [RE]', ytitle = 'Ratio: Expand/Convex'
;	;oplot, xbincntrs, double(n_expand_arr)/double(n_expand_arr+n_contract_arr), psym = 4
;	;makepng, pic_folder+'/ratio_expand'
;	;;;; radius of DFB
;	if strcmp(things(i_plot), 'radius') or strcmp(things(i_plot), 'width') then begin
;		if strcmp(things(i_plot), 'radius') then begin
;			string_q = 'Radius'
;			factor = 1.
;		endif
;		if strcmp(things(i_plot), 'width') then begin
;			string_q = 'Width'
;			factor = 2.
;		endif
;		usersym, 0.8*[-1,-1,1,1], 0.8*[1,-1,-1,1], /fill ;; square
;		plot, xbincntrs, factor*r_stat(*,1), xrange = xrange_plot, xtitle = xtitle, yrange = thing.yrange, ystyle = 1, ytitle = 'DFB '+string_q+' [R!dE!n]', position = positions(*,-(i_plot+1)), /noerase, xtickname = xticknames, xstyle = 1, /nodata
;		oplot, xbincntrs, factor*r_stat(*,1), color = 6, thick = l_thick
;		oplot, xbincntrs, factor*r_stat(*,1), psym = 8, color = 6
;		oplot, xbincntrs, factor*r_stat(*,0), color = 2
;		oplot, xbincntrs, factor*r_stat(*,0), psym = 8, color = 2, symsize = 0.7
;		oplot, xbincntrs, factor*r_stat(*,2), color = 2
;		oplot, xbincntrs, factor*r_stat(*,2), psym = 8, color = 2, symsize = 0.7
;		;makepng, pic_folder+'/r_dfb'
;		print, n_convex_arr
;		print, total(n_convex_arr)
;	endif
;	;;; satellite separation in y
;	;plot, xbincntrs, dfront_stat(*,1), xrange = [-6., -16.], xtitle = 'X [RE]', ytitle = 'Probes |dY| [RE]', color = 6
;	;oplot, xbincntrs, dfront_stat(*,1), psym = 4, color = 6
;	;oplot, xbincntrs, dfront_stat(*,0), color = 2
;	;oplot, xbincntrs, dfront_stat(*,0), psym = 4, color = 2
;	;oplot, xbincntrs, dfront_stat(*,2), color = 2
;	;oplot, xbincntrs, dfront_stat(*,2), psym = 4, color = 2
;	;makepng, pic_folder+'/dy_probes'
;	;;;; plot the bin boundaries
;	for i_bound = 0, n_elements(xbin_bounds)-1 do begin 
;		if strcmp(things(i_plot), 'radius') or strcmp(things(i_plot), 'width') then y_bars = [max(!y.crange)-1.2*rate_bar*max(!y.crange), max(!y.crange)]
;		if strcmp(things(i_plot), 'convex_ratio') then y_bars = [0, 1.2*rate_bar*max(!y.crange)]
;		oplot, [xbin_bounds(i_bound), xbin_bounds(i_bound)], y_bars, thick = 0.1
;	endfor
;	;;;; plot label
;	xyouts, !x.crange(0)+thing.x_abc*(!x.crange(1)-!x.crange(0)), !y.crange(0)+thing.y_abc*(!y.crange(1)-!y.crange(0)), abc(i_plot)
;endfor
;pclose
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; check satellite orbit; plot dY vs Y ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;y_range = [-6., 9.]
;ybin_size = 2.5
;dy_range = [0., 3.]
;y_mid = 0.5*(y[0,*]+y[1,*])
;stat_plot, transpose(y_mid), abs(dfront), k_c = 3, bin_range = y_range, binsize = ybin_size, qtt_2_range = dy_range, qtt_range = y_range, qtt_2_title = '|'+cap_delta_letter+'Y!n| [R!dE!n]', qtt_title = 'mean Y [R!dE!n]', kinbin = kinbin, bin_boundaries = binbounds, /no_mean, color_med = 6, color_quar = 2, type_med = 'square', /no_write, bar_thick = l_thick
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; for upward/dawnward dawn/dusk morning/evening paper: whether the shape of the bubble is squeezed circle ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
dawn_loc = 2.
dusk_loc = 2.
dawn_mlt = -1.
dusk_mlt = -1.
ny_crit = 0.0
phi_crit = 0.0
dy_sep_small = 0.5 ;; in RE, to use large seperation events only

dy_convex = y_convex[0,*]-y_convex[1,*]

;;;; use large seperation events only
i_small = where(abs(dy_convex) gt dy_sep_small, n_small)
if n_small gt 1 then begin
	mlt_convex = mlt_convex[*,i_small]
	nfront_convex = nfront_convex[*,i_small]
	r_convex = r_convex[*,i_small]
	dy_convex = dy_convex[*,i_small]
endif

;;dusk = (y_convex[0,*] gt dusk_loc) and (y_convex[1,*] gt dusk_loc)
;;dawn = (y_convex[0,*] lt dawn_loc) and (y_convex[1,*] lt dawn_loc)
dusk = (mlt_convex[0,*] lt dusk_mlt) and (mlt_convex[1,*] lt dusk_mlt)
dawn = (mlt_convex[0,*] gt dawn_mlt) and (mlt_convex[1,*] gt dawn_mlt)
;morning = (nfront_convex[0,*] lt -ny_crit) and (nfront_convex[1,*] lt -ny_crit)
;evening = (nfront_convex[0,*] gt ny_crit) and (nfront_convex[1,*] gt ny_crit)

;;dusk = 0.5*(y_convex[0,*]+y_convex[1,*]) gt dusk_loc
;;dawn = 0.5*(y_convex[0,*]+y_convex[1,*]) lt dawn_loc
;dusk = 0.5*(mlt_convex[0,*]+mlt_convex[1,*]) lt dusk_mlt
;dawn = 0.5*(mlt_convex[0,*]+mlt_convex[1,*]) gt dawn_mlt
;morning = 0.5*(nfront_convex[0,*]+nfront_convex[1,*]) lt -ny_crit
;evening = 0.5*(nfront_convex[0,*]+nfront_convex[1,*]) gt ny_crit
morning = 0.5*(asin(nfront_convex[0,*])+asin(nfront_convex[1,*])) lt -phi_crit
evening = 0.5*(asin(nfront_convex[0,*])+asin(nfront_convex[1,*])) gt phi_crit

i_dusk = where(dusk, n_dusk)
i_dawn = where(dawn, n_dawn)
i_dusk_morning = where(dusk and morning, n_dusk_morning)
i_dusk_evening = where(dusk and evening, n_dusk_evening)
i_dawn_morning = where(dawn and morning, n_dawn_morning)
i_dawn_evening = where(dawn and evening, n_dawn_evening)

;; get the numbers
rstat, r_convex[i_dusk_morning], med_dusk_morning, lq_dusk_morning, uq_dusk_morning
rstat, r_convex[i_dusk_evening], med_dusk_evening, lq_dusk_evening, uq_dusk_evening
rstat, r_convex[i_dawn_morning], med_dawn_morning, lq_dawn_morning, uq_dawn_morning
rstat, r_convex[i_dawn_evening], med_dawn_evening, lq_dawn_evening, uq_dawn_evening

print, 'Of all convex events:'
print, 'Median R: '+strcompress(string(median(r_convex, /even)))+' RE.'
print, 'Dusk Median Y seperation: '+strcompress(string(median(abs(dy_convex[i_dusk]), /even)))
print, 'Dusk, Morning: '+strcompress(string(n_dusk_morning))+' events.'
if n_dusk_morning gt 0 then begin
	print, 'Median R: '+strcompress(string(med_dusk_morning))+' RE. (Median Y separation: '+strcompress(string(median(abs(dy_convex[i_dusk_morning]), /even)))+')'
	print, 'LQ R: '+strcompress(string(lq_dusk_morning))+' RE.'
	print, 'UQ R: '+strcompress(string(uq_dusk_morning))+' RE.'
endif
print, 'Dusk, Evening: '+strcompress(string(n_dusk_evening))+' events.'
if n_dusk_evening gt 0 then begin
	print, 'Median R: '+strcompress(string(med_dusk_evening))+' RE. (Median Y separation: '+strcompress(string(median(abs(dy_convex[i_dusk_evening]), /even)))+')'
	print, 'LQ R: '+strcompress(string(lq_dusk_evening))+' RE.'
	print, 'UQ R: '+strcompress(string(uq_dusk_evening))+' RE.'
endif
print, 'Dawn Median Y seperation: '+strcompress(string(median(abs(dy_convex[i_dawn]), /even)))
print, 'Dawn, Morning: '+strcompress(string(n_dawn_morning))+' events.'
if n_dawn_morning gt 0 then begin
	print, 'Median R: '+strcompress(string(med_dawn_morning))+' RE. (Median Y separation: '+strcompress(string(median(abs(dy_convex[i_dawn_morning]), /even)))+')'
	print, 'LQ R: '+strcompress(string(lq_dawn_morning))+' RE.'
	print, 'UQ R: '+strcompress(string(uq_dawn_morning))+' RE.'
endif
print, 'Dawn, Evening: '+strcompress(string(n_dawn_evening))+' events.'
if n_dawn_evening gt 0 then begin
	print, 'Median R: '+strcompress(string(med_dawn_evening))+' RE. (Median Y separation: '+strcompress(string(median(abs(dy_convex[i_dawn_evening]), /even)))+')'
	print, 'LQ R: '+strcompress(string(lq_dawn_evening))+' RE.'
	print, 'UQ R: '+strcompress(string(uq_dawn_evening))+' RE.'
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

stop
end
