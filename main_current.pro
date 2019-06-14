pro main_current
;;; save the data for main_current_plot
computer = 'I:'

thm_init
@folders

seconds_check_out = 15.

list_suf = ''
;list_suf = '_earthward_df_binxl_nocare'
;list_suf = '_earthward_df_bavexl_nocare'
;list_suf = '_earthward_df_bdiravexl_nocare'
listname = 'dfb_list_lead_tail'+list_suf+'.txt'

events = load_list(listname, folder = list_folder)
;;; for test
;events = events(*,100:110)

;pos_list = event_status_instant(events, 'pos', time_length=3, datafolder=pos_folder)
;pos = pos_list(2:4,*)/6371.
;
;h = scale_height(events, b_folder = fgs_folder, p_folder = Pall_folder, Blobe = Bl)

;; DF thickness
thick_vec = df_thick(events, v_use = 'efs', care = 'no', seconds_check_in = 15., seconds_check_out = seconds_check_out, b_type = 'fgl', fgs_folder = fgs_folder, fgl_folder = fgl_folder, e_folder = efs_folder, maxgap = 10, b_in_list = b_in_list, b_out_list = b_out_list, /vec, t_in = t_max_arr, t_out = t_min_arr)

dataout_simple, save_folder+'/df_thick_vec'+list_suf, thick_vec
stop

b_in = b_in_list(2:4, *)
b_out = b_out_list(2:4, *)

if strcmp_or(list_suf, ['', '_earthward_df_bavexl_nocare', '_earthward_df_bdiravexl_nocare']) then begin
	b_ave_list = event_status_instant(events, 'fgl', vtrange = [t_min_arr(1,*),t_max_arr(1,*)], datafolder=fgl_folder)
	b_ave = b_ave_list(2:4, *)
endif
if strcmp_or(list_suf, ['', '_earthward_df_bdiravexl_nocare']) then begin
	b_dir_ave_list = event_status_instant(events, 'fgl', vtrange = [t_min_arr(1,*),t_max_arr(1,*)], datafolder=fgl_folder, /dir)
	b_dir_ave = b_dir_ave_list(2:4, *)
endif

;; L direction
dfb_mva, events, bfolder = fgl_folder, l_list = l, m_list = m_, n_list = n_, lambda_list = lambda, datatype = 'fgl', ratio_lm = 3 ;; for L to be good
;;;; re-adjust L direction: take l direction to be the direction that BL always jump
Bin_l = dotp_long(b_in, l)
Bout_l = dotp_long(b_out, l)
i_reverse = where(Bin_l-Bout_l lt 0, j_reverse)
if j_reverse gt 0 then l(*, i_reverse)=-l(*, i_reverse)

;;;; save
dataout_simple, save_folder+'/pos'+list_suf, pos
dataout_simple, save_folder+'/h'+list_suf, h
dataout_simple, save_folder+'/B_in'+list_suf, b_in
dataout_simple, save_folder+'/B_out'+list_suf, b_out
if strcmp_or(list_suf, ['', '_earthward_df_bavexl_nocare', '_earthward_df_bdiravexl_nocare']) then dataout_simple, save_folder+'/B_ave'+list_suf, b_ave
if strcmp_or(list_suf, ['', '_earthward_df_bdiravexl_nocare']) then dataout_simple, save_folder+'/B_dirave'+list_suf, b_dir_ave
dataout_simple, save_folder+'/l'+list_suf, l
dataout_simple, save_folder+'/lambda'+list_suf, lambda
dataout_simple, save_folder+'/df_thick_vec'+list_suf, thick_vec

stop

end
