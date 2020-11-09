

lam2 = Inp.lambda2;

if isfield(Inp, 'lambda3')
    lam3 = Inp.lambda3;
else
    lam3 = lam2;
end



Res = Big_Results.Res;
Mill_Res = Big_Results.Mill_Res;
Ell_Res = Big_Results.ell_Res;
rot_Res = Big_Results.rot_Res;
squi_Res = Big_Results.squi_Res;
squi_rot_Res = Big_Results.squi_rot_Res;

cell_Res = {Res, Mill_Res, Ell_Res, rot_Res, squi_Res, squi_rot_Res};

str_intro = 'Lin_sel_';
dyn_str = 'Dyn_';
mill_str = 'Mill_';
ell_str = 'ell_';
rot_str = 'rot_';
squi_str = 'squi_';
rot_squi_str = [rot_str, squi_str];
str_middle = {dyn_str, mill_str, ell_str, rot_str, squi_str, rot_squi_str};

cmax = 0;
for ii = 1:length(str_middle)
    
    switch str_middle{ii}
        case mill_str
            if mod(lam3, 1) ~= 0
                lam_str = num2str(lam3);
                lam_str = split(lam_str, '.');
                str_beg = ['lambda3_', lam_str{1},'_',lam_str{2}, '_'];
            else
                str_beg = ['lambda3_', num2str(lam3), '_'];
            end
            
        otherwise
            if mod(lam2, 1) ~= 0
                lam_str = num2str(lam2);
                lam_str = split(lam_str, '.');
                str_beg = ['lambda2_', lam_str{1},'_',lam_str{2}, '_'];
            else
                str_beg = ['lambda2_', num2str(lam2), '_'];
            end
    end
    
    lin_fname = [str_beg, str_intro, str_middle{ii}, ];
    lin_fname_RF = [str_beg, 'RF_', str_intro, str_middle{ii}];
    
    linearFeatureMaps(cell_Res{ii}, 0, 13, lin_fname);
    linearFeatureMaps(cell_Res{ii}, 1, 13, lin_fname_RF);
    
    cmax_temp = max(cell_Res{ii}.zThetaSel(:));
    if cmax_temp > cmax
        cmax = cmax_temp;
    end
    
end

OriPref_comp(Big_Results, Inp, 17, 1);
ori_fname = [str_beg, 'Ori_Pref']
saveas(99, ori_fname, 'epsc');

OriSel_comp(Big_Results, Inp, cmax)
ori_sel_fname = [str_beg, 'Ori_sel'];
saveas(999, ori_sel_fname, 'epsc');

comp_Results(Big_Results, Inp);
saveas(2122, [str_beg, 'RF_examples'], 'epsc')




