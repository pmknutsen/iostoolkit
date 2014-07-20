% the coordinates of the line on the former image is Xi_S,Yi_S
% CX_T=CX_B+box_lim(3); CY_T=CY_B+box_lim(1);
%

function [STR_LINE_SET]=CM_IOS_SET_LINE_MULTIPLE_SAND(AVG_im)
STR_LINE_SET=[];
for Int_Count=1:1:50
    [NEXT_OR_NOT]=Check_for_next();
    if (isempty(NEXT_OR_NOT)<1) % case where you did not press cancel   
        [Xi_S,Yi_S,box_lim,AVG_im]=CM_IOS_SET_LINE_SAND(AVG_im,Int_Count);
        STR_LINE_SET(Int_Count).Xi_S=Xi_S;
        STR_LINE_SET(Int_Count).Yi_S=Yi_S;
        STR_LINE_SET(Int_Count).box_lim=box_lim;
        STR_LINE_SET(Int_Count).AVG_im=AVG_im;
    else
        return
    end
end
end

function [NEXT_OR_NOT]=Check_for_next()
prompt = {'Press cancel to stop'};
dlg_title = 'Input';
num_lines = 1;
def = {'NEW LINE'};
NEXT_OR_NOT = inputdlg(prompt,dlg_title,num_lines,def);
end

