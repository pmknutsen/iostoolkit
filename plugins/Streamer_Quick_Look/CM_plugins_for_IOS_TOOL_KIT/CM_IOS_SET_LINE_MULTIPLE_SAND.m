% the coordinates of the line on the former image is Xi_S,Yi_S
% CX_T=CX_B+box_lim(3); CY_T=CY_B+box_lim(1);
%

function [STR_LINE_SET]=CM_IOS_SET_LINE_MULTIPLE_SAND(AVG_im,thickness)
STR_LINE_SET=[];

mastercounter=0;
for Int_Count=1:1:50
    mastercounter=mastercounter+1;
    [NEXT_OR_NOT]=Check_for_next();
    if (isempty(NEXT_OR_NOT)<1) % case where you did not press cancel
        if (findstr(NEXT_OR_NOT{1},'delete'))
            mastercounter=mastercounter-1;% move the counter one line back
            STR_LINE_SET(mastercounter)=[];% deletes the previous line
            [NEXT_OR_NOT]=Check_for_next();% allows to draw a new line instead if wanted
        end
        if (isempty(NEXT_OR_NOT)<1) % case where you did not press cancel
            if (findstr(NEXT_OR_NOT{1},'NEW LINE'))
                [Xi_S,Yi_S,box_lim,AVG_im,to_delete]=CM_IOS_SET_LINE_SAND(AVG_im,thickness,134,mastercounter);
                to_delete
                if (to_delete~=1)
                    
                    STR_LINE_SET(mastercounter).Xi_S=Xi_S;
                    STR_LINE_SET(mastercounter).Yi_S=Yi_S;
                    STR_LINE_SET(mastercounter).box_lim=box_lim;
                    STR_LINE_SET(mastercounter).AVG_im=AVG_im;
                else 
                     mastercounter=mastercounter-1;
                end
            end
        end
        
    else
        return
    end
    CM_IOS_PLOT_LINES_SAND(STR_LINE_SET,0);
end
end

function [NEXT_OR_NOT]=Check_for_next()
prompt = {'Press cancel to stop or type delete to erase previous line'};
dlg_title = 'Input';
num_lines = 1;
def = {'NEW LINE'};
NEXT_OR_NOT = inputdlg(prompt,dlg_title,num_lines,def);
end

