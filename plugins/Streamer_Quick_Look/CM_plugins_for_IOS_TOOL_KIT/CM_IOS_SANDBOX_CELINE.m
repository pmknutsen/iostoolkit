function [T]=CM_IOS_SANDBOX_CELINE(T)
 thickness=10;
 if ~isfield(T,'AA')
     [CORRECTED_BASELINE]=CM_SUBTRACT_BACKGROUND_SAND(T.mBaseline(:,:,2),200,400,0);
     [T.AA]=CM_IOS_SET_LINE_MULTIPLE_SAND(CORRECTED_BASELINE,thickness);
      CM_IOS_PLOT_LINES_SAND(T.AA,0)
 end
 
 [T]=CM_IOS_CALC_DIAM_SAND (T,thickness); % calculates the diam and initialize at the first instance
end

