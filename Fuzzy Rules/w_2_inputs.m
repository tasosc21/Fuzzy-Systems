%-----------------------------------------------------------------
%       function w=w_2_inputs(A,A_bar,B,B_bar)
%       This function provides the w coefficient of a 2-input
%       fuzzy rules with Rc and max-min
%-----------------------------------------------------------------

function w=w_2_inputs(A,A_bar,B,B_bar)
w_input1=max(min(A,A_bar))
w_input2=max(min(B,B_bar))
w=min(w_input1,w_input2)




