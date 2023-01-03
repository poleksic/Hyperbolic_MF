function new_x = h_vector_check_adjust(x)
    % check if x belongs to hyperbolic space (hyperboloid)
%     if abs(h_inner(x,x) + 1) > 0.001
%         fprintf('<x,x>=%.10f\n',h_inner(x,x));
%         warning('not on hyperboloid, will readjust')
%     end     
    new_x = h_vector_adjust(x);
end