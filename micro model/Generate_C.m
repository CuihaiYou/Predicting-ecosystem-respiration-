function y=Generate_C(C_op_0,CMin_0,CMax_0)
L_P =length(C_op_0); %length of the parameter vector
while(true)
  C_new=C_op_0+(rand(1,L_P)-0.5).*(CMax_0-CMin_0)/(L_P+1);
   Logic=C_new>CMin_0&C_new<CMax_0;
         if sum(Logic)==L_P
             break;
         end
end
   y=C_new;