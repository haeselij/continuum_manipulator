classdef segment
     properties (Constant)
%        r_ss;
       q_0 = 0.3014;
       w_g = [0; 0; 9.81;0];
       m = 0.256;
       A = 0.0183*0.0183*pi;  
       k=  0.0111;
       pressures_size = 10000;   
   end
   
   properties
        q_1;
        q_2;
        q_3;
        
        v_1; % = q_dot_1
        v_2;
        v_3;
        
        p_1;
        p_2;
        p_3;
        
        w_H_1h = zeros(4,4); % coordinatetransform matrix head to worldframe
        
        d_w_H_1h__d_q_1 = zeros(4,4); %first derrivative of H_1h after seg.q_1
        d_w_H_1h__d_q_2 = zeros(4,4);
        d_w_H_1h__d_q_3 = zeros(4,4);
% 
%         d2_w_H_1h__d_q_1_d_q_1 = zeros(4,1); %second derrivative of H_1h after seg.q_1
%         d2_w_H_1h__d_q_2_d_q_2 = zeros(4,1);
%         d2_w_H_1h__d_q_3_d_q_3 = zeros(4,1);
%         d2_w_H_1h__d_q_2_d_q_3 = zeros(4,1);
%         d2_w_H_1h__d_q_1_d_q_2 = zeros(4,1);
%         d2_w_H_1h__d_q_1_d_q_3 = zeros(4,1);
%         
        w_H_1h_2dot = zeros(4,1);
        
        q_container = zeros(3, 10000);  
        p_container = zeros(3, 10000);
   end
   
   methods
       function seg  = segment(p_bag)
           if(size(p_bag) == size(seg.p_container))
            seg.p_container= p_bag;
           else
            B = 'Size error'
           end
       end
      
       function GetKineticEnergy
       end
       function GetGeneralizedForce
       end
       function GetPotentialEnergy
  
   end
   
   end
end

