%% Cordinate transformation to get the tip position
j = 1;
tfTime = transform_time(1);
w_H_b_t = struct('tf',zeros(4,4));
coordinates_measured = zeros(3,transform_size-2);
quaternion = zeros(1,4);

% transform vicon measurement eliminate the markers
base_offset = 0.017;
tip_offset = 0.0075;
w_H_b = [1 0 0 0;
            0 1 0 0;
            0 0 1 -base_offset;
            0 0 0 1];
w_H_t = [1 0 0 0;
            0 1 0 0;
            0 0 1 -tip_offset;
            0 0 0 1];

for i = 2:transform_size - 1
    if (canTransform(bag,'segBase','tip',tfTime))
      tf_base_to_tip = getTransform(transform_select,'segBase', 'tip',tfTime);

      quaternion(1) = tf_base_to_tip.Transform.Rotation.X;
      quaternion(2) = tf_base_to_tip.Transform.Rotation.Y;
      quaternion(3) = tf_base_to_tip.Transform.Rotation.Z;
      quaternion(4) = tf_base_to_tip.Transform.Rotation.W;
      
      w_H_b_t(j).tf(:,:) =   quat2tform(quaternion);
      w_H_b_t(j).tf(1,4) =   w_H_b_t(j).tf(1,4) + tf_base_to_tip.Transform.Translation.X;
      w_H_b_t(j).tf(2,4) =   w_H_b_t(j).tf(2,4) + tf_base_to_tip.Transform.Translation.Y;
      w_H_b_t(j).tf(3,4) =   w_H_b_t(j).tf(3,4) + tf_base_to_tip.Transform.Translation.Z;

      w_H_b_t(j).tf = w_H_b * w_H_b_t(j).tf* w_H_t;
      
      coordinates_measured(1,j) = w_H_b_t(j).tf(1,4);
      coordinates_measured(2,j) = w_H_b_t(j).tf(2,4);
      coordinates_measured(3,j) = w_H_b_t(j).tf(3,4);
      
      j = j+1; 
    end
     
    tfTime = transform_time(i+1);
end

coordinate_size = j-1;
transform_time = transform_time - transform_time(1);



