#include <ros/ros.h>
#include <Eigen3/Eigen/Dense>
#include <string.h>
#include "continuum_manipulator_visualization.hpp"
#include <visualization_msgs/Marker.h>
#include "visualization_msgs/MarkerArray.h"


namespace muscle_ns{

  void Muscle::lengthCallback(const std_msgs::Float32MultiArray& m_length)  {

    for (int i = 0; i < 3; i++) {   //reading length of all 3 muscles
      muscle_length[i] = m_length.data[i];
      l_bar += muscle_length[i];
    }

    l_bar = l_bar/ 3.0;

    theta = (2.0/3.0)*sqrt(muscle_length[0]*muscle_length[0] + muscle_length[1]*muscle_length[1] + muscle_length[2]*muscle_length[2] -
            (muscle_length[0]*muscle_length[1] + muscle_length[0]*muscle_length[2] + muscle_length[1]*muscle_length[2]))/r_mm;
    phi = arctan2(sqrt(3)*(muscle_length[2] - muscle_length[1]), (muscle_length[1] + muscle_length[2] -2*muscle_length[0]));

    //draw base cylinder
    visualization_msgs::Marker cylinder;
    cylinder.type = visualization_msgs::Marker::CYLINDER;
    cylinder.header.frame_id  = "/my_frame";
    cylinder.header.stamp = ros::Time::now();
    cylinder.ns =  "Muscle";

    cylinder.id = 1 ;
    cylinder.scale.x =.10;
    cylinder.scale.y =.10;
    cylinder.scale.z =.02;

    cylinder.color.a =1.0;
    cylinder.color.r =1;
    cylinder.color.g =0;
    cylinder.color.b = 0;

    cylinder.pose.position.x = 0;
    cylinder.pose.position.y = 0;
    cylinder.pose.position.z = -0.010;

    cylinder.pose.orientation.x = 0;
    cylinder.pose.orientation.y = 0;
    cylinder.pose.orientation.z = 0;
    cylinder.pose.orientation.w = 0;

    cylinder.action = visualization_msgs::Marker::ADD;
    cylinder_pub.publish(cylinder);

    for (int i = 0; i < 3; i++) { // draw all 3 muscles

      Eigen::Matrix4f H_def = H_base[i];

      visualization_msgs::Marker sphere_list;
      sphere_list.type = visualization_msgs::Marker::SPHERE_LIST;
      sphere_list.header.frame_id = "/my_frame";
      sphere_list.header.stamp = ros::Time::now();
      sphere_list.ns = "spheres";
      sphere_list.pose.orientation.w = 1.0;
      sphere_list.id = 2 + i;

      sphere_list.color.a = 1.0;
      sphere_list.color.r = 1;
      sphere_list.color.g = 1;
      sphere_list.color.b = 1;

      sphere_list.scale.x = 0.035;
      sphere_list.scale.y = 0.035;
      sphere_list.scale.z = 0.035;

      geometry_msgs::Point p;


      for (int j = 0; j < resolution; j++)  {

        createTransform(H_muscle[i][j], theta/(resolution-1)*j, phi, muscle_length[i]/(resolution-1)*j);
        Eigen::Vector4f pos = H_def*H_muscle[i][j]*ZERO_ONE;
        sphere_list.header.stamp = ros::Time::now();
        sphere_list.ns = "Muscle";

        p.x = pos(0);
        p.y = pos(1);
        p.z = pos(2);

        sphere_list.action=visualization_msgs::Marker::ADD;

        if ( j < resolution - 14)  {

          sphere_list.points.push_back(p);
          }

        pos_tip[3*i + 0] = pos(0);
        pos_tip[3*i + 1] = pos(1);
        pos_tip[3*i + 2] = pos(2);
        }

    muscle_pub.publish(sphere_list);
    }
  // calculate position of the tip-midpoint
  pos_tip[0] = (pos_tip[0] + pos_tip[3] +pos_tip[6])/3;
  pos_tip[1] = (pos_tip[1] + pos_tip[4] +pos_tip[7])/3;
  pos_tip[2] = (pos_tip[2] + pos_tip[5] +pos_tip[8])/3;

 //draw tip cylinder
  cylinder.pose.position.x = pos_tip[0];
  cylinder.pose.position.y = pos_tip[1];
  cylinder.pose.position.z = pos_tip[2];

  Eigen::Vector4f ori;
  toQuaternion(ori, H_muscle[1][resolution -1 ]);
  cylinder.pose.orientation.x = ori(0);
  cylinder.pose.orientation.y = ori(1);
  cylinder.pose.orientation.z = ori(2);
  cylinder.pose.orientation.w = ori(3);
  cylinder.action = visualization_msgs::Marker::ADD;
  muscle_pub.publish(cylinder);
  }


  void Muscle::createTransform(Eigen::Matrix4f& H, float theta, float phi, float l)  {
  if(theta == 0)  { // avoid divide by zero

    H(0,0) = 1;
    H(0,1) = 0;
    H(0,2) = 0;
    H(0,3) = 0;

    H(1,0) = 0;
    H(1,1) = 1;
    H(1,2) = 0;
    H(1,3) = 0;

    H(2,0) = 0;
    H(2,1) = 0;
    H(2,2) = 1;
    H(2,3) = l;

    H(3,0) = 0;
    H(3,1) = 0;
    H(3,2) = 0;
    H(3,3) = 1;
  }
  else  {

    H(0,0) = cos(phi)*cos(phi)*(cos(theta)-1) + 1;
    H(0,1) = cos(phi)*sin(phi)*(cos(theta)-1);
    H(0,2) = cos(phi)*sin(theta);
    H(0,3) = l*cos(phi)*(1-cos(theta))/theta;

    H(1,0) = sin(phi)*cos(phi)*(cos(theta)-1);
    H(1,1) = cos(phi)*cos(phi)*(1-cos(theta)) + cos(theta);
    H(1,2) = sin(phi)*sin(theta);
    H(1,3) = l*sin(phi)*(1-cos(theta))/theta;

    H(2,0) = - cos(phi) * sin(theta);
    H(2,1) = - sin(phi) * sin(theta);
    H(2,2) = cos(theta);
    H(2,3) = l*sin(theta)/theta;

    H(3,0) = 0;
    H(3,1) = 0;
    H(3,2) = 0;
    H(3,3) = 1;

  }
  }


void Muscle::initializeMatrices() {

  ZERO_ONE = Eigen::Vector4f();
  ZERO_ONE << 0,0,0,1;
  X_AXIS   << 1,0,0,0;
  Y_AXIS   << 0,1,0,0;
  Z_AXIS   << 0,0,1,0;

  for (int i = 0; i < 3; i ++ )  {

    H_base[i] = Eigen::MatrixXf(4,4);
    H_base[i] << 1, 0, 0, r_mm * cos(2*M_PI/3*i),
                  0, 1, 0, r_mm * sin(2*M_PI/3*i),
                  0, 0, 1, 0,
                  0, 0, 0, 1;
  }

  theta = 0;
  phi = 0;
  for (int i = 0; i < 3; i++) {

    for (int j = 0; j < resolution; j++) {

      H_muscle[i][j] = Eigen::MatrixXf(4,4);
      createTransform(H_muscle[i][j], theta/(resolution-1)*j, phi, muscle_length[i]/(resolution-1)*j);
    }
  }
  }

void Muscle::toQuaternion(Eigen::Vector4f& quat, Eigen::Matrix4f H) {

  quat(3) = sqrt(1.0 + H(0,0) + H(1,1) + H(2,2)) / 2.0;
  double w4 = (4.0 * quat(3));
  quat(0) = (H(2,1) - H(1,2)) / w4 ;
  quat(1) = (H(0,2) - H(2,0)) / w4 ;
  quat(2) = (H(1,0) - H(0,1)) / w4 ;
}

float Muscle::arctan2(float y, float x)  {

  float ratio = y/x;

  if(x > 0)  {
	  return atan(ratio);
  }
  else if(x < 0 && y > 0)  {
 	  return atan(ratio) + M_PI;
  }
  else if(x < 0 && y == 0)  {
    return M_PI;
  }
  else if(x < 0 && y < 0)  {
	  return atan(ratio) - M_PI;
  }
  else if(x == 0 && y > 0)  {
	  return M_PI/2.0;
  }
  else return -M_PI/2.0;
  }
} /* namespace*/
