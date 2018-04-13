#include <ros/ros.h>
#include <Eigen3/Eigen/Dense>
#include <Eigen3/Eigen/Core>
#include <probo_msgs/distance.h>
#include <math.h>
#include <continuum_manipulator.hpp>
#include <probo_msgs/tip_pos.h>
#include <probo_msgs/angle_curr.h>

using namespace Eigen;

Segment::Segment(ros::NodeHandle& nodehandle, int pos)
{
  r_ss = 10;
  nodehandle_ = nodehandle;
  pos_in_trunk = pos;

  sub_l = nodehandle.subscribe("/laengen", 1, &Segment::lengthCallback, this);
  pub_pos = nodehandle.advertise<probo_msgs::tip_pos>("/tip_position", 1);
  pub_angle = nodehandle.advertise<probo_msgs::angle_curr>("/angle_curr", 1);
}


float Segment::sinX(float alpha)
{
  if (sin(alpha) == 0.0 && abs(alpha) < 3.1415)
  {
    return 1;
  }

  else return sin(alpha)/alpha;
}

float Segment::arctan2(float y, float x)
{
  float ratio = y/x;

  if(x > 0)
  {
	  return atan(ratio);
  }
  else if(x < 0 && y > 0)
  {
 	  return atan(ratio) + M_PI;
  }
  else if(x < 0 && y == 0)
  {
    return M_PI;
  }
  else if(x < 0 && y < 0)
  {
	  return atan(ratio) - M_PI;
  }
  else if(x == 0 && y > 0)
  {
	  return M_PI/2.0;
  }
  else return -M_PI/2.0;
}

void Segment::trafo(){
  Matrix3f R_0;
  R_0(0,0) = cos(phi_prev) * cos(theta_prev) * cos (phi_prev) + sin(phi_prev) * sin(phi_prev);
  R_0(0,1) = cos(phi_prev) * cos(theta_prev) * sin(phi_prev) - sin(phi_prev) * cos(phi_prev);
  R_0(0,2) = - cos(phi_prev) * sin(theta_prev);

  R_0(1,0) = sin(phi_prev) * cos(theta_prev) * cos (phi_prev) - cos(phi_prev) * sin(phi_prev);
  R_0(1,1) = sin(phi_prev) * cos(theta_prev) * sin(phi_prev) + cos(phi_prev) * cos(phi_prev);
  R_0(1,2) = - sin(phi_prev) * sin(theta_prev);

  R_0(2,0) = sin(theta_prev) * cos(phi_prev);
  R_0(2,1) = sin(theta_prev) * sin(phi_prev);
  R_0(2,2) = cos(theta_prev);

  tip_pos_0 = R_0*tip_pos + r_offset;
}

void Segment::lengthCallback(const probo_msgs::distance & l)
{
  float l_bar = 0;

  for(int i = 0; i < 3; i++)
  {
    length[i] = l.length[3*pos_in_trunk + i] - 15 + 300;
    l_bar += length[i];
  }

  l_bar /= 3.0;

  theta = (2.0/3.0) * sqrt(length[0]*length[0] + length[1]*length[1] + length[2]*length[2] -
                (length[0]*length[1] + length[0]*length[2] + length[1]*length[2]))/r_ss;
  phi = arctan2(sqrt(3)*(length[2] - length[1]),(length[1] + length[2] -2*length[0]));

  tip_pos(0) = -l_bar*sinX(theta/2.0)*cos((M_PI - theta)/2.0)*sin(phi - M_PI/2.0);
  tip_pos(1) = l_bar*sinX(theta/2.0)*cos((M_PI - theta)/2.0)*cos(phi - M_PI/2.0);
  tip_pos(2) = l_bar*sinX(theta/2.0)*sin((M_PI - theta)/2.0);

  if (pos_in_trunk == 0)
  {
      theta_prev = 0;
      phi_prev = 0;
  }
  else
  {
    theta_prev = (2.0/3.0) * sqrt(l.length[3*(pos_in_trunk-1)]*l.length[3*(pos_in_trunk-1) ] + l.length[3*(pos_in_trunk-1) + 1]*l.length[3*(pos_in_trunk-1) + 1] + l.length[3*(pos_in_trunk-1) + 2]*l.length[3*(pos_in_trunk-1) + 2] -
                  (l.length[3*(pos_in_trunk-1)]*l.length[3*(pos_in_trunk-1) + 1] + l.length[(pos_in_trunk-1)]*l.length[(pos_in_trunk-1) + 2] + l.length[(pos_in_trunk-1) + 1]*l.length[(pos_in_trunk-1) + 2]))/r_ss;
    phi_prev = arctan2(sqrt(3)*(l.length[3] - l.length[2]),(l.length[2] + l.length[3] -2*l.length[1]));

    r_offset (0) = -l_bar*sinX(theta_prev/2.0)*cos((M_PI - theta_prev)/2.0)*sin(phi_prev - M_PI/2.0);
    r_offset (1) = l_bar*sinX(theta_prev/2.0)*cos((M_PI - theta_prev)/2.0)*cos(phi_prev - M_PI/2.0);
    r_offset (2) = l_bar*sinX(theta_prev/2.0)*sin((M_PI - theta_prev)/2.0);
  }
  trafo();

  probo_msgs::angle_curr ac;
  probo_msgs::tip_pos tp;

  for(int i = 0; i<3; i++)
  {
    tp.tip_pos[i] = tip_pos_0(i);
  }

  ac.theta_curr[pos_in_trunk] = theta;
  ac.phi_curr[pos_in_trunk] = phi;

  pub_angle.publish(ac);
  pub_pos.publish(tp);
}
