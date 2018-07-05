#include <ros/ros.h>
#include <Eigen3/Eigen/Dense>
#include <Eigen3/Eigen/Core>
#include <probo_msgs/distance.h>
#include <math.h>
#include "continuum_manipulator.hpp"
#include "continuum_manipulator.cpp"



int main(int argc, char** argv)
{
  ros::init(argc, argv, "sub");
  ros::NodeHandle nodeHandle("~");
  Matrix4f baseFrame;
  baseFrame.setIdentity();
	Segment s0(nodeHandle, 0, baseFrame, "/winkel_seg_0", "/head_pos_0",0);
	Segment s1(nodeHandle, 1, s0.getH_prev(), "/winkel_seg_1", "/tip_pos_1",1);

  ros::spin();
  return 0;
}
