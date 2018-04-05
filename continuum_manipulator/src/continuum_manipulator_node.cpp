#include <ros/ros.h>
#include <Eigen3/Eigen/Dense>
#include <Eigen3/Eigen/Core>
#include <probo_msgs/distance.h>
#include <math.h>
#include "continuum_manipulator.hpp"
#include "continuum_manipulator.cpp"


int main(int argc, char** argv)
{
	float world_frame[3]={0,0,0};
	float *test_ptr=NULL;
  ros::init(argc, argv, "sub");
  ros::NodeHandle nodeHandle("~");
	Segment s1(10, test_ptr, nodeHandle);

  ros::spin();
  return 0;
}
