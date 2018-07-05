#include <ros/ros.h>
#include <visualization_msgs/Marker.h>
#include "std_msgs/String.h"
#include "continuum_manipulator_visualization.cpp"
#include "continuum_manipulator_visualization.hpp"
#include <Eigen3/Eigen/Dense>
#include <Eigen3/Eigen/Core>

namespace muscle_ns{
Muscle::Muscle(ros::NodeHandle& nodeHandle) : nodeHandle(nodeHandle)
{
  muscle_pub=nodeHandle.advertise<visualization_msgs::Marker>("visualization_marker", 10);
  cylinder_pub=nodeHandle.advertise<visualization_msgs::Marker>("visualization_marker_cylinder", 10);

  m_length_sub = nodeHandle.subscribe("/laengen", 1 , &Muscle::lengthCallback, this);
  initializeMatrices();



}
} /*namespace*/

int main( int argc, char** argv )
{
  ros::init(argc, argv, "continuum_manipulator_visualization");
  ros::NodeHandle nodeHandle("~");

muscle_ns:: Muscle muscle0(nodeHandle);
  ros::Rate r(1000000);

  while(ros::ok())  {
    ros::spinOnce();
  }

  return 0;
}
