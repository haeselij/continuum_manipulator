#include <ros/ros.h>
#include <probo_msgs/distance.h>
#include <math.h>
#include <Eigen3/Eigen/Dense>
#include <probo_msgs/distance.h>
#pragma once

using namespace std;
using namespace Eigen;




class Segment
{
  public:
    // constructor
    Segment(ros::NodeHandle& nodehandle, int pos);


  private:
    void lengthCallback(const probo_msgs::distance & l);
    float sinX(float alpha);
    float arctan2(float y, float x);
    void trafo();

    float r_ss;
    float length[3];
    Vector3f tip_pos;
    Vector3f tip_pos_0;
    Vector3f r_offset;
    float theta;
    float phi;
    int pos_in_trunk;
    float theta_prev;
    float phi_prev;

    ros::Subscriber sub_l;
    ros::Publisher pub_pos;
    ros::Publisher pub_angle;
    ros::NodeHandle nodehandle_;
};
