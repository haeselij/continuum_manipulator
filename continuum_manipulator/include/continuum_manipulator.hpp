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
    Segment(ros::NodeHandle& nodehandle, int pos, Matrix4f & H_prev, string str_angle, string str_head, bool end);
    // method to get transformation matrix of previous Segment
    Matrix4f & getH_prev();


  private:
    void lengthCallback(const probo_msgs::distance & l);
    //float sinX(float alpha);
    float arctan2(float y, float x);
    void build_trafo();


    float r_ss;
    float r_muscle;
    float length[3];
    float l_bar;
    Vector3f tip_pos_0;
    //Vector3f r_offset;
    float theta;
    float phi;
    int pos_in_trunk;
    //float theta_prev;
    //float phi_prev;

    Matrix4f  &H_prev_;
    Matrix4f H;
    Vector4f V_if;

    ros::Subscriber sub_l;
    ros::Publisher pub_pos;
    ros::Publisher pub_angle;
    ros::NodeHandle nodehandle_;
    bool end;
    void calculate_q();
    float q[3];
};
