#pragma once

#include <string.h>
#include <ros/ros.h>
#include <math.h>
#include <Eigen3/Eigen/Dense>
#include <Eigen3/Eigen/Core>
#include <probo_msgs/angle_ref.h>
#include <probo_msgs/distance.h>
#include <std_msgs/Int32.h>
#include <std_msgs/Float32MultiArray.h>

namespace muscle_ns{

class Muscle {
  public:
    Muscle (ros::NodeHandle& nodeHandle);

  private:
    void lengthCallback(const std_msgs::Float32MultiArray& m_length);
    float arctan2(float y, float x);

    void toQuaternion(Eigen::Vector4f& quat, Eigen::Matrix4f H);
    void initializeMatrices();
    void createTransform(Eigen::Matrix4f& H, float theta, float phi, float l) ;

    ros::NodeHandle nodeHandle;
    ros::Publisher muscle_pub;
    ros::Publisher cylinder_pub;
    ros::Subscriber muscle_sub;
    ros::Subscriber m_length_sub;

    visualization_msgs::Marker marker_muscle;
    const int resolution = 400;
    float muscle_length[3];
    const float muscle_dia = 0.045;
    const float r_mm = 0.03;
    float pos_tip[9];
    float theta;
    float phi;
    float l_bar;

    Eigen::Matrix4f H_base[3];
    Eigen::Matrix4f H_orient;
    Eigen::Matrix4f H_muscle[4][400];

    Eigen::Vector4f ZERO_ONE;
    Eigen::Vector4f X_AXIS;
    Eigen::Vector4f Y_AXIS;
    Eigen::Vector4f Z_AXIS;
};
}/*namespace*/
