#include <ros/ros.h>
#include <probo_msgs/distance.h>
#include <math.h>
#include <Eigen3/Eigen/Dense>
#include <probo_msgs/angle_curr.h>

using namespace std;



class Segment
{
  public:
    // constructor
    Segment(float r_ss, float * ptr_origin, ros::NodeHandle& nodehandle){

      
    }


  private:
    void distanceCallback(const probo_msgs::distance & msg_l);
    float *getTip();
    ros::Subscriber sub_l;
    ros::Publisher pub_pos;
    ros::Publisher pub_angle;
    ros::NodeHandle nodehandle;
    float r_ss;
    float *ptr_origin;
};
