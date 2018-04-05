#include <ros/ros.h>
#include <Eigen3/Eigen/Dense>
#include <Eigen3/Eigen/Core>
#include <probo_msgs/distance.h>
#include <math.h>

using namespace std;
using namespace Eigen;
float sinX(float alpha)
{
	if (sin(alpha)==0.0 && abs(alpha)<3.1415)
	{
	return 1;
	}
	else return sin(alpha)/alpha;

}

float arctan2(float y, float x)
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

void Callback(const probo_msgs::distance msg)
{
  float r_ss=10.0;
  float l1=msg.length[1]-15 +300;
  float l2=msg.length[0]-15 + 300;
  float l3=msg.length[2]-15 + 300;
  float l_bar=(l1 + l2 + l3 )/3;
  float theta= (2.0/3.0) * sqrt(l1*l1 +l2*l2 + l3*l3 - (l1*l2 + l1*l3 + l2*l3))/r_ss;
  float phi=arctan2(sqrt(3)*(l3-l2),l2 + l3 -2*l1);
  float x= - l_bar *sinX(theta/2.0)*cos((M_PI-theta)/2.0)*sin(phi - M_PI/2.0);
  float y= l_bar *sinX(theta/2.0)*cos((M_PI-theta)/2.0)*cos(phi - M_PI/2.0);
  float z=l_bar *sinX(theta/2.0)*sin((M_PI-theta)/2.0);

  ROS_INFO_STREAM("x=" << x<< endl);
  ROS_INFO_STREAM("y=" << y<< endl);
	ROS_INFO_STREAM("z=" << z<< endl);
	ROS_INFO_STREAM("theta=" << theta<< endl);
	ROS_INFO_STREAM("phi=" << phi<< endl);


}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "sub");
  ros::NodeHandle nodeHandle("~");
  ros::Subscriber subscriber = nodeHandle.subscribe("/laengen", 10, Callback);
  ros::spin();



}
