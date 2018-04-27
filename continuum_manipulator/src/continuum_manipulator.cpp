#include <ros/ros.h>
#include <Eigen3/Eigen/Dense>
#include <Eigen3/Eigen/Core>
#include <probo_msgs/distance.h>
#include <math.h>
#include <continuum_manipulator.hpp>
#include <probo_msgs/tip_pos.h>
#include <probo_msgs/angle_curr.h>

using namespace Eigen;
using namespace std;

Segment::Segment(ros::NodeHandle& nodehandle, int pos, Matrix4f& H_prev, string str_angle, string str_head, bool end) : H_prev_(H_prev)
{
  r_ss = 10;
  l_bar = 0;
  r_muscle = 26;
  nodehandle_ = nodehandle;
  pos_in_trunk = pos;
  // Translation from the interface
  V_if(0)=0;
  V_if(1)=0;
  V_if(2)=10;
  V_if(0)=0;

  /*
  // Check if segment is the first one. If this is the case there is no H_prev needed.
  if(pos_in_trunk == 0)
  {
    Matrix4f identity;
    identity.setIdentity(4, 4);

    for(int i = 0; i < 4; i++)
    {
      for(int j = 0; j < 4; j++)
      {
        H_prev_(i,j) = identity.coeff(i,j);
      }
    }
  //cout << H_prev_(0,0) << endl;
  }
  else
  {
  H_prev_ = H_prev;

  cout<<"H_prev_="<<H_prev_<< endl;
}*/

  sub_l = nodehandle.subscribe("/laengen", 1, &Segment::lengthCallback, this);
  pub_pos = nodehandle.advertise<probo_msgs::tip_pos>(str_head, 1);
  pub_angle = nodehandle.advertise<probo_msgs::angle_curr>(str_angle, 1);
}

Matrix4f & Segment::getH_prev()
{
  //return reference of H
  return H;
}

/*
float Segment::sinX(float alpha)
{
  if (sin(alpha) == 0.0 && abs(alpha) < 3.1415)
  {
    return 1;
  }

  else return sin(alpha)/alpha;
}
*/

// arctan2 function
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

// calculation of transformation matrix H
void Segment::build_trafo()
{
  // limit analysis for theta -> 0 in order to avoid devided by 0 error
  if(theta == 0)
  {
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
    H(2,3) = l_bar;

    H(3,0) = 0;
    H(3,1) = 0;
    H(3,2) = 0;
    H(3,3) = 1;
  }
  else
  {
    H(0,0) = cos(phi)*cos(phi)*(cos(theta)-1) + 1;
    H(0,1) = cos(phi)*sin(phi)*(cos(theta)-1);
    H(0,2) = cos(phi)*sin(theta);
    H(0,3) = l_bar*cos(phi)*(1-cos(theta))/theta;

    H(1,0) = sin(phi)*cos(phi)*(cos(theta)-1);
    H(1,1) = cos(phi)*cos(phi)*(1-cos(theta)) + cos(theta);
    H(1,2) = sin(phi)*sin(theta);
    H(1,3) = l_bar*sin(phi)*(1-cos(theta))/theta;

    H(2,0) = - cos(phi) * sin(theta);
    H(2,1) = - sin(phi) * sin(theta);
    H(2,2) = cos(theta);
    H(2,3) = l_bar*sin(theta)/theta;

    H(3,0) = 0;
    H(3,1) = 0;
    H(3,2) = 0;
    H(3,3) = 1;
  }

  if ( end == false)
  {
  Matrix4f H_if;

  H_if.setIdentity(4, 4);
  H_if(2,3) = 13.9;

  H = H*H_if;
  }
}

void Segment::lengthCallback(const probo_msgs::distance & l)
{
  probo_msgs::angle_curr ac;
  probo_msgs::tip_pos tp;
  l_bar = 0;
  for(int i = 0; i < 3; i++)
  {
    length[i] = l.length[3*pos_in_trunk + i] - 15 + 300;
    l_bar += length[i];
  }

  l_bar = l_bar/ 3.0;



  theta = (2.0/3.0)*sqrt(length[0]*length[0] + length[1]*length[1] + length[2]*length[2] -
          (length[0]*length[1] + length[0]*length[2] + length[1]*length[2]))/r_ss;
  phi = arctan2(sqrt(3)*(length[2] - length[1]), (length[1] + length[2] -2*length[0]));


  //tip_pos(0) = -l_bar*sinX(theta/2.0)*cos((M_PI - theta)/2.0)*sin(phi - M_PI/2.0);
  //tip_pos(1) = l_bar*sinX(theta/2.0)*cos((M_PI - theta)/2.0)*cos(phi - M_PI/2.0);
  //tip_pos(2) = l_bar*sinX(theta/2.0)*sin((M_PI - theta)/2.0);

  /*if (pos_in_trunk == 0)
  {
      theta_prev = 0;
      phi_prev = 0;
  }
  else
  {
    theta_prev = (2.0/3.0) * sqrt(l.length[3*(pos_in_trunk-1)]*l.length[3*(pos_in_trunk-1) ] + l.length[3*(pos_in_trunk-1) + 1]*l.length[3*(pos_in_trunk-1) + 1] + l.length[3*(pos_in_trunk-1) + 2]*l.length[3*(pos_in_trunk-1) + 2] -
                  (l.length[3*(pos_in_trunk-1)]*l.length[3*(pos_in_trunk-1) + 1] + l.length[(pos_in_trunk-1)]*l.length[(pos_in_trunk-1) + 2] + l.length[(pos_in_trunk-1) + 1]*l.length[(pos_in_trunk-1) + 2]))/r_ss;
    phi_prev = arctan2(sqrt(3)*(l.length[3] - l.length[2]),(l.length[2] + l.length[3] -2*l.length[1]));

    //r_offset (0) = -l_bar*sinX(theta_prev/2.0)*cos((M_PI - theta_prev)/2.0)*sin(phi_prev - M_PI/2.0);
    //(r_offset (1) = l_bar*sinX(theta_prev/2.0)*cos((M_PI - theta_prev)/2.0)*cos(phi_prev - M_PI/2.0);
    //r_offset (2) = l_bar*sinX(theta_prev/2.0)*sin((M_PI - theta_prev)/2.0);
  }*/

  build_trafo();
  calculate_q();
  /*if(pos_in_trunk==0)
  {
    cout << 0 << endl;
  }
  else cout << 1 << endl;*/
 //cout << "H="<< H<< endl;
 //cout << "H_prev=" << H_prev_ << endl;

  H = (H_prev_)*H;

 // adding the interface offset


  for(int i = 0; i<3; i++)
  {
    tp.tip_pos[3*pos_in_trunk+i] = H.coeff(i, 3);
  }

//cout << "theta=" << theta << endl;
//cout << "phi=" << phi << endl;
  ac.theta_curr[0] = theta;
  ac.phi_curr[0] = phi;

     pub_angle.publish(ac);
     pub_pos.publish(tp);


}

void Segment::calculate_q()
{
  // calculate the radii (curvature)
  float r_curv[3];
  r_curv[0] =  length[0] / theta;
  r_curv[1] =  length[1] / theta;
  r_curv[2] =  length[2] / theta;
  float l_im[3];
  l_im[0] = ( length[0] + length[1]) / 2;
  l_im[1] = ( length[1] + length [2]) / 2;
  l_im[2] = ( length[0] + length[2]) / 2;

  for (int i = 0; i < 3; i++)
  {
    q[i] = l_im[i] * ( 1 - r_muscle / r_curv[i] * cos( 2 * M_PI / 3 * i - phi +M_PI/6) );
    cout << pos_in_trunk << "Length " << i <<" = "<< q[i] << endl;
  }
}
