 #include <ros/ros.h>
#include <probo_msgs/pressure.h>
#include <math.h>
#include <Eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

Vector3f pressurefunction(float u_phi,float u_theta)
{
  u_phi = -u_phi;
  Matrix2f P1_inv;
  Matrix2f P2_inv;
  Matrix2f P3_inv;
  MatrixXf x = MatrixXf(2,1);
  MatrixXf b = MatrixXf(2,1);

  Matrix3f p;
  Matrix2f P1_;
  Matrix2f P2_;
  Matrix2f P3_;

  float r=1.0;  //normal 35 mm
  p(0,0)=0;
  p(0,1)= -r*sqrt(3.0)/2.0;
  p(0,2)= r*sqrt(3.0)/2.0;
  p(1,0)= -r;
  p(1,1)= r/2.0;
  p(1,2)= r/2.0;

  P1_ << p(0,1), p(0,2), //stimmt
      p(1,1), p(1,2);

  P1_inv << P1_(1,1), -P1_(0,1),
         -P1_(1,0), P1_(0,0);

 const float P1_det= P1_(0,0)*P1_(1,1)-P1_(0,1)*P1_(1,0) ;

 P1_inv=P1_inv*1/P1_det; //stimmt

  if(fabs(P1_det) < 1e-4){
   cout << "error" << P1_det << endl;
  }

  P2_ << p(0,0), p(0,2), //stimmt
     p(1,0), p(1,2);

  P2_inv << P2_(1,1), -P2_(0,1), //stimmt
        -P2_(1,0), P2_(0,0);

  const float P2_det= P2_(0,0)*P2_(1,1)-P2_(0,1)*P2_(1,0);
  P2_inv = P2_inv*1/P2_det;

  P3_<< p(0,0), p(0,1), //stimmt
         p(1,0), p(1,1);

  P3_inv << P3_(1,1), -P3_(0,1),
         -P3_(1,0), P3_(0,0);
  const float P3_det= P3_(0,0)*P3_(1,1)-P3_(0,1)*P3_(1,0);

  P3_inv = P3_inv *1/P3_det; //stimmt

  //solving Ax=b
  b <<  u_phi,
        u_theta;

  x = P1_inv*b;
  //  x=P_.ldlt().solve(b);
 Vector3f p_123;
 p_123 <<  0,
            x(0,0),
            x(1,0);

  //case p2 = 0
  // all P_123 have to be positiv, because there are no negative pressures
  if (p_123.minCoeff() < 0)
  {
    x =P2_inv*b;
    p_123 <<  x(0,0),
              0,
              x(1,0);

    //case p3 = 0
    if (p_123.minCoeff() < 0)
    {
      x=P3_inv*b;
      p_123 <<  x(0,0),
                x(1,0),
                0;
    }
  }

return p_123;

}


int main(int argc, char** argv)
{
  ros::init(argc, argv, "publisher");
  ros::NodeHandle nodeHandle("~");
  ros::Publisher publisher = nodeHandle.advertise<probo_msgs::pressure>("/druecke",1);
  ros::Rate loop_rate(100);
  float i = 0.0 ;
  float time = 500.0;
  float max_press = 0.8;
  int muscle_0 = 1;
  int muscle_1 = 6;
  int muscle_2 = 8;

  while (ros::ok())
  {
    probo_msgs::pressure msg_pressure;
    for (int j = 0; i < 12; i++)
    {
      msg_pressure.press[j] = 0;
    }
    // pressure skript single points one chamber


    if (i < time )
    {
      msg_pressure.press[muscle_0] = 0.25*max_press;

      publisher.publish(msg_pressure);
    }


     else if (i < 2 * time)
    {
      msg_pressure.press[muscle_0] = 0.5*max_press;
      publisher.publish(msg_pressure);
    }


      else if (i < 3 * time)
    {
      msg_pressure.press[muscle_0] = 0.75 * max_press;
      publisher.publish(msg_pressure);
    }


    else if (i < 4* time)
    {
      msg_pressure.press[muscle_0] = max_press;
      publisher.publish(msg_pressure);
    }

    // pressure skript points 2 chambers
    else if (i < 5*time)
    {
      publisher.publish(msg_pressure);
    }
      else if (i < 6 * time )
      {
        msg_pressure.press[muscle_0] = 0.25 * max_press;
        msg_pressure.press[muscle_1] = 0.25 * max_press;
        publisher.publish(msg_pressure);
      }

      else if (i < 7 * time)
      {
        msg_pressure.press[muscle_0] = 0.5*max_press;
        msg_pressure.press[muscle_1] = 0.5 * max_press;
        publisher.publish(msg_pressure);
    }
    else if (i < 8 * time)
    {
      msg_pressure.press[muscle_0] = 0.75 * max_press;
        msg_pressure.press[muscle_1] = 0.75 * max_press;
        publisher.publish(msg_pressure);
    }
    else if (i < 9 * time)
    {
      msg_pressure.press[muscle_0] =  max_press;
      msg_pressure.press[muscle_1] =  max_press;
      publisher.publish(msg_pressure);
    }


    // pressure skript points asynchron

    else if (i < 10*time)
    {
      publisher.publish(msg_pressure);
    }
    else if (i < 11 * time )
    {
      msg_pressure.press[muscle_0] = 0.25 * max_press;
      msg_pressure.press[muscle_1] = 0.5 * max_press;
      publisher.publish(msg_pressure);
    }

    else if (i < 12 * time)
    {
      msg_pressure.press[muscle_0] = 0.1*max_press;
      msg_pressure.press[muscle_1] = 0.5 * max_press;
      publisher.publish(msg_pressure);
    }
    else if (i < 13 * time)
    {
      msg_pressure.press[muscle_0] = 0.5 * max_press;
        msg_pressure.press[muscle_1] = 0.75 * max_press;
        publisher.publish(msg_pressure);
    }
    else if (i < 14 * time)
    {
      msg_pressure.press[muscle_0] =  0.3*max_press;
      msg_pressure.press[muscle_1] =  0.4*max_press;
      publisher.publish(msg_pressure);
    }

    // continous pressure change
    else if (i < 15 * time)
    {
      msg_pressure.press[muscle_0] = max_press / time*(i - 14*time);
      msg_pressure.press[muscle_1] = 0.0;
      publisher.publish(msg_pressure);
    }

    else if (i < 16 * time)
    {
      msg_pressure.press[muscle_0] = max_press*(1 - 1/time*(i - 15*time));
      msg_pressure.press[muscle_1] = 0.0;
      publisher.publish(msg_pressure);
    }

    else if ( i < 17*time)
    {
      msg_pressure.press[muscle_0] = max_press / time*(i - 16*time);
      msg_pressure.press[muscle_2] =  max_press / time*(i - 16*time);
      publisher.publish(msg_pressure);
    }

    else if ( i < 18*time)
    {
      msg_pressure.press[muscle_0] = max_press*(1 - 1/time*(i - 17.0*time));
      msg_pressure.press[muscle_2] = max_press*(1 - 1/time*(i - 17.0*time));
      publisher.publish(msg_pressure);
    }

  //circular movement
    else if (i < 20 * time)
    {
    float p_phi = 0.7*cos(M_PI/time*i);
    float p_theta = 0.7*sin(M_PI/time*i);
    Vector3f p_123 = pressurefunction(p_phi, p_theta);

    msg_pressure.press[muscle_0] = p_123[0];
    msg_pressure.press[muscle_1] = p_123[1];
    msg_pressure.press[muscle_2] = p_123[2];
    publisher.publish(msg_pressure);
    }



    //"left right movement"
    else if (i < 21*time)
    {
      publisher.publish(msg_pressure);
    }
    else if (i < 22*time)
   {
      msg_pressure.press[muscle_0] = max_press/time*(i-21*time);
      publisher.publish(msg_pressure);
    }
    else if (i < time * 23)
    {
      msg_pressure.press[muscle_0] = max_press-  max_press/time*(i - 22*time);
      publisher.publish(msg_pressure);
    }
    else if (i < time * 24)
    {
      msg_pressure.press[muscle_1] = max_press/time*(i - time*23);
      msg_pressure.press[muscle_2] = max_press/time*(i - time*23);
      //msg_pressure.press[4] = max/time*(i - time*23);
      //msg_pressure.press[5] = 1.5/time*(i - time*23);
      publisher.publish(msg_pressure);
    }
    else if ( i < time * 25)
      {
      msg_pressure.press[muscle_1] = max_press - max_press/time*(i - time*24);
      msg_pressure.press[muscle_2] = max_press - max_press/time*(i - time*24);
      //msg_pressure.press[4] = 1.5 - 1.5/time*(i - time*24);
      //msg_pressure.press[5] = 1.5 - 1.5/time*(i - time*24);
      publisher.publish(msg_pressure);
    }

    // "up down movement"

  /* if (i < time)
    {
      msg_pressure.press[muscle_0] = 1.0/2.0 *max_press/time*i;
      msg_pressure.press[1] = sqrt(3.0)/2.0*max_press/time*i;

      msg_pressure.press[3] = 1.0/2.0 *max_press/time*i;
      msg_pressure.press[4] = sqrt(3.0)/2.0*max_press/time*i;

      publisher.publish(msg_pressure);

    }

    else if (i < time*2)
    {
      msg_pressure.press[muscle_0] = 1.0/2.0 *max_press - 1.0/2.0 *max_press/time*(i-time);
      msg_pressure.press[1] = sqrt(3.0)/2.0*( 1.0 - 1.0/time*(i-time));

      msg_pressure.press[3] = 1.0/2.0 *max_press*( 1.0 - 1.0/time*(i-time));
      msg_pressure.press[4] = sqrt(3.0)/2.0*max_press*( 1.0 - 1.0/time*(i-time));

      publisher.publish(msg_pressure);

    }

    else if (i < time*3)
    {
      msg_pressure.press[muscle_0] = 1.0/2.0 *max_press/time*(i - time*2.0);
      msg_pressure.press[2] = sqrt(3.0)/2.0*max_press/time*(i - time*2.0);

      msg_pressure.press[3] = 1.0/2.0 *max_press/time*(i - time*2.0);
      msg_pressure.press[5] = sqrt(3.0)/2.0*max_press/time*(i - time*2.0);

      publisher.publish(msg_pressure);
    }

    else if ( i < time * 4)
      {
      msg_pressure.press[muscle_0] = 1.0/2.0 *max_press*( 1.0 - 1.0/time*(i - 3.0*time));
      msg_pressure.press[2] = sqrt(3.0)/2.0*( 1.0 - 1.0/time*(i - 3.0*time));

      msg_pressure.press[3] = 1.0/2.0 *max_press*( 1.0 - 1.0/time*(i - 3.0*time));
      msg_pressure.press[5] = sqrt(3.0)/2.0*max_press*( 1.0 - 1.0/time*(i - 3.0*time));

      publisher.publish(msg_pressure);
    }*/

    // fixed position1
    /*msg_pressure.press[muscle_0] = 1;
    msg_pressure.press[4] = 1;
    msg_pressure.press[5] = 1;
    publisher.publish(msg_pressure);*/

    //fixed position2
    /*msg_pressure.press[muscle_0] = max_press/2;
    msg_pressure.press[2] = max_press/4;
    msg_pressure.press[4] = max_press/2;
    msg_pressure.press[3] = max_press/5;
    publisher.publish(msg_pressure);*/

    ros::spinOnce();
    loop_rate.sleep();
    i = i +1.0;
  }

  return 0;
}
