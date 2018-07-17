
#include <math.h>
#include <ros.h>
#include <probo_msgs/distance.h>


ros::NodeHandle nh;
probo_msgs::distance sensor;
ros::Publisher laengen_pub("/laengen", &sensor);

void setup() {
  nh.initNode();
  nh.advertise(laengen_pub);
}



void loop() {
  sensor.length[0] = analogRead(A2);
  sensor.length[1] = analogRead(A1);
  sensor.length[2] = analogRead(A0);

  laengen_pub.publish(&sensor);

  nh.spinOnce();
}
