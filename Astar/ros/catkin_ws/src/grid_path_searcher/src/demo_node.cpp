#include <fstream>
#include <iostream>
#include <math.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>
#include <ros/console.h>
#include <ros/ros.h>
#include <sensor_msgs/PointCloud2.h>

#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>

#include "Astar_searcher.h"
#include "JPS_searcher.h"
#include "backward.hpp"
#include <memory>

using namespace std;
using namespace Eigen;

namespace backward {
backward::SignalHandling sh;
}

// simulation param from launch file
double _resolution, _inv_resolution, _cloud_margin;
double _x_size, _y_size, _z_size;

// useful global variables
bool _has_map = false;

Vector3d _start_pt;
Vector3d _map_lower, _map_upper;
int _max_x_id, _max_y_id, _max_z_id;

// ros related
ros::Subscriber _map_sub, _pts_sub;
ros::Publisher _grid_path_vis_pub, _visited_nodes_vis_pub, _grid_map_vis_pub;

std::shared_ptr<AstarPathFinder> _astar_path_finder (new AstarPathFinder);
std::shared_ptr<JPSPathFinder> _jps_path_finder (new JPSPathFinder);

void rcvWaypointsCallback(const nav_msgs::Path &wp);
void rcvPointCloudCallBack(const sensor_msgs::PointCloud2 &pointcloud_map);

void visGridPath(vector<Vector3d> nodes, bool is_use_jps);
void visVisitedNode(vector<Vector3d> nodes);
void innerCycle(std::shared_ptr<AstarPathFinder> path_finder, const Vector3d start_pt, const Vector3d target_pt);
void innerCycle(std::shared_ptr<JPSPathFinder> path_finder, const Vector3d start_pt, const Vector3d target_pt);
void pathFinding(const Vector3d start_pt, const Vector3d target_pt);

void rcvWaypointsCallback(const nav_msgs::Path &wp) {
  if (wp.poses[0].pose.position.z < 0.0 || _has_map == false)
    return;

  Vector3d target_pt;
  target_pt << wp.poses[0].pose.position.x, wp.poses[0].pose.position.y,
      wp.poses[0].pose.position.z;

  ROS_INFO("[node] receive the planning target");
  pathFinding(_start_pt, target_pt);
}

void rcvPointCloudCallBack(const sensor_msgs::PointCloud2 &pointcloud_map) {
  if (_has_map)
    return;

  pcl::PointCloud<pcl::PointXYZ> cloud;
  pcl::PointCloud<pcl::PointXYZ> cloud_vis;
  sensor_msgs::PointCloud2 map_vis;

  pcl::fromROSMsg(pointcloud_map, cloud);

  if ((int)cloud.points.size() == 0)
    return;

  pcl::PointXYZ pt;
  for (int idx = 0; idx < (int)cloud.points.size(); idx++) {
    pt = cloud.points[idx];

    // set obstalces into grid map for path planning
    _astar_path_finder->setObs(pt.x, pt.y, pt.z);
    _jps_path_finder->setObs(pt.x, pt.y, pt.z);

    // for visualize only
    Vector3d cor_round =
        _astar_path_finder->coordRounding(Vector3d(pt.x, pt.y, pt.z));
    pt.x = cor_round(0);
    pt.y = cor_round(1);
    pt.z = cor_round(2);
    cloud_vis.points.push_back(pt);
  }

  cloud_vis.width = cloud_vis.points.size();
  cloud_vis.height = 1;
  cloud_vis.is_dense = true;

  pcl::toROSMsg(cloud_vis, map_vis);

  map_vis.header.frame_id = "/world";
  _grid_map_vis_pub.publish(map_vis);

  _has_map = true;
}

void innerCycle(std::shared_ptr<AstarPathFinder> path_finder, const Vector3d start_pt, const Vector3d target_pt) {
  // Call A*/JPS to search for a path
  path_finder->AstarGraphSearch(start_pt, target_pt);

  // Retrieve the path
  auto grid_path = path_finder->getPath();
  auto visited_nodes = path_finder->getVisitedNodes();

  // Visualize the result
  visGridPath(grid_path, false);
  visVisitedNode(visited_nodes);

  // Reset map for next call
  path_finder->resetUsedGrids();
}

void innerCycle(std::shared_ptr<JPSPathFinder> path_finder, const Vector3d start_pt, const Vector3d target_pt) {
  // Call A*/JPS to search for a path
  path_finder->JPSGraphSearch(start_pt, target_pt);

  // Retrieve the path
  auto grid_path = path_finder->getPath();
  auto visited_nodes = path_finder->getVisitedNodes();

  // Visualize the result
  visGridPath(grid_path, true);
  visVisitedNode(visited_nodes);

  // Reset map for next call
  path_finder->resetUsedGrids();
}

void pathFinding(const Vector3d start_pt, const Vector3d target_pt) {
  // --------------------- 1. dijkstra ---------------------
  _astar_path_finder->set_h_type(4);
  _astar_path_finder->set_tie_breaker(0);
  innerCycle(_astar_path_finder, start_pt, target_pt);

  // --------------------- 2. A* with Manhatton heuristic ---
  _astar_path_finder->set_h_type(1);
  _astar_path_finder->set_tie_breaker(0);
  innerCycle(_astar_path_finder, start_pt, target_pt);

  // --------------------- 3. A* with Euclidean heuristic ---
  _astar_path_finder->set_h_type(0);
  _astar_path_finder->set_tie_breaker(0);
  innerCycle(_astar_path_finder, start_pt, target_pt);

  // --------------------- 4. A* with L(inf) norm heuristic ---
  _astar_path_finder->set_h_type(2);
  _astar_path_finder->set_tie_breaker(0);
  innerCycle(_astar_path_finder, start_pt, target_pt);

  // --------------------- 5. A* with diagonal heuristic ---
  _astar_path_finder->set_h_type(3);
  _astar_path_finder->set_tie_breaker(0);
  innerCycle(_astar_path_finder, start_pt, target_pt);

  // --------------------- 6. A* with diagonal & tie breaker ---
  _astar_path_finder->set_h_type(3);
  _astar_path_finder->set_tie_breaker(1);
  innerCycle(_astar_path_finder, start_pt, target_pt);

//_use_jps = 0 -> Do not use JPS
//_use_jps = 1 -> Use JPS
// you just need to change the #define value of _use_jps
#define _use_jps 1
#if _use_jps
  {
    // --------------------- 7. JPS with diagonal ---------------------
    _jps_path_finder->set_h_type(3);
    _jps_path_finder->set_tie_breaker(0);
    innerCycle(_jps_path_finder, start_pt, target_pt);

    // --------------------- 7. JPS with diagonal & tie breaker ---
    _jps_path_finder->set_h_type(3);
    _jps_path_finder->set_tie_breaker(1);
    innerCycle(_jps_path_finder, start_pt, target_pt);
  }
#endif
}

int main(int argc, char **argv) {
  ros::init(argc, argv, "demo_node");
  ros::NodeHandle nh("~");

  _map_sub = nh.subscribe("map", 1, rcvPointCloudCallBack);
  _pts_sub = nh.subscribe("waypoints", 1, rcvWaypointsCallback);

  _grid_map_vis_pub = nh.advertise<sensor_msgs::PointCloud2>("grid_map_vis", 1);
  _grid_path_vis_pub =
      nh.advertise<visualization_msgs::Marker>("grid_path_vis", 1);
  _visited_nodes_vis_pub =
      nh.advertise<visualization_msgs::Marker>("visited_nodes_vis", 1);

  nh.param("map/cloud_margin", _cloud_margin, 0.0);
  nh.param("map/resolution", _resolution, 0.2);

  nh.param("map/x_size", _x_size, 50.0);
  nh.param("map/y_size", _y_size, 50.0);
  nh.param("map/z_size", _z_size, 5.0);

  nh.param("planning/start_x", _start_pt(0), 0.0);
  nh.param("planning/start_y", _start_pt(1), 0.0);
  nh.param("planning/start_z", _start_pt(2), 0.0);

  _map_lower << -_x_size / 2.0, -_y_size / 2.0, 0.0;
  _map_upper << +_x_size / 2.0, +_y_size / 2.0, _z_size;

  _inv_resolution = 1.0 / _resolution;

  _max_x_id = (int)(_x_size * _inv_resolution);
  _max_y_id = (int)(_y_size * _inv_resolution);
  _max_z_id = (int)(_z_size * _inv_resolution);

  // _astar_path_finder = new AstarPathFinder();
  _astar_path_finder->initGridMap(_resolution, _map_lower, _map_upper,
                                  _max_x_id, _max_y_id, _max_z_id);

  // _jps_path_finder = new JPSPathFinder();
  _jps_path_finder->initGridMap(_resolution, _map_lower, _map_upper, _max_x_id,
                                _max_y_id, _max_z_id);

  ros::Rate rate(100);
  bool status = ros::ok();
  while (status) {
    ros::spinOnce();
    status = ros::ok();
    rate.sleep();
  }

  // delete _astar_path_finder;
  // delete _jps_path_finder;
  return 0;
}

void visGridPath(vector<Vector3d> nodes, bool is_use_jps) {
  visualization_msgs::Marker node_vis;
  node_vis.header.frame_id = "world";
  node_vis.header.stamp = ros::Time::now();

  if (is_use_jps)
    node_vis.ns = "demo_node/jps_path";
  else
    node_vis.ns = "demo_node/astar_path";

  node_vis.type = visualization_msgs::Marker::CUBE_LIST;
  node_vis.action = visualization_msgs::Marker::ADD;
  node_vis.id = 0;

  node_vis.pose.orientation.x = 0.0;
  node_vis.pose.orientation.y = 0.0;
  node_vis.pose.orientation.z = 0.0;
  node_vis.pose.orientation.w = 1.0;

  if (is_use_jps) {
    node_vis.color.a = 1.0;
    node_vis.color.r = 0.0;
    node_vis.color.g = 1.0;
    node_vis.color.b = 1.0;
  } else {
    node_vis.color.a = 1.0;
    node_vis.color.r = 0.0;
    node_vis.color.g = 0.1;
    node_vis.color.b = 0.0;
  }

  node_vis.scale.x = _resolution;
  node_vis.scale.y = _resolution;
  node_vis.scale.z = _resolution;

  geometry_msgs::Point pt;
  for (int i = 0; i < int(nodes.size()); i++) {
    Vector3d coord = nodes[i];
    pt.x = coord(0);
    pt.y = coord(1);
    pt.z = coord(2);

    node_vis.points.push_back(pt);
    // we are filling blanks btwn two straight nodes for JPS
    if (i < int(nodes.size() - 1)) {
      geometry_msgs::Point pt1;
      coord = nodes[i+1];
      pt1.x = coord(0);
      pt1.y = coord(1);
      pt1.z = coord(2);
      // at least a one node blank >= _resolution * 2 for (x/y/z) direction
      // if (fabs(pt1.x - pt.x) > _resolution * 1.9) {
      //   int num = ceil(fabs((pt1.x - pt.x) / _resolution));
      //   double delta = (pt1.x - pt.x) / num;
      //   for (int j = 0; j < num; ++j) {
      //     geometry_msgs::Point pt2;
      //     pt1.x = pt.x + delta * j;
      //     pt1.y = pt.y;
      //     pt1.z = pt.z;
      //     node_vis.points.push_back(pt1);
      //   }
      // } else if (fabs(pt1.y - pt.y) > _resolution * 1.9) {
      //   int num = ceil(fabs((pt1.y - pt.y) / _resolution));
      //   double delta = (pt1.y - pt.y) / num;
      //   for (int j = 0; j < num; ++j) {
      //     geometry_msgs::Point pt2;
      //     pt1.x = pt.x;
      //     pt1.y = pt.y + delta * j;
      //     pt1.z = pt.z;
      //     node_vis.points.push_back(pt1);
      //   }
      // } else if (fabs(pt1.z - pt.z) > _resolution * 1.9) {
      //   int num = ceil(fabs((pt1.z - pt.z) / _resolution));
      //   double delta = (pt1.z - pt.z) / num;
      //   for (int j = 0; j < num; ++j) {
      //     geometry_msgs::Point pt2;
      //     pt1.x = pt.x;
      //     pt1.y = pt.y;
      //     pt1.z = pt.z + delta * j;
      //     node_vis.points.push_back(pt1);
      //   }
      // }
    }
  }

  _grid_path_vis_pub.publish(node_vis);
}

void visVisitedNode(vector<Vector3d> nodes) {
  visualization_msgs::Marker node_vis;
  node_vis.header.frame_id = "world";
  node_vis.header.stamp = ros::Time::now();
  node_vis.ns = "demo_node/expanded_nodes";
  node_vis.type = visualization_msgs::Marker::CUBE_LIST;
  node_vis.action = visualization_msgs::Marker::ADD;
  node_vis.id = 0;

  node_vis.pose.orientation.x = 0.0;
  node_vis.pose.orientation.y = 0.0;
  node_vis.pose.orientation.z = 0.0;
  node_vis.pose.orientation.w = 1.0;
  node_vis.color.a = 0.5;
  node_vis.color.r = 0.0;
  node_vis.color.g = 0.0;
  node_vis.color.b = 1.0;

  node_vis.scale.x = _resolution;
  node_vis.scale.y = _resolution;
  node_vis.scale.z = _resolution;

  geometry_msgs::Point pt;
  for (int i = 0; i < int(nodes.size()); i++) {
    Vector3d coord = nodes[i];
    pt.x = coord(0);
    pt.y = coord(1);
    pt.z = coord(2);

    node_vis.points.push_back(pt);
  }

  _visited_nodes_vis_pub.publish(node_vis);
}