#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  // Leftmost lane = 0. Start in lane 1.
  int lane = 1;

  // Reference velocity in m/s
  double ref_speed = 0;

  // Reference speed increment in m/s
  double ref_speed_delta = 0.1;

  // Min. safe distance in m, from other cars
  int safe_dist = 32;

  // Speed limit in m/s
  auto speed_limit = 22.1;

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,
    &map_waypoints_dx,&map_waypoints_dy,&lane,&ref_speed,&ref_speed_delta,
    &safe_dist,&speed_limit]
    (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

            // Define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
            int prev_size = previous_path_x.size();
            if (prev_size > 0) {
              car_s = end_path_s;
            }
            // Check other cars around us using sensor fusion data,
            // to determine if lane / speed should be changed
            // Variables to check for closeness of cars
            // in current lane, left lane, and right lane
            bool too_close_in_current_lane = false;
            bool too_close_in_left_lane = false;
            bool too_close_in_right_lane = false;
            for (int i = 0; i < sensor_fusion.size(); i++) {
              // Get the other car's location and speed
              float d = sensor_fusion[i][6];
              double vx = sensor_fusion[i][3];
              double vy = sensor_fusion[i][4];
              double other_car_speed = distance(0, 0, vx, vy);
              // Projected s value for the other car - where it will be in the future
              double other_car_s = sensor_fusion[i][5];
              other_car_s += prev_size * 0.02 * other_car_speed;
              double s_diff = other_car_s - car_s;

              // Check if the other car is at an unsafe distance around us
              if (s_diff < safe_dist) {
                // Check if the other car is ahead of us in our lane at an unsafe distance
                if (s_diff > 0 && d > 4 * lane && d < 4 * (lane + 1)) {
                  too_close_in_current_lane = true;
                }
                // Check if the other car is behind us at an unsafe distance
                else if (s_diff > -safe_dist) {
                  // Check if the other car is in left lane
                  if (d > 4 * (lane - 1) && d < 4 * lane) {
                    too_close_in_left_lane = true;
                  }
                  // Check if the other car is in right lane
                  else if (d > 4 * (lane + 1) && d < 4 * (lane + 2)) {
                    too_close_in_right_lane = true;
                  }
                }
              }
            }

            // Variable to check if we need to reduce speed
            bool reduce_speed = false;
            // Check if we need to change lane
            if (too_close_in_current_lane) {
              // Need to change lane
              // Check if changing lane to left lane is feasible and safe.
              if (!too_close_in_left_lane && lane > 0) {
                lane--;
              }
              // Left lane is not feasible or is unsafe
              // Check if changing lane to right lane is feasible and safe.
              else if (!too_close_in_right_lane && lane < 2) {
                lane++;
              }
              // Cannot change lane. Just reduce speed (later)
              else {
                reduce_speed = true;
              }
            }

            // Create a list of widely spaced (x, y) waypoints
            vector<double> pts_x;
            vector<double> pts_y;

            // Reference x, y, yaw states
            double ref_x = car_x;
            double ref_y = car_y;
            double ref_yaw = deg2rad(car_yaw);

            if (prev_size < 2) {
              // Use two points that make the path tangent to the car
              pts_x.push_back(car_x - cos(car_yaw));
              pts_x.push_back(car_x);
              pts_y.push_back(car_y - sin(car_yaw));
              pts_y.push_back(car_y);
            }
            // Use the previous path's endpoint as starting reference
            else {
              // Redefine reference state as previous path endpoint
              ref_x = previous_path_x[prev_size - 1];
              ref_y = previous_path_y[prev_size - 1];

              double ref_x_prev = previous_path_x[prev_size - 2];
              double ref_y_prev = previous_path_y[prev_size - 2];

              ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

              // Use two points that make the path tangent to the previous path's endpoint
              pts_x.push_back(ref_x_prev);
              pts_x.push_back(ref_x);
              pts_y.push_back(ref_y_prev);
              pts_y.push_back(ref_y);
            }

            for (int i = 1; i < 4; i++) {
              vector<double> next_waypoint = getXY(car_s + safe_dist * i, 4 * lane + 2, map_waypoints_s, map_waypoints_x, map_waypoints_y);
              pts_x.push_back(next_waypoint[0]);
              pts_y.push_back(next_waypoint[1]);
            }

            // Transform waypoints to local car coordinates
            for (int i = 0; i < pts_x.size(); i++) {
              double x = pts_x[i] - ref_x;
              double y = pts_y[i] - ref_y;

              pts_x[i] = x * cos(-ref_yaw) - y * sin(-ref_yaw);
              pts_y[i] = x * sin(-ref_yaw) + y * cos(-ref_yaw);
            }

            // Create spline
            tk::spline spline;
            // Set anchor points
            spline.set_points(pts_x, pts_y);

            // Calculate how to space spline points so that velocity is ref_speed
            double target_x = safe_dist;
            double target_y = spline(target_x);
            double target_dist = distance(0, 0, target_x, target_y);

            // Future path
            vector<double> next_x_vals;
          	vector<double> next_y_vals;

          	// Start with all previous path points to get the future path
            for(int i = 0; i < prev_size; i++)
            {
                  next_x_vals.push_back(previous_path_x[i]);
                  next_y_vals.push_back(previous_path_y[i]);
            }

            double x_addon = 0;

            for (int i = 1; i <= 50 - prev_size; i++) {
              // Accelerate or decelrate avoiding jerks
              if (reduce_speed) {
                ref_speed -= ref_speed_delta;
              }
              else if (ref_speed < speed_limit) {
                ref_speed += ref_speed_delta;
              }

              double N = target_dist / (0.02 * ref_speed);
              double x = x_addon + target_x / N;
              double y = spline(x);

              x_addon = x;

              double x_ref = x;
              double y_ref = y;

              // Transform waypoints back to map coordinates
              x = ref_x + x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw);
              y = ref_y + x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw);

              next_x_vals.push_back(x);
              next_y_vals.push_back(y);
            }

            json msgJson;
            msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
