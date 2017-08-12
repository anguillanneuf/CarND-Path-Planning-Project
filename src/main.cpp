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

double distance(double x1, double y1, double x2, double y2) {
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y) {

    double closestLen = 100000; //large number
    int closestWaypoint = 0;

    for (int i = 0; i < maps_x.size(); i++) {
        double map_x = maps_x[i];
        double map_y = maps_y[i];
        double dist = distance(x, y, map_x, map_y);
        if (dist < closestLen) {
            closestLen = dist;
            closestWaypoint = i;
        }

    }

    return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y) {

    int closestWaypoint = ClosestWaypoint(x, y, maps_x, maps_y);

    double map_x = maps_x[closestWaypoint];
    double map_y = maps_y[closestWaypoint];

    double heading = atan2((map_y - y), (map_x - x)); // atan2 in [-PI, PI], where wp is w.r.t. car

    double angle = abs(theta - heading); // difference in car yaw and heading

    if (angle > pi() / 4) // if point is not in car's line of sight, consider it behind
        closestWaypoint++;

    return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y) {
    int next_wp = NextWaypoint(x, y, theta, maps_x, maps_y);

    int prev_wp;
    prev_wp = next_wp - 1;
    if (next_wp == 0) {
        prev_wp = maps_x.size() - 1;
    }

    double n_x = maps_x[next_wp] - maps_x[prev_wp];
    double n_y = maps_y[next_wp] - maps_y[prev_wp];
    double x_x = x - maps_x[prev_wp];
    double x_y = y - maps_y[prev_wp];

    // find the projection of x onto n
    double proj_norm = (x_x * n_x + x_y * n_y) / (n_x * n_x + n_y * n_y);
    double proj_x = proj_norm * n_x;
    double proj_y = proj_norm * n_y;

    double frenet_d = distance(x_x, x_y, proj_x, proj_y);

    //see if d value is positive or negative by comparing it to a center point

    double center_x = 1000 - maps_x[prev_wp];
    double center_y = 2000 - maps_y[prev_wp];
    double centerToPos = distance(center_x, center_y, x_x, x_y);
    double centerToRef = distance(center_x, center_y, proj_x, proj_y);

    if (centerToPos <= centerToRef) {
        frenet_d *= -1;
    }

    // calculate s value
    double frenet_s = 0;
    for (int i = 0; i < prev_wp; i++) {
        frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
    }

    frenet_s += distance(0, 0, proj_x, proj_y);

    return {frenet_s, frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y) {
    int prev_wp = -1;

    while (s > maps_s[prev_wp + 1] && (prev_wp < (int) (maps_s.size() - 1))) {
        prev_wp++;
    }

    int wp2 = (prev_wp + 1) % maps_x.size();

    double heading = atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));
    // the x,y,s along the segment
    double seg_s = (s - maps_s[prev_wp]);

    double seg_x = maps_x[prev_wp] + seg_s * cos(heading);
    double seg_y = maps_y[prev_wp] + seg_s * sin(heading);

    double perp_heading = heading - pi() / 2;

    double x = seg_x + d * cos(perp_heading);
    double y = seg_y + d * sin(perp_heading);

    return {x, y};

}

// Get approximate Ego vs, vd, as, ad
vector<double> getEgoReadings(vector<double> previous_path_x, vector<double> previous_path_y, double theta, vector<double> maps_x, vector<double> maps_y){
    double vs=0.0, vd=0.0, as=0.0, ad=0.0;

    if (previous_path_x.size()>3){
        vector<double> vs_, vd_, theta_;

        for (int i = 0; i < previous_path_x.size() -1; i ++){
            double temp_theta = atan2((previous_path_y[i+1] - previous_path_y[i]), (previous_path_x[i+1] - previous_path_x[i]));
            theta_.push_back(temp_theta);
            vector<double> temp = getFrenet(previous_path_x[i+1], previous_path_y[i+1], temp_theta, maps_x, maps_y);
            vs_.push_back(temp[0]);
            vd_.push_back(temp[1]);
            vs += temp[0];
            vd += temp[1];
        }
        // s / t
        vs = vs/((previous_path_x.size()-1)*0.02);
        vd = vd/((previous_path_y.size()-1)*0.02);

        for (int j = 0; j < vs_.size()-1; j ++){
            as += (vs_[j+1]-vs_[j]);
            ad += (vd_[j+1]-vd_[j]);
        }
        // vs / t
        as = as/((vs_.size()-1)*0.02);
        ad = ad/((vd_.size()-1)*0.02);
    }

    return {vs, vd, as, ad};
}

// Fit polynomial
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals, int order) {
    assert(xvals.size() == yvals.size());
    assert(order >= 1 && order <= xvals.size() - 1);
    Eigen::MatrixXd A(xvals.size(), order + 1);

    for (int i = 0; i < xvals.size(); i++) {
        A(i, 0) = 1.0;
    }

    for (int j = 0; j < xvals.size(); j++) {
        for (int i = 0; i < order; i++) {
            A(j, i + 1) = A(j, i) * xvals(j);
        }
    }

    auto Q = A.householderQr();
    auto result = Q.solve(yvals);
    return result;
}

// Find road curvature
double getRoadCurvature(double car_s, int lane, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y){

    Eigen::VectorXd x(10), y(10);

    for (int i = 0; i < 10; i ++){
        vector<double> temp = getXY(car_s + 10 * i, (2 + 4 * lane), maps_s, maps_x, maps_y);
        x[i] = temp[0];
        y[i] = temp[1];
    }

    auto coeffs = polyfit(x, y, 3);

    double dfdx = coeffs[1] + 2 * coeffs[2] * x[0] + 3 * coeffs[3] * pow(x[0], 2);
    double dfdx2 = 2 * coeffs[2] + 2 * 3 * coeffs[3] * x[0];
    return pow(1 + dfdx * dfdx, 1.5) / abs(dfdx2);
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


    int lane = 1;
    double ref_v = 0.0;


    h.onMessage([&ref_v, &map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy, &lane](
            uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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
                    double car_yaw = j[1]["yaw"]; // in degrees
                    double car_speed = j[1]["speed"];


                    // Previous path data given to the Planner
                    auto previous_path_x = j[1]["previous_path_x"];
                    auto previous_path_y = j[1]["previous_path_y"];
                    // Previous path's end s and d values
                    double end_path_s = j[1]["end_path_s"];
                    double end_path_d = j[1]["end_path_d"];

                    // Sensor Fusion Data, a list of all other cars on the same side of the road.
                    auto sensor_fusion = j[1]["sensor_fusion"];

                    json msgJson;

                    vector<double> next_x_vals;
                    vector<double> next_y_vals;

                    // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds

                    int prev_size = previous_path_x.size();

                    if (prev_size > 0)
                        car_s = end_path_s;

                    // pick dd threshold
                    // generate anchor points
                    // generate trajectories from anchor points
                    // use trajectories, map, sensor_fusion to calculate cost
                    // choose trajectory with minimum cost
                    // special cases handling: when trajectory meets changes (sunday)


                    bool too_close = false;

                    //find ref_v to use
                    for (int i = 0; i < sensor_fusion.size(); i ++){
                        float d = sensor_fusion[i][6];
                        if (d > 4 * lane & d < 4 * (lane + 1)){
                            double vx = sensor_fusion[i][3];
                            double vy = sensor_fusion[i][4];
                            double check_speed = sqrt(vx * vx + vy * vy);
                            double check_car_s = sensor_fusion[i][5];

                            check_car_s += ((double)prev_size*0.02*check_speed); // when ego gets to the end of
                            // the previous trajectory, where would the other car be

                            if ((check_car_s > car_s) & (check_car_s - car_s < 30)){
                                //ref_v = 30;
                                too_close = true;
                                // time to change lanes, lower speed
                                if (lane > 0)
                                    lane = 0;
                            }
                        }
                    }

                    if (too_close)
                        ref_v -= 0.25;
                    else if (ref_v < 49.5){
                        ref_v += 0.25; // more efficient if done in below
                    }


                    vector<double> ptsx;
                    vector<double> ptsy;
                    double ref_x = car_x;
                    double ref_y = car_y;
                    double ref_yaw = deg2rad(car_yaw);

                    if(prev_size < 2){
                        // create points tangent to the car
                        double prev_car_x = car_x - cos(deg2rad(car_yaw)) ;
                        double prev_car_y = car_y - sin(deg2rad(car_yaw));

                        ptsx.push_back(prev_car_x);
                        ptsx.push_back(car_x);

                        ptsy.push_back(prev_car_y);
                        ptsy.push_back(car_y);

                    } else {
                        // redefine reference state using previous path
                        ref_x = previous_path_x[prev_size-1];
                        ref_y = previous_path_y[prev_size-1];

                        double ref_prev_x = previous_path_x[prev_size - 2];
                        double ref_prev_y = previous_path_y[prev_size - 2];
                        ref_yaw = atan2((ref_y - ref_prev_y), (ref_x - ref_prev_x));

                        ptsx.push_back(ref_prev_x);
                        ptsx.push_back(ref_x);
                        ptsy.push_back(ref_prev_y);
                        ptsy.push_back(ref_y);
                    }

                    // TODO: find ego_readings vd [vs, as, vd, ad]
                    vector<double> ego_readings = getEgoReadings(previous_path_x, previous_path_y, car_yaw, map_waypoints_x, map_waypoints_y);
                    double ad = ego_readings[3];
                    cout << "d_dot: " << ad << endl; // abs(ad0 goes above 1 when changing lanes.

                    // TODO: use map xy to find road curvature
                    double curvature = getRoadCurvature(car_s, lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
                    cout << "Road curvature: " << curvature << endl;// R_curve < 50.0

                    // add more points to generate trajectory
                    vector<double> next_wp0 = getXY(car_s+30, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
                    vector<double> next_wp1 = getXY(car_s+60, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
                    vector<double> next_wp2 = getXY(car_s+90, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

                    ptsx.push_back(next_wp0[0]);
                    ptsx.push_back(next_wp1[0]);
                    ptsx.push_back(next_wp2[0]);

                    ptsy.push_back(next_wp0[1]);
                    ptsy.push_back(next_wp1[1]);
                    ptsy.push_back(next_wp2[1]);

                    for(int i = 0; i < ptsx.size(); i++)
                    {
                        double shift_x = ptsx[i]-ref_x;
                        double shift_y = ptsy[i]-ref_y;

                        ptsx[i] = shift_x*cos(0 - ref_yaw) - shift_y*sin(0 - ref_yaw); // yaw is +ve when turning right
                        ptsy[i] = shift_x*sin(0 - ref_yaw) + shift_y*cos(0 - ref_yaw);
                    }

                    // adding the anchor points to ptsx and ptsy.
                    tk::spline s;
                    s.set_points(ptsx, ptsy);


                    for(int i = 0; i < previous_path_x.size(); i++){
                        next_x_vals.push_back(previous_path_x[i]);
                        next_y_vals.push_back(previous_path_y[i]);
                    }

                    double target_x = 30.0;// plan 30m ahead along car's x direction
                    double target_y = s(target_x);
                    double target_distance = sqrt(target_x*target_x+target_y*target_y);//distance(0,0,target_x,target_y);
                    double x_add_on = 0.0;


                    for (int i = 1; i <= 50 - previous_path_x.size(); i ++){

                        double N = target_distance/(ref_v/2.24*0.02);

                        double x_point = x_add_on + target_x/N;
                        double y_point = s(x_point);

                        x_add_on = x_point;

                        double x_ref = x_point;
                        double y_ref = y_point;

                        x_point = ref_x + x_ref*cos(ref_yaw)-y_ref*sin(ref_yaw);
                        y_point = ref_y + x_ref*sin(ref_yaw)+y_ref*cos(ref_yaw);

                        next_x_vals.push_back(x_point);
                        next_y_vals.push_back(y_point);

                    }



                    // TODO: end

//                    for (auto i: sensor_fusion)
//                        cout << i << " ";
//                    cout << endl;


                    msgJson["next_x"] = next_x_vals;
                    msgJson["next_y"] = next_y_vals;

                    auto msg = "42[\"control\"," + msgJson.dump() + "]";

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
















































































