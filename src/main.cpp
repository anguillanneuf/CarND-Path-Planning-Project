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

struct Ego{
    string state;
    int goal_lane;
};

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }

double deg2rad(double x) { return x * pi() / 180; }

double rad2deg(double x) { return x * 180 / pi(); }

double logistic(double x){
    return 2.0 / (1 + exp(-x)) - 1.0;
}

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

// Calculate lane number for d
int calculateLane(double d){
    int lane;
    if (d > 0.0 & d <= 4.0){
        lane = 0;
    } else if (d >4.0 & d <= 8.0) {
        lane = 1;
    } else{
        lane = 2;
    }
    return lane;
}

// Get Frenet for the end of a trajectory
vector<double> getFrenetFromTrajectory(vector<vector<double>> trajectory, vector<double> maps_x, vector<double> maps_y){
    int n = trajectory.size();

    double end_yaw = atan2((trajectory[n-1][1] - trajectory[n-2][1]), (trajectory[n-1][0] - trajectory[n-2][0]));

    vector<double> end_sd = getFrenet(trajectory[n-1][0], trajectory[n-1][1], end_yaw, maps_x, maps_y);

    return end_sd;
}

// Get approximate Ego readings given a path or trajectory
vector<double> getEgoReadings(vector<vector<double>> trajectory, vector<double> maps_x, vector<double> maps_y){
    double vs=0.0, vd=0.0, as=0.0, ad=0.0, js = 0.0, jd = 0.0, vx=0.0, vy=0.0, ax=0.0, ay=0.0, jx=0.0, jy=0.0;
    double vs_max=0.0, vd_max=0.0, as_max=0.0, ad_max=0.0, js_max = 0.0, jd_max = 0.0,
            vx_max=0.0, vy_max=0.0, ax_max=0.0, ay_max=0.0, jx_max = 0.0, jy_max = 0.0;
    double v_max = 0.0, a_max = 0.0, j_max = 0.0, v = 0.0;

    vector<double> previous_path_x, previous_path_y;
    for(int i = 0; i < trajectory.size(); i++){
        previous_path_x.push_back(trajectory[i][0]);
        previous_path_y.push_back(trajectory[i][1]);
    }

    if (previous_path_x.size()>5){
        vector<double> s_, d_, vs_, vd_, as_, ad_, theta_, vx_, vy_, ax_, ay_, v_;

        for (int i = 0; i < previous_path_x.size() -1; i ++){
            double temp_theta = atan2((previous_path_y[i+1] - previous_path_y[i]), (previous_path_x[i+1] - previous_path_x[i]));
            theta_.push_back(temp_theta);
            vector<double> temp = getFrenet(previous_path_x[i+1], previous_path_y[i+1], temp_theta, maps_x, maps_y);
            s_.push_back(temp[0]);
            d_.push_back(temp[1]);
        }
        // s / t
        for (int j = 0; j < s_.size()-1; j ++){
            vs_.push_back(s_[j+1]-s_[j]);
            vd_.push_back(d_[j+1]-d_[j]);
            vs += (s_[j+1]-s_[j]);
            vd += (d_[j+1]-d_[j]);
            if (vs_max < vs)
                vs_max = vs;
            if (vd_max < vd)
                vd_max = vd;
            vx_.push_back(previous_path_x[j+1]-previous_path_x[j]);
            vy_.push_back(previous_path_y[j+1]-previous_path_y[j]);
            vx += previous_path_x[j+1]-previous_path_x[j];
            vy += previous_path_y[j+1]-previous_path_y[j];
            if (vx_max < vx)
                vx_max = vx;
            if (vy_max < vy)
                vy_max = vy;
            if (sqrt(vx*vx+vy*vy) > v_max)
                v_max = sqrt(vx*vx+vy*vy);
            v_.push_back(sqrt(vx*vx+vy*vy));
            v += sqrt(vx*vx+vy*vy);
        }

        vs = vs/((s_.size()-1)*0.02);
        vd = vd/((d_.size()-1)*0.02);
        vx = vx/((previous_path_x.size()-1)*0.02);
        vy = vy/((previous_path_y.size()-1)*0.02);
        v = v/((v_.size()-1)*0.02);

        for (int j = 0; j < vs_.size()-1; j ++){
            as_.push_back(vs_[j+1]-vs_[j]);
            ad_.push_back(vd_[j+1]-vd_[j]);
            as += (vs_[j+1]-vs_[j]);
            ad += (vd_[j+1]-vd_[j]);
            if (as_max < as)
                as_max = as;
            if (as_max < ad)
                ad_max = ad;
            ax_.push_back(vx_[j+1]-vx_[j]);
            ay_.push_back(vy_[j+1]-vy_[j]);
            ax += vx_[j+1]-vx_[j];
            ay += vy_[j+1]-vy_[j];
            if (ax_max < ax)
                ax_max = ax;
            if (ay_max < ay)
                ay_max = ay;
            if (sqrt(ax*ax+ay*ay) > a_max)
                a_max = sqrt(ax*ax+ay*ay);
        }
        // vs / t
        as = as/((vs_.size()-1)*0.02);
        ad = ad/((vd_.size()-1)*0.02);
        ax = ax/((vx_.size()-1)*0.02);
        ay = ay/((vy_.size()-1)*0.02);

        for (int j = 0; j < as_.size()-1; j ++){
            js += (as_[j+1]-as_[j]);
            jd += (ad_[j+1]-ad_[j]);
            if (js_max < js)
                js_max = js;
            if (jd_max < jd)
                jd_max = jd;
            jx += ax_[j+1]-ax_[j];
            jy += ay_[j+1]-ay_[j];
            if (jx_max < jx)
                jx_max = jx;
            if (jy_max < jy)
                jy_max = jy;
            if (sqrt(jx*jx+jy*jy) > j_max)
                j_max = sqrt(jx*jx+jy*jy);
        }

        // vs / t
        js = js/((as_.size()-1)*0.02);
        jd = jd/((ad_.size()-1)*0.02);
        jx = jx/((ax_.size()-1)*0.02);
        jy = jy/((ay_.size()-1)*0.02);

    }

    return {vs, vd, as, ad, js, jd,
            vs_max, vd_max, as_max, ad_max, js_max, jd_max,
            vx, vy, ax, ay, jx, jy,
            vx_max, vd_max, as_max, ad_max, js_max, jd_max,
            v_max, a_max, j_max, v};
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

// Generate anchor sd points 4 seconds into the future (after the end of the previous path)
vector<vector<double>> generateAnchors(double car_s, int prev_size, vector<vector<double>> sensor_fusion, int lane){

    vector<int> lanes;
    lanes.push_back(lane-1); lanes.push_back(lane+1);

    vector<vector<double>> anchors;

    for(auto l: lanes){
        if(l > -1 | l < 3) {
            vector<vector<double>> cars_in_lane_sf;

            // I hope auto takes care of empty sensor_fusion
            for (auto sf: sensor_fusion) {
                if (sf[6] < l * 4 & sf[6] > (l + 1) * 4) {
                    cars_in_lane_sf.push_back(sf);
                }
            }

            vector<double> marks;
            marks.push_back(car_s + 90);

            // Look ahead 4 seconds (= 400 intervals) ahead of the end of the previous path
            for(auto sf: cars_in_lane_sf){
                double vx = sf[3];
                double vy = sf[4];
                double check_speed = sqrt(vx * vx + vy * vy);
                double check_car_s = sf[5];
                check_car_s += ((double)(prev_size + 200) * 0.02 * check_speed);
                marks.push_back(check_car_s);
            }

            // Drop anchors between marks
            // 1). sort marks;
            // 2). look 60m ahead of car_s, drop anchors every 2m apart, giving other cars a buffer zone of 15m
            // sort(marks.begin(), marks.end());
            for (int j = 1; j < 60; j += 2){
                bool ok_to_drop = true;
                for (auto m: marks){
                    if (abs(car_s + j - m) < 15)
                        ok_to_drop = false;
                }
                if (ok_to_drop)
                    anchors.push_back({car_s + j, (double)2 + l * 4});
            }
        }
    }
    return anchors;
}

// Generate trajectory from goal_lane
vector<vector<double>> generateTrajectoryFromGoalLane(vector<double> sd, vector<double> prev_path_x, vector<double> prev_path_y, double ref_v, int goal_lane,
                                                     vector<double> maps_s, vector<double> maps_x, vector<double> maps_y){

    vector<vector<double>> trajectory;
    vector<double> pts_x, pts_y;

    // add two points to pts_x and pts_y
    for (int i = 2; i > 0; i --){
        pts_x.push_back(prev_path_x[prev_path_x.size()-i]);
        pts_y.push_back(prev_path_y[prev_path_y.size()-i]);
    }

    double ref_x = pts_x[1];
    double ref_y = pts_y[1];
    double ref_prev_x = pts_x[0];
    double ref_prev_y = pts_y[0];
    double ref_yaw = atan2((ref_y - ref_prev_y), (ref_x - ref_prev_x));

    // add more points to pts_x and pts_y
    for (int i = 1; i <= 6; i ++){
        vector<double> xy = getXY(sd[0]+30*i, (2+4*goal_lane), maps_s, maps_x, maps_y);
        pts_x.push_back(xy[0]);
        pts_y.push_back(xy[1]);
    }

    // transform from global to local coordinates
    for(int i = 0; i < pts_x.size(); i++)
    {
        double shift_x = pts_x[i]-ref_x;
        double shift_y = pts_y[i]-ref_y;

        pts_x[i] = shift_x*cos(0 - ref_yaw) - shift_y*sin(0 - ref_yaw);
        pts_y[i] = shift_x*sin(0 - ref_yaw) + shift_y*cos(0 - ref_yaw);
    }

    // fit spline
    tk::spline s;
    s.set_points(pts_x, pts_y);

    // populate trajectory with previous path first
    for(int i = 0; i < prev_path_x.size(); i++)
        trajectory.push_back({prev_path_x[i], prev_path_y[i]});

    // plan 30m ahead along car's x direction
    double target_x = 30.0;
    double target_y = s(target_x);
    double target_distance = sqrt(target_x*target_x+target_y*target_y);
    double x_add_on = 0.0;

    // add on to previous path using points from spline
    for (int i = 1; i <= 50 - prev_path_x.size(); i ++){

        double N = target_distance/(ref_v/2.24*0.02);

        double x_point = x_add_on + target_x/N;
        double y_point = s(x_point);

        x_add_on = x_point;

        double x_ref = x_point;
        double y_ref = y_point;

        // transform from local to global coordinates
        x_point = ref_x + x_ref*cos(ref_yaw)-y_ref*sin(ref_yaw);
        y_point = ref_y + x_ref*sin(ref_yaw)+y_ref*cos(ref_yaw);

        trajectory.push_back({x_point, y_point});
    }

    return trajectory;
}

// Generate trajectory (x,y) from anchor (s, d)
pair<int, vector<vector<double>>> generateTrajectory(vector<double> sd, vector<double> prev_path_x, vector<double> prev_path_y,
                                          double ref_v,
                                          vector<double> maps_s, vector<double> maps_x, vector<double> maps_y){

    pair<int, vector<vector<double>>> laneAndTrajectory;
    vector<vector<double>> trajectory;
    vector<double> pts_x, pts_y;
    int curLane, endLane;
    curLane = calculateLane(sd[1]);

    // add two points to pts_x and pts_y
    for (int i = 2; i > 0; i --){
        pts_x.push_back(prev_path_x[prev_path_x.size()-i]);
        pts_y.push_back(prev_path_y[prev_path_y.size()-i]);
    }

    double ref_x = pts_x[1];
    double ref_y = pts_y[1];
    double ref_prev_x = pts_x[0];
    double ref_prev_y = pts_y[0];
    double ref_yaw = atan2((ref_y - ref_prev_y), (ref_x - ref_prev_x));

    vector<vector<double>> faraway;
    // add another four points to pts_x and pts_y
    for (int i = 1; i <= 4; i ++){
        vector<double> xy = getXY(sd[0]+30*i, (2+4*curLane), maps_s, maps_x, maps_y);

        pts_x.push_back(xy[0]);
        pts_y.push_back(xy[1]);
        faraway.push_back({xy[0], xy[1]});
    }


    // transform from global to local coordinates
    for(int i = 0; i < pts_x.size(); i++)
    {
        double shift_x = pts_x[i]-ref_x;
        double shift_y = pts_y[i]-ref_y;

        pts_x[i] = shift_x*cos(0 - ref_yaw) - shift_y*sin(0 - ref_yaw);
        pts_y[i] = shift_x*sin(0 - ref_yaw) + shift_y*cos(0 - ref_yaw);
    }

    // fit spline
    tk::spline s;
    s.set_points(pts_x, pts_y);

    // populate trajectory with previous path first
    for(int i = 0; i < prev_path_x.size(); i++)
        trajectory.push_back({prev_path_x[i], prev_path_y[i]});

    // plan 30m ahead along car's x direction
    double target_x = 30.0;
    double target_y = s(target_x);
    double target_distance = sqrt(target_x*target_x+target_y*target_y);
    double x_add_on = 0.0;

    // add on to previous path using points from spline
    for (int i = 1; i <= 60 - prev_path_x.size(); i ++){

        double N = target_distance/(ref_v/2.24*0.02);

        double x_point = x_add_on + target_x/N;
        double y_point = s(x_point);

        x_add_on = x_point;

        double x_ref = x_point;
        double y_ref = y_point;

        // transform from local to global coordinates
        x_point = ref_x + x_ref*cos(ref_yaw)-y_ref*sin(ref_yaw);
        y_point = ref_y + x_ref*sin(ref_yaw)+y_ref*cos(ref_yaw);

        trajectory.push_back({x_point, y_point});
    }

    vector<double> endsd = getFrenetFromTrajectory(faraway, maps_x, maps_y);
    endLane = calculateLane(endsd[1]);

    // trajectory has 50 elements now
    laneAndTrajectory.first = endLane;
    trajectory.resize(50);
    laneAndTrajectory.second = trajectory;
    return laneAndTrajectory;
}

// calculate cost
double calculateCost(vector<vector<double>> trajectory, vector<double> ego_begin_sd, vector<double> ego_end_sd, vector<double> ego_readings,
                     vector<vector<double>> sensor_fusion){

    double cost;
    double ego_end_s = ego_end_sd[0], ego_end_d = ego_end_sd[1];
    int ego_end_lane = calculateLane(ego_end_d);
    double ego_v_max = ego_readings[24], ego_a_max = ego_readings[25], ego_j_max = ego_readings[26], ego_v = ego_readings[27];

    double colli = 0.0, buffer = 0.0, v_lim = 0.0, a_lim = 0.0, j_lim = 0.0, effi;
    double COLLISION_COST = 10.0;
    double BUFFER_COST = 1.0;
    double EFFICIENCY_COST = 1.0;

    for (auto ego_xy: trajectory){

        for(auto sf: sensor_fusion){

            double check_car_speed = sqrt(sf[3] * sf[3] + sf[4] * sf[4]);
            double check_car_s = sf[5];
            check_car_s += ((double)trajectory.size() * 0.02 * check_car_speed);
            double check_car_d = sf[6];
            int check_car_lane = calculateLane(check_car_d);

            if (ego_end_lane == check_car_lane){
                if ((ego_end_s > check_car_s & ego_begin_sd[0] < sf[5]) | (ego_end_s < check_car_s & ego_begin_sd[0] > sf[5])){
                    colli += (1.0 * COLLISION_COST);
                } else {
                    buffer += (1/exp(abs(ego_end_s - check_car_s) - 15.0) * BUFFER_COST);
                }
            }
        }
    }

    if(ego_v_max * 2.24 > 50.0)
        v_lim = 1.0;
    if(ego_a_max > 10.0)
        a_lim = 1.0;
    if(ego_j_max > 50.0)
        j_lim = 1.0;

    effi = logistic(50.0 - ego_v * 2.24) * EFFICIENCY_COST;

    cost = colli + buffer + v_lim + a_lim + j_lim + effi;

    cout << cost << ": "<< colli << ", " << buffer << ", " << v_lim << ", " << a_lim << ", " << j_lim << ", " << effi <<endl;

    return cost;
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

    double ref_v = 0.0;
    struct Ego ego;
    ego.state = "KL";
    ego.goal_lane = 1;

    h.onMessage([&ref_v, &map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy, &ego](
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

                    struct Ego ego_ = ego;

                    int prev_size = previous_path_x.size();
                    int cur_lane = calculateLane(car_d);

                    bool too_close_ahead = false;
                    bool check_car_ahead = false;
                    bool too_close_behind = false;
                    double check_car_ahead_vs = 60.0;

                    if (prev_size > 0){
                        car_d = end_path_d;
                        car_s = end_path_s;
                    }

                    double car_s0 = j[1]["s"];
                    double adjacent_speed = 100.0;

                    //find ref_v to use
                    for (int i = 0; i < sensor_fusion.size(); i ++){
                        float d = sensor_fusion[i][6];
                        double closest_distance = 100.0;

                        if (d > 4 * cur_lane & d < 4 * (cur_lane + 1)){
                            check_car_ahead = true;
                            double vx = sensor_fusion[i][3];
                            double vy = sensor_fusion[i][4];
                            double check_speed = sqrt(vx * vx + vy * vy);
                            double check_car_s0 = sensor_fusion[i][5];
                            double check_car_s;

                            check_car_s = check_car_s0 + ((double)prev_size*0.02*check_speed); // when ego gets to the end of
                            // the previous trajectory, where would the other car be

                            if (((check_car_s > car_s) && (check_car_s - car_s < 30) && (check_car_s0 > car_s0))){

                                too_close_ahead = true;
                                // determine the speed of the first car ahead
                                if ((check_car_s - car_s) < closest_distance){
                                    closest_distance = check_car_s - car_s;
                                    check_car_ahead_vs = check_speed; // m/s
                                }
                            }

                            if(((check_car_s < car_s) && (car_s - check_car_s < 30) && (check_car_s0 < car_s0))){
                                too_close_behind = true;
                            }
                        }
                    }

                    if ((too_close_ahead && !too_close_behind))
                        ref_v -= 0.25;
                    else if (ref_v < 49.5 || (!too_close_ahead && too_close_behind)){
                        ref_v += 0.25; // more efficient if done in below
                        ref_v = min(ref_v, 49.5);
                    } else if (too_close_ahead || too_close_behind){
                        ref_v = check_car_ahead_vs;
                    }

                    vector<double> ptsx;
                    vector<double> ptsy;

                    // makes sure that previous_path has 2 points at least
                    if(prev_size < 2){
                        // create points tangent to the car
                        double prev_car_x = car_x - cos(deg2rad(car_yaw)) ;
                        double prev_car_y = car_y - sin(deg2rad(car_yaw));

                        ptsx.push_back(prev_car_x);
                        ptsx.push_back(car_x);

                        ptsy.push_back(prev_car_y);
                        ptsy.push_back(car_y);

                        previous_path_x = ptsx;
                        previous_path_y = ptsy;

                    }

                    // TODO: use map xy to find road curvature
//                    double curvature = getRoadCurvature(car_s, lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
//                    cout << "d_dot: " << vd << " Road curvature: " << curvature << endl;

                    vector<vector<double>> trajectory;

                    struct Ego ego_prev = ego;

                    if (ego.state == "LCL") {

                        if (abs(ego.goal_lane * 4 + 2 - car_d) < 1.0){
                            ego.state = "KL";
                            goto KL;
                        }
                        else{
                            trajectory = generateTrajectoryFromGoalLane({car_s, car_d}, previous_path_x, previous_path_y, ref_v, ego.goal_lane,
                                                                        map_waypoints_s, map_waypoints_x, map_waypoints_y);
                        }
                    } else if (ego.state == "LCR") {

                        if (abs(ego.goal_lane * 4 + 2 - car_d) < 1.0){
                            ego.state = "KL";
                            goto KL;
                        }
                        else{
                            trajectory = generateTrajectoryFromGoalLane({car_s, car_d}, previous_path_x, previous_path_y, ref_v, ego.goal_lane,
                                                                        map_waypoints_s, map_waypoints_x, map_waypoints_y);
                        }
                    } else if ((car_speed < 45.0) && (check_car_ahead) && (check_car_ahead_vs < 45.0/2.24)) {

                        cout << "Choose PLCL or PLCR" << endl;

                        vector<vector<double>> anchors;
                        int anchor_lane;
                        vector<double> costs;
                        double cost = 9999;

                        anchors = generateAnchors(car_s, prev_size, sensor_fusion, cur_lane);

                        for (auto anchor: anchors) {

                            pair<int, vector<vector<double>>> temp_laneAndTrajectory;
                            vector<vector<double>> temp_trajectory;
                            vector<double> ego_end_sd;
                            int temp_lane = calculateLane(anchor[1]);

                            if(temp_lane < cur_lane){
                                ego.state = "PLCL";
                                ego.goal_lane = cur_lane - 1;
                                if (ego.goal_lane < 0)
                                    continue;
                            } else {
                                ego.state = "PLCR";
                                ego.goal_lane = cur_lane + 1;
                                if (ego.goal_lane > 2)
                                    continue;
                            }

                            temp_laneAndTrajectory = generateTrajectory(anchor, previous_path_x, previous_path_y, ref_v,
                                                                        map_waypoints_s, map_waypoints_x, map_waypoints_y);

                            temp_trajectory = temp_laneAndTrajectory.second;

                            ego_end_sd = getFrenetFromTrajectory(temp_trajectory, map_waypoints_x, map_waypoints_y);

                            vector<double> ego_readings;
                            ego_readings = getEgoReadings(temp_trajectory, map_waypoints_x, map_waypoints_y);

                            double temp_cost = calculateCost(temp_trajectory, {car_s, car_d}, ego_end_sd, ego_readings, sensor_fusion);

                            costs.push_back(temp_cost);

                            if (temp_cost < cost){
                                cost = temp_cost;
                                trajectory = temp_trajectory;
                                anchor_lane = temp_lane;
                            }
                        }

                        if (cost > 10) {
                            goto KL;
                        }else if (anchor_lane < cur_lane){
                            cout << ego.state << endl;
                            ego.state = "LCL";
                            ego.goal_lane = cur_lane - 1;
                        } else {
                            cout << ego.state << endl;
                            ego.state = "LCR";
                            ego.goal_lane = cur_lane + 1;
                        }

                    } else {
                        KL:
                        ego.state = "KL";
                        trajectory = generateTrajectoryFromGoalLane({car_s, car_d}, previous_path_x, previous_path_y, ref_v, ego.goal_lane,
                                                                    map_waypoints_s, map_waypoints_x, map_waypoints_y);
                    }

                    if (ego_.state != ego.state)
                        cout << ego.state << " from " << cur_lane <<" to " << ego.goal_lane << endl;

                    // TODO: end

                    for (int i = 0; i < 50; i ++){
                        next_x_vals.push_back(trajectory[i][0]);
                        next_y_vals.push_back(trajectory[i][1]);
                    }

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