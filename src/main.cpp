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
    double goal_s;
};

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }

double deg2rad(double x) { return x * pi() / 180; }

double rad2deg(double x) { return x * 180 / pi(); }

double logistic(double x){
    return 2.0 / (1 + exp(-x)) - 1.0;
}

double meanOfVector(vector<double> x){
    double tot;
    for(auto i: x)
        tot += i;
    return tot/x.size();
}

double maxofVector(vector<double> x){
    double biggest = -9999;
    for(auto i: x) {
        if (biggest < i)
            biggest = i;
    }
    return biggest;
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
    } else if (d > 8.0 & d < 12.0){
        lane = 2;
    } else{
        lane = -1;
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

// Get trajectory readings s[], d[], v[], a[], j[], given a trajectory (x[],y[])
vector<vector<double>> getTrajectoryReadings(vector<vector<double>> trajectory, vector<double> maps_x, vector<double> maps_y){

    vector<double> x, y, theta, s, d, vx, vy, vxy, vs, vd, vsd, ax, ay, axy, as, ad, asd, jx, jy, jxy, js, jd, jsd;

    for(int i = 0; i < trajectory.size(); i++){
        x.push_back(trajectory[i][0]);
        y.push_back(trajectory[i][1]);
    }

    if (x.size()>5){

        for (int i = 0; i < x.size() -1; i ++){
            double temp_theta = atan2((y[i+1] - y[i]), (x[i+1] - x[i]));
            theta.push_back(temp_theta);
            vector<double> temp_sd = getFrenet(x[i+1], y[i+1], temp_theta, maps_x, maps_y);
            s.push_back(temp_sd[0]);
            d.push_back(temp_sd[1]);
        }

        for (int j = 0; j < s.size()-1; j ++){
            double temp_vs = (s[j+1]-s[j])/0.02;
            double temp_vd = (d[j+1]-d[j])/0.02;
            double temp_vx = (x[j+1]-x[j])/0.02;
            double temp_vy = (y[j+1]-y[j])/0.02;

            vs.push_back(temp_vs);
            vd.push_back(temp_vd);
            vx.push_back(temp_vx);
            vy.push_back(temp_vy);

            vsd.push_back(sqrt(temp_vs*temp_vs + temp_vd*temp_vd));
            vxy.push_back(sqrt(temp_vx*temp_vx + temp_vy*temp_vy));
        }

        for (int j = 0; j < vs.size()-1; j ++){
            double temp_as = (vs[j+1]-vs[j])/0.02;
            double temp_ad = (vd[j+1]-vd[j])/0.02;
            double temp_ax = (vx[j+1]-vx[j])/0.02;
            double temp_ay = (vy[j+1]-vy[j])/0.02;

            as.push_back(temp_as);
            ad.push_back(temp_ad);
            ax.push_back(temp_ax);
            ay.push_back(temp_ay);

            asd.push_back(sqrt(temp_as*temp_as + temp_ad*temp_ad));
            axy.push_back(sqrt(temp_ax*temp_ax + temp_ay*temp_ay));
        }

        for (int j = 0; j < as.size()-1; j ++){
            double temp_js = (as[j+1]-as[j])/0.02;
            double temp_jd = (ad[j+1]-ad[j])/0.02;
            double temp_jx = (ax[j+1]-ax[j])/0.02;
            double temp_jy = (ay[j+1]-ay[j])/0.02;

            js.push_back(temp_js);
            jd.push_back(temp_jd);
            jx.push_back(temp_jx);
            jy.push_back(temp_jy);

            jsd.push_back(sqrt(temp_js*temp_js + temp_jd*temp_jd));
            jxy.push_back(sqrt(temp_jx*temp_jx + temp_jy*temp_jy));
        }

    }

    return {s, d, vsd, vxy, asd, axy, jsd, jxy};
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
double getRoadCurvature(double ix, double iy, double itheta, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y){

    Eigen::VectorXd x(10), y(10);

    int next_wp = NextWaypoint(ix, iy, itheta, maps_x, maps_y);

    int prev_wp;
    prev_wp = next_wp - 1;

    if (next_wp == 0) {
        prev_wp = maps_x.size() - 1;
    }

    for (int i = 0; i < 10; i ++){
        x[i] = maps_x[prev_wp+i];
        y[i] = maps_y[prev_wp+i];
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
pair<int, vector<vector<double>>> generateTrajectory(vector<double> sd, vector<double> prev_path_x, vector<double> prev_path_y, double ref_v,
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
double calculateCost(vector<vector<double>> trajectory, vector<vector<double>> ego_readings,
                     vector<vector<double>> sensor_fusion, int ego_goal_lane, double car_ahead_speed){

    double cost = 0.0;

    vector<double> s, d, vsd, vxy, asd, axy, jsd, jxy;
    s = ego_readings[0];
    d = ego_readings[1];
    vsd = ego_readings[2];
    vxy = ego_readings[3];
    asd = ego_readings[4];
    axy = ego_readings[5];
    jsd = ego_readings[6];
    jxy = ego_readings[7];

    double ego_end_s = s[s.size()-1];
    double ego_end_d = d[d.size()-1];
    double ego_vsd_max = maxofVector(vsd), ego_vxy_max = maxofVector(vxy);
    double ego_asd_max = maxofVector(asd), ego_axy_max = maxofVector(asd);
    double ego_jsd_max = maxofVector(jsd), ego_jxy_max = maxofVector(jsd);
    double ego_vsd_mean = meanOfVector(vsd), ego_vxy_mean = meanOfVector(vxy);
    double ego_asd_mean = meanOfVector(asd), ego_axy_mean = meanOfVector(axy);
    double ego_jsd_mean = meanOfVector(jsd), ego_jxy_mean = meanOfVector(jxy);

    double colli = 0.0, buffer = 0.0, v_lim = 0.0, a_lim = 0.0, j_lim = 0.0, effi = 0.0, road_lim = 0.0;

    double COLLISION_COST = 10.0;
    double BUFFER_COST = 1.0;
    double EFFICIENCY_COST = 1.0;

    int ego_cur_lane = calculateLane(d[0]);

//    for (int i = 0; i < s.size(); i++){

    for(auto sf: sensor_fusion){

        double check_car_speed = sqrt(sf[3] * sf[3] + sf[4] * sf[4]);
        double check_car_s0 = sf[5];
        double check_car_s = check_car_s0 + ((double)trajectory.size() * 0.02 * check_car_speed);
        double check_car_d = sf[6];
        int check_car_lane = calculateLane(check_car_d);

        // check both cur_lane and goal_lane for collision and buffer using s
        if (((ego_goal_lane == check_car_lane))){
            if ((ego_end_s > check_car_s && ego_end_s - check_car_s < 15) | (ego_end_s < check_car_s && check_car_s - ego_end_s < 15)){
                colli += (1.0 * COLLISION_COST);
                goto Out;
            } else if (abs(check_car_d - d[0]) < 3.0) {
                buffer += (1-logistic(abs(check_car_s - ego_end_s)) * BUFFER_COST);
            }
        }

//        for (int i = 0; i < s.size(); i ++){
//            check_car_s = check_car_s0 + ((double)i * 0.02 * check_car_speed);
//            if (((ego_goal_lane == check_car_lane) || (ego_cur_lane == check_car_lane))){
//                if ((s[i] > check_car_s && s[i] - check_car_s < 15) | (s[i] < check_car_s && check_car_s - s[i] < 15)){
//                    colli += (1.0 * COLLISION_COST);
//                    goto Out;
//                }
//            }
//        }

    }
//    }

    if(ego_vxy_max * 2.24 > 49.5||ego_vsd_max * 2.24 > 49.5)
        v_lim = 1.0;
    if(ego_axy_max > 10.0 || ego_asd_max > 10.0)
        a_lim = 1.0;
    if(ego_jxy_max > 50.0 || ego_jsd_max > 50.0)
        j_lim = 1.0;

    for(auto i: d){
        int temp_lane = calculateLane(i);
        if (temp_lane < -1){
            road_lim = 10.0;
            goto Out;
        }
    }

    effi = logistic(abs(49.5 - ego_vxy_mean * 2.24)/50.0) * EFFICIENCY_COST;

    Out:
    cost = colli + buffer + v_lim + a_lim + j_lim + effi + road_lim;

    if (cost < 10)
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
    ego.state = "START";
    ego.goal_lane = 1;
    ego.goal_s = 0.0;

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
                    double car_s0 = j[1]["s"];
                    double car_d0 = j[1]["d"];
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
                    int cur_lane = calculateLane(car_d0);

                    bool too_close_ahead = false;
                    bool check_car_ahead = false;
                    bool too_close_behind = false;
                    double check_car_ahead_vs = 60.0;

                    if (prev_size > 0){
                        car_d = end_path_d;
                        car_s = end_path_s;
                    }

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

                            if ((ego.state == "KL") && ((check_car_s > car_s) && (check_car_s - car_s < 30) && (check_car_s0 > car_s0))){

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
                    double curvature = getRoadCurvature(car_x, car_y, car_yaw, map_waypoints_s, map_waypoints_x, map_waypoints_y);
//                    cout << "Road curvature: " << curvature << endl;

                    vector<vector<double>> trajectory;

                    struct Ego ego_prev = ego;

                    if (ego.state == "LCL") {

                        if (abs(ego.goal_lane * 4 + 2 - car_d) < 1.0 && car_s0 - ego.goal_s > 30.0){
                            cout << "LCL completed" << endl;
                            ego.state = "KL";
                            goto KL;
                        }
                        else{
                            trajectory = generateTrajectoryFromGoalLane({car_s, car_d}, previous_path_x, previous_path_y, ref_v, ego.goal_lane,
                                                                        map_waypoints_s, map_waypoints_x, map_waypoints_y);
                        }
                    } else if (ego.state == "LCR") {

                        if (abs(ego.goal_lane * 4 + 2 - car_d) < 1.0 && car_s0 - ego.goal_s > 30.0){
                            cout << "LCR completed" << endl;
                            ego.state = "KL";
                            goto KL;
                        }
                        else{
                            trajectory = generateTrajectoryFromGoalLane({car_s, car_d}, previous_path_x, previous_path_y, ref_v, ego.goal_lane,
                                                                        map_waypoints_s, map_waypoints_x, map_waypoints_y);
                        }
                    } else if ((car_speed < 45.0) && (check_car_ahead) && (check_car_ahead_vs < 45.0/2.24) && curvature > 800.0) {

                        vector<vector<double>> anchors;
                        int anchor_lane;
                        vector<double> costs;
                        double cost = 9999;

                        anchors = generateAnchors(car_s, prev_size, sensor_fusion, cur_lane);

                        for (auto anchor: anchors) {

                            pair<int, vector<vector<double>>> temp_laneAndTrajectory;
                            vector<vector<double>> temp_trajectory;
                            int temp_lane = calculateLane(anchor[1]);
                            int goal_lane;

                            temp_laneAndTrajectory = generateTrajectory(anchor, previous_path_x, previous_path_y, ref_v,
                                                                        map_waypoints_s, map_waypoints_x, map_waypoints_y);

                            temp_trajectory = temp_laneAndTrajectory.second;

                            vector<vector<double>> ego_readings;
                            ego_readings = getTrajectoryReadings(temp_trajectory, map_waypoints_x, map_waypoints_y);

                            double temp_cost = calculateCost(temp_trajectory, ego_readings, sensor_fusion, goal_lane, check_car_ahead_vs);

                            costs.push_back(temp_cost);

                            if (temp_cost < cost){
                                cost = temp_cost;
                                trajectory = temp_trajectory;
                                anchor_lane = temp_lane;
                                ego.goal_s = ego_readings[0][ego_readings[0].size()-1];
                            }
                        }

                        if (cost > 10 | anchor_lane < 0 || anchor_lane > 2) {
                            goto KL;
                        }else if (anchor_lane < cur_lane){
                            ego.state = "PLCL";
                            cout << ego.state << endl;
                            cout << "LCL started" << endl;
                            ego.state = "LCL";
                            ego.goal_lane = cur_lane - 1;
                        } else {
                            ego.state = "PLCR";
                            cout << ego.state << endl;
                            cout << "LCR started" << endl;
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