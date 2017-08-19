# Model Documentation  

The model in `main.cpp` makes use of finite state machines and cost functions to help Ego execute highway driving. There are five states in this finite state machine: Keep Lane (KL), Prepare for Lane Change Left (PLCL), Prepare for Lane Change Right (PCLR), Lane Change Left (LCL), and Lane Change Right (LCR). These describe the motion that Ego needs to undertake. Cost functions are only activated in PLCL or PLCR states to help Ego determine the safest and smoothest trajectory to follow next.  You may click on the image below to view a YouTube clip that shows Ego go past the 20 mile mark without any road accident.  

[![img](/Screen%20Shot%202017-08-18%20at%2011.09.02%20PM.png)](https://youtu.be/kIvFfGi8Wew)  

## Basic Structure  

The `TODO` section of `main.cpp` has my code for generating the (x, y) points every 0.02 seconds. Line 700-780 contains mostly code taken from the project walkthrough [video](https://youtu.be/7sI3VHFPP0w). The main purpose there is to modify Ego's speed `ref_v` by checking the speed of traffic ahead of it, if any. In addition to that, speed of the traffic ahead, the road curvature, plus Ego's own speed inform me later if I need to trigger the PLCL/PCLR states. The pseudo-code below describes my logic.  

```C++
// what's Ego's current state?
// is there any car closely in front of Ego?
// what is this car's speed?

if (state is LCL) {
  if (LCL is complete) {
    // goto KL
  } else {
    // continue with its planning using its previous path in (x,y), end of its previous path in (s,d), 
    // its goal lane, ref_v, map
  }
} else if (state is LCR) {
  if (LCR is complete) {
    // goto KL
  } else {
    // continue with its planning
  }
} else if (there is a car in front & the car is slow & Ego is slow too & road is straight) {
    // generate some anchor points in prospective lanes
    for (each of these anchor points) {
      // generate a trajectory
      // calculate the cost of this trajectory and pick one with the smallest cost
    }
    
    if (even the smallest cost is too big) {
      // goto KL
    } else {
      // update Ego state to either LCL or LCR
    }
} else {
    // continue with its planning
}     
```
As shown, the procedure to follow in LCL/LCR/KL is quite stragihtforward. During PLCL/PLCR, the only complication is to generate anchor points first, then evaluate them based on a cost function, the rest is the same as LCL/LCR/KL.  

## Anchor Points  

The code for `generateAnchors` is in Line 318-398. Given the end (s,d) of Ego's previous path, sensor fusion data on other cars on the road, Ego's current lane, and Ego's current speed, this is how anchor points are created and selected.  

```C++
// consider the lane to the left and right of Ego's current lane

if (lane is valid) {
  // find all the cars in this target lane
  
  // find out which ones can overtake Ego under 1.5 seconds provided that they are only 0-60m behind
  // gather their IDs, speed, time to overtake Ego, distance along the road when the overtake happens
  
  if (no cars exist) {
    // drop some anchors 0-30m beyond the end of the previous path
  } else {
    for (consider every point 0-30m from the end of Ego's previous path) {
      if (this point doesn't come too close into contact with other cars in the same lane){
        // drop an anchor
      }
    }
  }
}

```

After the anchors are in place, I make use of `spline.h` to fit a line that goes through the last two points of Ego's previous path, the anchor point, as well as a few more points beyond the anchor point in the anchor lane. This code is very similar to the code given in the project walkthrough and it's in Line 404-477 in `main.cpp`.  

## Cost Function  

Lastly, after the trajectories are generated, I use a cost function to decide which ones are the smoothest and safest. Collision, going beyond the road limit, and being baffled by the traffic situations immediately cause the cost function to exit and report a large cost. Here is how they are triggered and how other costs are calculated in my code (Line 480-606).  

1. **Collision**  
  When Ego approaches the center line between its current lane and its goal lane, any car is within 15m in the goal lane.  
  When Ego reaches the end of its planned trajectory, any car is within 15m in its current lane. 
2. **Road Limit**  
  For all the points (x,y) along the trajectory, does any of them fall outside of the road?  
3. **Baffled**  
  If the traffic slightly behind and ahead of the slow car in the goal lane is not much faster, Ego is baffled about whether to change lane or not. 
4. **Buffer**  
  The closer the other cars are to Ego at the end of the planned trajectory, the bigger the buffer cost. Conversely, the buffer goes down.  
5. **Maximum Speed/Acceleration/Jerk**  
  As specified in the project `README.md`, maximum speed is 50 mph, maxinum acceleration is 10 m/s<sup>2</sup>, and maximum jerk is 50 m/s<sup>3</sup>. 
6. **Efficiency**  
  The average speed for the considered trajectory shall be as fast as it can be.  

# Reflection  

For someone who doesn't own a driver license and doesn't know how to drive, I approached this project with the simplest road manuvre logic there can be, which is to only change lanes when the traffic ahead is too slow. I image that there are better and smarter ways to navigate highway traffic, and I imagine that other advanced search algorithms or reinforcement learning can realize them. I would be excited to explore them in the future.  

In the beginning of doing this project, I only thought of looking ahead of where Ego is going. After I got that to work, I realized that it could be more efficient if I looked for gaps behind Ego. That led me to change my `generateAnchors` function to first calculate the time it takes for vehicles behind Ego to catch up with it, and then drop anchors right behind where that happens. I really liked the fact that I came to this realization again while doing the exercise, although I have been taught about it in the lecture video.  

My model is able to drive in the simulator for rounds most of the time. However, I noticed that sometimes it had a difficult time when dealing with cars suddenly changing lanes as well as road accidents caused by other cars. I believe that doing predictions on each car on the road using Unscented Kalman Filters would help in those situations. Besides that, I haven't touched the part of the code that could vary acceleration in the same cycle to allow Ego to accelerate or deccelerate more efficiently. That will be the first thing I will tweak to enter for the Bosch Challenge. 
