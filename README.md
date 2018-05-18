# FCND - Estimation Project

![C++ trajectory](images/Cover.png?raw=true)

## Overview

The goal of the project is to develop an estimator to be used by the controller from the previous project to successfully fly a desired flight path using realistic sensors.
<br><br><br>

## Associated files

Program and configuration:<br>
**C++ program:** `src/QuadEstimatorEKF.cpp`<br>
**C++ program:** `src/QuadControl.cpp`<br>
**C++ header:** `src/QuadControl.h`<br>
**Configuration:** `config/QuadEstimatorEKF.txt`<br>
**Configuration:** `config/QuadControlParams.txt`
<br><br><br>


## Contents

The following report consists of 3 sections:

**01 Implementation**<br>
    01.1 Sensor noise<br>
    01.2 Attitude estimation<br>
    01.3 Prediction<br>
    01.4 Magnetometer update<br>
    01.5 GPS update

**02 Flight evaluation**<br>
    02.1 Sensor noise<br>
    02.2 Attitude estimation<br>
    02.3 Prediction<br>
    02.4 Magnetometer update<br>
    02.5 GPS update

**03 References**<br>
    Books, research papers, and tools
<br><br><br>

## 01 Implementation

### 01.1 Sensor noise

Produced 9.27 seconds of data. Read into Excel. Use STDEV() function to calculate standard deviation of measurements. For GPS, the value is 0.723312 and for IMU, 0.490335. With some adjustments and a capture of 32.57s of data, then GPS: 0.6782  IMU: 0.4955. These are close to the simulated settings of 0.7 for GPU and 0.5 for IMU.

[see lines 3-4 in `06_SensorNoise.txt`]

### 01.2 Attitude estimation

According to Quan, "nonlinear complementary filters are based on a nonlinear relationship between
the angular velocity and the angle of rotation."

To reduce attitude estimation error, first the attitude estimate is converted to a quaternion, then integrated with the Quaternion class function `IntegrateBodyRate()`. Finally, it is converted back to Euler angles using the convenience functions `Roll()`, `Pitch()`, and `Yaw()`.

[see lines 100-109 in `QuadEstimatorEKF.cpp`]

### 01.3 Prediction

Implemented integrate forward, calculate the partial derivative of the body-to-global rotation matrix, and predict the state covariance forward.

[see lines 172-181 in `QuadEstimatorEKF.cpp`]<br>
[see lines 209-219 in `QuadEstimatorEKF.cpp`]<br>
[see lines 266-276 in `QuadEstimatorEKF.cpp`]<br>

Note: in reference paper [2], there is an error in formula 52 on page 9. (I reported this as issue #320 in Waffle.) The code in QuadEstimatorEKF::GetRbgPrime is implemented with the correction (at line 215) to match Diebel[3] equation 71.

### 01.4 Magnetometer update

Implemented `UpdateFromMag` to include the magnetometer data into the state. The solution normalizes the difference between the measured and estimated yaw.

[see lines 328-335 in `QuadEstimatorEKF.cpp`]<br>

### 01.5 GPS update

Implemented GPS update. `hprime` is assigned as an identity matrix as in Diebel[3] equation 55.

[see lines 302-305 in `QuadEstimatorEKF.cpp`]<br>
<br><br><br>

## 02 Flight evaluation

### 02.1 Sensor noise

The calculated standard deviation correctly captures ~68% of the sensor measurements.

Unless otherwise noted, results are from tests made with the supplied controller (described as "relaxed to work with an estimated state").

**RESULTS:**

    Simulation #2 (../config/06_SensorNoise.txt)
    PASS: ABS(Quad.GPS.X-Quad.Pos.X) was less than MeasuredStdDev_GPSPosXY for 72% of the time
    PASS: ABS(Quad.IMU.AX-0.000000) was less than MeasuredStdDev_AccelXY for 67% of the time

![Scenario 6](images/Figure1.png?raw=true)<br>
Figure 1: Scenario 6 - measuring standard deviation

## 

### 02.2 Attitude estimation

**RESULTS:**

    Simulation #5 (../config/07_AttitudeEstimation.txt)
    PASS: ABS(Quad.Est.E.MaxEuler) was less than 0.100000 for at least 3.000000 seconds

![Scenario 6](images/Figure2.png?raw=true)<br>
Figure 2: Scenario 7 - reducing estimated attitude error

## 

### 02.3 Prediction

Estimator integrates forward to track the actual state.

![Scenario 8](images/Figure3.png?raw=true)<br>
Figure 3: Scenario 8 - prediction with reasonably slow drift

## 

![Scenario 9](images/Figure4.png?raw=true)<br>
Figure 4: Scenario 9 - before tuning

## 

Tuning the InitStdDevs and the process noise model brings the errors in range.

![Scenario 9](images/Figure5.png?raw=true)<br>
Figure 5: Scenario 9 - after tuning in X dimension

## 

![Scenario 9](images/Figure6.png?raw=true)<br>
Figure 6: Scenario 9 - after tuning in Y dimension

## 

![Scenario 9](images/Figure7.png?raw=true)<br>
Figure 7: Scenario 9 - after tuning in Z dimension

## 

### 02.4 Magnetometer update

Maintain an error of less than 0.1 radians in heading for at least 10 seconds.

**RESULTS:**

    Simulation #3 (../config/10_MagUpdate.txt)
    PASS: ABS(Quad.Est.E.Yaw) was less than 0.100000 for at least 10.000000 seconds
    PASS: ABS(Quad.Est.E.Yaw-0.000000) was less than Quad.Est.S.Yaw for 73% of the time

![Scenario 10](images/Figure8.png?raw=true)<br>
Figure 8: Scenario 10 - magnetometer update

## 

### 02.5 GPS update

Tuned the GPSPosXYStd and GPSPosZStd values in QuadEstimatorEKF.txt

**RESULTS:**

    Simulation #6 (../config/11_GPSUpdate.txt)
    PASS: ABS(Quad.Est.E.Pos) was less than 1.000000 for at least 20.000000 seconds

![Scenario 11](images/Figure9.png?raw=true)<br>
Figure 9: Scenario 11 - GPS update

## 

**NOTE:** Results below are from tests made with the controller from the prior project. (Located in: https://github.com/jwdunn1/FCND-Controls-CPP)<br>

Initial testing using the prior project controller (before tuning) was metrically successful (estimated position error of < 1m), however the aggressive tuning causes the path to overshoot the corners of the box trajectory.

**RESULTS:**

    Simulation #4 (../config/11_GPSUpdate.txt)
    PASS: ABS(Quad.Est.E.Pos) was less than 1.000000 for at least 20.000000 seconds

![Scenario 11](images/Figure10.png?raw=true)<br>
Figure 10: Scenario 11 - GPS update pretuning

![Scenario 11](images/Figure11.png?raw=true)<br>
Figure 11: Scenario 11 - GPS update pretuning (overhead)

## 

To mitigate the overshoots, the position and velocity gains were reduced from the more aggressive initial settings until a box shape is achieved.

![Scenario 11](images/Figure12.png?raw=true)<br>
Figure 12: Scenario 11 - GPS update tuned

![Scenario 11](images/Figure13.png?raw=true)<br>
Figure 13: Scenario 11 - GPS update tuned (overhead)

## 

**NOTE:** The original (more aggressive) controller gains are required for scenarios 3, 4, and 5. Successful passing of these scenarios is not required for this project.

<br><br><br>

## 03 References

[1] ***Introduction to Multicopter Design and Control***<br>
Quan Quan<br>
Springer Singapore, 2017

[2] ***Estimation for Quadrotots***<br>
Tellex, Brown, and Lupashin<br>
https://www.overleaf.com/read/vymfngphcccj

[3] ***Representing Attitude: Euler Angles, Unit Quaternions, and Rotation Vectors***<br>
James Diebel<br>
https://www.astro.rug.nl/software/kapteyn-beta/_downloads/attitude.pdf

[4] ***Kalman Tutorial***<br>
Simon Levy<br>
https://home.wlu.edu/~levys/kalman_tutorial/

[5] ***Kalman Filters on Youtube***<br>
Michel van Biezen<br>
https://www.youtube.com/playlist?list=PLX2gX-ftPVXU3oUFNATxGXY90AULiqnWT

[6] ***A New Extension of the Kalman Filter to Nonlinear
Systems***<br>
Julier and Uhlmann<br>
https://people.eecs.berkeley.edu/~pabbeel/cs287-fa15/optreadings/JulierUhlmann-UKF.pdf

[7] ***The Unscented Kalman Filter for Nonlinear Estimation***<br>
Wan and van der Merwe<br>
https://www.seas.harvard.edu/courses/cs281/papers/unscented.pdf

[8] ***A New Method for the Nonlinear Transformation of Means and Covariances in Filters and Estimators***<br>
Julier, Uhlmann, and Durrant-Whyte<br>
https://pdfs.semanticscholar.org/76c5/2206888eedef8d8dead3007992e53e3c4ae8.pdf

[9] ***Python Kalman filter***<br>
Daniel Duckworth<br>
https://pykalman.github.io/

[10] ***Quadcopter Dynamics, Simulation, and Control***<br>
Andrew Gibiansky<br>
http://andrew.gibiansky.com/downloads/pdf/Quadcopter%20Dynamics,%20Simulation,%20and%20Control.pdf

[11] ***Computing Euler angles from a rotation matrix***<br>
Gregory Slabaugh<br>
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.371.6578&rep=rep1&type=pdf