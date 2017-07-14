# Unscented_Kalman_Filter
SDC Term 2 Project 2

The code from lessons were helpful.
* RMSE obtained for Dataset 1 (0.0661, 0.0856, 0.2887, 0.2182).
* RMSE obtained for Dataset 2 (0.0762, 0.0622, 0.5363, 0.2080).

For tuning process noise parameters, I calculated the maximum logitudinal and angular acceleration from the dataset 1.
* Max Longitudinal Acceleration - 1.00 m/s^2.
* Max Angular Acceleration - 0.55 rad/s^2.

So I tunes and arrived at values close to half of the max accelerations.
* stda = 0.44.
* stdyawdd = 0.275.

With these parameters, RMSE was reducing but was not less than the required.
I tuned the intial P_covariance matrix values t obtain the required RMSE values.
P = 0.15,0,0,0,0,
0,0.15,0,0,0,
0,0,0.5,0,0,
0,0,0,1.1,0,
0,0,0,0,1;

The NIS was calculated with each updation step.
Dataset 1.
* For Radar 94 % of the cases NIS is lesser than 7.18.
* For Laser 98 % of the cases NIS is lesser than 5.99.

Dataset 2.
* For Radar 93 % of the cases NIS is lesser than 7.18.
* For Laser 95 % of the cases NIS is lesser than 5.99.

Also evaluated RMSE for Radar only and Laser only cases.
* Laser Only Dataset 1 (0.1774,0.1465, 0.6417,0.2469).
* Laser Only Dataset 2 (0.1718,0.1378,0.6334,0.2680).
* Radar Only Dataset 1 (0.2116, 0.2186, 0.4601, 0.2421).
* Radar Only Dataset 2 (0.5211,0.2315,0.8942,0.3165).

With both sensors combined the RMSE has definitely reduced.

EKF RMSE from previous project.
RMSE obtained for Dataset 1. (0.0973,0.0855, 0.4513, 0.4399). 
RMSE obtained for Dataset 2. (0.0726,0.0967, 0.4579, 0.4966).

RMSE with UKF has reduced in most of the aspects especially in case of velocity magnitude.




