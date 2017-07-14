#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a = 0.44;//0.44

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd = 0.275; //0.275

  // Laser measurement noise standard deviation position1 in m
  std_laspx = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.
  
  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized=false;
  n_x = 5;
  n_aug = 7;
  lambda = 3 - n_aug;
  Xsig_pred_ = MatrixXd(n_x, 2*n_aug+1);
  weights = VectorXd(2*n_aug+1);
  x_ << 0,0,0,0,0;
  P_ << 0.15,0,0,0,0,
        0,0.15,0,0,0,
        0,0,0.5,0,0,
        0,0,0,1.1,0,
        0,0,0,0,1;
  laser_meas_cnt = 0;
  radar_meas_cnt = 0;
  laser_chi_thresh = 0;
  radar_chi_thresh = 0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */

int laser_cnt = 0;
int radar_cnt = 0;
int laser_95_cnt = 0;
int radar_95_cnt = 0;
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized)
  {
    cout << "First" << endl;
    prev_timestamp_ = meas_package.timestamp_;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ == true)
    {
      x_[0] = meas_package.raw_measurements_[0]*cos(meas_package.raw_measurements_[1]);
      x_[1] = meas_package.raw_measurements_[0]*sin(meas_package.raw_measurements_[1]); 
      is_initialized = true;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == true)
    {
      x_[0] = meas_package.raw_measurements_[0];
      x_[1] = meas_package.raw_measurements_[1];
      is_initialized = true;
    }

    cout << "Initialized" << endl;
    // set weights
    double weight_0 = lambda/(lambda+n_aug);
    weights(0) = weight_0;
    for (int i=1; i<2*n_aug+1; i++) 
    {  //2n+1 weights
      double weight = 0.5/(n_aug+lambda);
      weights(i) = weight;
    }
    return;
  }

  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ == true) ||
      (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == true))
  { 
    /*Calcualte the time difference betwen the measurements*/
    double delta_t =  (meas_package.timestamp_ - prev_timestamp_)/1000000.0;
    prev_timestamp_ = meas_package.timestamp_;
    /*Predict the state mean and covariance*/
    Prediction(delta_t);
    cout << "Predicted P" << P_ << endl;
  }

  /*Update the state mean and covariance based on the measurement*/
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ == true)
  {
    UpdateRadar(meas_package);
    cout << "Radar Updated P " << P_ << endl;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == true)
  {
    UpdateLidar(meas_package);
    cout << "Laser Updated P " << P_ << endl;
  }
  return;
}

/*Generate the Augmented Sigma points*/
void UKF::GenerateAugSigmaPoints(MatrixXd* Xsig_out)
{

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug, n_aug);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

  //create augmented mean state
  x_aug << 0,0,0,0,0,0,0;
  x_aug.head(n_x) = x_;
  //create augmented covariance matrix
  int n_nu = n_aug - n_x;
  MatrixXd Q = MatrixXd(n_nu,n_nu);
  Q << (std_a*std_a),0,
       0, (std_yawdd*std_yawdd);    
  P_aug = MatrixXd::Zero(n_aug,n_aug);
  P_aug.topLeftCorner(n_x, n_x) = P_;
  P_aug.bottomRightCorner(n_nu, n_nu) = Q;
  //create square root matrix
  
  MatrixXd A = P_aug.llt().matrixL();
  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug; i++)
  {
      Xsig_aug.col(i+1) = x_aug + sqrt(lambda+n_aug)*A.col(i);
      Xsig_aug.col(i+n_aug+1) = x_aug - sqrt(lambda+n_aug)*A.col(i);
  }

  *Xsig_out = Xsig_aug;

}

/*Pregict the Sigma points based on the time difference*/
void UKF::PredictSigmaPoints(MatrixXd Xsig_aug, double delta_t)
{

  for (int i=0; i<(2*n_aug + 1);i++)
  {
    float p_x = Xsig_aug(0,i);
    float p_y = Xsig_aug(1,i);
    float v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
}

/*From predicted sigma points calculate the state mean and covariance*/
void UKF::PredictMeanCovariance(VectorXd *x_out, MatrixXd *P_out)
{

   //create vector for predicted state
  VectorXd x = VectorXd(n_x);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x, n_x);

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) 
  {  //iterate over sigma points
    x = x+ weights(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) 
  {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    x_diff(3)=normalizeAngle(x_diff(3));

    P = P + weights(i) * x_diff * x_diff.transpose() ;
  }
  *x_out = x;
  *P_out = P;
}
/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
  GenerateAugSigmaPoints(&Xsig_aug);
  PredictSigmaPoints(Xsig_aug, delta_t);
  VectorXd x_out = VectorXd(n_x);
  MatrixXd P_out = MatrixXd(n_x,n_x);
  PredictMeanCovariance(&x_out, &P_out);
  x_ = x_out;
  P_ = P_out;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  double NIS;
  VectorXd z_pred = VectorXd(2);
  MatrixXd S = MatrixXd(2,2);
  MatrixXd Z_sig = MatrixXd(2, 2*n_aug+1);
  /*Convert from measurement space to state space*/
  PredictLaserMeasurements(&z_pred, &S, &Z_sig);
  /*Update space mean and covariance*/
  UpdateState(meas_package.sensor_type_, z_pred, S, Z_sig, meas_package.raw_measurements_);
  /*Calculate normalized innovative square*/
  NIS=CalculateNIS(z_pred, meas_package.raw_measurements_,S);
  cout << "Laser Measurement NIS " << NIS << endl;
  laser_meas_cnt++;
  if (NIS > 5.99)
    laser_chi_thresh++;
  cout << "Laser count=" << laser_meas_cnt << " NIS above 95%=" << laser_chi_thresh << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  double NIS;
  VectorXd z_pred = VectorXd(3);
  MatrixXd S = MatrixXd(3,3);
  MatrixXd Z_sig = MatrixXd(3, 2*n_aug+1);
  /*Convert from measurement space to state space*/
  PredictRadarMeasurements(&z_pred, &S, &Z_sig);
  /*Update state mean and coveriance*/
  UpdateState(meas_package.sensor_type_, z_pred, S, Z_sig, meas_package.raw_measurements_);
  /*Calculate normalize innovatice square*/
  NIS=CalculateNIS(z_pred, meas_package.raw_measurements_,S);
  cout << "Radar Measurement NIS " << NIS << endl;
  radar_meas_cnt++;
  if (NIS > 7.18)
    radar_chi_thresh++;
  cout << "Radar count=" << radar_meas_cnt << " NIS above 95%=" << radar_chi_thresh << endl;
}

/*Update the State mean and covariance based on the measurements*/
void UKF::UpdateState(MeasurementPackage::SensorType sensor, VectorXd z_pred, MatrixXd S, MatrixXd Zsig, VectorXd z)
{
  //create matrix for cross correlation Tc
  int n_z;
  if (sensor == MeasurementPackage::RADAR)
  {
    n_z = 3;
  }
  else if (sensor == MeasurementPackage::LASER)
  {
    n_z = 2;
  }
  MatrixXd Tc = MatrixXd(n_x, n_z);


  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    if (sensor == MeasurementPackage::RADAR)
    {
      //angle normalization
      z_diff(1)=normalizeAngle(z_diff(1));
    }

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3)=normalizeAngle(x_diff(3));

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  if (sensor == MeasurementPackage::RADAR)
  {
    //angle normalization
    z_diff(1)=normalizeAngle(z_diff(1));
  }

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}

/*Convert from Radar measurement space to state space*/
void UKF::PredictRadarMeasurements(VectorXd* z_out, MatrixXd* S_out, MatrixXd *Zsig_out)
{
  //create matrix for sigma points in measurement space
  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);


  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    if (Zsig(0,i) < 0.0001)
    {
      cout << "PredictRadarMeasurements: Divide by zero" << endl;
      Zsig(0,i) = 0.0001; 
    }
    Zsig(1,i) = atan2(p_y,p_x);                    //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / Zsig(0,i);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug+1; i++) {
      z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    z_diff(1)=normalizeAngle(z_diff(1));

    S = S + weights(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr*std_radr, 0, 0,
          0, std_radphi*std_radphi, 0,
          0, 0,std_radrd*std_radrd;
  S = S + R;

  *z_out = z_pred;
  *S_out = S;
  *Zsig_out = Zsig;
}

/*Predict from  Laset measurement space to state space*/
void UKF::PredictLaserMeasurements(VectorXd* z_out, MatrixXd* S_out, MatrixXd *Zsig_out)
{
  //create matrix for sigma points in measurement space
  int n_z = 2;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);


  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);


    // measurement model
    Zsig(0,i) = p_x;//px
    Zsig(1,i) = p_y;//py
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug+1; i++) {
      z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S = S + weights(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx*std_laspx, 0,
          0, std_laspy*std_laspy;
  S = S + R;

  *z_out = z_pred;
  *S_out = S;
  *Zsig_out = Zsig;
}


double UKF::normalizeAngle(double x) {

  while (x > M_PI) x -=2.*M_PI;
  while (x < -M_PI) x +=2.*M_PI;
  return x;
}


double UKF::CalculateNIS(VectorXd z_pred, VectorXd z, MatrixXd S)
{
  double NIS;

  VectorXd z_diff = z - z_pred;
  NIS = z_diff.transpose()*S.inverse()*z_diff;
  return NIS;
}
