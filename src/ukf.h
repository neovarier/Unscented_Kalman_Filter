#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;


  ///* previous timestamp
  long long prev_timestamp_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd;

  ///* Weights of sigma points
  VectorXd weights;

  ///* State dimension
  int n_x;

  ///* Augmented state dimension
  int n_aug;

  ///* Sigma point spreading parameter
  double lambda;

  /*Total number of laser measurements*/
  int laser_meas_cnt = 0;
  /*Total number of radar measurements*/
  int radar_meas_cnt = 0;
  /*Number of cases where Laser NIS is above 95%*/
  int laser_chi_thresh = 0;
  /*Number of cases where Radar NIS is above 95%*/
  int radar_chi_thresh = 0;



  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  void GenerateAugSigmaPoints(MatrixXd* Xsig_out);

  void PredictSigmaPoints(MatrixXd Xsig_aug, double delta_t);

  void PredictMeanCovariance(VectorXd *x_out, MatrixXd *P_out);

  void UpdateState(MeasurementPackage::SensorType sensor, VectorXd z_pred, MatrixXd S, MatrixXd Zsig, VectorXd z);
  
  void PredictRadarMeasurements(VectorXd* z_out, MatrixXd* S_out, MatrixXd *Zsig_out);

  void PredictLaserMeasurements(VectorXd* z_out, MatrixXd* S_out, MatrixXd *Zsig_out);

  double CalculateNIS(VectorXd z_pred, VectorXd z, MatrixXd S);

  double normalizeAngle(double x);
};

#endif /* UKF_H */
