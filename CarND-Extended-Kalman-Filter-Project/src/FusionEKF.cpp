#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
   R_laser_ = Eigen::MatrixXd::Zero(2, 2);
   R_radar_ = Eigen::MatrixXd::Zero(3, 3);
   H_laser_ = Eigen::MatrixXd::Zero(2, 4);
   Hj_ = Eigen::MatrixXd::Zero(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ << 1,0,0,0,
              0,1,0,0;
  Hj_ << 1,1,0,0,
         1,1,0,0,
         1,1,1,1;
  
  //Initializing the state transition matrix 
  MatrixXd F = Eigen::MatrixXd::Zero(4, 4);
  F<< 1,0,1,0,
      0,1,0,1,
      0,0,1,0,
      0,0,0,1;

  //Initialiazing the process covariance matrix
  MatrixXd P = Eigen::MatrixXd::Zero(4, 4);
  P << 1,0,0,0,
       0,1,0,0,
       0,0,1000,0,
       0,0,0,1000;
  // process covariance matrix
  MatrixXd Q = Eigen::MatrixXd(4, 4);
  Q << 0, 0, 0, 0,
       0, 0, 0, 0,
       0, 0, 0, 0,
       0, 0, 0, 0;
  cout<<"Function Calling"<<endl;

  ekf_.Init(P, F, H_laser_, R_laser_, Q);

  cout<<"Done Init"<<endl;



  //Initializing the state vector x

  float noise_ax = 9;
  float noise_ay = 9;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
     cout << "EKF: " << endl;
     ekf_.x_=Eigen::VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
          /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
      if (fabs(ekf_.x_(0)) < 0.0001 and fabs(ekf_.x_(1)) < 0.0001)
      {
        ekf_.x_(0) = 0.0001;
        ekf_.x_(1) = 0.0001;
    }
    }

    else if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

      /**
      Convert radar from polar to cartesian coordinates and initialize state.x  
      */
  

      float rho=measurement_pack.raw_measurements_[0];
      float phi=measurement_pack.raw_measurements_[1];
      float rho_dot=measurement_pack.raw_measurements_[2];

      float x = rho * cos(phi); 
	    float y = rho * sin(phi);
	    float vx = rho_dot * cos(phi);
	    float vy = rho_dot * sin(phi);
      ekf_.x_ << x, y, vx , vy;
    }
    // done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ = measurement_pack.timestamp_;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt   * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  //set the acceleration noise componentsx
  float noise_ax = 9;
  float noise_ay = 9;

  //set the process covariance matrix Q
  ekf_.Q_ = Eigen::MatrixXd(4, 4);
  ekf_.Q_ << dt_4/4 *noise_ax, 0, dt_3/2 *noise_ax,0,
             0, dt_4 / 4 * noise_ay, 0, dt_3/2*noise_ay,
             dt_3 / 2 * noise_ax, 0, dt_2*noise_ax, 0,
             0, dt_3 / 2 * noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    Tools tools;
    ekf_.H_=tools.CalculateJacobian(ekf_.x_);
    ekf_.R_=R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_=H_laser_;
    ekf_.R_=R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
