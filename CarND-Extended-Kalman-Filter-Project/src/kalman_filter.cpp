#include "kalman_filter.h"
#include <math.h> //for math operations


// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init( Eigen::MatrixXd &P_in, Eigen::MatrixXd &F_in,
                        Eigen::MatrixXd &H_in, Eigen::MatrixXd &R_in, Eigen::MatrixXd &Q_in) {
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_=F_*x_;
  Eigen::MatrixXd Ft = F_.transpose();
  P_=F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  Eigen::VectorXd y=z-H_*x_;
  Eigen::MatrixXd Ht = H_.transpose();
  Eigen::MatrixXd S=H_*P_*Ht + R_;
  Eigen::MatrixXd S_inv = S.inverse();
  Eigen::MatrixXd PHt = P_ * Ht;
  Eigen::MatrixXd K = PHt *S_inv;
  /** New Measurement update**/
  x_=x_ + (K*y);
  long x_size = x_.size();
	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations

  */
  float rho= sqrt(pow(x_(0),2) + pow(x_(1),2));
  float phi= atan2(x_(1),x_(0));
  float rho_dot;
  if(fabs(rho)<=0.0001){
    rho_dot=0;
  }else{
    rho_dot=(x_(0)*x_(2) + x_(1)*x_(3))/rho;
  }

  //Using the above equations to to calulate z prediction
  Eigen::VectorXd z_pred(3);
  z_pred<<rho,phi,rho_dot;

  Eigen::VectorXd y= z-z_pred;
  Eigen::MatrixXd Ht = H_.transpose();
  Eigen::MatrixXd S=H_*P_*Ht + R_;
  Eigen::MatrixXd S_inv = S.inverse();
  Eigen::MatrixXd PHt = P_ * Ht;
  Eigen::MatrixXd K = PHt *S_inv;
  /** New Measurement update**/
  x_=x_ + (K*y);
  long x_size = x_.size();
	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

}
