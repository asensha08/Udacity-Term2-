#include <iostream>
#include "tools.h"
#include <math.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  //check the valid inputs
  VectorXd rmse;
  if(estimations.size()!=ground_truth.size() || estimations.size()==0){
    cout<<"Check size of estimations. Estimation and Ground Truth Vector doesn't match in size";
    return rmse;
  }
  for(unsigned int i=0; i<estimations.size(); i++){
    VectorXd residual= estimations[i]- ground_truth[i];
    residual=residual.array()*residual.array();
    rmse+=residual;

    
	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;

  }
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);

  float px= x_state(0);
  float py= x_state(1);
  float vx= x_state(2);
  float vy= x_state(3);

  float c1= pow(px,2) + pow(py,2);
  float c2= sqrt(c1);
  float c3=pow(c1,1.5);

  if(fabs(c1) < 0.001){
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return Hj;
	}

  Hj<<px/c2,py/c2,0,0,
     -py/c1,px/c1,0,0,
      py*((vx*py)-(vy*px))/c3,px*((vy*px)-(vx*py))/c3,px/c2,py/c2;

  return Hj;

}
