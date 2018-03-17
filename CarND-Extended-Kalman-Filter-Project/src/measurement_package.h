#ifndef MEASUREMENT_PACKAGE_H_
#define MEASUREMENT_PACKAGE_H_

#include "Eigen/Dense"

using Eigen::VectorXd;

class MeasurementPackage {
public:
  long long timestamp_;

  enum SensorType{LASER,RADAR} sensor_type_;
  VectorXd raw_measurements_;
};

#endif /* MEASUREMENT_PACKAGE_H_ */
