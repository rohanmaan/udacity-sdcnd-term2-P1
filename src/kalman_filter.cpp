#include "kalman_filter.h"
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
	//cout<<"Starting Prediction Step"<<endl;
	x_=F_*x_;
	MatrixXd Ft_=F_.transpose();
	P_=F_*P_*Ft_ + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
	//cout<<"Updating the KF Equations for LIDAR....."<< endl;
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	cout<<"     Error y = "<< y <<endl;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;
	//new estimate
	x_ = x_ + (K * y);
	MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z , const VectorXd &z_pred) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
	//cout<<"Updating the EKF Equations for RADAR..."<< endl;
	VectorXd y = z - z_pred;
	float pi =3.141592;

	while (y[1] > pi)
	    {
	        y[1] -= 2*pi;
	    }

	    while (y[1] < -pi)
	    {
	        y[1] += 2*pi;
	    }

	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
	P_ = (I - K * H_) * P_;
}
