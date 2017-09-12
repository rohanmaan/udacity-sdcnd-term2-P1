#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	//cout<<"Calculating the RMSE for the filter....."<<endl;
	VectorXd rmse(4);
	rmse<<0,0,0,0;
	if (estimations.size()!=ground_truth.size() || estimations.size()==0){
		cout<<"estimations and ground truth data is not valid";
		return rmse;
	}


	for(unsigned int i=0;i<estimations.size();++i){
		VectorXd residual = estimations[i]-ground_truth[i];
		residual=residual.array()*residual.array();
		rmse = residual +rmse;

	}
	rmse=rmse/estimations.size();
	rmse=rmse.array().sqrt();
	cout<<"RMSE = "<<rmse <<endl;
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
	//cout<<"Calculating Jacobian Matrix for EKF radar....."<<endl;
	MatrixXd Hj(3,4);
	float p_x=x_state(0);
	float p_y=x_state(1);
	float vx=x_state(2);
	float vy=x_state(3);

	float c1=p_x*p_x + p_y*p_y;
	float c2=pow(c1,0.5);
	float c3=pow(c1,1.5);

	if(fabs(c1) < 0.0001){
			cout << "CalculateJacobian () - Error - Division by Zero" << endl;
			return Hj;
		}
	Hj <<	p_x/c2,p_y/c2,0,0,
			-p_y/c1,p_x/c1,0,0,
			p_y*(vx*p_y - vy*p_x)/c3,p_x*(vy*p_x-vx*p_y)/c3,p_x/c2,p_y/c2;

	return Hj;

}
