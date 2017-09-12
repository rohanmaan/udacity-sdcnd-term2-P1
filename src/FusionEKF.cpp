#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>



using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  ekf_.x_ = VectorXd(4);
  ekf_.Q_ = MatrixXd(4,4);
  ekf_.F_=MatrixXd(4,4);

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
        ekf_.P_=MatrixXd(4,4);


  
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
    ekf_.F_<< 1,0,1,0,
            0,1,0,1,
          0,0,1,0,
          0,0,0,1;
 
    ekf_.P_<<1,0,0,0,
             0,1,0,0,
             0,0,1000,0,
             0,0,0,1000;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      //"EKF: First State Due to Radar Initialization Started.........."
    	float Ro = measurement_pack.raw_measurements_(0);
    	float Phi =measurement_pack.raw_measurements_(1);
    	float Ro_Dot=measurement_pack.raw_measurements_(2);
      float Px=cos(Phi)*Ro;
      float Py=tan(Phi)*Px;
      ekf_.x_<<Px, Py,0,0;
     // "     x_  Initial State if First Measurement is from a Radar= "<<ekf_.x_
      previous_timestamp_= measurement_pack.timestamp_;


    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      //cout << "EKF: First State Due to Lidar Initialization Started.........." << endl;
      ekf_.x_<<measurement_pack.raw_measurements_(0),measurement_pack.raw_measurements_(1),0,0;
      //cout << "     x_  Initial State if First Measurement is from a Lidar= "<<ekf_.x_ << endl;
      previous_timestamp_= measurement_pack.timestamp_;
    }

    // done initializing, no need to predict or update
    //cout << "EKF: First State Initialization Completed.........." << endl;
    //cout<<"     x_ = "<< ekf_.x_<< endl;
    is_initialized_ = true;
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

    long long  measured_time_stamp=measurement_pack.timestamp_;
    float dt=(measured_time_stamp - previous_timestamp_)/1000000.0;
    previous_timestamp_=measured_time_stamp;
    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;
    //cout<<"Before Starting Prediction Step......."<<endl;
    //cout<<"P_= "<<ekf_.P_<<endl;
    // Process acceleration noise componenets ie the sigma2 square ax and sigma square ay
   float noise_ax=10.0;
   float noise_ay=10.0; 
   ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
               0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
               dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
               0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

//Predict 
   
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
      ekf_.H_=tool.CalculateJacobian(ekf_.x_);
      //cout<<"  Jacobian Matrix Hj_ = "<<ekf_.H_<<endl;
      ekf_.R_=R_radar_;
      float Pred_px=ekf_.x_(0);
      float Pred_py=ekf_.x_(1);
      float Pred_vx=ekf_.x_(2);
      float Pred_vy=ekf_.x_(3);
      float Pred_Ro=sqrt((Pred_px*Pred_px) + (Pred_py*Pred_py) );
      float Pred_Phi=atan2(Pred_py,Pred_px);
      float Pred_RoDot=(Pred_px*Pred_vx + Pred_py * Pred_vy)/Pred_Ro;
      VectorXd Z_pred=VectorXd(3);
      Z_pred<<Pred_Ro,Pred_Phi,Pred_RoDot;
      VectorXd Z =VectorXd(3);
      Z <<measurement_pack.raw_measurements_(0),measurement_pack.raw_measurements_(1),measurement_pack.raw_measurements_(2);
      ekf_.UpdateEKF(Z ,Z_pred);
  } else {
    // Laser updates
       H_laser_<<1,0,0,0,
                 0,1,0,0;
      ekf_.H_=H_laser_;
      ekf_.R_=R_laser_;
      VectorXd Z =VectorXd(2);
      Z <<measurement_pack.raw_measurements_(0),measurement_pack.raw_measurements_(1);
      ekf_.Update(Z);
  }


  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;


}
