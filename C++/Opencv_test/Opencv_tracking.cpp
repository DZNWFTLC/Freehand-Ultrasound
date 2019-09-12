#include <opencv2/opencv.hpp>
#include <opencv2/aruco.hpp>
#include <opencv2/core/matx.hpp>
#include <vector>
#include <opencv2/opencv.hpp>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "iostream"
#include <vector>
#include <stdio.h>
#include <opencv2/imgcodecs.hpp>
#include <string>
#include <math.h>


using namespace cv;
using namespace std;
using namespace aruco;
const double PI = 3.141592653589793238463;


int main(int argc, char** argv)
{
int wait_time = 10;


Mat distortion_coefficients1;
Mat camera_matrix1;
Mat distortion_coefficients2;
Mat camera_matrix2;

Ptr<DetectorParameters> parameter=DetectorParameters::create();
parameter->errorCorrectionRate=0.8f;
parameter->adaptiveThreshWinSizeStep = 3;
parameter->adaptiveThreshWinSizeMax = 16;
parameter->maxMarkerPerimeterRate = 0.5;
parameter->minMarkerPerimeterRate = 0.01;
parameter->cornerRefinementMaxIterations = 50;
parameter->cornerRefinementWinSize = 10;
parameter->cornerRefinementMinAccuracy = 0.05;
FileStorage calibration1("cameraL.yaml", FileStorage::READ);
calibration1["camera_matrix"] >> camera_matrix1;
calibration1["cistortion_coefficients"] >> distortion_coefficients1;
calibration1.release();

FileStorage calibration2("cameraR.yaml", FileStorage::READ);
calibration2["camera_matrix"] >> camera_matrix2;
calibration2["cistortion_coefficients"] >> distortion_coefficients2;
calibration2.release();

VideoCapture video1("..\\Data\\SyncedWithTimestamps\\EndoscopeImageMemory_0_de_test.avi");
VideoCapture video2("..\\Data\\SyncedWithTimestamps\\EndoscopeImageMemory_1_de_test.avi");
int turn = 0;
while (video1.grab() && video2.grab())
{
	Vec3d rvec_sum(0, 0, 0);
	Vec3d tvec_sum(0, 0, 0);
	int num_markers = 0;

	Mat image1, image2, image_distortion, image_copy;
	video1.retrieve(image1);
	video2.retrieve(image2);

	//undistort(image, image_distortion, camera_matrix, distortion_coefficients);

	image1.copyTo(image_copy);
	vector<int> ids1, ids2;
	vector<vector<Point2f> > corners1, corners2;
	Ptr<aruco::Dictionary> dictionary =aruco::generateCustomDictionary(6, 3);
	//Ptr<Dictionary> dictionary = getPredefinedDictionary(DICT_4X4_50);

	Mat T_m = Mat::eye(4, 4, CV_64F); 
	Mat T1_m = Mat::eye(4, 4, CV_64F); 
	Mat Rot_m, Rot_m2;
	Vec3d t4(-0.025, 0.025, -0.056);
	Vec3d t0(-0.275, 0.025, -0.056);
	Mat T1_2 = (Mat_<double>(4, 4) << 1,-0.0003,-0.0055,0.0542,0.0003,1,0.0095,0.0023,0.0055,-0.0095,0.9999,0,0,0,0,1);
	//cout << "T1_2  = " << endl << T1_2 << endl << endl;

	detectMarkers(image1, dictionary, corners1, ids1, parameter);
	//detectMarkers(image2, dictionary, corners2, ids2, parameter);

	// if at least one marker detected
	vector<Vec3d> rvecs, tvecs;
	if (ids1.size() > 0) {
		drawDetectedMarkers(image_copy, corners1, ids1);
		cv::aruco::estimatePoseSingleMarkers(corners1, 0.04, camera_matrix1, distortion_coefficients1, rvecs, tvecs);
		for (int i = 0; i < ids1.size(); i++)
			if (ids1[i] == 0 || ids1[i] == 4) {
				//drawDetectedMarkers(image_copy, corners, ids);

				switch (ids1[i]) {
				case 4: {
					Vec3d tvec(tvecs[i] + t4);
					//cout << "ID  = " << endl << ids1[i] << endl << endl;
					//cout << "new  = " << endl << tvec << endl << endl;
					//cout << "old  = " << endl << tvecs[i] << endl << endl;
					cout << "L1  = " << endl << rvecs[i] << endl << endl;

					cv::aruco::drawAxis(image_copy, camera_matrix1, distortion_coefficients1, rvecs[i], tvec, 0.1);
					rvec_sum = rvec_sum + rvecs[i];
					tvec_sum = tvec_sum + tvec;
					num_markers++;
					break;
				}
				case 5: {
					Vec3d tvec(tvecs[i] + t0);
					//cout << "ID  = " << endl << ids1[i] << endl << endl;
					//cout << "new  = " << endl << tvec << endl << endl;
					//cout << "old  = " << endl << tvecs[i] << endl << endl;
					cout << "L2  = " << endl << rvecs[i] << endl << endl;

					cv::aruco::drawAxis(image_copy, camera_matrix1, distortion_coefficients1, rvecs[i], tvec, 0.1);
					rvec_sum = rvec_sum + rvecs[i];
					tvec_sum = tvec_sum + tvec;
					num_markers++;
					break;
				}
				default: drawAxis(image_copy, camera_matrix1, distortion_coefficients1, rvecs[i], tvecs[i], 0.1);
				}
				/*if (ids[i] == 4) {
					Rodrigues(rvecs[i], Rot_m1);
					Vec3d tvec(tvecs[i] + t4);
					for (int j = 0; j < 3; j++) {
						for (int k = 0; k < 4; k++) {
							if (k < 3)
								T_m1.at<double>(j, k) = Rot_m1.at<double>(j, k);
							else
								T_m1.at<double>(j, k) = tvec[j];
						}
					}
				*/	
			}
	}/*
	if (ids2.size() > 0) {
		//drawDetectedMarkers(image_copy, corners1, ids1);
		cv::aruco::estimatePoseSingleMarkers(corners2, 0.04, camera_matrix2, distortion_coefficients2, rvecs, tvecs);
		for (int i = 0; i < ids2.size(); i++)
			if (ids2[i] == 0 || ids2[i] == 1) {
				Vec3d tvec;
				Rodrigues(rvecs[i], Rot_m);
				for (int j = 0; j < 3; j++) {
					for (int k = 0; k < 4; k++) {
						if (k < 3)
							T_m.at<double>(j, k) = Rot_m.at<double>(j, k);
						else
							T_m.at<double>(j, k) = tvecs[i][j];
					}
				}
				T1_m = T1_2 * T_m;
				for (int j = 0; j < 3; j++) {
					for (int k = 0; k < 4; k++) {
						if (k < 3)
							Rot_m.at<double>(j, k) = T1_m.at<double>(j, k);
						else
							tvec[j] = T1_m.at<double>(j, k);
					}
				}
				switch (ids2[i]) {
				case 0: {
					Vec3d tvec(tvecs[i] + t4);
					//cout << "ID  = " << endl << ids2[i] << endl << endl;
					//cout << "new  = " << endl << tvec << endl << endl;
					//cout << "old  = " << endl << tvecs[i] << endl << endl;
					cout << "R1  = " << endl << rvecs[i] << endl << endl;

					cv::aruco::drawAxis(image_copy, camera_matrix1, distortion_coefficients1, rvecs[i], tvec, 0.1);
					rvec_sum = rvec_sum + rvecs[i];
					tvec_sum = tvec_sum + tvec;
					num_markers++;
					break;
				}
				case 1: {
					Vec3d tvec(tvecs[i] + t0);
					//cout << "ID  = " << endl << ids2[i] << endl << endl;
					//cout << "new  = " << endl << tvec << endl << endl;
					//cout << "old  = " << endl << tvecs[i] << endl << endl;
					cout << "R2  = " << endl << rvecs[i] << endl << endl;

					cv::aruco::drawAxis(image_copy, camera_matrix1, distortion_coefficients1, rvecs[i], tvec, 0.1);
					rvec_sum = rvec_sum + rvecs[i];
					tvec_sum = tvec_sum + tvec;
					num_markers++;
					break;
				}
				default: drawAxis(image_copy, camera_matrix1, distortion_coefficients1, rvecs[i], tvecs[i], 0.1);
				}
			}
	}
	
	if (num_markers ==2) {
		rvec_sum = rvec_sum / num_markers;
		tvec_sum = tvec_sum / num_markers;
		cout << "rvec_sum  = " << endl << rvec_sum << endl << endl;
		drawAxis(image_copy, camera_matrix1, distortion_coefficients1, rvec_sum, tvec_sum, 0.2);
	}*/
	imshow("Detected markers", image_copy);
	char key = (char)waitKey(wait_time);
	if (key == 27)
		break;
}
video1.release();
video2.release();

}